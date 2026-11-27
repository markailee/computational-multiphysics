#include "neutron/NeutronSolver.hpp"

#include "fvm/TecplotWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace mp {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

std::vector<std::string> strings(const std::string& text) {
    std::istringstream input(text);
    std::vector<std::string> result;
    std::string value;
    while (input >> value) result.push_back(value);
    return result;
}

std::vector<double> numbers(const std::string& text, int count, const std::string& key) {
    std::istringstream input(text);
    std::vector<double> result;
    double value = 0.0;
    while (input >> value) result.push_back(value);
    if (!input.eof()) throw std::runtime_error("Invalid numeric list for setting: " + key);
    if (result.size() == 1 && count > 1) result.resize(count, result.front());
    if (static_cast<int>(result.size()) != count) {
        throw std::runtime_error("Setting " + key + " requires one value or " + std::to_string(count) + " values");
    }
    for (double entry : result) {
        if (!std::isfinite(entry)) throw std::runtime_error("Non-finite value in setting: " + key);
    }
    return result;
}

std::vector<double> variableNumbers(const std::string& text, const std::string& key) {
    std::istringstream input(text);
    std::vector<double> result;
    double value = 0.0;
    while (input >> value) result.push_back(value);
    if (!input.eof() || result.empty()) throw std::runtime_error("Invalid numeric list for setting: " + key);
    for (double entry : result) {
        if (!std::isfinite(entry)) throw std::runtime_error("Non-finite value in setting: " + key);
    }
    return result;
}

std::vector<double> optionalNumbers(const Config& config, const std::string& key, int count, double fallback) {
    return config.has(key) ? numbers(config.require(key), count, key) : std::vector<double>(count, fallback);
}

std::vector<double> optionalMatrix(const Config& config, const std::string& key, int groups) {
    if (!config.has(key)) return std::vector<double>(groups * groups, 0.0);
    return numbers(config.require(key), groups * groups, key);
}

void requireNonnegative(const std::vector<double>& values, const std::string& key) {
    for (double value : values) {
        if (value < 0.0) throw std::runtime_error("Setting " + key + " must be non-negative");
    }
}

double robinGamma(const NeutronBoundaryCondition& condition, int group) {
    if (condition.kind == NeutronBoundaryKind::Vacuum) return 0.5;
    if (condition.kind == NeutronBoundaryKind::Marshak) return condition.value(group);
    if (condition.kind == NeutronBoundaryKind::Albedo) {
        const double beta = condition.value(group);
        return 0.5 * (1.0 - beta) / (1.0 + beta);
    }
    return 0.0;
}

std::string csvNumber(double value) {
    std::ostringstream output;
    output << std::setprecision(16) << value;
    return output.str();
}

} // namespace

NeutronSolver::NeutronSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory)
    : mesh_(mesh), config_(config), caseDirectory_(std::move(caseDirectory)),
      cellCount_(static_cast<int>(mesh.cells.size())), boundaries_(mesh.faces.size()) {}

void NeutronSolver::initialize() {
    if (initialized_) return;
    readConfiguration();
    if (model_ == NeutronModel::DiscreteOrdinates) {
        buildQuadrature();
        buildSnSweepOrder();
        system_ = std::make_unique<PetscLinearSystem>(cellCount_);
    } else {
        system_ = std::make_unique<PetscLinearSystem>(groupCount_ * cellCount_);
    }
    buildBoundaries();
    if (model_ == NeutronModel::DiscreteOrdinates) {
        buildReflectionMaps();
        outgoingBoundaryAngular_.assign(
            static_cast<std::size_t>(groupCount_) * directionCount_ * mesh_.boundaryFaceCount(), 0.0);
        if (snAcceleration_ == NeutronSnAcceleration::DSA) buildDsaOperators();
    } else {
        buildDiffusionOperator();
    }
    scalarFlux_.assign(groupCount_ * cellCount_, initialFlux_);
    currentX_.assign(groupCount_ * cellCount_, 0.0);
    currentY_.assign(groupCount_ * cellCount_, 0.0);
    if (model_ == NeutronModel::DiscreteOrdinates) {
        angularFlux_.assign(groupCount_ * directionCount_ * cellCount_, initialFlux_ / (4.0 * pi));
    }
    initialized_ = true;
}

void NeutronSolver::readConfiguration() {
    const std::string model = lowercase(config_.get("neutron_model", "sn"));
    if (model == "sn" || model == "transport" || model == "discrete_ordinates") {
        model_ = NeutronModel::DiscreteOrdinates;
    } else if (model == "diffusion") {
        model_ = NeutronModel::Diffusion;
    } else {
        throw std::runtime_error("Unknown neutron_model: " + model);
    }

    const std::string problem = lowercase(config_.get("problem_type", "eigenvalue"));
    if (problem == "fixed_source" || problem == "fixed" || problem == "source") {
        problem_ = NeutronProblem::FixedSource;
    } else if (problem == "eigenvalue" || problem == "eigen" || problem == "criticality") {
        problem_ = NeutronProblem::Eigenvalue;
    } else if (problem == "transient" || problem == "kinetics" || problem == "time_dependent") {
        problem_ = NeutronProblem::Transient;
    } else {
        throw std::runtime_error("Unknown neutron problem_type: " + problem);
    }

    groupCount_ = config_.integer("energy_groups", 1);
    if (groupCount_ <= 0) throw std::runtime_error("energy_groups must be positive");

    snOrder_ = config_.integer("sn_order", 0);
    polarCount_ = config_.integer("sn_polar", 2);
    azimuthCount_ = config_.integer("sn_azimuth", 4);
    const std::string defaultQuad = (snOrder_ > 0) ? "level_symmetric" : "product_midpoint";
    const std::string quad = lowercase(config_.get("sn_quadrature", defaultQuad));
    if (quad == "product_midpoint" || quad == "midpoint" || quad == "chebyshev_legendre") {
        quadratureKind_ = NeutronQuadratureKind::ProductMidpoint;
    } else if (quad == "product_axis" || quad == "axis") {
        quadratureKind_ = NeutronQuadratureKind::ProductAxis;
    } else if (quad == "level_symmetric" || quad == "lq" || quad == "ls") {
        quadratureKind_ = NeutronQuadratureKind::LevelSymmetric;
        if (snOrder_ <= 0) snOrder_ = 4;
    } else {
        throw std::runtime_error("Unknown sn_quadrature: " + quad);
    }

    if (model_ == NeutronModel::DiscreteOrdinates && quadratureKind_ != NeutronQuadratureKind::LevelSymmetric &&
        (polarCount_ < 2 || polarCount_ % 2 != 0 || azimuthCount_ < 4 || azimuthCount_ % 4 != 0)) {
        throw std::runtime_error("sn_polar must be even and at least 2; sn_azimuth must be a multiple of 4");
    }

    const std::string scheme = lowercase(config_.get("sn_scheme", "upwind"));
    if (scheme == "upwind" || scheme == "step" || scheme == "first_order") {
        snScheme_ = NeutronSnScheme::Upwind;
    } else if (scheme == "second_order" || scheme == "muscl" || scheme == "linear") {
        snScheme_ = NeutronSnScheme::SecondOrder;
    } else {
        throw std::runtime_error("Unknown sn_scheme: " + scheme);
    }
    snLimiter_ = lowercase(config_.get("sn_limiter", "barth_jespersen"));

    const std::string accel = lowercase(config_.get("sn_acceleration", "none"));
    if (accel == "none" || accel == "off" || accel == "false") {
        snAcceleration_ = NeutronSnAcceleration::None;
    } else if (accel == "dsa" || accel == "diffusion") {
        snAcceleration_ = NeutronSnAcceleration::DSA;
    } else {
        throw std::runtime_error("Unknown sn_acceleration: " + accel);
    }

    maximumIterations_ = config_.integer("maximum_iterations", 100);
    maximumSourceIterations_ = config_.integer("maximum_source_iterations", 100);
    maximumLinearIterations_ = config_.integer("maximum_linear_iterations", 1000);
    eigenModesRequested_ = config_.integer("eigen_modes", 1);
    maximumModeIterations_ = config_.integer("maximum_mode_iterations", 150);
    logInterval_ = config_.integer("log_interval", 1);
    transientOutputInterval_ = config_.integer("transient_output_interval", 1);

    tolerance_ = config_.number("tolerance", 1.0e-8);
    eigenTolerance_ = config_.number("eigen_tolerance", 1.0e-10);
    sourceTolerance_ = config_.number("source_tolerance", 1.0e-9);
    balanceTolerance_ = config_.number("balance_tolerance", 1.0e-8);
    modeTolerance_ = config_.number("mode_tolerance", 1.0e-8);
    linearTolerance_ = config_.number("linear_tolerance", 1.0e-10);
    linearAbsoluteTolerance_ = config_.number("linear_absolute_tolerance", 1.0e-12);
    initialFlux_ = config_.number("initial_flux", 1.0);
    targetPower_ = config_.number("target_power", 0.0);
    axialBuckling_ = config_.number("axial_buckling", 0.0);
    computeAdjoint_ = config_.boolean("compute_adjoint", false);
    timeStep_ = config_.number("time_step", 0.01);
    endTime_ = config_.number("end_time", 0.5);
    criticalInitialState_ = config_.boolean("critical_initial_state", true);

    const std::string timeScheme = lowercase(config_.get("time_scheme", "backward_euler"));
    if (timeScheme == "backward_euler" || timeScheme == "implicit_euler" || timeScheme == "be") {
        timeScheme_ = NeutronTimeScheme::BackwardEuler;
    } else if (timeScheme == "crank_nicolson" || timeScheme == "cn") {
        timeScheme_ = NeutronTimeScheme::CrankNicolson;
    } else {
        throw std::runtime_error("Unknown time_scheme: " + timeScheme);
    }

    linearSolver_ = config_.get("linear_solver", "gmres");
    preconditioner_ = config_.get("preconditioner", "ilu");
    if (maximumIterations_ <= 0 || maximumSourceIterations_ <= 0 || maximumLinearIterations_ <= 0 ||
        eigenModesRequested_ <= 0 || maximumModeIterations_ <= 0 || logInterval_ <= 0 ||
        transientOutputInterval_ <= 0 || tolerance_ <= 0.0 || eigenTolerance_ <= 0.0 || sourceTolerance_ <= 0.0 ||
        balanceTolerance_ <= 0.0 || modeTolerance_ <= 0.0 || linearTolerance_ <= 0.0 ||
        linearAbsoluteTolerance_ <= 0.0 || initialFlux_ < 0.0 || targetPower_ < 0.0 || axialBuckling_ < 0.0 ||
        timeStep_ <= 0.0 || endTime_ <= 0.0) {
        throw std::runtime_error("Invalid neutron iteration, tolerance, or normalization setting");
    }
    if ((problem_ == NeutronProblem::Eigenvalue || problem_ == NeutronProblem::Transient) && initialFlux_ <= 0.0) {
        throw std::runtime_error("A positive initial_flux is required for neutron eigenvalue and transient problems");
    }
    if (problem_ == NeutronProblem::FixedSource && targetPower_ > 0.0) {
        throw std::runtime_error("target_power normalization is only available for eigenvalue or transient calculations");
    }

    // Kinetics / delayed neutron parameters (optional)
    if (config_.has("neutron_velocity")) {
        neutronVelocity_ = numbers(config_.require("neutron_velocity"), groupCount_, "neutron_velocity");
        for (double v : neutronVelocity_) {
            if (v <= 0.0) throw std::runtime_error("neutron_velocity must be positive");
        }
    }
    if (config_.has("delayed_beta") || config_.has("delayed_lambda")) {
        delayedBeta_ = variableNumbers(config_.require("delayed_beta"), "delayed_beta");
        delayedLambda_ = numbers(config_.require("delayed_lambda"), static_cast<int>(delayedBeta_.size()), "delayed_lambda");
        requireNonnegative(delayedBeta_, "delayed_beta");
        const double betaSum = std::accumulate(delayedBeta_.begin(), delayedBeta_.end(), 0.0);
        if (betaSum >= 1.0) throw std::runtime_error("Sum of delayed_beta must be less than 1");
        for (double lam : delayedLambda_) {
            if (lam <= 0.0) throw std::runtime_error("delayed_lambda must be positive");
        }
        const int delayedCount = static_cast<int>(delayedBeta_.size());
        if (config_.has("delayed_chi")) {
            const std::vector<double> raw = variableNumbers(config_.require("delayed_chi"), "delayed_chi");
            if (static_cast<int>(raw.size()) == groupCount_) {
                delayedChi_.resize(delayedCount * groupCount_);
                for (int i = 0; i < delayedCount; ++i) {
                    for (int g = 0; g < groupCount_; ++g) delayedChi_[i * groupCount_ + g] = raw[g];
                }
            } else if (static_cast<int>(raw.size()) == delayedCount * groupCount_) {
                delayedChi_ = raw;
            } else {
                throw std::runtime_error("delayed_chi must have G or I*G entries");
            }
            for (int i = 0; i < delayedCount; ++i) {
                double sum = 0.0;
                for (int g = 0; g < groupCount_; ++g) {
                    if (delayedChi_[i * groupCount_ + g] < 0.0) throw std::runtime_error("delayed_chi must be non-negative");
                    sum += delayedChi_[i * groupCount_ + g];
                }
                if (sum <= 0.0) throw std::runtime_error("Each delayed_chi spectrum must have positive sum");
                for (int g = 0; g < groupCount_; ++g) delayedChi_[i * groupCount_ + g] /= sum;
            }
        }
    }

    // Materials
    const std::vector<std::string> materialNames = strings(config_.get("materials", "default"));
    if (materialNames.empty()) throw std::runtime_error("Neutron materials cannot be empty");
    std::set<std::string> uniqueNames;
    for (const std::string& name : materialNames) {
        if (!uniqueNames.insert(name).second) throw std::runtime_error("Duplicate neutron material: " + name);
    }

    hasP1Scattering_ = false;
    materials_.clear();
    for (const std::string& name : materialNames) {
        const std::string prefix = "material." + name + ".";
        NeutronMaterial material;
        material.name = name;
        material.scattering = optionalMatrix(config_, prefix + "scattering", groupCount_);
        material.scatteringP1 = optionalMatrix(config_, prefix + "scattering_p1", groupCount_);
        requireNonnegative(material.scattering, prefix + "scattering");
        for (double v : material.scatteringP1) {
            if (std::abs(v) > 0.0) hasP1Scattering_ = true;
        }

        std::vector<double> outgoing(groupCount_, 0.0);
        std::vector<double> outOfGroup(groupCount_, 0.0);
        for (int source = 0; source < groupCount_; ++source) {
            for (int destination = 0; destination < groupCount_; ++destination) {
                const double s = material.scattering[destination * groupCount_ + source];
                outgoing[source] += s;
                if (destination != source) outOfGroup[source] += s;
            }
        }

        const bool hasTotal = config_.has(prefix + "sigma_total");
        const bool hasAbs = config_.has(prefix + "sigma_absorption");
        const bool hasRem = config_.has(prefix + "sigma_removal");

        if (hasTotal) {
            material.total = numbers(config_.require(prefix + "sigma_total"), groupCount_, prefix + "sigma_total");
            requireNonnegative(material.total, prefix + "sigma_total");
            if (hasAbs) {
                material.absorption = numbers(config_.require(prefix + "sigma_absorption"), groupCount_, prefix + "sigma_absorption");
            } else {
                material.absorption.resize(groupCount_);
                for (int g = 0; g < groupCount_; ++g) material.absorption[g] = material.total[g] - outgoing[g];
            }
        } else if (hasAbs) {
            material.absorption = numbers(config_.require(prefix + "sigma_absorption"), groupCount_, prefix + "sigma_absorption");
            material.total.resize(groupCount_);
            for (int g = 0; g < groupCount_; ++g) material.total[g] = material.absorption[g] + outgoing[g];
        } else if (hasRem) {
            const std::vector<double> rem = numbers(config_.require(prefix + "sigma_removal"), groupCount_, prefix + "sigma_removal");
            requireNonnegative(rem, prefix + "sigma_removal");
            material.absorption.resize(groupCount_);
            material.total.resize(groupCount_);
            for (int g = 0; g < groupCount_; ++g) {
                material.absorption[g] = rem[g] - outOfGroup[g];
                material.total[g] = rem[g] + material.scattering[g * groupCount_ + g];
            }
        } else {
            throw std::runtime_error("Material " + name + " requires sigma_total, sigma_absorption, or sigma_removal");
        }

        requireNonnegative(material.absorption, prefix + "sigma_absorption");
        for (int g = 0; g < groupCount_; ++g) {
            if (model_ == NeutronModel::Diffusion && !config_.has(prefix + "diffusion_coefficient") && material.total[g] <= 0.0) {
                throw std::runtime_error("Diffusion material " + name + " requires positive sigma_total or explicit diffusion_coefficient");
            }
            const double expected = material.absorption[g] + outgoing[g];
            if (std::abs(expected - material.total[g]) > 1.0e-8 * std::max(1.0, material.total[g])) {
                throw std::runtime_error("Neutron total cross section must equal absorption plus outgoing scattering in material " + name);
            }
        }

        const bool hasNuSigmaF = config_.has(prefix + "nu_sigma_fission");
        const bool hasKappaSigmaF = config_.has(prefix + "kappa_sigma_fission");
        if (hasNuSigmaF) {
            material.nuSigmaFission = numbers(config_.require(prefix + "nu_sigma_fission"), groupCount_, prefix + "nu_sigma_fission");
            requireNonnegative(material.nuSigmaFission, prefix + "nu_sigma_fission");
            if (config_.has(prefix + "sigma_fission")) {
                material.fission = numbers(config_.require(prefix + "sigma_fission"), groupCount_, prefix + "sigma_fission");
                requireNonnegative(material.fission, prefix + "sigma_fission");
            } else {
                const std::vector<double> nu = optionalNumbers(config_, prefix + "nu", groupCount_, 2.43);
                material.fission.resize(groupCount_);
                for (int g = 0; g < groupCount_; ++g) {
                    material.fission[g] = (nu[g] > 0.0) ? material.nuSigmaFission[g] / nu[g] : 0.0;
                }
            }
            if (hasKappaSigmaF) {
                material.kappaSigmaFission = numbers(config_.require(prefix + "kappa_sigma_fission"), groupCount_, prefix + "kappa_sigma_fission");
                requireNonnegative(material.kappaSigmaFission, prefix + "kappa_sigma_fission");
            } else {
                const std::vector<double> energyRelease = optionalNumbers(config_, prefix + "energy_release", groupCount_, 1.0);
                requireNonnegative(energyRelease, prefix + "energy_release");
                material.kappaSigmaFission.resize(groupCount_);
                for (int g = 0; g < groupCount_; ++g) {
                    material.kappaSigmaFission[g] = config_.has(prefix + "energy_release")
                        ? energyRelease[g] * material.fission[g]
                        : material.nuSigmaFission[g];
                }
            }
        } else {
            material.fission = optionalNumbers(config_, prefix + "sigma_fission", groupCount_, 0.0);
            requireNonnegative(material.fission, prefix + "sigma_fission");
            for (int g = 0; g < groupCount_; ++g) {
                if (material.fission[g] > material.absorption[g] + 1.0e-12) {
                    throw std::runtime_error("Neutron sigma_fission cannot exceed sigma_absorption");
                }
            }
            const bool fissile = std::any_of(material.fission.begin(), material.fission.end(), [](double v) { return v > 0.0; });
            std::vector<double> nu = fissile
                ? numbers(config_.require(prefix + "nu"), groupCount_, prefix + "nu")
                : optionalNumbers(config_, prefix + "nu", groupCount_, 0.0);
            requireNonnegative(nu, prefix + "nu");
            const std::vector<double> energyRelease = optionalNumbers(config_, prefix + "energy_release", groupCount_, 1.0);
            requireNonnegative(energyRelease, prefix + "energy_release");
            material.nuSigmaFission.resize(groupCount_);
            material.kappaSigmaFission.resize(groupCount_);
            for (int g = 0; g < groupCount_; ++g) {
                material.nuSigmaFission[g] = nu[g] * material.fission[g];
                material.kappaSigmaFission[g] = hasKappaSigmaF
                    ? numbers(config_.require(prefix + "kappa_sigma_fission"), groupCount_, prefix + "kappa_sigma_fission")[g]
                    : energyRelease[g] * material.fission[g];
            }
        }

        const bool fissile = std::any_of(material.nuSigmaFission.begin(), material.nuSigmaFission.end(),
                                         [](double v) { return v > 0.0; });
        material.chi = fissile ? optionalNumbers(config_, prefix + "chi", groupCount_, 0.0)
                               : std::vector<double>(groupCount_, 0.0);
        if (fissile && !config_.has(prefix + "chi")) material.chi.front() = 1.0;
        requireNonnegative(material.chi, prefix + "chi");
        const double chiSum = std::accumulate(material.chi.begin(), material.chi.end(), 0.0);
        if (fissile) {
            if (chiSum <= 0.0) throw std::runtime_error("Fissile neutron material requires a nonzero chi spectrum");
            for (double& value : material.chi) value /= chiSum;
        } else if (chiSum > 0.0) {
            throw std::runtime_error("Non-fissile neutron material cannot have a nonzero chi spectrum");
        }

        material.diffusion = config_.has(prefix + "diffusion_coefficient")
            ? numbers(config_.require(prefix + "diffusion_coefficient"), groupCount_, prefix + "diffusion_coefficient")
            : std::vector<double>(groupCount_, 0.0);
        for (int g = 0; g < groupCount_; ++g) {
            if (!config_.has(prefix + "diffusion_coefficient")) {
                const double sigmaTr = material.total[g] - material.scatteringP1[g * groupCount_ + g];
                material.diffusion[g] = (sigmaTr > 1.0e-12) ? (1.0 / (3.0 * sigmaTr)) : 0.0;
            }
            if (model_ == NeutronModel::Diffusion && material.diffusion[g] <= 0.0) {
                throw std::runtime_error("Neutron diffusion coefficients must be positive");
            }
        }
        material.externalSource = optionalNumbers(config_, prefix + "external_source", groupCount_, 0.0);
        requireNonnegative(material.externalSource, prefix + "external_source");
        materials_.push_back(std::move(material));
    }

    // Regions (rectangular bounds or circular pins)
    const std::string defaultName = config_.get("default_material", materialNames.front());
    const auto defaultIterator = std::find(materialNames.begin(), materialNames.end(), defaultName);
    if (defaultIterator == materialNames.end()) throw std::runtime_error("Unknown default_material: " + defaultName);
    cellMaterials_.assign(cellCount_, static_cast<int>(defaultIterator - materialNames.begin()));
    cellRegions_.assign(cellCount_, -1);

    regions_.clear();
    for (const std::string& name : strings(config_.get("regions", ""))) {
        const std::string prefix = "region." + name + ".";
        const std::string materialName = config_.require(prefix + "material");
        const auto iterator = std::find(materialNames.begin(), materialNames.end(), materialName);
        if (iterator == materialNames.end()) throw std::runtime_error("Unknown material in region " + name + ": " + materialName);
        NeutronRegion region;
        region.name = name;
        region.material = static_cast<int>(iterator - materialNames.begin());
        if (config_.has(prefix + "circle")) {
            const std::vector<double> circle = numbers(config_.require(prefix + "circle"), 3, prefix + "circle");
            if (circle[2] <= 0.0) throw std::runtime_error("Circle radius must be positive in region " + name);
            region.isCircle = true;
            region.centre = {circle[0], circle[1]};
            region.radius = circle[2];
        } else {
            const std::vector<double> bounds = numbers(config_.require(prefix + "bounds"), 4, prefix + "bounds");
            if (bounds[1] <= bounds[0] || bounds[3] <= bounds[2]) throw std::runtime_error("Invalid bounds in region " + name);
            region.isCircle = false;
            region.lower = {bounds[0], bounds[2]};
            region.upper = {bounds[1], bounds[3]};
        }
        regions_.push_back(region);
    }

    for (int cell = 0; cell < cellCount_; ++cell) {
        const Vector2& centre = mesh_.cells[cell].centre;
        for (int rIdx = 0; rIdx < static_cast<int>(regions_.size()); ++rIdx) {
            const NeutronRegion& region = regions_[rIdx];
            bool inside = false;
            if (region.isCircle) {
                const double dx = centre.x - region.centre.x;
                const double dy = centre.y - region.centre.y;
                inside = (dx * dx + dy * dy <= region.radius * region.radius);
            } else {
                inside = (centre.x >= region.lower.x && centre.x <= region.upper.x &&
                          centre.y >= region.lower.y && centre.y <= region.upper.y);
            }
            if (inside) {
                cellMaterials_[cell] = region.material;
                cellRegions_[cell] = rIdx;
            }
        }
    }

    // First-order perturbation configuration (optional)
    if (config_.has("perturbation_delta_absorption") || config_.has("perturbation_delta_nu_sigma_fission") ||
        config_.has("perturbation_delta_diffusion") || config_.has("perturbation_delta_scattering")) {
        perturbation_.enabled = true;
        perturbation_.regionName = config_.get("perturbation_region", "");
        perturbation_.materialName = config_.get("perturbation_material", "");
        perturbation_.deltaAbsorption = optionalNumbers(config_, "perturbation_delta_absorption", groupCount_, 0.0);
        perturbation_.deltaNuSigmaFission = optionalNumbers(config_, "perturbation_delta_nu_sigma_fission", groupCount_, 0.0);
        perturbation_.deltaDiffusion = optionalNumbers(config_, "perturbation_delta_diffusion", groupCount_, 0.0);
        perturbation_.deltaScattering = optionalMatrix(config_, "perturbation_delta_scattering", groupCount_);
        computeAdjoint_ = true;
    }

    // Transient perturbation events (optional)
    transientEvents_.clear();
    for (const std::string& eventName : strings(config_.get("transient_events", ""))) {
        const std::string prefix = "event." + eventName + ".";
        NeutronTransientEvent ev;
        ev.name = eventName;
        ev.regionName = config_.get(prefix + "region", "");
        ev.materialName = config_.get(prefix + "material", "");
        const std::string type = lowercase(config_.get(prefix + "type", "step"));
        if (type == "step") {
            ev.isRamp = false;
        } else if (type == "ramp") {
            ev.isRamp = true;
        } else {
            throw std::runtime_error("Unknown transient event type: " + type);
        }
        ev.startTime = config_.number(prefix + "start_time", 0.0);
        ev.endTime = config_.number(prefix + "end_time", ev.startTime);
        if (ev.isRamp && ev.endTime <= ev.startTime) {
            throw std::runtime_error("Ramp event " + eventName + " requires end_time > start_time");
        }
        ev.deltaAbsorption = optionalNumbers(config_, prefix + "delta_absorption", groupCount_, 0.0);
        ev.deltaNuSigmaFission = optionalNumbers(config_, prefix + "delta_nu_sigma_fission", groupCount_, 0.0);
        ev.deltaDiffusion = optionalNumbers(config_, prefix + "delta_diffusion", groupCount_, 0.0);
        ev.deltaScattering = optionalMatrix(config_, prefix + "delta_scattering", groupCount_);
        transientEvents_.push_back(std::move(ev));
    }
}

void NeutronSolver::buildBoundaries() {
    for (const BoundaryPatch& patch : mesh_.patches) {
        const std::string key = "boundary." + patch.name + ".neutron";
        const std::vector<std::string> words = strings(config_.require(key));
        if (words.empty()) throw std::runtime_error("Empty neutron boundary on patch " + patch.name);
        const std::string kind = lowercase(words.front());
        NeutronBoundaryCondition condition;

        auto parseValues = [&](double fallback, bool required) {
            if (words.size() == 1) {
                if (required) throw std::runtime_error("Neutron boundary requires value(s) on patch " + patch.name);
                return std::vector<double>(groupCount_, fallback);
            }
            std::ostringstream joined;
            for (std::size_t i = 1; i < words.size(); ++i) {
                if (i > 1) joined << ' ';
                joined << words[i];
            }
            return numbers(joined.str(), groupCount_, key);
        };

        const std::string whiteKey = "boundary." + patch.name + ".white_fraction";
        if (config_.has(whiteKey)) {
            condition.whiteFraction = config_.number(whiteKey, 0.0);
            if (condition.whiteFraction < 0.0 || condition.whiteFraction > 1.0) {
                throw std::runtime_error("white_fraction must be in [0,1] on patch " + patch.name);
            }
        }

        if (kind == "vacuum") {
            condition.kind = NeutronBoundaryKind::Vacuum;
            condition.values = parseValues(0.0, false);
        } else if (kind == "reflective" || kind == "reflection") {
            condition.kind = NeutronBoundaryKind::Reflective;
            condition.values = parseValues(1.0, false);
        } else if (kind == "white_reflective" || kind == "white") {
            condition.kind = NeutronBoundaryKind::WhiteReflective;
            condition.values = parseValues(1.0, false);
            condition.whiteFraction = 1.0;
            for (double v : condition.values) {
                if (v < 0.0 || v > 1.0) throw std::runtime_error("white_reflective albedo must be in [0,1]");
            }
        } else if (kind == "albedo") {
            condition.kind = NeutronBoundaryKind::Albedo;
            condition.values = parseValues(1.0, true);
            for (double v : condition.values) {
                if (v < 0.0 || v > 1.0) throw std::runtime_error("Neutron albedo must be in [0,1]");
            }
        } else if (kind == "incident" || kind == "incident_flux") {
            condition.kind = NeutronBoundaryKind::IncidentFlux;
            condition.values = parseValues(0.0, true);
            if (model_ != NeutronModel::DiscreteOrdinates) {
                throw std::runtime_error("incident_flux requires the SN model");
            }
            requireNonnegative(condition.values, key);
        } else if (kind == "incident_current") {
            condition.kind = NeutronBoundaryKind::IncidentCurrent;
            condition.values = parseValues(0.0, true);
            if (model_ != NeutronModel::Diffusion) {
                throw std::runtime_error("incident_current requires the diffusion model");
            }
            requireNonnegative(condition.values, key);
        } else if (kind == "fixed_flux" || kind == "fixedvalue") {
            condition.kind = NeutronBoundaryKind::FixedScalarFlux;
            condition.values = parseValues(0.0, true);
            if (model_ != NeutronModel::Diffusion) throw std::runtime_error("fixed_flux is a diffusion-only boundary");
        } else if (kind == "marshak" || kind == "robin") {
            condition.kind = NeutronBoundaryKind::Marshak;
            condition.values = parseValues(0.5, false);
            if (model_ != NeutronModel::Diffusion) throw std::runtime_error("marshak/robin is a diffusion-only boundary");
            for (double v : condition.values) {
                if (v < 0.0) throw std::runtime_error("marshak/robin coefficient must be non-negative");
            }
        } else if (kind == "zero_current" || kind == "zerogradient") {
            condition.kind = NeutronBoundaryKind::ZeroCurrent;
            condition.values = parseValues(0.0, false);
            if (model_ != NeutronModel::Diffusion) throw std::runtime_error("Use reflective, not zero_current, for SN boundaries");
        } else {
            throw std::runtime_error("Unknown neutron boundary '" + kind + "' on patch " + patch.name);
        }

        for (int face : patch.faces) boundaries_[face] = condition;
    }
}

void NeutronSolver::buildDiffusionOperatorFromMaterials(const std::vector<NeutronMaterial>& materials,
                                                        double invDtScale,
                                                        SparseMatrix& matrix,
                                                        std::vector<double>& boundaryRhs) const {
    const int size = groupCount_ * cellCount_;
    matrix = SparseMatrix(size);
    boundaryRhs.assign(size, 0.0);

    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials[cellMaterials_[cell]];
        const double volume = mesh_.cells[cell].volume;
        for (int destination = 0; destination < groupCount_; ++destination) {
            const int row = scalarIndex(destination, cell);
            const double selfScattering = material.scattering[destination * groupCount_ + destination];
            const double removal = material.total[destination] - selfScattering +
                                   material.diffusion[destination] * axialBuckling_;
            const double timeDiag = (invDtScale > 0.0 && !neutronVelocity_.empty())
                ? (invDtScale / neutronVelocity_[destination])
                : 0.0;
            matrix.add(row, row, volume * (removal + timeDiag));
            for (int source = 0; source < groupCount_; ++source) {
                if (source == destination) continue;
                const double s = material.scattering[destination * groupCount_ + source];
                if (s != 0.0) matrix.add(row, scalarIndex(source, cell), -volume * s);
            }
        }
    }

    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        const double ownerDistance = norm(face.centre - mesh_.cells[owner].centre);
        const double neighbourDistance = norm(face.centre - mesh_.cells[neighbour].centre);
        for (int group = 0; group < groupCount_; ++group) {
            const double ownerD = materials[cellMaterials_[owner]].diffusion[group];
            const double neighbourD = materials[cellMaterials_[neighbour]].diffusion[group];
            const double effectiveD = (ownerDistance + neighbourDistance) /
                (ownerDistance / ownerD + neighbourDistance / neighbourD);
            const double conductance = effectiveD * face.area / face.delta;
            const int ownerRow = scalarIndex(group, owner);
            const int neighbourRow = scalarIndex(group, neighbour);
            matrix.add(ownerRow, ownerRow, conductance);
            matrix.add(ownerRow, neighbourRow, -conductance);
            matrix.add(neighbourRow, neighbourRow, conductance);
            matrix.add(neighbourRow, ownerRow, -conductance);
        }
    }

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
        for (int group = 0; group < groupCount_; ++group) {
            const int row = scalarIndex(group, owner);
            const double diffusion = materials[cellMaterials_[owner]].diffusion[group];
            if (condition.kind == NeutronBoundaryKind::FixedScalarFlux) {
                const double conductance = diffusion * face.area / face.delta;
                matrix.add(row, row, conductance);
                boundaryRhs[row] += conductance * condition.value(group);
            } else if (condition.kind == NeutronBoundaryKind::IncidentCurrent) {
                // Exact Marshak incident partial current (Eq. 5.297):
                // J^- = phi_b/4 - J_n/2  ==>  J_n = phi_b/2 - 2*J^-_in
                // Eliminating phi_b with J_n = (D/d)*(phi_P - phi_b):
                // J_n = [D/(2D + d)] * phi_P - [4D/(2D + d)] * J^-_in
                const double diagonal = diffusion * face.area / (2.0 * diffusion + face.delta);
                const double source = 4.0 * diffusion * face.area / (2.0 * diffusion + face.delta);
                matrix.add(row, row, diagonal);
                boundaryRhs[row] += source * condition.value(group);
            } else if (condition.kind == NeutronBoundaryKind::Vacuum ||
                       condition.kind == NeutronBoundaryKind::Marshak ||
                       condition.kind == NeutronBoundaryKind::Albedo) {
                const double gamma = robinGamma(condition, group);
                const double conductance = diffusion * gamma * face.area / (diffusion + gamma * face.delta);
                matrix.add(row, row, conductance);
            }
        }
    }
}

void NeutronSolver::buildDiffusionOperator() {
    buildDiffusionOperatorFromMaterials(materials_, 0.0, diffusionMatrix_, diffusionBoundaryRhs_);
    system_->setMatrix(diffusionMatrix_);
}

void NeutronSolver::buildDsaOperators() {
    dsaSystems_.clear();
    dsaSystems_.resize(groupCount_);
    for (int group = 0; group < groupCount_; ++group) {
        bool hasSelfScatter = false;
        for (const NeutronMaterial& mat : materials_) {
            if (mat.scattering[group * groupCount_ + group] > 0.0) hasSelfScatter = true;
        }
        if (!hasSelfScatter) continue;

        SparseMatrix dsaMatrix(cellCount_);
        double totalRemovalAndBoundary = 0.0;
        std::vector<double> dsaD(cellCount_, 0.0);

        for (int cell = 0; cell < cellCount_; ++cell) {
            const NeutronMaterial& mat = materials_[cellMaterials_[cell]];
            const double sigmaTr = std::max(mat.total[group] - mat.scatteringP1[group * groupCount_ + group], 1.0e-3);
            dsaD[cell] = 1.0 / (3.0 * sigmaTr);
            const double selfScat = mat.scattering[group * groupCount_ + group];
            const double rem = std::max(mat.total[group] - selfScat, 0.0) + mat.diffusion[group] * axialBuckling_;
            const double diag = mesh_.cells[cell].volume * rem;
            dsaMatrix.add(cell, cell, diag);
            totalRemovalAndBoundary += diag;
        }

        for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const int neighbour = face.neighbour;
            const double dO = norm(face.centre - mesh_.cells[owner].centre);
            const double dN = norm(face.centre - mesh_.cells[neighbour].centre);
            const double effD = (dO + dN) / (dO / dsaD[owner] + dN / dsaD[neighbour]);
            const double cond = effD * face.area / face.delta;
            dsaMatrix.add(owner, owner, cond);
            dsaMatrix.add(owner, neighbour, -cond);
            dsaMatrix.add(neighbour, neighbour, cond);
            dsaMatrix.add(neighbour, owner, -cond);
        }

        for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            double gamma = 0.0;
            if (condition.kind == NeutronBoundaryKind::Vacuum || condition.kind == NeutronBoundaryKind::IncidentFlux) {
                gamma = 0.5;
            } else if (condition.kind == NeutronBoundaryKind::Albedo) {
                const double beta = condition.value(group);
                gamma = 0.5 * (1.0 - beta) / (1.0 + beta);
            } else if (condition.kind == NeutronBoundaryKind::WhiteReflective) {
                const double beta = condition.value(group);
                gamma = 0.5 * (1.0 - beta) / (1.0 + beta);
            }
            if (gamma > 0.0) {
                const double cond = dsaD[owner] * gamma * face.area / (dsaD[owner] + gamma * face.delta);
                dsaMatrix.add(owner, owner, cond);
                totalRemovalAndBoundary += cond;
            }
        }

        if (totalRemovalAndBoundary <= 1.0e-14) {
            // Purely scattering + all-reflective group: matrix is singular, skip DSA for this group.
            continue;
        }

        dsaSystems_[group] = std::make_unique<PetscLinearSystem>(cellCount_);
        dsaSystems_[group]->setMatrix(dsaMatrix);
    }
}

std::vector<double> NeutronSolver::diffusionRightHandSide(double eigenvalue,
                                                           const std::vector<double>& fissionState,
                                                           bool includeExternal) const {
    std::vector<double> rhs = includeExternal ? diffusionBoundaryRhs_
                                              : std::vector<double>(groupCount_ * cellCount_, 0.0);
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        double fission = 0.0;
        for (int source = 0; source < groupCount_; ++source) {
            fission += material.nuSigmaFission[source] * fissionState[scalarIndex(source, cell)];
        }
        for (int group = 0; group < groupCount_; ++group) {
            const double ext = includeExternal ? material.externalSource[group] : 0.0;
            rhs[scalarIndex(group, cell)] += mesh_.cells[cell].volume *
                (ext + material.chi[group] * fission / eigenvalue);
        }
    }
    return rhs;
}

void NeutronSolver::run() {
    initialize();
    if (!solved_) {
        if (problem_ == NeutronProblem::Eigenvalue) {
            solveEigenvalue();
            if (targetPower_ > 0.0) {
                const double power = balance(eigenvalue_, true).power;
                if (power <= 0.0) throw std::runtime_error("Cannot normalize a zero-power neutron solution");
                normalize(targetPower_ / power);
                normalization_ = "target_power";
            } else {
                double integral = 0.0;
                double volume = 0.0;
                for (int cell = 0; cell < cellCount_; ++cell) {
                    volume += mesh_.cells[cell].volume;
                    for (int group = 0; group < groupCount_; ++group) {
                        integral += mesh_.cells[cell].volume * scalarFlux_[scalarIndex(group, cell)];
                    }
                }
                if (integral <= 0.0) throw std::runtime_error("Cannot normalize a zero-flux neutron solution");
                normalize(volume / integral);
                normalization_ = "mean_total_flux_1";
            }
            if (eigenModesRequested_ > 1) solveHigherModes();
            if (computeAdjoint_) {
                solveAdjoint();
                computeKineticsAndPerturbation();
            }
        } else if (problem_ == NeutronProblem::FixedSource) {
            solveFixedSource();
        } else {
            solveTransient();
        }
        solved_ = true;
    }

    updateOutputFields();
    const std::filesystem::path output = outputPath();
    write(output);
    writeCells(outputPath("_cells.csv"));
    writeSummary(outputPath("_summary.csv"));
    writeHistory(outputPath("_history.csv"));
    if (!regions_.empty() || materials_.size() > 1) writeRegions(outputPath("_regions.csv"));
    if (model_ == NeutronModel::DiscreteOrdinates) writeDirections(outputPath("_directions.csv"));
    if (!modes_.empty()) writeModes(outputPath("_modes.csv"));
    if (!transientHistory_.empty()) writeTransient(outputPath("_transient.csv"));

    if (problem_ != NeutronProblem::Transient) {
        const Balance finalBalance = balance(eigenvalue_, problem_ == NeutronProblem::Eigenvalue);
        std::cout << "Neutron " << name() << " completed: " << history_.size()
                  << " outer iterations, balance " << finalBalance.relativeResidual;
        if (problem_ == NeutronProblem::Eigenvalue) {
            std::cout << ", k_eff " << std::setprecision(12) << eigenvalue_
                      << ", eigen_res " << lastEigenResidual_;
            if (modes_.size() > 1) {
                std::cout << ", DR " << modes_[1].dominanceRatio;
            }
            if (kinetics_.adjointComputed) {
                std::cout << ", k_adj " << kinetics_.adjointEigenvalue;
            }
        }
        std::cout << "\nResults: " << output << '\n';
    } else {
        std::cout << "Neutron " << name() << " transient completed: " << transientHistory_.size() - 1
                  << " steps, final relative power " << std::setprecision(10)
                  << transientHistory_.back().relativePower << "\nResults: " << output << '\n';
    }
}

void NeutronSolver::advance(double dt) {
    (void)dt;
    run();
}

void NeutronSolver::solveFixedSource() {
    history_.clear();
    totalLinearIterations_ = 0;
    lastLinearResidual_ = 0.0;
    eigenvalue_ = 1.0;
    for (int iteration = 1; iteration <= maximumIterations_; ++iteration) {
        const std::vector<double> previous = scalarFlux_;
        const int previousLinearIterations = totalLinearIterations_;
        if (model_ == NeutronModel::DiscreteOrdinates) {
            solveSnInner(1.0, nullptr, true);
        } else {
            lastSourceIterations_ = 1;
            solveDiffusionStep(1.0, previous, true);
        }
        const double change = fluxChange(previous);
        const Balance current = balance(1.0, false);
        history_.push_back({iteration, 1.0, change, 0.0, current.relativeResidual, lastSourceIterations_,
                            totalLinearIterations_ - previousLinearIterations, lastLinearResidual_});
        if (iteration == 1 || iteration % logInterval_ == 0) {
            std::cout << "Neutron fixed-source iteration " << iteration << ": flux change " << change
                      << ", balance " << current.relativeResidual << '\n';
        }
        const bool hasFission = std::any_of(materials_.begin(), materials_.end(), [](const NeutronMaterial& m) {
            return std::any_of(m.nuSigmaFission.begin(), m.nuSigmaFission.end(), [](double v) { return v > 0.0; });
        });
        if ((change < tolerance_ || (!hasFission && iteration == 1 && model_ == NeutronModel::DiscreteOrdinates)) &&
            current.relativeResidual < balanceTolerance_) {
            return;
        }
    }
    throw std::runtime_error("Neutron fixed-source calculation did not converge within maximum_iterations");
}

void NeutronSolver::solveEigenvalue() {
    history_.clear();
    totalLinearIterations_ = 0;
    lastLinearResidual_ = 0.0;
    eigenvalue_ = 1.0;
    lastEigenResidual_ = 0.0;
    for (int iteration = 1; iteration <= maximumIterations_; ++iteration) {
        const double oldEigenvalue = eigenvalue_;
        const std::vector<double> previous = scalarFlux_;
        const double previousProduction = fissionProduction(previous);
        if (!(previousProduction > 0.0)) {
            throw std::runtime_error("Neutron eigenvalue iteration has zero fission production");
        }
        const int previousLinearIterations = totalLinearIterations_;
        if (model_ == NeutronModel::DiscreteOrdinates) {
            solveSnInner(eigenvalue_, &previous, false);
        } else {
            lastSourceIterations_ = 1;
            solveDiffusionStep(eigenvalue_, previous, false);
        }
        const double newProduction = fissionProduction(scalarFlux_);
        if (!(newProduction > 0.0)) throw std::runtime_error("Neutron eigenvalue update has zero fission production");
        eigenvalue_ = oldEigenvalue * newProduction / previousProduction;
        normalize(1.0 / newProduction);
        std::vector<double> normalizedPrevious = previous;
        for (double& value : normalizedPrevious) value /= previousProduction;
        const double change = fluxChange(normalizedPrevious);
        const double eigenChange = std::abs(eigenvalue_ - oldEigenvalue) / std::abs(eigenvalue_);
        const Balance current = balance(eigenvalue_, true);
        lastEigenResidual_ = (model_ == NeutronModel::Diffusion) ? eigenResidual(eigenvalue_) : change;
        history_.push_back({iteration, eigenvalue_, change, lastEigenResidual_, current.relativeResidual,
                            lastSourceIterations_, totalLinearIterations_ - previousLinearIterations,
                            lastLinearResidual_});
        if (iteration == 1 || iteration % logInterval_ == 0) {
            std::cout << "Neutron power iteration " << iteration << ": k_eff " << std::setprecision(12) << eigenvalue_
                      << ", flux change " << change << ", eigen_res " << lastEigenResidual_
                      << ", balance " << current.relativeResidual << '\n';
        }
        if (change < tolerance_ && eigenChange < eigenTolerance_ && current.relativeResidual < balanceTolerance_ &&
            std::abs(current.eigenvalue - eigenvalue_) <= balanceTolerance_ * std::abs(eigenvalue_)) {
            return;
        }
    }
    throw std::runtime_error("Neutron eigenvalue calculation did not converge within maximum_iterations");
}

void NeutronSolver::computeSnHighOrderCorrection(int group, int direction,
                                                 const std::vector<double>& angularState,
                                                 std::vector<double>& boundaryDelta,
                                                 std::vector<double>& cellRhsCorrection) const {
    const std::size_t bCount = mesh_.boundaryFaceCount();
    boundaryDelta.assign(bCount, 0.0);
    cellRhsCorrection.assign(cellCount_, 0.0);

    const NeutronDirection& omega = directions_[direction];
    std::vector<double> psi(cellCount_, 0.0);
    for (int c = 0; c < cellCount_; ++c) psi[c] = angularState[angularIndex(group, direction, c)];

    // Green-Gauss cell gradient of psi_{g,m}
    std::vector<Vector2> grad(cellCount_);
    std::vector<double> minNbr = psi;
    std::vector<double> maxNbr = psi;

    for (std::size_t fIdx = 0; fIdx < mesh_.internalFaceCount(); ++fIdx) {
        const Face& face = mesh_.faces[fIdx];
        const int o = face.owner;
        const int n = face.neighbour;
        const double fVal = face.ownerWeight * psi[o] + (1.0 - face.ownerWeight) * psi[n];
        const Vector2 s = face.normal * face.area;
        grad[o] += s * fVal;
        grad[n] -= s * fVal;
        minNbr[o] = std::min(minNbr[o], psi[n]);
        maxNbr[o] = std::max(maxNbr[o], psi[n]);
        minNbr[n] = std::min(minNbr[n], psi[o]);
        maxNbr[n] = std::max(maxNbr[n], psi[o]);
    }
    for (std::size_t fIdx = mesh_.internalFaceCount(); fIdx < mesh_.faces.size(); ++fIdx) {
        const Face& face = mesh_.faces[fIdx];
        const int o = face.owner;
        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
        double fVal = psi[o];
        if (flux < 0.0) {
            const NeutronBoundaryCondition& cond = boundaries_[fIdx];
            if (cond.kind == NeutronBoundaryKind::Vacuum) {
                fVal = 0.0;
            } else if (cond.kind == NeutronBoundaryKind::IncidentFlux) {
                fVal = cond.value(group);
            } else if (cond.kind == NeutronBoundaryKind::Reflective ||
                       (cond.kind == NeutronBoundaryKind::Albedo && cond.whiteFraction < 1.0e-12)) {
                const int refl = reflectionDirections_[reflectionIndex(static_cast<int>(fIdx), direction)];
                fVal = (refl >= 0) ? cond.value(group) * angularState[angularIndex(group, refl, o)] : psi[o];
            }
            minNbr[o] = std::min(minNbr[o], fVal);
            maxNbr[o] = std::max(maxNbr[o], fVal);
        }
        grad[o] += face.normal * (face.area * fVal);
    }
    for (int c = 0; c < cellCount_; ++c) grad[c] /= mesh_.cells[c].volume;

    // Barth-Jespersen / Venkatakrishnan limiter per cell
    std::vector<double> alpha(cellCount_, 1.0);
    if (snLimiter_ != "none") {
        const bool venkat = (snLimiter_ == "venkatakrishnan");
        for (int c = 0; c < cellCount_; ++c) {
            const double dMax = maxNbr[c] - psi[c];
            const double dMin = minNbr[c] - psi[c];
            double lim = 1.0;
            const double eps2 = venkat ? std::pow(0.3 * std::sqrt(mesh_.cells[c].volume), 3.0) : 0.0;
            for (int fIdx : mesh_.cells[c].faces) {
                const Vector2 dr = mesh_.faces[fIdx].centre - mesh_.cells[c].centre;
                const double delta = dot(grad[c], dr);
                if (delta > 1.0e-15) {
                    if (venkat) {
                        const double num = (dMax * dMax + eps2) * delta + 2.0 * delta * delta * dMax;
                        const double den = delta * (dMax * dMax + 2.0 * delta * delta + dMax * delta + eps2);
                        lim = std::min(lim, den > 0.0 ? num / den : 1.0);
                    } else {
                        lim = std::min(lim, dMax > 0.0 ? std::min(1.0, dMax / delta) : 0.0);
                    }
                } else if (delta < -1.0e-15) {
                    if (venkat) {
                        const double num = (dMin * dMin + eps2) * delta + 2.0 * delta * delta * dMin;
                        const double den = delta * (dMin * dMin + 2.0 * delta * delta + dMin * delta + eps2);
                        lim = std::min(lim, den > 0.0 ? num / den : 1.0);
                    } else {
                        lim = std::min(lim, dMin < 0.0 ? std::min(1.0, dMin / delta) : 0.0);
                    }
                }
            }
            alpha[c] = std::clamp(lim, 0.0, 1.0);
        }
    }

    // Internal face conservative deferred correction: flux_HO = flux_upwind + F * delta_upwind
    for (std::size_t fIdx = 0; fIdx < mesh_.internalFaceCount(); ++fIdx) {
        const Face& face = mesh_.faces[fIdx];
        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
        if (std::abs(flux) < 1.0e-15) continue;
        const int upwind = (flux >= 0.0) ? face.owner : face.neighbour;
        const int downwind = (flux >= 0.0) ? face.neighbour : face.owner;
        double delta = alpha[upwind] * dot(grad[upwind], face.centre - mesh_.cells[upwind].centre);
        if (psi[upwind] + delta < 0.0) delta = -psi[upwind];
        const double corr = std::abs(flux) * delta;
        cellRhsCorrection[upwind] -= corr;
        cellRhsCorrection[downwind] += corr;
    }

    // Outgoing boundary faces: consistent extrapolation for conservation and 2nd-order accuracy
    const std::size_t nInt = mesh_.internalFaceCount();
    for (std::size_t fIdx = nInt; fIdx < mesh_.faces.size(); ++fIdx) {
        const Face& face = mesh_.faces[fIdx];
        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
        if (flux <= 1.0e-15) continue;
        const int owner = face.owner;
        double delta = alpha[owner] * dot(grad[owner], face.centre - mesh_.cells[owner].centre);
        if (psi[owner] + delta < 0.0) delta = -psi[owner];
        boundaryDelta[fIdx - nInt] = delta;
        cellRhsCorrection[owner] -= flux * delta;
    }
}

void NeutronSolver::sweepSnDirection(int group, int direction, double eigenvalue,
                                     const std::vector<double>& sourceScalar,
                                     const std::vector<double>& sourceCurrentX,
                                     const std::vector<double>& sourceCurrentY,
                                     const std::vector<double>* fissionState,
                                     const std::vector<double>& previousAngular,
                                     const std::vector<double>& whiteIncoming,
                                     const std::vector<double>& hoCorrection,
                                     bool includeExternal,
                                     std::vector<double>& angularOut) {
    if (!(eigenvalue > 0.0)) throw std::runtime_error("Neutron eigenvalue must be positive");
    const NeutronDirection& omega = directions_[direction];
    const std::vector<double>& fissionScalar = fissionState ? *fissionState : sourceScalar;
    const std::size_t nInt = mesh_.internalFaceCount();
    const std::size_t bCount = mesh_.boundaryFaceCount();

    if (!directionAcyclic_[direction]) {
        SparseMatrix matrix(cellCount_);
        std::vector<double> rhs(cellCount_, 0.0);
        assembleSnEquation(group, direction, eigenvalue, sourceScalar, sourceCurrentX, sourceCurrentY,
                           fissionState, previousAngular, whiteIncoming, hoCorrection, includeExternal,
                           matrix, rhs);
        std::vector<double> solution(cellCount_, 0.0);
        for (int c = 0; c < cellCount_; ++c) solution[c] = angularOut[angularIndex(group, direction, c)];
        const bool zeroRhs = std::all_of(rhs.begin(), rhs.end(), [](double v) { return v == 0.0; });
        const bool zeroGuess = std::all_of(solution.begin(), solution.end(), [](double v) { return v == 0.0; });
        if (zeroRhs && zeroGuess) return;
        system_->setMatrix(matrix);
        const LinearSolveInfo info = system_->solve(rhs, solution, linearSolver_, preconditioner_,
            linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
        totalLinearIterations_ += info.iterations;
        lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);
        for (int c = 0; c < cellCount_; ++c) angularOut[angularIndex(group, direction, c)] = solution[c];
        return;
    }

    // Fast O(N) direct topological transport sweep
    const int baseOffset = (group * directionCount_ + direction) * cellCount_;
    for (int cell : sweepOrder_[direction]) {
        const Cell& cellGeom = mesh_.cells[cell];
        const NeutronMaterial& mat = materials_[cellMaterials_[cell]];

        double scattering = 0.0;
        double fission = 0.0;
        for (int src = 0; src < groupCount_; ++src) {
            const int sIdx = scalarIndex(src, cell);
            scattering += mat.scattering[group * groupCount_ + src] * sourceScalar[sIdx];
            if (hasP1Scattering_) {
                const double s1 = mat.scatteringP1[group * groupCount_ + src];
                if (s1 != 0.0) {
                    scattering += 3.0 * s1 * (omega.x * sourceCurrentX[sIdx] + omega.y * sourceCurrentY[sIdx]);
                }
            }
            fission += mat.nuSigmaFission[src] * fissionScalar[sIdx];
        }
        const double ext = includeExternal ? mat.externalSource[group] : 0.0;
        double rhs = cellGeom.volume / (4.0 * pi) * (scattering + mat.chi[group] * fission / eigenvalue + ext);
        if (!hoCorrection.empty()) rhs += hoCorrection[cell];

        double diag = cellGeom.volume * (mat.total[group] + mat.diffusion[group] * axialBuckling_);

        for (std::size_t local = 0; local < cellGeom.faces.size(); ++local) {
            const int fIdx = cellGeom.faces[local];
            const int sign = cellGeom.faceSigns[local];
            const Face& face = mesh_.faces[fIdx];
            const double flux = sign * face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
            if (flux >= 0.0) {
                diag += flux;
            } else if (!face.isBoundary()) {
                const int nbr = (face.owner == cell) ? face.neighbour : face.owner;
                rhs -= flux * angularOut[baseOffset + nbr];
            } else {
                const NeutronBoundaryCondition& cond = boundaries_[fIdx];
                if (cond.kind == NeutronBoundaryKind::Vacuum) continue;
                if (cond.kind == NeutronBoundaryKind::IncidentFlux) {
                    if (includeExternal) rhs -= flux * cond.value(group);
                } else if (cond.kind == NeutronBoundaryKind::WhiteReflective) {
                    rhs -= flux * whiteIncoming[fIdx - nInt];
                } else if (cond.kind == NeutronBoundaryKind::Reflective ||
                           cond.kind == NeutronBoundaryKind::Albedo) {
                    const double beta = cond.value(group);
                    const double aw = cond.whiteFraction;
                    double inc = 0.0;
                    if (aw < 1.0 - 1.0e-12) {
                        const int refl = reflectionDirections_[reflectionIndex(fIdx, direction)];
                        if (refl < 0) throw std::runtime_error("Missing reflected ordinate in SN sweep");
                        const std::size_t outIdx = (static_cast<std::size_t>(group) * directionCount_ + refl) * bCount + (fIdx - nInt);
                        const double reflVal = (snScheme_ == NeutronSnScheme::SecondOrder)
                            ? outgoingBoundaryAngular_[outIdx]
                            : previousAngular[angularIndex(group, refl, cell)];
                        inc += beta * (1.0 - aw) * reflVal;
                    }
                    if (aw > 1.0e-12) {
                        inc += aw * whiteIncoming[fIdx - nInt];
                    }
                    rhs -= flux * inc;
                } else {
                    throw std::runtime_error("Invalid boundary kind in SN sweep");
                }
            }
        }
        angularOut[baseOffset + cell] = rhs / diag;
    }
}

void NeutronSolver::applyDsaCorrection(int group, const std::vector<double>& oldScalarGroup) {
    if (group < 0 || group >= static_cast<int>(dsaSystems_.size()) || !dsaSystems_[group]) return;

    std::vector<double> rhs(cellCount_, 0.0);
    double diffNormSq = 0.0;
    for (int c = 0; c < cellCount_; ++c) {
        const double diff = scalarFlux_[scalarIndex(group, c)] - oldScalarGroup[c];
        diffNormSq += mesh_.cells[c].volume * diff * diff;
        const double selfScat = materials_[cellMaterials_[c]].scattering[group * groupCount_ + group];
        rhs[c] = mesh_.cells[c].volume * selfScat * diff;
    }
    if (diffNormSq <= 1.0e-30) return;

    std::vector<double> deltaPhi(cellCount_, 0.0);
    const LinearSolveInfo info = dsaSystems_[group]->solve(rhs, deltaPhi, linearSolver_, preconditioner_,
        linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
    totalLinearIterations_ += info.iterations;
    lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);

    double corrNormSq = 0.0;
    for (int c = 0; c < cellCount_; ++c) {
        if (!std::isfinite(deltaPhi[c])) return;
        corrNormSq += mesh_.cells[c].volume * deltaPhi[c] * deltaPhi[c];
    }
    const double diffNorm = std::sqrt(diffNormSq);
    const double corrNorm = std::sqrt(corrNormSq);
    // Safeguard damping for optically thick coarse cells where continuous DSA can overshoot
    double relax = 1.0;
    if (corrNorm > 20.0 * diffNorm) relax = 20.0 * diffNorm / corrNorm;

    for (int c = 0; c < cellCount_; ++c) {
        const double corr = relax * deltaPhi[c];
        scalarFlux_[scalarIndex(group, c)] += corr;
        for (int m = 0; m < directionCount_; ++m) {
            angularFlux_[angularIndex(group, m, c)] += corr / (4.0 * pi);
        }
    }
}

void NeutronSolver::solveSnInner(double eigenvalue, const std::vector<double>* fissionState, bool includeExternal) {
    const std::size_t nInt = mesh_.internalFaceCount();
    const std::size_t bCount = mesh_.boundaryFaceCount();
    std::vector<double> whiteIncoming(bCount, 0.0);
    std::vector<double> boundaryDelta;
    std::vector<double> hoCorrection;

    for (int iteration = 1; iteration <= maximumSourceIterations_; ++iteration) {
        const std::vector<double> previousScalar = scalarFlux_;
        const std::vector<double> previousAngular = angularFlux_;
        double maxAngularChange = 0.0;
        double maxAngularScale = 0.0;

        for (int group = 0; group < groupCount_; ++group) {
            std::vector<double> oldGroupScalar(cellCount_, 0.0);
            for (int c = 0; c < cellCount_; ++c) oldGroupScalar[c] = scalarFlux_[scalarIndex(group, c)];

            // Compute white-reflective incoming isotropic flux on boundary faces for this group
            std::fill(whiteIncoming.begin(), whiteIncoming.end(), 0.0);
            for (std::size_t fIdx = nInt; fIdx < mesh_.faces.size(); ++fIdx) {
                const NeutronBoundaryCondition& cond = boundaries_[fIdx];
                const bool usesWhite = (cond.kind == NeutronBoundaryKind::WhiteReflective) ||
                                       (cond.kind == NeutronBoundaryKind::Albedo && cond.whiteFraction > 1.0e-12);
                if (!usesWhite) continue;
                const Face& face = mesh_.faces[fIdx];
                const int owner = face.owner;
                double jOut = 0.0;
                for (int m = 0; m < directionCount_; ++m) {
                    const NeutronDirection& omega = directions_[m];
                    const double proj = omega.x * face.normal.x + omega.y * face.normal.y;
                    if (proj > 1.0e-14) {
                        const std::size_t outIdx = (static_cast<std::size_t>(group) * directionCount_ + m) * bCount + (fIdx - nInt);
                        const double psiOut = (snScheme_ == NeutronSnScheme::SecondOrder && iteration > 1)
                            ? outgoingBoundaryAngular_[outIdx]
                            : angularFlux_[angularIndex(group, m, owner)];
                        jOut += omega.weight * proj * psiOut;
                    }
                }
                const double wIn = faceIncomingWeightSum_[fIdx];
                if (wIn > 1.0e-15) {
                    whiteIncoming[fIdx - nInt] = cond.value(group) * jOut / wIn;
                }
            }

            // Sweep all directions in this group (using Gauss-Seidel updated scalar/current moments)
            for (int direction = 0; direction < directionCount_; ++direction) {
                if (snScheme_ == NeutronSnScheme::SecondOrder) {
                    computeSnHighOrderCorrection(group, direction, angularFlux_, boundaryDelta, hoCorrection);
                } else {
                    hoCorrection.clear();
                }

                sweepSnDirection(group, direction, eigenvalue, scalarFlux_, currentX_, currentY_,
                                 fissionState, previousAngular, whiteIncoming, hoCorrection,
                                 includeExternal, angularFlux_);

                if (snScheme_ == NeutronSnScheme::SecondOrder) {
                    const std::size_t outBase = (static_cast<std::size_t>(group) * directionCount_ + direction) * bCount;
                    for (std::size_t b = 0; b < bCount; ++b) {
                        const int owner = mesh_.faces[nInt + b].owner;
                        outgoingBoundaryAngular_[outBase + b] =
                            angularFlux_[angularIndex(group, direction, owner)] + boundaryDelta[b];
                    }
                }
            }

            // Update group scalar flux and current moments immediately (Gauss-Seidel across groups)
            for (int c = 0; c < cellCount_; ++c) {
                double phi = 0.0;
                double jx = 0.0;
                double jy = 0.0;
                for (int m = 0; m < directionCount_; ++m) {
                    const NeutronDirection& omega = directions_[m];
                    const double psi = angularFlux_[angularIndex(group, m, c)];
                    const double wp = omega.weight * psi;
                    phi += wp;
                    jx += wp * omega.x;
                    jy += wp * omega.y;
                }
                const int sIdx = scalarIndex(group, c);
                scalarFlux_[sIdx] = phi;
                currentX_[sIdx] = jx;
                currentY_[sIdx] = jy;
            }

            if (snAcceleration_ == NeutronSnAcceleration::DSA && snScheme_ == NeutronSnScheme::Upwind) {
                applyDsaCorrection(group, oldGroupScalar);
            }
        }

        for (std::size_t i = 0; i < angularFlux_.size(); ++i) {
            maxAngularChange = std::max(maxAngularChange, std::abs(angularFlux_[i] - previousAngular[i]));
            maxAngularScale = std::max(maxAngularScale, std::abs(angularFlux_[i]));
        }

        const double change = fluxChange(previousScalar);
        const double angularRel = maxAngularChange / std::max(maxAngularScale, 1.0e-30);
        lastSourceIterations_ = iteration;
        if (change < sourceTolerance_ && angularRel < std::max(sourceTolerance_ * 10.0, 1.0e-12)) return;
    }
    throw std::runtime_error("Neutron SN source iteration did not converge within maximum_source_iterations");
}

void NeutronSolver::assembleSnEquation(int group, int direction, double eigenvalue,
                                       const std::vector<double>& sourceScalar,
                                       const std::vector<double>& sourceCurrentX,
                                       const std::vector<double>& sourceCurrentY,
                                       const std::vector<double>* fissionState,
                                       const std::vector<double>& previousAngular,
                                       const std::vector<double>& whiteIncoming,
                                       const std::vector<double>& hoCorrection,
                                       bool includeExternal,
                                       SparseMatrix& matrix, std::vector<double>& rhs) const {
    if (!(eigenvalue > 0.0)) throw std::runtime_error("Neutron eigenvalue must be positive");
    const NeutronDirection& omega = directions_[direction];
    const std::vector<double>& fissionScalar = fissionState ? *fissionState : sourceScalar;
    const std::size_t nInt = mesh_.internalFaceCount();
    const std::size_t bCount = mesh_.boundaryFaceCount();

    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        matrix.add(cell, cell, mesh_.cells[cell].volume * (material.total[group] + material.diffusion[group] * axialBuckling_));
        double scattering = 0.0;
        double fission = 0.0;
        for (int source = 0; source < groupCount_; ++source) {
            const int sIdx = scalarIndex(source, cell);
            scattering += material.scattering[group * groupCount_ + source] * sourceScalar[sIdx];
            if (hasP1Scattering_) {
                const double s1 = material.scatteringP1[group * groupCount_ + source];
                if (s1 != 0.0) {
                    scattering += 3.0 * s1 * (omega.x * sourceCurrentX[sIdx] + omega.y * sourceCurrentY[sIdx]);
                }
            }
            fission += material.nuSigmaFission[source] * fissionScalar[sIdx];
        }
        const double ext = includeExternal ? material.externalSource[group] : 0.0;
        rhs[cell] = mesh_.cells[cell].volume / (4.0 * pi) *
            (scattering + material.chi[group] * fission / eigenvalue + ext);
        if (!hoCorrection.empty()) rhs[cell] += hoCorrection[cell];
    }

    for (std::size_t faceIndex = 0; faceIndex < nInt; ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        matrix.add(owner, owner, std::max(flux, 0.0));
        matrix.add(owner, neighbour, std::min(flux, 0.0));
        matrix.add(neighbour, neighbour, std::max(-flux, 0.0));
        matrix.add(neighbour, owner, -std::max(flux, 0.0));
    }

    for (std::size_t faceIndex = nInt; faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
        if (flux >= 0.0) {
            matrix.add(owner, owner, flux);
            continue;
        }
        const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
        if (condition.kind == NeutronBoundaryKind::Vacuum) continue;
        if (condition.kind == NeutronBoundaryKind::IncidentFlux) {
            if (includeExternal) rhs[owner] -= flux * condition.value(group);
        } else if (condition.kind == NeutronBoundaryKind::WhiteReflective) {
            rhs[owner] -= flux * whiteIncoming[faceIndex - nInt];
        } else if (condition.kind == NeutronBoundaryKind::Reflective || condition.kind == NeutronBoundaryKind::Albedo) {
            const double beta = condition.value(group);
            const double aw = condition.whiteFraction;
            double inc = 0.0;
            if (aw < 1.0 - 1.0e-12) {
                const int reflected = reflectionDirections_[reflectionIndex(static_cast<int>(faceIndex), direction)];
                if (reflected < 0) throw std::runtime_error("Missing reflected ordinate in SN boundary assembly");
                const std::size_t outIdx = (static_cast<std::size_t>(group) * directionCount_ + reflected) * bCount + (faceIndex - nInt);
                const double reflVal = (snScheme_ == NeutronSnScheme::SecondOrder)
                    ? outgoingBoundaryAngular_[outIdx]
                    : previousAngular[angularIndex(group, reflected, owner)];
                inc += beta * (1.0 - aw) * reflVal;
            }
            if (aw > 1.0e-12) inc += aw * whiteIncoming[faceIndex - nInt];
            rhs[owner] -= flux * inc;
        } else {
            throw std::runtime_error("Invalid boundary kind in SN equation");
        }
    }
}

void NeutronSolver::solveDiffusionStep(double eigenvalue, const std::vector<double>& fissionState, bool includeExternal) {
    std::vector<double> rhs = diffusionRightHandSide(eigenvalue, fissionState, includeExternal);
    const LinearSolveInfo info = system_->solve(rhs, scalarFlux_, linearSolver_, preconditioner_,
        linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
    totalLinearIterations_ += info.iterations;
    lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);
}

void NeutronSolver::updateScalarAndCurrentMoments() {
    std::fill(scalarFlux_.begin(), scalarFlux_.end(), 0.0);
    std::fill(currentX_.begin(), currentX_.end(), 0.0);
    std::fill(currentY_.begin(), currentY_.end(), 0.0);
    for (int group = 0; group < groupCount_; ++group) {
        for (int direction = 0; direction < directionCount_; ++direction) {
            const NeutronDirection& omega = directions_[direction];
            for (int cell = 0; cell < cellCount_; ++cell) {
                const int sIdx = scalarIndex(group, cell);
                const double wp = omega.weight * angularFlux_[angularIndex(group, direction, cell)];
                scalarFlux_[sIdx] += wp;
                currentX_[sIdx] += wp * omega.x;
                currentY_[sIdx] += wp * omega.y;
            }
        }
    }
}

double NeutronSolver::fluxChange(const std::vector<double>& previous) const {
    if (previous.size() != scalarFlux_.size()) throw std::runtime_error("Neutron flux-history size mismatch");
    // Volume-weighted dimensionless L2 relative norm across ALL groups and cells (no inter-group cancellation)
    double numerator = 0.0;
    double denominator = 0.0;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const double vol = mesh_.cells[cell].volume;
        for (int group = 0; group < groupCount_; ++group) {
            const int index = scalarIndex(group, cell);
            const double diff = scalarFlux_[index] - previous[index];
            numerator += vol * diff * diff;
            denominator += vol * scalarFlux_[index] * scalarFlux_[index];
        }
    }
    return std::sqrt(numerator / std::max(denominator, 1.0e-30));
}

double NeutronSolver::eigenResidual(double eigenvalue) const {
    if (model_ != NeutronModel::Diffusion || !(eigenvalue > 0.0)) return 0.0;
    const int n = groupCount_ * cellCount_;
    const std::vector<double> fissionRhs = diffusionRightHandSide(eigenvalue, scalarFlux_, false);
    double numSq = 0.0;
    double denSq = 0.0;
    for (int row = 0; row < n; ++row) {
        double ax = 0.0;
        for (const auto& [col, val] : diffusionMatrix_.row(row)) {
            ax += val * scalarFlux_[col];
        }
        const double r = ax - fissionRhs[row];
        numSq += r * r;
        denSq += fissionRhs[row] * fissionRhs[row];
    }
    return std::sqrt(numSq) / (std::sqrt(denSq) + 1.0e-30);
}

double NeutronSolver::fissionProduction(const std::vector<double>& scalar) const {
    double result = 0.0;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        for (int group = 0; group < groupCount_; ++group) {
            result += mesh_.cells[cell].volume * material.nuSigmaFission[group] * scalar[scalarIndex(group, cell)];
        }
    }
    return result;
}

void NeutronSolver::normalize(double factor) {
    if (!(factor > 0.0) || !std::isfinite(factor)) throw std::runtime_error("Invalid neutron normalization factor");
    for (double& value : scalarFlux_) value *= factor;
    for (double& value : currentX_) value *= factor;
    for (double& value : currentY_) value *= factor;
    for (double& value : angularFlux_) value *= factor;
    for (double& value : outgoingBoundaryAngular_) value *= factor;
}

NeutronSolver::Balance NeutronSolver::balance(double eigenvalue, bool eigenproblem) const {
    Balance result;
    result.production = fissionProduction(scalarFlux_);
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        for (int group = 0; group < groupCount_; ++group) {
            const double flux = scalarFlux_[scalarIndex(group, cell)];
            result.absorption += mesh_.cells[cell].volume * material.absorption[group] * flux;
            result.axialLeakage += mesh_.cells[cell].volume * material.diffusion[group] * axialBuckling_ * flux;
            result.externalSource += mesh_.cells[cell].volume * material.externalSource[group];
            result.power += mesh_.cells[cell].volume * material.kappaSigmaFission[group] * flux;
        }
    }

    const std::size_t nInt = mesh_.internalFaceCount();
    const std::size_t bCount = mesh_.boundaryFaceCount();

    if (model_ == NeutronModel::DiscreteOrdinates) {
        for (std::size_t faceIndex = nInt; faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            for (int group = 0; group < groupCount_; ++group) {
                // Precompute outgoing partial current if white reflection is active on this face
                double whiteVal = 0.0;
                const bool usesWhite = (condition.kind == NeutronBoundaryKind::WhiteReflective) ||
                                       (condition.kind == NeutronBoundaryKind::Albedo && condition.whiteFraction > 1.0e-12);
                if (usesWhite && faceIncomingWeightSum_[faceIndex] > 1.0e-15) {
                    double jOut = 0.0;
                    for (int m = 0; m < directionCount_; ++m) {
                        const NeutronDirection& omega = directions_[m];
                        const double proj = omega.x * face.normal.x + omega.y * face.normal.y;
                        if (proj > 1.0e-14) {
                            const std::size_t outIdx = (static_cast<std::size_t>(group) * directionCount_ + m) * bCount + (faceIndex - nInt);
                            const double psiOut = (snScheme_ == NeutronSnScheme::SecondOrder)
                                ? outgoingBoundaryAngular_[outIdx]
                                : angularFlux_[angularIndex(group, m, owner)];
                            jOut += omega.weight * proj * psiOut;
                        }
                    }
                    whiteVal = condition.value(group) * jOut / faceIncomingWeightSum_[faceIndex];
                }

                for (int direction = 0; direction < directionCount_; ++direction) {
                    const NeutronDirection& omega = directions_[direction];
                    const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
                    double faceValue = 0.0;
                    if (flux >= 0.0) {
                        const std::size_t outIdx = (static_cast<std::size_t>(group) * directionCount_ + direction) * bCount + (faceIndex - nInt);
                        faceValue = (snScheme_ == NeutronSnScheme::SecondOrder)
                            ? outgoingBoundaryAngular_[outIdx]
                            : angularFlux_[angularIndex(group, direction, owner)];
                    } else if (condition.kind == NeutronBoundaryKind::IncidentFlux) {
                        if (!eigenproblem) faceValue = condition.value(group);
                    } else if (condition.kind == NeutronBoundaryKind::WhiteReflective) {
                        faceValue = whiteVal;
                    } else if (condition.kind == NeutronBoundaryKind::Reflective || condition.kind == NeutronBoundaryKind::Albedo) {
                        const double beta = condition.value(group);
                        const double aw = condition.whiteFraction;
                        if (aw < 1.0 - 1.0e-12) {
                            const int reflected = reflectionDirections_[reflectionIndex(static_cast<int>(faceIndex), direction)];
                            const std::size_t outIdx = (static_cast<std::size_t>(group) * directionCount_ + reflected) * bCount + (faceIndex - nInt);
                            const double reflVal = (snScheme_ == NeutronSnScheme::SecondOrder)
                                ? outgoingBoundaryAngular_[outIdx]
                                : angularFlux_[angularIndex(group, reflected, owner)];
                            faceValue += beta * (1.0 - aw) * reflVal;
                        }
                        if (aw > 1.0e-12) faceValue += aw * whiteVal;
                    }
                    const double partial = omega.weight * flux * faceValue;
                    result.leakage += partial;
                    if (partial < 0.0) result.incoming -= partial;
                }
            }
        }
    } else {
        for (std::size_t faceIndex = nInt; faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            for (int group = 0; group < groupCount_; ++group) {
                const double diffusion = materials_[cellMaterials_[owner]].diffusion[group];
                const double scalar = scalarFlux_[scalarIndex(group, owner)];
                double current = 0.0;
                if (condition.kind == NeutronBoundaryKind::FixedScalarFlux) {
                    const double bVal = eigenproblem ? 0.0 : condition.value(group);
                    current = diffusion * face.area / face.delta * (scalar - bVal);
                } else if (condition.kind == NeutronBoundaryKind::IncidentCurrent) {
                    const double diagonal = diffusion * face.area / (2.0 * diffusion + face.delta);
                    const double source = 4.0 * diffusion * face.area / (2.0 * diffusion + face.delta);
                    const double jIn = eigenproblem ? 0.0 : condition.value(group);
                    current = diagonal * scalar - source * jIn;
                } else if (condition.kind == NeutronBoundaryKind::Vacuum ||
                           condition.kind == NeutronBoundaryKind::Marshak ||
                           condition.kind == NeutronBoundaryKind::Albedo) {
                    const double gamma = robinGamma(condition, group);
                    current = diffusion * gamma * face.area / (diffusion + gamma * face.delta) * scalar;
                }
                result.leakage += current;
                if (current < 0.0) result.incoming -= current;
            }
        }
    }

    const double totalLeakage = result.leakage + result.axialLeakage;
    const double expected = eigenproblem ? result.production / eigenvalue : result.production + result.externalSource;
    result.residual = expected - result.absorption - totalLeakage;
    result.relativeResidual = std::abs(result.residual) /
        std::max({std::abs(expected), std::abs(result.absorption), std::abs(totalLeakage), result.incoming, 1.0e-30});
    const double losses = result.absorption + totalLeakage;
    if (eigenproblem && losses > 0.0) result.eigenvalue = result.production / losses;
    return result;
}

void NeutronSolver::updateOutputFields() const {
    outputScalars_.clear();
    outputCurrent_ = std::make_unique<VectorField>("current", mesh_);
    auto field = [this](const std::string& name) {
        auto result = std::make_unique<ScalarField>(name, mesh_);
        ScalarField* pointer = result.get();
        outputScalars_[name] = std::move(result);
        return pointer;
    };
    std::vector<ScalarField*> groupFields(groupCount_);
    for (int group = 0; group < groupCount_; ++group) {
        groupFields[group] = field("phi_g" + std::to_string(group + 1));
    }
    ScalarField* total = field("phi_total");
    ScalarField* absorption = field("absorption_rate");
    ScalarField* fission = field("fission_rate");
    ScalarField* production = field("neutron_production");
    ScalarField* power = field("power_density");
    ScalarField* materialIndex = field("material_index");

    std::vector<ScalarField*> adjointFields;
    ScalarField* adjointTotal = nullptr;
    if (kinetics_.adjointComputed && !adjointScalarFlux_.empty()) {
        adjointFields.resize(groupCount_);
        for (int g = 0; g < groupCount_; ++g) {
            adjointFields[g] = field("phi_adj_g" + std::to_string(g + 1));
        }
        adjointTotal = field("phi_adj_total");
    }

    std::vector<ScalarField*> modeFields;
    for (std::size_t m = 1; m < modes_.size(); ++m) {
        modeFields.push_back(field("mode_" + std::to_string(m) + "_phi_total"));
    }

    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        (*materialIndex)[cell] = cellMaterials_[cell];
        for (int group = 0; group < groupCount_; ++group) {
            const double flux = scalarFlux_[scalarIndex(group, cell)];
            (*groupFields[group])[cell] = flux;
            (*total)[cell] += flux;
            (*absorption)[cell] += material.absorption[group] * flux;
            (*fission)[cell] += material.fission[group] * flux;
            (*production)[cell] += material.nuSigmaFission[group] * flux;
            (*power)[cell] += material.kappaSigmaFission[group] * flux;
            if (adjointTotal) {
                const double adj = adjointScalarFlux_[scalarIndex(group, cell)];
                (*adjointFields[group])[cell] = adj;
                (*adjointTotal)[cell] += adj;
            }
        }
        for (std::size_t m = 1; m < modes_.size(); ++m) {
            double sum = 0.0;
            for (int g = 0; g < groupCount_; ++g) sum += modes_[m].scalarFlux[scalarIndex(g, cell)];
            (*modeFields[m - 1])[cell] = sum;
        }
    }

    if (model_ == NeutronModel::DiscreteOrdinates) {
        for (int cell = 0; cell < cellCount_; ++cell) {
            Vector2 current;
            for (int group = 0; group < groupCount_; ++group) {
                for (int direction = 0; direction < directionCount_; ++direction) {
                    const NeutronDirection& omega = directions_[direction];
                    const double value = omega.weight * angularFlux_[angularIndex(group, direction, cell)];
                    current.x += value * omega.x;
                    current.y += value * omega.y;
                }
            }
            (*outputCurrent_)[cell] = current;
        }
        return;
    }

    // Diffusion cell-centered current J = -D * grad(phi) using harmonic flux-continuous interface interpolation
    for (int group = 0; group < groupCount_; ++group) {
        std::vector<Vector2> gradient(cellCount_);
        for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int o = face.owner;
            const int n = face.neighbour;
            const double dO = norm(face.centre - mesh_.cells[o].centre);
            const double dN = norm(face.centre - mesh_.cells[n].centre);
            const double gO = materials_[cellMaterials_[o]].diffusion[group] / std::max(dO, 1.0e-30);
            const double gN = materials_[cellMaterials_[n]].diffusion[group] / std::max(dN, 1.0e-30);
            const double value = (gO * scalarFlux_[scalarIndex(group, o)] +
                                  gN * scalarFlux_[scalarIndex(group, n)]) / (gO + gN);
            const Vector2 area = face.normal * face.area;
            gradient[o] += area * value;
            gradient[n] -= area * value;
        }
        for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            const double diffusion = materials_[cellMaterials_[owner]].diffusion[group];
            const double scalar = scalarFlux_[scalarIndex(group, owner)];
            double boundaryValue = scalar;
            if (condition.kind == NeutronBoundaryKind::FixedScalarFlux) {
                boundaryValue = (problem_ == NeutronProblem::Eigenvalue) ? 0.0 : condition.value(group);
            } else if (condition.kind == NeutronBoundaryKind::IncidentCurrent) {
                const double jIn = (problem_ == NeutronProblem::Eigenvalue) ? 0.0 : condition.value(group);
                boundaryValue = (diffusion / face.delta * scalar + 2.0 * jIn) / (diffusion / face.delta + 0.5);
            } else if (condition.kind == NeutronBoundaryKind::Vacuum ||
                       condition.kind == NeutronBoundaryKind::Marshak ||
                       condition.kind == NeutronBoundaryKind::Albedo) {
                const double gamma = robinGamma(condition, group);
                boundaryValue = diffusion * scalar / (diffusion + gamma * face.delta);
            }
            gradient[owner] += face.normal * (face.area * boundaryValue);
        }
        for (int cell = 0; cell < cellCount_; ++cell) {
            const double diffusion = materials_[cellMaterials_[cell]].diffusion[group];
            (*outputCurrent_)[cell] -= gradient[cell] * (diffusion / mesh_.cells[cell].volume);
        }
    }
}

void NeutronSolver::write(const std::filesystem::path& file) const {
    updateOutputFields();
    std::map<std::string, const ScalarField*> scalars;
    for (const auto& [name, values] : outputScalars_) scalars[name] = values.get();
    std::map<std::string, const VectorField*> vectors{{"current", outputCurrent_.get()}};
    writeTecplot(file, mesh_, scalars, vectors, "neutron");
}

void NeutronSolver::writeCells(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron cell output: " + file.string());
    output << "cell,x,y,material,region";
    for (int group = 0; group < groupCount_; ++group) output << ",phi_g" << group + 1;
    output << ",phi_total,current_x,current_y,absorption_rate,fission_rate,neutron_production,power_density";
    if (kinetics_.adjointComputed && !adjointScalarFlux_.empty()) {
        for (int group = 0; group < groupCount_; ++group) output << ",phi_adj_g" << group + 1;
        output << ",phi_adj_total";
    }
    for (std::size_t m = 1; m < modes_.size(); ++m) {
        output << ",mode_" << m << "_phi_total";
    }
    output << '\n';
    output << std::setprecision(16) << std::scientific;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        const std::string regName = (cellRegions_[cell] >= 0) ? regions_[cellRegions_[cell]].name : "default";
        output << cell << ',' << mesh_.cells[cell].centre.x << ',' << mesh_.cells[cell].centre.y << ','
               << material.name << ',' << regName;
        double total = 0.0;
        double absorption = 0.0;
        double fission = 0.0;
        double production = 0.0;
        double power = 0.0;
        for (int group = 0; group < groupCount_; ++group) {
            const double flux = scalarFlux_[scalarIndex(group, cell)];
            output << ',' << flux;
            total += flux;
            absorption += material.absorption[group] * flux;
            fission += material.fission[group] * flux;
            production += material.nuSigmaFission[group] * flux;
            power += material.kappaSigmaFission[group] * flux;
        }
        const Vector2& current = (*outputCurrent_)[cell];
        output << ',' << total << ',' << current.x << ',' << current.y << ','
               << absorption << ',' << fission << ',' << production << ',' << power;
        if (kinetics_.adjointComputed && !adjointScalarFlux_.empty()) {
            double adjTotal = 0.0;
            for (int group = 0; group < groupCount_; ++group) {
                const double adj = adjointScalarFlux_[scalarIndex(group, cell)];
                output << ',' << adj;
                adjTotal += adj;
            }
            output << ',' << adjTotal;
        }
        for (std::size_t m = 1; m < modes_.size(); ++m) {
            double sum = 0.0;
            for (int g = 0; g < groupCount_; ++g) sum += modes_[m].scalarFlux[scalarIndex(g, cell)];
            output << ',' << sum;
        }
        output << '\n';
    }
}

void NeutronSolver::writeSummary(const std::filesystem::path& file) const {
    const bool eigenproblem = (problem_ == NeutronProblem::Eigenvalue);
    const Balance current = balance(eigenvalue_, eigenproblem);
    double minimum = std::numeric_limits<double>::max();
    double maximum = std::numeric_limits<double>::lowest();
    double volume = 0.0;
    double fluxIntegral = 0.0;
    double maxPowerDensity = 0.0;
    double fissileVolume = 0.0;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const double vol = mesh_.cells[cell].volume;
        volume += vol;
        const NeutronMaterial& mat = materials_[cellMaterials_[cell]];
        double cellPower = 0.0;
        bool cellFissile = false;
        for (int group = 0; group < groupCount_; ++group) {
            const double value = scalarFlux_[scalarIndex(group, cell)];
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            fluxIntegral += vol * value;
            cellPower += mat.kappaSigmaFission[group] * value;
            if (mat.kappaSigmaFission[group] > 0.0) cellFissile = true;
        }
        if (cellFissile) {
            fissileVolume += vol;
            maxPowerDensity = std::max(maxPowerDensity, cellPower);
        }
    }
    const double meanFissilePower = (fissileVolume > 0.0) ? (current.power / fissileVolume) : 0.0;
    const double peakingFactor = (meanFissilePower > 0.0) ? (maxPowerDensity / meanFissilePower) : 0.0;

    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron summary: " + file.string());
    output << "problem_type,model,energy_groups,directions,cells,outer_iterations,k_eff,k_balance,reactivity,"
              "eigen_residual,dominance_ratio,k_adjoint,beta_eff,generation_time,prompt_lifetime,"
              "perturbation_dk_over_k,perturbation_drho,power_peaking,"
              "normalization,volume,flux_integral,flux_min,flux_max,fission_production,absorption,"
              "axial_leakage,net_leakage,incoming_current,external_source,fission_power,"
              "balance_absolute,balance_relative,total_linear_iterations,linear_residual\n";

    const std::string pType = (problem_ == NeutronProblem::Eigenvalue) ? "eigenvalue"
                            : (problem_ == NeutronProblem::FixedSource) ? "fixed_source" : "transient";
    output << pType << ',' << name() << ',' << groupCount_ << ','
           << directionCount_ << ',' << cellCount_ << ',' << history_.size() << ',';
    if (eigenproblem || problem_ == NeutronProblem::Transient) {
        output << csvNumber(eigenvalue_) << ',' << csvNumber(current.eigenvalue) << ','
               << csvNumber((eigenvalue_ - 1.0) / eigenvalue_) << ',' << csvNumber(lastEigenResidual_) << ',';
    } else {
        output << ",,,,";
    }
    if (modes_.size() > 1) output << csvNumber(modes_[1].dominanceRatio);
    output << ',';
    if (kinetics_.adjointComputed) output << csvNumber(kinetics_.adjointEigenvalue);
    output << ',';
    if (kinetics_.kineticsComputed) {
        output << csvNumber(kinetics_.betaEffective) << ',' << csvNumber(kinetics_.generationTime) << ','
               << csvNumber(kinetics_.promptLifetime);
    } else {
        output << ",,";
    }
    output << ',';
    if (kinetics_.perturbationComputed) {
        output << csvNumber(kinetics_.perturbationDeltaKOverK) << ',' << csvNumber(kinetics_.perturbationDeltaRho);
    } else {
        output << ",";
    }
    output << ',' << csvNumber(peakingFactor) << ',' << normalization_ << ',' << csvNumber(volume) << ','
           << csvNumber(fluxIntegral) << ',' << csvNumber(minimum) << ',' << csvNumber(maximum) << ','
           << csvNumber(current.production) << ',' << csvNumber(current.absorption) << ','
           << csvNumber(current.axialLeakage) << ',' << csvNumber(current.leakage) << ','
           << csvNumber(current.incoming) << ',' << csvNumber(current.externalSource) << ','
           << csvNumber(current.power) << ',' << csvNumber(current.residual) << ','
           << csvNumber(current.relativeResidual) << ',' << totalLinearIterations_ << ','
           << csvNumber(lastLinearResidual_) << '\n';
}

void NeutronSolver::writeHistory(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron iteration history: " + file.string());
    output << "iteration,k_eff,flux_change,eigen_residual,balance_relative,source_iterations,linear_iterations,linear_residual\n";
    for (const IterationRecord& record : history_) {
        output << record.iteration << ',';
        if (problem_ == NeutronProblem::Eigenvalue || problem_ == NeutronProblem::Transient) {
            output << csvNumber(record.eigenvalue);
        }
        output << ',' << csvNumber(record.fluxChange) << ',' << csvNumber(record.eigenResidual) << ','
               << csvNumber(record.balance) << ',' << record.sourceIterations << ','
               << record.linearIterations << ',' << csvNumber(record.linearResidual) << '\n';
    }
}

void NeutronSolver::writeDirections(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron direction output: " + file.string());
    output << "direction,omega_x,omega_y,omega_z,weight\n" << std::setprecision(16) << std::scientific;
    for (int direction = 0; direction < directionCount_; ++direction) {
        const NeutronDirection& omega = directions_[direction];
        output << direction << ',' << omega.x << ',' << omega.y << ',' << omega.z << ',' << omega.weight << '\n';
    }
}

void NeutronSolver::writeRegions(const std::filesystem::path& file) const {
    struct RegionAccum {
        std::string name;
        std::string material;
        double volume = 0.0;
        std::vector<double> fluxIntegral;
        double absorption = 0.0;
        double fission = 0.0;
        double production = 0.0;
        double power = 0.0;
    };
    std::map<std::string, RegionAccum> stats;
    double totalPower = 0.0;
    double totalFissileVol = 0.0;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& mat = materials_[cellMaterials_[cell]];
        const std::string key = (cellRegions_[cell] >= 0) ? regions_[cellRegions_[cell]].name : ("material:" + mat.name);
        auto& entry = stats[key];
        if (entry.fluxIntegral.empty()) {
            entry.name = key;
            entry.material = mat.name;
            entry.fluxIntegral.assign(groupCount_, 0.0);
        }
        const double vol = mesh_.cells[cell].volume;
        entry.volume += vol;
        double cellPower = 0.0;
        bool fissile = false;
        for (int g = 0; g < groupCount_; ++g) {
            const double phi = scalarFlux_[scalarIndex(g, cell)];
            entry.fluxIntegral[g] += vol * phi;
            entry.absorption += vol * mat.absorption[g] * phi;
            entry.fission += vol * mat.fission[g] * phi;
            entry.production += vol * mat.nuSigmaFission[g] * phi;
            cellPower += vol * mat.kappaSigmaFission[g] * phi;
            if (mat.kappaSigmaFission[g] > 0.0) fissile = true;
        }
        entry.power += cellPower;
        totalPower += cellPower;
        if (fissile) totalFissileVol += vol;
    }
    const double meanPowerDensity = (totalFissileVol > 0.0) ? (totalPower / totalFissileVol) : 0.0;

    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron region output: " + file.string());
    output << "region,material,volume";
    for (int g = 0; g < groupCount_; ++g) output << ",mean_phi_g" << g + 1;
    output << ",absorption,fission,production,power,normalized_power_density\n";
    output << std::setprecision(16) << std::scientific;
    for (const auto& [key, entry] : stats) {
        output << entry.name << ',' << entry.material << ',' << entry.volume;
        for (int g = 0; g < groupCount_; ++g) output << ',' << entry.fluxIntegral[g] / std::max(entry.volume, 1.0e-30);
        const double normP = (meanPowerDensity > 0.0 && entry.volume > 0.0)
            ? ((entry.power / entry.volume) / meanPowerDensity) : 0.0;
        output << ',' << entry.absorption << ',' << entry.fission << ',' << entry.production << ','
               << entry.power << ',' << normP << '\n';
    }
}

void NeutronSolver::writeModes(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron modes output: " + file.string());
    output << "mode,k_eff,k_imag,reactivity,dominance_ratio,residual\n";
    for (const NeutronEigenmode& mode : modes_) {
        const double rho = (std::abs(mode.eigenvalue) > 1.0e-30) ? (mode.eigenvalue - 1.0) / mode.eigenvalue : 0.0;
        output << mode.index << ',' << csvNumber(mode.eigenvalue) << ',' << csvNumber(mode.imaginary) << ','
               << csvNumber(rho) << ',' << csvNumber(mode.dominanceRatio) << ',' << csvNumber(mode.residual) << '\n';
    }
}

void NeutronSolver::writeTransient(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron transient output: " + file.string());
    output << "step,time,total_power,relative_power,peak_power_density,flux_integral,precursor_integral,linear_iterations\n";
    for (const TransientRecord& rec : transientHistory_) {
        output << rec.step << ',' << csvNumber(rec.time) << ',' << csvNumber(rec.totalPower) << ','
               << csvNumber(rec.relativePower) << ',' << csvNumber(rec.peakPowerDensity) << ','
               << csvNumber(rec.fluxIntegral) << ',' << csvNumber(rec.precursorIntegral) << ','
               << rec.linearIterations << '\n';
    }
}

std::filesystem::path NeutronSolver::outputPath(const std::string& suffix) const {
    std::filesystem::path output = config_.require("output");
    if (output.is_relative()) output = caseDirectory_ / output;
    if (suffix.empty()) return output;
    return output.parent_path() / (output.stem().string() + suffix);
}

} // namespace mp
