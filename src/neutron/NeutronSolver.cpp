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

std::vector<double> optionalNumbers(const Config& config, const std::string& key, int count, double fallback) {
    return config.has(key) ? numbers(config.require(key), count, key) : std::vector<double>(count, fallback);
}

std::vector<double> optionalScattering(const Config& config, const std::string& key, int groups) {
    if (!config.has(key)) return std::vector<double>(groups * groups, 0.0);
    return numbers(config.require(key), groups * groups, key);
}

void requireNonnegative(const std::vector<double>& values, const std::string& key) {
    for (double value : values) {
        if (value < 0.0) throw std::runtime_error("Setting " + key + " must be non-negative");
    }
}

double robinGamma(const NeutronBoundaryCondition& condition) {
    if (condition.kind == NeutronBoundaryKind::Vacuum) return 0.5;
    if (condition.kind == NeutronBoundaryKind::Marshak) return condition.value;
    if (condition.kind == NeutronBoundaryKind::Albedo) {
        return 0.5 * (1.0 - condition.value) / (1.0 + condition.value);
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
        system_ = std::make_unique<PetscLinearSystem>(cellCount_);
    } else {
        system_ = std::make_unique<PetscLinearSystem>(groupCount_ * cellCount_);
    }
    buildBoundaries();
    if (model_ == NeutronModel::Diffusion) buildDiffusionOperator();
    scalarFlux_.assign(groupCount_ * cellCount_, initialFlux_);
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
    } else {
        throw std::runtime_error("Unknown neutron problem_type: " + problem);
    }

    groupCount_ = config_.integer("energy_groups", 1);
    polarCount_ = config_.integer("sn_polar", 2);
    azimuthCount_ = config_.integer("sn_azimuth", 4);
    if (groupCount_ <= 0) throw std::runtime_error("energy_groups must be positive");
    if (model_ == NeutronModel::DiscreteOrdinates &&
        (polarCount_ < 2 || polarCount_ % 2 != 0 || azimuthCount_ < 4 || azimuthCount_ % 4 != 0)) {
        throw std::runtime_error("sn_polar must be even and at least 2; sn_azimuth must be a multiple of 4");
    }

    maximumIterations_ = config_.integer("maximum_iterations", 100);
    maximumSourceIterations_ = config_.integer("maximum_source_iterations", 100);
    maximumLinearIterations_ = config_.integer("maximum_linear_iterations", 1000);
    logInterval_ = config_.integer("log_interval", 1);
    tolerance_ = config_.number("tolerance", 1.0e-8);
    eigenTolerance_ = config_.number("eigen_tolerance", 1.0e-10);
    sourceTolerance_ = config_.number("source_tolerance", 1.0e-9);
    balanceTolerance_ = config_.number("balance_tolerance", 1.0e-8);
    linearTolerance_ = config_.number("linear_tolerance", 1.0e-10);
    linearAbsoluteTolerance_ = config_.number("linear_absolute_tolerance", 1.0e-12);
    initialFlux_ = config_.number("initial_flux", 1.0);
    targetPower_ = config_.number("target_power", 0.0);
    linearSolver_ = config_.get("linear_solver", "gmres");
    preconditioner_ = config_.get("preconditioner", "ilu");
    if (maximumIterations_ <= 0 || maximumSourceIterations_ <= 0 || maximumLinearIterations_ <= 0 || logInterval_ <= 0 ||
        tolerance_ <= 0.0 || eigenTolerance_ <= 0.0 || sourceTolerance_ <= 0.0 || balanceTolerance_ <= 0.0 ||
        linearTolerance_ <= 0.0 || linearAbsoluteTolerance_ <= 0.0 || initialFlux_ < 0.0 || targetPower_ < 0.0) {
        throw std::runtime_error("Invalid neutron iteration, tolerance, or normalization setting");
    }
    if (problem_ == NeutronProblem::Eigenvalue && initialFlux_ <= 0.0) {
        throw std::runtime_error("A positive initial_flux is required for the neutron eigenvalue problem");
    }
    if (problem_ == NeutronProblem::FixedSource && targetPower_ > 0.0) {
        throw std::runtime_error("target_power normalization is only available for neutron eigenvalue calculations");
    }

    const std::vector<std::string> materialNames = strings(config_.get("materials", "default"));
    if (materialNames.empty()) throw std::runtime_error("Neutron materials cannot be empty");
    std::set<std::string> uniqueNames;
    for (const std::string& name : materialNames) {
        if (!uniqueNames.insert(name).second) throw std::runtime_error("Duplicate neutron material: " + name);
    }

    for (const std::string& name : materialNames) {
        const std::string prefix = "material." + name + ".";
        NeutronMaterial material;
        material.name = name;
        material.total = numbers(config_.require(prefix + "sigma_total"), groupCount_, prefix + "sigma_total");
        material.scattering = optionalScattering(config_, prefix + "scattering", groupCount_);
        requireNonnegative(material.total, prefix + "sigma_total");
        requireNonnegative(material.scattering, prefix + "scattering");

        std::vector<double> outgoing(groupCount_, 0.0);
        for (int source = 0; source < groupCount_; ++source) {
            for (int destination = 0; destination < groupCount_; ++destination) {
                outgoing[source] += material.scattering[destination * groupCount_ + source];
            }
        }
        if (config_.has(prefix + "sigma_absorption")) {
            material.absorption = numbers(config_.require(prefix + "sigma_absorption"), groupCount_, prefix + "sigma_absorption");
        } else {
            material.absorption.resize(groupCount_);
            for (int group = 0; group < groupCount_; ++group) {
                material.absorption[group] = material.total[group] - outgoing[group];
            }
        }
        requireNonnegative(material.absorption, prefix + "sigma_absorption");
        for (int group = 0; group < groupCount_; ++group) {
            if (material.total[group] <= 0.0) throw std::runtime_error("Neutron sigma_total must be positive");
            const double expected = material.absorption[group] + outgoing[group];
            if (std::abs(expected - material.total[group]) > 1.0e-8 * std::max(1.0, material.total[group])) {
                throw std::runtime_error("Neutron total cross section must equal absorption plus outgoing scattering in material " + name);
            }
        }

        material.fission = optionalNumbers(config_, prefix + "sigma_fission", groupCount_, 0.0);
        requireNonnegative(material.fission, prefix + "sigma_fission");
        for (int group = 0; group < groupCount_; ++group) {
            if (material.fission[group] > material.absorption[group] + 1.0e-12) {
                throw std::runtime_error("Neutron sigma_fission cannot exceed sigma_absorption");
            }
        }
        const bool fissile = std::any_of(material.fission.begin(), material.fission.end(), [](double value) { return value > 0.0; });
        std::vector<double> nu = fissile
            ? numbers(config_.require(prefix + "nu"), groupCount_, prefix + "nu")
            : optionalNumbers(config_, prefix + "nu", groupCount_, 0.0);
        requireNonnegative(nu, prefix + "nu");
        const std::vector<double> energyRelease = optionalNumbers(config_, prefix + "energy_release", groupCount_, 1.0);
        requireNonnegative(energyRelease, prefix + "energy_release");
        material.nuSigmaFission.resize(groupCount_);
        material.kappaSigmaFission.resize(groupCount_);
        for (int group = 0; group < groupCount_; ++group) {
            material.nuSigmaFission[group] = nu[group] * material.fission[group];
            material.kappaSigmaFission[group] = energyRelease[group] * material.fission[group];
        }
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
        for (int group = 0; group < groupCount_; ++group) {
            if (!config_.has(prefix + "diffusion_coefficient")) {
                material.diffusion[group] = 1.0 / (3.0 * material.total[group]);
            }
            if (material.diffusion[group] <= 0.0) throw std::runtime_error("Neutron diffusion coefficients must be positive");
        }
        material.externalSource = optionalNumbers(config_, prefix + "external_source", groupCount_, 0.0);
        requireNonnegative(material.externalSource, prefix + "external_source");
        materials_.push_back(std::move(material));
    }

    const std::string defaultName = config_.get("default_material", materialNames.front());
    const auto defaultIterator = std::find(materialNames.begin(), materialNames.end(), defaultName);
    if (defaultIterator == materialNames.end()) throw std::runtime_error("Unknown default_material: " + defaultName);
    cellMaterials_.assign(cellCount_, static_cast<int>(defaultIterator - materialNames.begin()));
    struct Region { Vector2 lower; Vector2 upper; int material = 0; };
    std::vector<Region> regions;
    for (const std::string& name : strings(config_.get("regions", ""))) {
        const std::string materialName = config_.require("region." + name + ".material");
        const auto iterator = std::find(materialNames.begin(), materialNames.end(), materialName);
        if (iterator == materialNames.end()) throw std::runtime_error("Unknown material in region " + name + ": " + materialName);
        const std::vector<double> bounds = numbers(config_.require("region." + name + ".bounds"), 4, "region." + name + ".bounds");
        if (bounds[1] <= bounds[0] || bounds[3] <= bounds[2]) throw std::runtime_error("Invalid bounds in region " + name);
        regions.push_back({{bounds[0], bounds[2]}, {bounds[1], bounds[3]}, static_cast<int>(iterator - materialNames.begin())});
    }
    for (int cell = 0; cell < cellCount_; ++cell) {
        const Vector2& centre = mesh_.cells[cell].centre;
        for (const Region& region : regions) {
            if (centre.x >= region.lower.x && centre.x <= region.upper.x &&
                centre.y >= region.lower.y && centre.y <= region.upper.y) {
                cellMaterials_[cell] = region.material;
            }
        }
    }
}

void NeutronSolver::buildQuadrature() {
    std::vector<double> nodes(polarCount_);
    std::vector<double> weights(polarCount_);
    for (int i = 0; i < (polarCount_ + 1) / 2; ++i) {
        double z = std::cos(pi * (i + 0.75) / (polarCount_ + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double p0 = 1.0;
            double p1 = z;
            for (int degree = 2; degree <= polarCount_; ++degree) {
                const double p2 = ((2.0 * degree - 1.0) * z * p1 - (degree - 1.0) * p0) / degree;
                p0 = p1;
                p1 = p2;
            }
            derivative = polarCount_ * (z * p1 - p0) / (z * z - 1.0);
            const double step = p1 / derivative;
            z -= step;
            if (std::abs(step) < 1.0e-15) break;
        }
        const double weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
        nodes[i] = -z;
        nodes[polarCount_ - 1 - i] = z;
        weights[i] = weight;
        weights[polarCount_ - 1 - i] = weight;
    }
    directions_.clear();
    for (int polar = 0; polar < polarCount_; ++polar) {
        const double transverse = std::sqrt(std::max(0.0, 1.0 - nodes[polar] * nodes[polar]));
        for (int azimuth = 0; azimuth < azimuthCount_; ++azimuth) {
            const double angle = 2.0 * pi * azimuth / azimuthCount_;
            directions_.push_back({transverse * std::cos(angle), transverse * std::sin(angle), nodes[polar],
                                   weights[polar] * 2.0 * pi / azimuthCount_});
        }
    }
    directionCount_ = static_cast<int>(directions_.size());
    const double weightSum = std::accumulate(directions_.begin(), directions_.end(), 0.0,
        [](double sum, const NeutronDirection& direction) { return sum + direction.weight; });
    for (NeutronDirection& direction : directions_) direction.weight *= 4.0 * pi / weightSum;
}

void NeutronSolver::buildBoundaries() {
    for (const BoundaryPatch& patch : mesh_.patches) {
        const std::vector<std::string> words = strings(config_.require("boundary." + patch.name + ".neutron"));
        if (words.empty()) throw std::runtime_error("Empty neutron boundary on patch " + patch.name);
        const std::string kind = lowercase(words.front());
        NeutronBoundaryCondition condition;
        auto value = [&](double fallback, bool required) {
            if (words.size() > 2) throw std::runtime_error("Unexpected neutron boundary data on patch " + patch.name);
            if (words.size() == 2) return std::stod(words[1]);
            if (required) throw std::runtime_error("Neutron boundary requires a value on patch " + patch.name);
            return fallback;
        };
        if (kind == "vacuum") {
            condition.kind = NeutronBoundaryKind::Vacuum;
            value(0.0, false);
        } else if (kind == "reflective" || kind == "reflection") {
            condition.kind = NeutronBoundaryKind::Reflective;
            condition.value = 1.0;
            value(1.0, false);
        } else if (kind == "albedo") {
            condition.kind = NeutronBoundaryKind::Albedo;
            condition.value = value(1.0, true);
            if (condition.value < 0.0 || condition.value > 1.0) throw std::runtime_error("Neutron albedo must be in [0,1]");
        } else if (kind == "incident" || kind == "incident_flux") {
            condition.kind = NeutronBoundaryKind::IncidentFlux;
            condition.value = value(0.0, true);
            if (condition.value < 0.0 || model_ != NeutronModel::DiscreteOrdinates) {
                throw std::runtime_error("incident_flux requires a non-negative value and the SN model");
            }
        } else if (kind == "incident_current") {
            condition.kind = NeutronBoundaryKind::IncidentCurrent;
            condition.value = value(0.0, true);
            if (condition.value < 0.0 || model_ != NeutronModel::Diffusion) {
                throw std::runtime_error("incident_current requires a non-negative value and the diffusion model");
            }
        } else if (kind == "fixed_flux" || kind == "fixedvalue") {
            condition.kind = NeutronBoundaryKind::FixedScalarFlux;
            condition.value = value(0.0, true);
            if (model_ != NeutronModel::Diffusion) throw std::runtime_error("fixed_flux is a diffusion-only boundary");
        } else if (kind == "marshak") {
            condition.kind = NeutronBoundaryKind::Marshak;
            condition.value = value(0.5, false);
            if (condition.value <= 0.0 || model_ != NeutronModel::Diffusion) {
                throw std::runtime_error("marshak requires a positive Robin coefficient and the diffusion model");
            }
        } else if (kind == "zero_current" || kind == "zerogradient") {
            condition.kind = NeutronBoundaryKind::ZeroCurrent;
            value(0.0, false);
            if (model_ != NeutronModel::Diffusion) throw std::runtime_error("Use reflective, not zero_current, for SN boundaries");
        } else {
            throw std::runtime_error("Unknown neutron boundary '" + kind + "' on patch " + patch.name);
        }
        for (int face : patch.faces) boundaries_[face] = condition;
    }

    if (model_ != NeutronModel::DiscreteOrdinates) return;
    reflectionDirections_.assign(mesh_.faces.size() * directionCount_, -1);
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
        if (condition.kind != NeutronBoundaryKind::Reflective && condition.kind != NeutronBoundaryKind::Albedo) continue;
        const bool xAligned = std::abs(std::abs(face.normal.x) - 1.0) < 1.0e-12 && std::abs(face.normal.y) < 1.0e-12;
        const bool yAligned = std::abs(std::abs(face.normal.y) - 1.0) < 1.0e-12 && std::abs(face.normal.x) < 1.0e-12;
        if (!xAligned && !yAligned) {
            throw std::runtime_error("The product SN quadrature currently supports reflective boundaries aligned with x or y");
        }
        for (int direction = 0; direction < directionCount_; ++direction) {
            const NeutronDirection& omega = directions_[direction];
            const double normalProjection = omega.x * face.normal.x + omega.y * face.normal.y;
            if (normalProjection >= 0.0) continue;
            const double reflectedX = omega.x - 2.0 * normalProjection * face.normal.x;
            const double reflectedY = omega.y - 2.0 * normalProjection * face.normal.y;
            int best = -1;
            double bestDistance = std::numeric_limits<double>::max();
            for (int candidate = 0; candidate < directionCount_; ++candidate) {
                const NeutronDirection& other = directions_[candidate];
                const double distance = (other.x - reflectedX) * (other.x - reflectedX) +
                    (other.y - reflectedY) * (other.y - reflectedY) + (other.z - omega.z) * (other.z - omega.z);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    best = candidate;
                }
            }
            if (best < 0 || bestDistance > 1.0e-20 ||
                std::abs(directions_[best].weight - omega.weight) > 1.0e-12 * omega.weight) {
                throw std::runtime_error("No exact reflected ordinate exists for this boundary and quadrature");
            }
            reflectionDirections_[reflectionIndex(static_cast<int>(faceIndex), direction)] = best;
        }
    }
}

void NeutronSolver::buildDiffusionOperator() {
    const int size = groupCount_ * cellCount_;
    diffusionMatrix_ = SparseMatrix(size);
    diffusionBoundaryRhs_.assign(size, 0.0);
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        const double volume = mesh_.cells[cell].volume;
        for (int destination = 0; destination < groupCount_; ++destination) {
            const int row = scalarIndex(destination, cell);
            const double selfScattering = material.scattering[destination * groupCount_ + destination];
            diffusionMatrix_.add(row, row, volume * (material.total[destination] - selfScattering));
            for (int source = 0; source < groupCount_; ++source) {
                if (source == destination) continue;
                diffusionMatrix_.add(row, scalarIndex(source, cell),
                    -volume * material.scattering[destination * groupCount_ + source]);
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
            const double ownerD = materials_[cellMaterials_[owner]].diffusion[group];
            const double neighbourD = materials_[cellMaterials_[neighbour]].diffusion[group];
            const double effectiveD = (ownerDistance + neighbourDistance) /
                (ownerDistance / ownerD + neighbourDistance / neighbourD);
            const double conductance = effectiveD * face.area / face.delta;
            const int ownerRow = scalarIndex(group, owner);
            const int neighbourRow = scalarIndex(group, neighbour);
            diffusionMatrix_.add(ownerRow, ownerRow, conductance);
            diffusionMatrix_.add(ownerRow, neighbourRow, -conductance);
            diffusionMatrix_.add(neighbourRow, neighbourRow, conductance);
            diffusionMatrix_.add(neighbourRow, ownerRow, -conductance);
        }
    }

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
        for (int group = 0; group < groupCount_; ++group) {
            const int row = scalarIndex(group, owner);
            const double diffusion = materials_[cellMaterials_[owner]].diffusion[group];
            if (condition.kind == NeutronBoundaryKind::FixedScalarFlux) {
                const double conductance = diffusion * face.area / face.delta;
                diffusionMatrix_.add(row, row, conductance);
                diffusionBoundaryRhs_[row] += conductance * condition.value;
            } else if (condition.kind == NeutronBoundaryKind::IncidentCurrent) {
                const double diagonal = diffusion * face.area / (2.0 * diffusion + face.delta);
                const double source = 2.0 * diffusion * face.area / (2.0 * diffusion + face.delta);
                diffusionMatrix_.add(row, row, diagonal);
                diffusionBoundaryRhs_[row] += source * condition.value;
            } else if (condition.kind == NeutronBoundaryKind::Vacuum ||
                       condition.kind == NeutronBoundaryKind::Marshak ||
                       condition.kind == NeutronBoundaryKind::Albedo) {
                const double gamma = robinGamma(condition);
                const double conductance = diffusion * gamma * face.area / (diffusion + gamma * face.delta);
                diffusionMatrix_.add(row, row, conductance);
            }
        }
    }
    system_->setMatrix(diffusionMatrix_);
}

std::vector<double> NeutronSolver::diffusionRightHandSide(double eigenvalue,
                                                           const std::vector<double>& fissionState) const {
    std::vector<double> rhs = diffusionBoundaryRhs_;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        double fission = 0.0;
        for (int source = 0; source < groupCount_; ++source) {
            fission += material.nuSigmaFission[source] * fissionState[scalarIndex(source, cell)];
        }
        for (int group = 0; group < groupCount_; ++group) {
            rhs[scalarIndex(group, cell)] += mesh_.cells[cell].volume *
                (material.externalSource[group] + material.chi[group] * fission / eigenvalue);
        }
    }
    return rhs;
}

void NeutronSolver::run() {
    initialize();
    if (!solved_) {
        if (problem_ == NeutronProblem::Eigenvalue) solveEigenvalue();
        else solveFixedSource();
        if (problem_ == NeutronProblem::Eigenvalue) {
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
        }
        solved_ = true;
    }
    updateOutputFields();
    const std::filesystem::path output = outputPath();
    write(output);
    writeCells(outputPath("_cells.csv"));
    writeSummary(outputPath("_summary.csv"));
    writeHistory(outputPath("_history.csv"));
    if (model_ == NeutronModel::DiscreteOrdinates) writeDirections(outputPath("_directions.csv"));
    const Balance finalBalance = balance(eigenvalue_, problem_ == NeutronProblem::Eigenvalue);
    std::cout << "Neutron " << name() << " completed: " << history_.size()
              << " outer iterations, balance " << finalBalance.relativeResidual;
    if (problem_ == NeutronProblem::Eigenvalue) std::cout << ", k_eff " << std::setprecision(12) << eigenvalue_;
    std::cout << "\nResults: " << output << '\n';
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
            solveSnInner(1.0, nullptr);
        } else {
            lastSourceIterations_ = 1;
            solveDiffusionStep(1.0, previous);
        }
        const double change = fluxChange(previous);
        const Balance current = balance(1.0, false);
        history_.push_back({iteration, 1.0, change, current.relativeResidual, lastSourceIterations_,
                            totalLinearIterations_ - previousLinearIterations, lastLinearResidual_});
        if (iteration == 1 || iteration % logInterval_ == 0) {
            std::cout << "Neutron fixed-source iteration " << iteration << ": flux change " << change
                      << ", balance " << current.relativeResidual << '\n';
        }
        if (change < tolerance_ && current.relativeResidual < balanceTolerance_) return;
    }
    throw std::runtime_error("Neutron fixed-source calculation did not converge within maximum_iterations");
}

void NeutronSolver::solveEigenvalue() {
    history_.clear();
    totalLinearIterations_ = 0;
    lastLinearResidual_ = 0.0;
    eigenvalue_ = 1.0;
    for (int iteration = 1; iteration <= maximumIterations_; ++iteration) {
        const double oldEigenvalue = eigenvalue_;
        const std::vector<double> previous = scalarFlux_;
        const double previousProduction = fissionProduction(previous);
        if (!(previousProduction > 0.0)) {
            throw std::runtime_error("Neutron eigenvalue iteration has zero fission production");
        }
        const int previousLinearIterations = totalLinearIterations_;
        if (model_ == NeutronModel::DiscreteOrdinates) {
            solveSnInner(eigenvalue_, &previous);
        } else {
            lastSourceIterations_ = 1;
            solveDiffusionStep(eigenvalue_, previous);
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
        history_.push_back({iteration, eigenvalue_, change, current.relativeResidual, lastSourceIterations_,
                            totalLinearIterations_ - previousLinearIterations, lastLinearResidual_});
        if (iteration == 1 || iteration % logInterval_ == 0) {
            std::cout << "Neutron power iteration " << iteration << ": k_eff " << std::setprecision(12) << eigenvalue_
                      << ", flux change " << change << ", balance " << current.relativeResidual << '\n';
        }
        if (change < tolerance_ && eigenChange < eigenTolerance_ && current.relativeResidual < balanceTolerance_ &&
            std::abs(current.eigenvalue - eigenvalue_) <= balanceTolerance_ * std::abs(eigenvalue_)) return;
    }
    throw std::runtime_error("Neutron eigenvalue calculation did not converge within maximum_iterations");
}

void NeutronSolver::solveSnInner(double eigenvalue, const std::vector<double>* fissionState) {
    for (int iteration = 1; iteration <= maximumSourceIterations_; ++iteration) {
        const std::vector<double> previousScalar = scalarFlux_;
        const std::vector<double> previousAngular = angularFlux_;
        for (int group = 0; group < groupCount_; ++group) {
            for (int direction = 0; direction < directionCount_; ++direction) {
                SparseMatrix matrix(cellCount_);
                std::vector<double> rhs(cellCount_, 0.0);
                assembleSnEquation(group, direction, eigenvalue, previousScalar, fissionState,
                                   previousAngular, matrix, rhs);
                std::vector<double> solution(cellCount_, 0.0);
                for (int cell = 0; cell < cellCount_; ++cell) {
                    solution[cell] = angularFlux_[angularIndex(group, direction, cell)];
                }
                const bool zeroRhs = std::all_of(rhs.begin(), rhs.end(), [](double value) { return value == 0.0; });
                const bool zeroGuess = std::all_of(solution.begin(), solution.end(), [](double value) { return value == 0.0; });
                if (zeroRhs && zeroGuess) {
                    for (int cell = 0; cell < cellCount_; ++cell) angularFlux_[angularIndex(group, direction, cell)] = 0.0;
                    continue;
                }
                system_->setMatrix(matrix);
                const LinearSolveInfo info = system_->solve(rhs, solution, linearSolver_, preconditioner_,
                    linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
                totalLinearIterations_ += info.iterations;
                lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);
                for (int cell = 0; cell < cellCount_; ++cell) {
                    angularFlux_[angularIndex(group, direction, cell)] = solution[cell];
                }
            }
        }
        updateScalarFlux();
        const double change = fluxChange(previousScalar);
        lastSourceIterations_ = iteration;
        if (change < sourceTolerance_) return;
    }
    throw std::runtime_error("Neutron SN source iteration did not converge within maximum_source_iterations");
}

void NeutronSolver::assembleSnEquation(int group, int direction, double eigenvalue,
                                       const std::vector<double>& sourceScalar,
                                       const std::vector<double>* fissionState,
                                       const std::vector<double>& previousAngular,
                                       SparseMatrix& matrix, std::vector<double>& rhs) const {
    if (!(eigenvalue > 0.0)) throw std::runtime_error("Neutron eigenvalue must be positive");
    const NeutronDirection& omega = directions_[direction];
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        matrix.add(cell, cell, mesh_.cells[cell].volume * material.total[group]);
        double scattering = 0.0;
        double fission = 0.0;
        const std::vector<double>& fissionScalar = fissionState ? *fissionState : sourceScalar;
        for (int source = 0; source < groupCount_; ++source) {
            scattering += material.scattering[group * groupCount_ + source] * sourceScalar[scalarIndex(source, cell)];
            fission += material.nuSigmaFission[source] * fissionScalar[scalarIndex(source, cell)];
        }
        rhs[cell] = mesh_.cells[cell].volume / (4.0 * pi) *
            (scattering + material.chi[group] * fission / eigenvalue + material.externalSource[group]);
    }

    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        matrix.add(owner, owner, std::max(flux, 0.0));
        matrix.add(owner, neighbour, std::min(flux, 0.0));
        matrix.add(neighbour, neighbour, std::max(-flux, 0.0));
        matrix.add(neighbour, owner, -std::max(flux, 0.0));
    }

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
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
            rhs[owner] -= flux * condition.value;
        } else if (condition.kind == NeutronBoundaryKind::Reflective || condition.kind == NeutronBoundaryKind::Albedo) {
            const int reflected = reflectionDirections_[reflectionIndex(static_cast<int>(faceIndex), direction)];
            if (reflected < 0) throw std::runtime_error("Missing reflected ordinate in SN boundary assembly");
            rhs[owner] -= flux * condition.value * previousAngular[angularIndex(group, reflected, owner)];
        } else {
            throw std::runtime_error("Invalid boundary kind in SN equation");
        }
    }
}

void NeutronSolver::solveDiffusionStep(double eigenvalue, const std::vector<double>& fissionState) {
    std::vector<double> rhs = diffusionRightHandSide(eigenvalue, fissionState);
    const LinearSolveInfo info = system_->solve(rhs, scalarFlux_, linearSolver_, preconditioner_,
        linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
    totalLinearIterations_ += info.iterations;
    lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);
}

void NeutronSolver::updateScalarFlux() {
    std::fill(scalarFlux_.begin(), scalarFlux_.end(), 0.0);
    for (int group = 0; group < groupCount_; ++group) {
        for (int direction = 0; direction < directionCount_; ++direction) {
            const double weight = directions_[direction].weight;
            for (int cell = 0; cell < cellCount_; ++cell) {
                scalarFlux_[scalarIndex(group, cell)] += weight * angularFlux_[angularIndex(group, direction, cell)];
            }
        }
    }
}

double NeutronSolver::fluxChange(const std::vector<double>& previous) const {
    if (previous.size() != scalarFlux_.size()) throw std::runtime_error("Neutron flux-history size mismatch");
    double numerator = 0.0;
    double scale = 0.0;
    for (int cell = 0; cell < cellCount_; ++cell) {
        double difference = 0.0;
        double total = 0.0;
        for (int group = 0; group < groupCount_; ++group) {
            const int index = scalarIndex(group, cell);
            difference += scalarFlux_[index] - previous[index];
            total += scalarFlux_[index];
        }
        numerator += mesh_.cells[cell].volume * difference * difference;
        scale += mesh_.cells[cell].volume * std::abs(total);
    }
    return std::sqrt(numerator) / std::max(scale, 1.0e-30);
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
    for (double& value : angularFlux_) value *= factor;
}

NeutronSolver::Balance NeutronSolver::balance(double eigenvalue, bool eigenproblem) const {
    Balance result;
    result.production = fissionProduction(scalarFlux_);
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        for (int group = 0; group < groupCount_; ++group) {
            const double flux = scalarFlux_[scalarIndex(group, cell)];
            result.absorption += mesh_.cells[cell].volume * material.absorption[group] * flux;
            result.externalSource += mesh_.cells[cell].volume * material.externalSource[group];
            result.power += mesh_.cells[cell].volume * material.kappaSigmaFission[group] * flux;
        }
    }

    if (model_ == NeutronModel::DiscreteOrdinates) {
        for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            for (int group = 0; group < groupCount_; ++group) {
                for (int direction = 0; direction < directionCount_; ++direction) {
                    const NeutronDirection& omega = directions_[direction];
                    const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
                    double faceValue = 0.0;
                    if (flux >= 0.0) {
                        faceValue = angularFlux_[angularIndex(group, direction, owner)];
                    } else if (condition.kind == NeutronBoundaryKind::IncidentFlux) {
                        faceValue = condition.value;
                    } else if (condition.kind == NeutronBoundaryKind::Reflective || condition.kind == NeutronBoundaryKind::Albedo) {
                        const int reflected = reflectionDirections_[reflectionIndex(static_cast<int>(faceIndex), direction)];
                        faceValue = condition.value * angularFlux_[angularIndex(group, reflected, owner)];
                    }
                    const double partial = omega.weight * flux * faceValue;
                    result.leakage += partial;
                    if (partial < 0.0) result.incoming -= partial;
                }
            }
        }
    } else {
        for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            for (int group = 0; group < groupCount_; ++group) {
                const double diffusion = materials_[cellMaterials_[owner]].diffusion[group];
                const double scalar = scalarFlux_[scalarIndex(group, owner)];
                double current = 0.0;
                if (condition.kind == NeutronBoundaryKind::FixedScalarFlux) {
                    current = diffusion * face.area / face.delta * (scalar - condition.value);
                } else if (condition.kind == NeutronBoundaryKind::IncidentCurrent) {
                    const double diagonal = diffusion * face.area / (2.0 * diffusion + face.delta);
                    const double source = 2.0 * diffusion * face.area / (2.0 * diffusion + face.delta);
                    current = diagonal * scalar - source * condition.value;
                } else if (condition.kind == NeutronBoundaryKind::Vacuum ||
                           condition.kind == NeutronBoundaryKind::Marshak ||
                           condition.kind == NeutronBoundaryKind::Albedo) {
                    const double gamma = robinGamma(condition);
                    current = diffusion * gamma * face.area / (diffusion + gamma * face.delta) * scalar;
                }
                result.leakage += current;
                if (current < 0.0) result.incoming -= current;
            }
        }
    }

    const double expected = eigenproblem ? result.production / eigenvalue : result.production + result.externalSource;
    result.residual = expected - result.absorption - result.leakage;
    result.relativeResidual = std::abs(result.residual) /
        std::max({std::abs(expected), std::abs(result.absorption), std::abs(result.leakage), result.incoming, 1.0e-30});
    const double losses = result.absorption + result.leakage;
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

    for (int group = 0; group < groupCount_; ++group) {
        std::vector<Vector2> gradient(cellCount_);
        for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const Vector2 area = face.normal * face.area;
            const double value = face.ownerWeight * scalarFlux_[scalarIndex(group, face.owner)] +
                (1.0 - face.ownerWeight) * scalarFlux_[scalarIndex(group, face.neighbour)];
            gradient[face.owner] += area * value;
            gradient[face.neighbour] -= area * value;
        }
        for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const int owner = face.owner;
            const NeutronBoundaryCondition& condition = boundaries_[faceIndex];
            const double diffusion = materials_[cellMaterials_[owner]].diffusion[group];
            const double scalar = scalarFlux_[scalarIndex(group, owner)];
            double boundaryValue = scalar;
            if (condition.kind == NeutronBoundaryKind::FixedScalarFlux) {
                boundaryValue = condition.value;
            } else if (condition.kind == NeutronBoundaryKind::IncidentCurrent) {
                boundaryValue = (diffusion / face.delta * scalar + condition.value) / (diffusion / face.delta + 0.5);
            } else if (condition.kind == NeutronBoundaryKind::Vacuum ||
                       condition.kind == NeutronBoundaryKind::Marshak ||
                       condition.kind == NeutronBoundaryKind::Albedo) {
                const double gamma = robinGamma(condition);
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
    output << "cell,x,y,material";
    for (int group = 0; group < groupCount_; ++group) output << ",phi_g" << group + 1;
    output << ",phi_total,current_x,current_y,absorption_rate,fission_rate,neutron_production,power_density\n";
    output << std::setprecision(16) << std::scientific;
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& material = materials_[cellMaterials_[cell]];
        output << cell << ',' << mesh_.cells[cell].centre.x << ',' << mesh_.cells[cell].centre.y << ',' << material.name;
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
               << absorption << ',' << fission << ',' << production << ',' << power << '\n';
    }
}

void NeutronSolver::writeSummary(const std::filesystem::path& file) const {
    const bool eigenproblem = problem_ == NeutronProblem::Eigenvalue;
    const Balance current = balance(eigenvalue_, eigenproblem);
    double minimum = std::numeric_limits<double>::max();
    double maximum = std::numeric_limits<double>::lowest();
    double volume = 0.0;
    double fluxIntegral = 0.0;
    for (int cell = 0; cell < cellCount_; ++cell) {
        volume += mesh_.cells[cell].volume;
        for (int group = 0; group < groupCount_; ++group) {
            const double value = scalarFlux_[scalarIndex(group, cell)];
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            fluxIntegral += mesh_.cells[cell].volume * value;
        }
    }
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron summary: " + file.string());
    output << "problem_type,model,energy_groups,directions,cells,outer_iterations,k_eff,k_balance,reactivity,"
              "normalization,volume,flux_integral,flux_min,flux_max,fission_production,absorption,net_leakage,"
              "incoming_current,external_source,fission_power,balance_absolute,balance_relative,"
              "total_linear_iterations,linear_residual\n";
    output << (eigenproblem ? "eigenvalue" : "fixed_source") << ',' << name() << ',' << groupCount_ << ','
           << directionCount_ << ',' << cellCount_ << ',' << history_.size() << ',';
    if (eigenproblem) output << csvNumber(eigenvalue_) << ',' << csvNumber(current.eigenvalue) << ','
                             << csvNumber((eigenvalue_ - 1.0) / eigenvalue_);
    else output << ",,";
    output << ',' << normalization_ << ',' << csvNumber(volume) << ',' << csvNumber(fluxIntegral) << ','
           << csvNumber(minimum) << ',' << csvNumber(maximum) << ',' << csvNumber(current.production) << ','
           << csvNumber(current.absorption) << ',' << csvNumber(current.leakage) << ',' << csvNumber(current.incoming) << ','
           << csvNumber(current.externalSource) << ',' << csvNumber(current.power) << ',' << csvNumber(current.residual) << ','
           << csvNumber(current.relativeResidual) << ',' << totalLinearIterations_ << ',' << csvNumber(lastLinearResidual_) << '\n';
}

void NeutronSolver::writeHistory(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create neutron iteration history: " + file.string());
    output << "iteration,k_eff,flux_change,balance_relative,source_iterations,linear_iterations,linear_residual\n";
    for (const IterationRecord& record : history_) {
        output << record.iteration << ',';
        if (problem_ == NeutronProblem::Eigenvalue) output << csvNumber(record.eigenvalue);
        output << ',' << csvNumber(record.fluxChange) << ',' << csvNumber(record.balance) << ','
               << record.sourceIterations << ',' << record.linearIterations << ',' << csvNumber(record.linearResidual) << '\n';
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

std::filesystem::path NeutronSolver::outputPath(const std::string& suffix) const {
    std::filesystem::path output = config_.require("output");
    if (output.is_relative()) output = caseDirectory_ / output;
    if (suffix.empty()) return output;
    return output.parent_path() / (output.stem().string() + suffix);
}

} // namespace mp
