#include "flow/TurbulenceModel.hpp"

#include "fvm/Operators.hpp"
#include "fvm/ScalarTransport.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace mp {
namespace {

constexpr double cMu = 0.09;
constexpr double cEpsilon1 = 1.44;
constexpr double cEpsilon2 = 1.92;
constexpr double sigmaKEpsilon = 1.0;
constexpr double sigmaEpsilon = 1.3;
constexpr double betaStar = 0.09;
constexpr double kappa = 0.41;
constexpr double wallE = 9.8;
constexpr double a1 = 0.31;
constexpr double sigmaK1 = 0.85;
constexpr double sigmaK2 = 1.0;
constexpr double sigmaOmega1 = 0.5;
constexpr double sigmaOmega2 = 0.856;
constexpr double beta1 = 0.075;
constexpr double beta2 = 0.0828;
constexpr double gamma1 = 5.0 / 9.0;
constexpr double gamma2 = 0.44;

double pointSegmentDistance(const Vector2& point, const Vector2& a, const Vector2& b) {
    const Vector2 segment = b - a;
    const double lengthSquared = normSquared(segment);
    if (lengthSquared <= 0.0) return norm(point - a);
    const double fraction = std::clamp(dot(point - a, segment) / lengthSquared, 0.0, 1.0);
    return norm(point - (a + segment * fraction));
}

double normalizedUpdate(const ScalarField& field, const ScalarField& previous, double floor) {
    double scale = floor;
    double update = 0.0;
    for (std::size_t cell = 0; cell < field.size(); ++cell) {
        scale = std::max(scale, std::abs(field[cell]));
        update = std::max(update, std::abs(field[cell] - previous[cell]));
    }
    return update / scale;
}

} // namespace

TurbulenceModel::TurbulenceModel(const Mesh& mesh, const Config& config,
                                 const FlowBoundaryConditions& boundaryConditions,
                                 double density, double viscosity, TurbulenceModelKind kind)
    : mesh_(mesh), config_(config), boundaryConditions_(boundaryConditions), density_(density),
      viscosity_(viscosity), kind_(kind), relaxation_(config.number("turbulence_relaxation", 0.5)),
      minimumK_(config.number("minimum_k", 1.0e-12)),
      minimumEpsilon_(config.number("minimum_epsilon", 1.0e-12)),
      minimumOmega_(config.number("minimum_omega", 1.0e-8)),
      maximumEddyViscosityRatio_(config.number("maximum_eddy_viscosity_ratio", 1000.0)),
      k_("k", mesh, config.number("initial_k", 1.0e-6)),
      epsilon_("epsilon", mesh, config.number("initial_epsilon", 1.0e-6)),
      omega_("omega", mesh, config.number("initial_omega", 1.0)),
      previousK_("previousK", mesh), previousEpsilon_("previousEpsilon", mesh),
      previousOmega_("previousOmega", mesh), oldK_("oldK", mesh), oldEpsilon_("oldEpsilon", mesh),
      oldOmega_("oldOmega", mesh), oldOldK_("oldOldK", mesh), oldOldEpsilon_("oldOldEpsilon", mesh),
      oldOldOmega_("oldOldOmega", mesh), eddyViscosity_("mu_t", mesh, 0.0),
      blendingF1_("F1", mesh, kind == TurbulenceModelKind::KEpsilon ? 0.0 : 1.0),
      blendingF2_("F2", mesh, 0.0), wallDistance_(mesh.cells.size(), 1.0e30),
      production_(mesh.cells.size(), 0.0), strainInvariant_(mesh.cells.size(), 0.0),
      sigmaK_(mesh.cells.size(), kind == TurbulenceModelKind::KEpsilon ? sigmaKEpsilon : sigmaK1),
      sigmaOmega_(mesh.cells.size(), kind == TurbulenceModelKind::KEpsilon ? sigmaEpsilon : sigmaOmega1),
      beta_(mesh.cells.size(), beta1), gamma_(mesh.cells.size(), gamma1),
      constrainedCells_(mesh.cells.size(), 0), constrainedValues_(mesh.cells.size(), 0.0),
      kSystem_(static_cast<int>(mesh.cells.size())), scaleSystem_(static_cast<int>(mesh.cells.size())) {
    if (relaxation_ <= 0.0 || relaxation_ > 1.0) {
        throw std::runtime_error("turbulence_relaxation must be in (0, 1]");
    }
    if (minimumK_ <= 0.0 || minimumEpsilon_ <= 0.0 || minimumOmega_ <= 0.0 ||
        maximumEddyViscosityRatio_ <= 0.0) {
        throw std::runtime_error("Turbulence bounds and viscosity-ratio limit must be positive");
    }
    if (k_[0] < minimumK_ || epsilon_[0] < minimumEpsilon_ || omega_[0] < minimumOmega_) {
        throw std::runtime_error("Initial turbulence values must exceed their lower bounds");
    }

    for (const BoundaryPatch& patch : mesh_.patches) {
        const bool noSlip = boundaryConditions_.noSlip(patch.faces.front());
        kBoundaries_.push_back(readBoundary(patch.name, "k", noSlip));
        if (kind_ == TurbulenceModelKind::KEpsilon) {
            epsilonBoundaries_.push_back(readBoundary(patch.name, "epsilon", noSlip));
        } else {
            omegaBoundaries_.push_back(readBoundary(patch.name, "omega", noSlip));
        }
    }

    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
            if (!boundaryConditions_.noSlip(face)) continue;
            const Face& wall = mesh_.faces[face];
            wallDistance_[cell] = std::min(wallDistance_[cell],
                pointSegmentDistance(mesh_.cells[cell].centre, mesh_.nodes[wall.node0], mesh_.nodes[wall.node1]));
        }
        wallDistance_[cell] = std::max(wallDistance_[cell], 1.0e-12);
        eddyViscosity_[cell] = kind_ == TurbulenceModelKind::KEpsilon
            ? std::min(cMu * density_ * k_[cell] * k_[cell] / epsilon_[cell],
                       maximumEddyViscosityRatio_ * viscosity_)
            : density_ * k_[cell] / omega_[cell];
    }
    previousK_.values() = k_.values();
    previousEpsilon_.values() = epsilon_.values();
    previousOmega_.values() = omega_.values();
    oldK_.values() = k_.values();
    oldEpsilon_.values() = epsilon_.values();
    oldOmega_.values() = omega_.values();
    oldOldK_.values() = k_.values();
    oldOldEpsilon_.values() = epsilon_.values();
    oldOldOmega_.values() = omega_.values();
    updateBoundaries();
}

std::string TurbulenceModel::name() const {
    return kind_ == TurbulenceModelKind::KEpsilon ? "standard-k-epsilon-wall-function" : "menter-sst-2003";
}

TurbulenceBoundaryCondition TurbulenceModel::readBoundary(const std::string& patch,
                                                           const std::string& variable,
                                                           bool noSlip) const {
    const std::string key = "boundary." + patch + "." + variable;
    std::istringstream input(config_.require(key));
    std::string typeName;
    input >> typeName;
    const std::string type = lowercase(typeName);
    TurbulenceBoundaryCondition condition;
    if (type == "fixedvalue") {
        condition.type = TurbulenceBoundaryType::FixedValue;
        input >> condition.value;
        if (!input || (variable == "k" ? condition.value < 0.0 : condition.value <= 0.0)) {
            throw std::runtime_error("Invalid turbulence boundary data for " + key);
        }
    } else if (type == "zerogradient") {
        condition.type = TurbulenceBoundaryType::ZeroGradient;
    } else if (type == "wallfunction") {
        condition.type = TurbulenceBoundaryType::WallFunction;
    } else {
        throw std::runtime_error("Unknown turbulence boundary condition for " + key);
    }

    if (condition.type == TurbulenceBoundaryType::WallFunction &&
        !(noSlip && ((kind_ == TurbulenceModelKind::KEpsilon && variable == "epsilon") ||
                     (kind_ == TurbulenceModelKind::SSTMenter2003 && variable == "omega")))) {
        throw std::runtime_error("wallFunction is not valid for " + key);
    }
    if (noSlip && kind_ == TurbulenceModelKind::KEpsilon && variable == "k" &&
        condition.type != TurbulenceBoundaryType::ZeroGradient) {
        throw std::runtime_error("Standard k-epsilon requires zeroGradient k at no-slip patch " + patch);
    }
    if (noSlip && kind_ == TurbulenceModelKind::KEpsilon && variable == "epsilon" &&
        condition.type != TurbulenceBoundaryType::WallFunction) {
        throw std::runtime_error("Standard k-epsilon requires wallFunction epsilon at no-slip patch " + patch);
    }
    if (noSlip && kind_ == TurbulenceModelKind::SSTMenter2003 && variable == "omega" &&
        condition.type != TurbulenceBoundaryType::WallFunction) {
        throw std::runtime_error("SST requires wallFunction omega at no-slip patch " + patch);
    }
    if (noSlip && kind_ == TurbulenceModelKind::SSTMenter2003 && variable == "k" &&
        (condition.type != TurbulenceBoundaryType::FixedValue || condition.value != 0.0)) {
        throw std::runtime_error("SST requires fixedValue 0 for k at no-slip patch " + patch);
    }
    return condition;
}

void TurbulenceModel::updateFieldBoundaries(
    ScalarField& field, const std::vector<TurbulenceBoundaryCondition>& conditions) {
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const TurbulenceBoundaryCondition& condition = conditions[face.patch];
        switch (condition.type) {
        case TurbulenceBoundaryType::FixedValue:
            field.boundary(faceIndex) = condition.value;
            field.boundaryDirichlet(faceIndex) = true;
            break;
        case TurbulenceBoundaryType::ZeroGradient:
            field.boundary(faceIndex) = field[face.owner];
            field.boundaryDirichlet(faceIndex) = false;
            break;
        case TurbulenceBoundaryType::WallFunction:
            if (kind_ == TurbulenceModelKind::SSTMenter2003) {
                const double y = face.delta;
                field.boundary(faceIndex) = 6.0 * viscosity_ / (density_ * beta1 * y * y);
                field.boundaryDirichlet(faceIndex) = true;
            } else {
                field.boundary(faceIndex) = field[face.owner];
                field.boundaryDirichlet(faceIndex) = false;
            }
            break;
        }
    }
}

void TurbulenceModel::updateBoundaries() {
    updateFieldBoundaries(k_, kBoundaries_);
    if (kind_ == TurbulenceModelKind::KEpsilon) updateFieldBoundaries(epsilon_, epsilonBoundaries_);
    else updateFieldBoundaries(omega_, omegaBoundaries_);
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const int owner = mesh_.faces[face].owner;
        eddyViscosity_.boundary(face) = kind_ == TurbulenceModelKind::SSTMenter2003 &&
            boundaryConditions_.noSlip(face) ? 0.0 : eddyViscosity_[owner];
        blendingF1_.boundary(face) = blendingF1_[owner];
        blendingF2_.boundary(face) = blendingF2_[owner];
    }
}

void TurbulenceModel::beginIteration() {
    previousK_.values() = k_.values();
    previousEpsilon_.values() = epsilon_.values();
    previousOmega_.values() = omega_.values();
}

void TurbulenceModel::computeProduction(const std::vector<VectorGradient>& velocityGradient) {
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        const double duDx = velocityGradient[cell][0].x;
        const double duDy = velocityGradient[cell][0].y;
        const double dvDx = velocityGradient[cell][1].x;
        const double dvDy = velocityGradient[cell][1].y;
        const double divergence = duDx + dvDy;
        const double strainSquared = 2.0 * duDx * duDx + 2.0 * dvDy * dvDy +
            (duDy + dvDx) * (duDy + dvDx);
        strainInvariant_[cell] = std::sqrt(std::max(strainSquared, 0.0));
        production_[cell] = eddyViscosity_[cell] * strainSquared -
            (2.0 / 3.0) * density_ * k_[cell] * divergence;
        production_[cell] = std::max(production_[cell], 0.0);
    }
}

void TurbulenceModel::updateSSTProperties(const GradientScheme& gradient,
                                          const std::vector<VectorGradient>& velocityGradient) {
    const auto kGradient = gradient.scalar(k_);
    const auto omegaGradient = gradient.scalar(omega_);
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        const double k = std::max(k_[cell], minimumK_);
        const double omega = std::max(omega_[cell], minimumOmega_);
        const double distance = wallDistance_[cell];
        const double kinematicViscosity = viscosity_ / density_;
        const double crossDiffusion = 2.0 * density_ * sigmaOmega2 *
            dot(kGradient[cell], omegaGradient[cell]) / omega;
        const double cd = std::max(crossDiffusion, 1.0e-10);
        const double first = std::max(std::sqrt(k) / (betaStar * omega * distance),
                                      500.0 * kinematicViscosity / (distance * distance * omega));
        const double second = 4.0 * density_ * sigmaOmega2 * k / (cd * distance * distance);
        const double argument1 = std::min(first, second);
        const double argument2 = std::max(2.0 * std::sqrt(k) / (betaStar * omega * distance),
                                          500.0 * kinematicViscosity / (distance * distance * omega));
        blendingF1_[cell] = std::tanh(std::pow(argument1, 4));
        blendingF2_[cell] = std::tanh(argument2 * argument2);
        const double f1 = blendingF1_[cell];
        sigmaK_[cell] = f1 * sigmaK1 + (1.0 - f1) * sigmaK2;
        sigmaOmega_[cell] = f1 * sigmaOmega1 + (1.0 - f1) * sigmaOmega2;
        beta_[cell] = f1 * beta1 + (1.0 - f1) * beta2;
        gamma_[cell] = f1 * gamma1 + (1.0 - f1) * gamma2;

        const double duDx = velocityGradient[cell][0].x;
        const double duDy = velocityGradient[cell][0].y;
        const double dvDx = velocityGradient[cell][1].x;
        const double dvDy = velocityGradient[cell][1].y;
        const double strainSquared = 2.0 * duDx * duDx + 2.0 * dvDy * dvDy +
            (duDy + dvDx) * (duDy + dvDx);
        strainInvariant_[cell] = std::sqrt(std::max(strainSquared, 0.0));
        eddyViscosity_[cell] = density_ * a1 * k /
            std::max(a1 * omega, strainInvariant_[cell] * blendingF2_[cell]);
    }
}

LinearSolveInfo TurbulenceModel::solveField(
    ScalarField& field, ScalarField& previous, ScalarField& old, ScalarField& oldOld,
    PetscLinearSystem& system, const std::vector<double>& massFlux, const GradientScheme& gradient,
    const std::vector<double>& diffusivity, const std::vector<double>& source,
    const std::vector<double>& destruction, const std::vector<char>& constrained,
    const std::vector<double>& constrainedValues, double dt, bool hasOldOld) {
    const TimeDerivative timeDerivative(config_.get("time_scheme", "steady"));
    ScalarTransportEquation equation{field, previous, old, oldOld, massFlux, diffusivity, source,
        destruction, constrained, constrainedValues, density_, TimeCoefficients{}, relaxation_,
        config_.number("linear_tolerance", 1.0e-9),
        config_.number("linear_absolute_tolerance", 1.0e-12),
        config_.integer("maximum_linear_iterations", 1000),
        config_.get("turbulence_ksp", "bcgs"), config_.get("turbulence_pc", "ilu")};
    equation.gradientLimiter = readGradientLimiter(config_);
    equation.time = timeDerivative.coefficients(dt, hasOldOld);
    return solveScalarTransport(mesh_, gradient, system, equation);
}

void TurbulenceModel::clip(ScalarField& field, double floor) {
    for (double& value : field.values()) {
        if (!std::isfinite(value)) throw std::runtime_error("Non-finite value in turbulence field " + field.name());
        value = std::max(value, floor);
    }
}

void TurbulenceModel::solveKEpsilon(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                                    const GradientScheme& gradient) {
    const std::size_t count = mesh_.cells.size();
    std::vector<double> diffusivityK(count), diffusivityEpsilon(count);
    std::vector<double> destructionK(count), destructionEpsilon(count), sourceEpsilon(count);
    std::fill(constrainedCells_.begin(), constrainedCells_.end(), 0);
    std::fill(constrainedValues_.begin(), constrainedValues_.end(), 0.0);

    for (std::size_t cell = 0; cell < count; ++cell) {
        diffusivityK[cell] = viscosity_ + eddyViscosity_[cell] / sigmaKEpsilon;
        diffusivityEpsilon[cell] = viscosity_ + eddyViscosity_[cell] / sigmaEpsilon;
        destructionK[cell] = density_ * epsilon_[cell] / std::max(k_[cell], minimumK_);
    }
    const LinearSolveInfo kInfo = solveField(k_, previousK_, oldK_, oldOldK_, kSystem_, massFlux,
        gradient, diffusivityK, production_, destructionK, constrainedCells_, constrainedValues_, dt, hasOldOld);
    clip(k_, minimumK_);
    updateBoundaries();

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        if (!boundaryConditions_.noSlip(faceIndex)) continue;
        const Face& face = mesh_.faces[faceIndex];
        const double epsilon = std::pow(cMu, 0.75) *
            std::pow(std::max(k_[face.owner], minimumK_), 1.5) / (kappa * face.delta);
        if (!constrainedCells_[face.owner] || epsilon > constrainedValues_[face.owner]) {
            constrainedCells_[face.owner] = 1;
            constrainedValues_[face.owner] = epsilon;
        }
    }
    for (std::size_t cell = 0; cell < count; ++cell) {
        const double k = std::max(k_[cell], minimumK_);
        sourceEpsilon[cell] = cEpsilon1 * epsilon_[cell] * production_[cell] / k;
        destructionEpsilon[cell] = cEpsilon2 * density_ * epsilon_[cell] / k;
    }
    const LinearSolveInfo epsilonInfo = solveField(epsilon_, previousEpsilon_, oldEpsilon_, oldOldEpsilon_,
        scaleSystem_, massFlux, gradient, diffusivityEpsilon, sourceEpsilon, destructionEpsilon,
        constrainedCells_, constrainedValues_, dt, hasOldOld);
    clip(epsilon_, minimumEpsilon_);
    for (std::size_t cell = 0; cell < count; ++cell) {
        eddyViscosity_[cell] = std::min(cMu * density_ * k_[cell] * k_[cell] / epsilon_[cell],
                                        maximumEddyViscosityRatio_ * viscosity_);
    }
    linearIterations_ = kInfo.iterations + epsilonInfo.iterations;
    linearResidual_ = std::max(kInfo.residual, epsilonInfo.residual);
    kLinearResidual_ = kInfo.residual;
    scaleLinearResidual_ = epsilonInfo.residual;
}

void TurbulenceModel::solveSST(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                               const GradientScheme& gradient) {
    const std::size_t count = mesh_.cells.size();
    const auto kGradient = gradient.scalar(k_);
    const auto omegaGradient = gradient.scalar(omega_);
    std::vector<double> diffusivityK(count), diffusivityOmega(count);
    std::vector<double> limitedProduction(count), destructionK(count);
    std::vector<double> sourceOmega(count), destructionOmega(count), crossDiffusion(count);
    std::fill(constrainedCells_.begin(), constrainedCells_.end(), 0);
    std::fill(constrainedValues_.begin(), constrainedValues_.end(), 0.0);

    for (std::size_t cell = 0; cell < count; ++cell) {
        const double k = std::max(k_[cell], minimumK_);
        const double omega = std::max(omega_[cell], minimumOmega_);
        limitedProduction[cell] = std::min(production_[cell], 10.0 * betaStar * density_ * k * omega);
        diffusivityK[cell] = viscosity_ + sigmaK_[cell] * eddyViscosity_[cell];
        diffusivityOmega[cell] = viscosity_ + sigmaOmega_[cell] * eddyViscosity_[cell];
        destructionK[cell] = betaStar * density_ * omega;
        crossDiffusion[cell] = 2.0 * (1.0 - blendingF1_[cell]) * density_ * sigmaOmega2 *
            dot(kGradient[cell], omegaGradient[cell]) / omega;
    }
    const LinearSolveInfo kInfo = solveField(k_, previousK_, oldK_, oldOldK_, kSystem_, massFlux,
        gradient, diffusivityK, limitedProduction, destructionK, constrainedCells_, constrainedValues_,
        dt, hasOldOld);
    clip(k_, minimumK_);
    updateBoundaries();

    for (std::size_t cell = 0; cell < count; ++cell) {
        sourceOmega[cell] = crossDiffusion[cell];
        if (eddyViscosity_[cell] > 1.0e-30) {
            sourceOmega[cell] += gamma_[cell] * density_ * limitedProduction[cell] / eddyViscosity_[cell];
        }
        destructionOmega[cell] = beta_[cell] * density_ * omega_[cell];
    }
    const LinearSolveInfo omegaInfo = solveField(omega_, previousOmega_, oldOmega_, oldOldOmega_,
        scaleSystem_, massFlux, gradient, diffusivityOmega, sourceOmega, destructionOmega,
        constrainedCells_, constrainedValues_, dt, hasOldOld);
    clip(omega_, minimumOmega_);
    linearIterations_ = kInfo.iterations + omegaInfo.iterations;
    linearResidual_ = std::max(kInfo.residual, omegaInfo.residual);
    kLinearResidual_ = kInfo.residual;
    scaleLinearResidual_ = omegaInfo.residual;
}

void TurbulenceModel::solve(double dt, bool hasOldOld, const VectorField& velocity,
                            const std::vector<double>& massFlux, const GradientScheme& gradient) {
    const auto velocityGradient = gradient.vector(velocity);
    if (kind_ == TurbulenceModelKind::SSTMenter2003) updateSSTProperties(gradient, velocityGradient);
    computeProduction(velocityGradient);
    if (kind_ == TurbulenceModelKind::KEpsilon) {
        for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
            if (!boundaryConditions_.noSlip(faceIndex)) continue;
            const Face& face = mesh_.faces[faceIndex];
            const double shear = norm(wallShearStress(faceIndex, velocity[face.owner]));
            const double frictionVelocity = std::sqrt(shear / density_);
            production_[face.owner] = shear * frictionVelocity / (kappa * face.delta);
        }
        solveKEpsilon(dt, hasOldOld, massFlux, gradient);
    } else {
        solveSST(dt, hasOldOld, massFlux, gradient);
        updateSSTProperties(gradient, velocityGradient);
    }
    updateBoundaries();
    updateChange();
}

void TurbulenceModel::updateChange() {
    change_ = normalizedUpdate(k_, previousK_, minimumK_);
    change_ = std::max(change_, kind_ == TurbulenceModelKind::KEpsilon
        ? normalizedUpdate(epsilon_, previousEpsilon_, minimumEpsilon_)
        : normalizedUpdate(omega_, previousOmega_, minimumOmega_));
}

double TurbulenceModel::wallFunctionViscosity(std::size_t faceIndex, const Vector2& cellVelocity) const {
    (void)cellVelocity;
    const Face& face = mesh_.faces[faceIndex];
    const double uTau = std::pow(cMu, 0.25) * std::sqrt(std::max(k_[face.owner], minimumK_));
    const double yPlus = density_ * uTau * face.delta / viscosity_;
    if (yPlus <= 11.53) return viscosity_;
    const double uPlus = std::log(wallE * yPlus) / kappa;
    return density_ * uTau * face.delta / std::max(uPlus, 1.0e-12);
}

Vector2 TurbulenceModel::wallShearStress(std::size_t faceIndex, const Vector2& cellVelocity) const {
    const Face& face = mesh_.faces[faceIndex];
    const Vector2 tangentialVelocity = cellVelocity - face.normal * dot(cellVelocity, face.normal);
    const double wallViscosity = usesStandardWallFunction()
        ? wallFunctionViscosity(faceIndex, cellVelocity) : viscosity_;
    return tangentialVelocity * (wallViscosity / face.delta);
}

TurbulenceWallDiagnostics TurbulenceModel::wallDiagnostics(std::size_t faceIndex,
                                                            const Vector2& cellVelocity) const {
    const Face& face = mesh_.faces[faceIndex];
    TurbulenceWallDiagnostics result;
    result.face = faceIndex;
    result.owner = face.owner;
    result.distance = face.delta;
    result.shear = wallShearStress(faceIndex, cellVelocity);
    result.frictionVelocity = std::sqrt(norm(result.shear) / density_);
    result.yPlus = density_ * result.frictionVelocity * face.delta / viscosity_;
    return result;
}

void TurbulenceModel::updateTimeHistory(bool& hasOldOld) {
    oldOldK_.values() = oldK_.values();
    oldK_.values() = k_.values();
    oldOldEpsilon_.values() = oldEpsilon_.values();
    oldEpsilon_.values() = epsilon_.values();
    oldOldOmega_.values() = oldOmega_.values();
    oldOmega_.values() = omega_.values();
    hasOldOld = true;
}

std::unique_ptr<TurbulenceModel> makeTurbulenceModel(const Mesh& mesh, const Config& config,
                                                     const FlowBoundaryConditions& boundaryConditions,
                                                     double density, double viscosity) {
    const std::string model = lowercase(config.get("turbulence_model", "laminar"));
    if (model == "laminar" || model == "none" || model == "off") return nullptr;
    if (model == "k_epsilon" || model == "k-epsilon" || model == "kepsilon" ||
        model == "standard_k_epsilon") {
        return std::make_unique<TurbulenceModel>(mesh, config, boundaryConditions, density, viscosity,
                                                 TurbulenceModelKind::KEpsilon);
    }
    if (model == "sst" || model == "sst_k_omega" || model == "sst-k-omega" ||
        model == "k_omega_sst" || model == "sst_2003") {
        return std::make_unique<TurbulenceModel>(mesh, config, boundaryConditions, density, viscosity,
                                                 TurbulenceModelKind::SSTMenter2003);
    }
    throw std::runtime_error("Unknown turbulence_model: " + model);
}

} // namespace mp
