#include "flow/EnergyModel.hpp"

#include "fvm/Operators.hpp"
#include "fvm/ScalarTransport.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace mp {

EnergyModel::EnergyModel(const Mesh& mesh, const Config& config, bool compressible)
    : mesh_(mesh), config_(config), compressible_(compressible),
      relaxation_(config.number("energy_relaxation", 0.7)),
      turbulentPrandtl_(config.number("turbulent_prandtl", 0.9)),
      heatSource_(config.number("heat_source", 0.0)),
      viscousDissipation_(config.boolean("viscous_dissipation", compressible)),
      linearRelativeTolerance_(config.number("linear_tolerance", 1.0e-9)),
      linearAbsoluteTolerance_(config.number("linear_absolute_tolerance", 1.0e-12)),
      maximumLinearIterations_(config.integer("maximum_linear_iterations", 1000)),
      timeDerivative_(config.get("time_scheme", "steady")),
      temperature_("T", mesh, config.number("initial_temperature", config.number("temperature", 300.0))),
      previousTemperature_("previousT", mesh, temperature_[0]),
      oldTemperature_("oldT", mesh, temperature_[0]),
      oldOldTemperature_("oldOldT", mesh, temperature_[0]),
      totalEnthalpy_("H0", mesh), previousEnthalpy_("previousH0", mesh),
      oldEnthalpy_("oldH0", mesh), oldOldEnthalpy_("oldOldH0", mesh),
      kineticEnergy_("K", mesh), diffusivity_(mesh.cells.size(), 0.0), source_(mesh.cells.size(), 0.0),
      additionalFaceFlux_(mesh.faces.size(), 0.0), destruction_(mesh.cells.size(), 0.0),
      constrainedCells_(mesh.cells.size(), 0), constrainedValues_(mesh.cells.size(), 0.0),
      energySystem_(static_cast<int>(mesh.cells.size())) {
    if (relaxation_ <= 0.0 || relaxation_ > 1.0) throw std::runtime_error("energy_relaxation must be in (0, 1]");
    if (!(turbulentPrandtl_ > 0.0)) throw std::runtime_error("turbulent_prandtl must be positive");
    for (double value : temperature_.values()) {
        if (!(value > 0.0) || !std::isfinite(value)) throw std::runtime_error("initial_temperature must be positive");
    }
    for (const BoundaryPatch& patch : mesh.patches) boundaries_.push_back(readBoundary(patch.name));
}

EnergyBoundaryCondition EnergyModel::readBoundary(const std::string& patch) const {
    const std::string key = "boundary." + patch + ".temperature";
    std::istringstream input(config_.require(key));
    std::string typeName;
    input >> typeName;
    const std::string type = lowercase(typeName);
    EnergyBoundaryCondition condition;
    if (type == "fixedvalue") {
        condition.type = EnergyBoundaryType::FixedValue;
        input >> condition.value;
        if (!(condition.value > 0.0)) throw std::runtime_error("Fixed temperature must be positive for " + key);
    } else if (type == "totaltemperature" || type == "stagnationtemperature" ||
               type == "total" || type == "stagnation") {
        if (!compressible_) {
            throw std::runtime_error("totalTemperature boundary requires compressible flow for " + key);
        }
        condition.type = EnergyBoundaryType::TotalTemperature;
        input >> condition.value;
        if (!(condition.value > 0.0)) throw std::runtime_error("Total temperature must be positive for " + key);
    } else if (type == "zerogradient" || type == "adiabatic") {
        condition.type = EnergyBoundaryType::ZeroGradient;
        condition.value = 0.0;
    } else if (type == "fixedheatflux" || type == "heatflux") {
        condition.type = EnergyBoundaryType::FixedHeatFlux;
        input >> condition.value;
    } else {
        throw std::runtime_error("Unknown temperature boundary condition for " + key);
    }
    if (!input || !(input >> std::ws).eof()) throw std::runtime_error("Invalid temperature boundary data for " + key);
    return condition;
}

void EnergyModel::initialize(const VectorField& velocity, const std::vector<double>& heatCapacity) {
    if (initialized_) return;
    updateKineticEnergy(velocity);
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        if (!(heatCapacity[cell] > 0.0)) throw std::runtime_error("Energy transport requires positive heat capacity");
        totalEnthalpy_[cell] = heatCapacity[cell] * temperature_[cell] + kineticEnergy_[cell];
    }
    updateBoundaries(velocity, heatCapacity);
    previousTemperature_.values() = temperature_.values();
    oldTemperature_.values() = temperature_.values();
    oldOldTemperature_.values() = temperature_.values();
    previousEnthalpy_.values() = totalEnthalpy_.values();
    oldEnthalpy_.values() = totalEnthalpy_.values();
    oldOldEnthalpy_.values() = totalEnthalpy_.values();
    initialized_ = true;
}

void EnergyModel::updateKineticEnergy(const VectorField& velocity) {
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        kineticEnergy_[cell] = 0.5 * normSquared(velocity[cell]);
    }
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        kineticEnergy_.boundary(face) = 0.5 * normSquared(velocity.boundary(face));
        kineticEnergy_.boundaryDirichlet(face) = true;
    }
}

bool EnergyModel::isTotalTemperature(std::size_t faceIndex) const {
    return boundaries_[mesh_.faces[faceIndex].patch].type == EnergyBoundaryType::TotalTemperature;
}

double EnergyModel::boundaryValue(std::size_t faceIndex) const {
    return boundaries_[mesh_.faces[faceIndex].patch].value;
}

void EnergyModel::updateBoundaries(const VectorField& velocity, const std::vector<double>& heatCapacity) {
    updateKineticEnergy(velocity);
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const int owner = mesh_.faces[face].owner;
        const EnergyBoundaryCondition& condition = boundaries_[mesh_.faces[face].patch];
        if (condition.type == EnergyBoundaryType::TotalTemperature) {
            const double staticTemperature = condition.value - kineticEnergy_.boundary(face) / heatCapacity[owner];
            temperature_.boundary(face) = std::max(1.0, staticTemperature);
            temperature_.boundaryDirichlet(face) = true;
            totalEnthalpy_.boundary(face) = heatCapacity[owner] * condition.value;
            totalEnthalpy_.boundaryDirichlet(face) = true;
            continue;
        }
        const bool fixedTemperature = condition.type == EnergyBoundaryType::FixedValue;
        temperature_.boundaryDirichlet(face) = fixedTemperature;
        temperature_.boundary(face) = fixedTemperature ? condition.value : temperature_[owner];
        if (compressible_) {
            totalEnthalpy_.boundary(face) = heatCapacity[owner] * temperature_.boundary(face) +
                kineticEnergy_.boundary(face);
            totalEnthalpy_.boundaryDirichlet(face) = fixedTemperature;
        }
    }
}

void EnergyModel::beginIteration() {
    previousTemperature_.values() = temperature_.values();
    previousEnthalpy_.values() = totalEnthalpy_.values();
}

LinearSolveInfo EnergyModel::solve(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                                   const GradientScheme& gradient, const std::vector<double>& density,
                                   const std::vector<double>& oldDensity, const std::vector<double>& oldOldDensity,
                                   const std::vector<double>& viscosity, const std::vector<double>& conductivity,
                                   const std::vector<double>& heatCapacity, const ScalarField* turbulentViscosity,
                                   const VectorField& velocity, const ScalarField& pressure,
                                   const ScalarField& oldPressure, const ScalarField& oldOldPressure,
                                   const Vector2& bodyForce) {
    initialize(velocity, heatCapacity);
    updateBoundaries(velocity, heatCapacity);
    const auto velocityGradient = compressible_ && viscousDissipation_
        ? gradient.vector(velocity) : std::vector<VectorGradient>{};
    const auto kineticGradient = compressible_
        ? gradient.scalar(kineticEnergy_) : std::vector<ScalarGradient>{};
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOld);
    std::fill(source_.begin(), source_.end(), 0.0);
    std::fill(additionalFaceFlux_.begin(), additionalFaceFlux_.end(), 0.0);
    std::vector<std::array<Vector2, 2>> viscousStress(mesh_.cells.size());
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        const double cp = heatCapacity[cell];
        if (!(cp > 0.0) || !(density[cell] > 0.0) || conductivity[cell] < 0.0) {
            throw std::runtime_error("Energy transport requires positive density and heat capacity, and non-negative conductivity");
        }
        const double eddy = turbulentViscosity ? (*turbulentViscosity)[cell] : 0.0;
        diffusivity_[cell] = conductivity[cell] / cp + eddy / turbulentPrandtl_;
        double energySource = heatSource_;
        if (compressible_) {
            const double pressureRate = time.current * pressure[cell] - time.previous * oldPressure[cell]
                - time.previousPrevious * oldOldPressure[cell];
            energySource += pressureRate + density[cell] * dot(bodyForce, velocity[cell]);
            if (viscousDissipation_) {
                const VectorGradient& derivatives = velocityGradient[cell];
                const double divergence = derivatives[0].x + derivatives[1].y;
                const double shear = derivatives[0].y + derivatives[1].x;
                viscousStress[cell][0] = viscosity[cell] *
                    Vector2{2.0 * derivatives[0].x - (2.0 / 3.0) * divergence, shear};
                viscousStress[cell][1] = viscosity[cell] *
                    Vector2{shear, 2.0 * derivatives[1].y - (2.0 / 3.0) * divergence};
            }
        }
        source_[cell] = compressible_ ? energySource : energySource / cp;
    }
    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount() && compressible_; ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        const double faceDiffusivity = face.ownerWeight * diffusivity_[owner] +
            (1.0 - face.ownerWeight) * diffusivity_[neighbour];
        const ScalarGradient faceKineticGradient = face.ownerWeight * kineticGradient[owner] +
            (1.0 - face.ownerWeight) * kineticGradient[neighbour];
        double flux = faceDiffusivity * (face.area / face.delta * (kineticEnergy_[neighbour] - kineticEnergy_[owner]) +
            dot(nonOrthogonalVector(mesh_, faceIndex), faceKineticGradient));
        if (viscousDissipation_) {
            std::array<Vector2, 2> faceStress;
            for (int row = 0; row < 2; ++row) {
                faceStress[row] = face.ownerWeight * viscousStress[owner][row] +
                    (1.0 - face.ownerWeight) * viscousStress[neighbour][row];
            }
            const Vector2 faceVelocity = interpolate(mesh_, velocity, faceIndex);
            const Vector2 work{dot(faceStress[0], faceVelocity), dot(faceStress[1], faceVelocity)};
            flux -= dot(work, faceVector(mesh_, faceIndex));
        }
        additionalFaceFlux_[faceIndex] = flux;
    }
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const EnergyBoundaryCondition& condition = boundaries_[mesh_.faces[face].patch];
        const int owner = mesh_.faces[face].owner;
        if (condition.type == EnergyBoundaryType::FixedHeatFlux) {
            source_[owner] += condition.value * mesh_.faces[face].area / mesh_.cells[owner].volume /
                (compressible_ ? 1.0 : heatCapacity[owner]);
        }
        if (!compressible_) continue;
        if (condition.type == EnergyBoundaryType::FixedValue ||
            condition.type == EnergyBoundaryType::TotalTemperature) {
            additionalFaceFlux_[face] = diffusivity_[owner] *
                (mesh_.faces[face].area / mesh_.faces[face].delta *
                    (kineticEnergy_.boundary(face) - kineticEnergy_[owner]) +
                 dot(nonOrthogonalVector(mesh_, face), kineticGradient[owner]));
        }
        if (viscousDissipation_ && velocity.boundaryDirichlet(face)) {
            const Vector2 boundaryVelocity = velocity.boundary(face);
            const Vector2 work{dot(viscousStress[owner][0], boundaryVelocity),
                dot(viscousStress[owner][1], boundaryVelocity)};
            additionalFaceFlux_[face] -= dot(work, faceVector(mesh_, face));
        }
    }

    ScalarField& field = compressible_ ? totalEnthalpy_ : temperature_;
    const ScalarField& previous = compressible_ ? previousEnthalpy_ : previousTemperature_;
    const ScalarField& old = compressible_ ? oldEnthalpy_ : oldTemperature_;
    const ScalarField& oldOld = compressible_ ? oldOldEnthalpy_ : oldOldTemperature_;
    ScalarTransportEquation equation{field, previous, old, oldOld,
        massFlux, diffusivity_, source_, destruction_, constrainedCells_, constrainedValues_, 1.0, time,
        compressible_ ? 1.0 : relaxation_, linearRelativeTolerance_, linearAbsoluteTolerance_,
        maximumLinearIterations_, config_.get("energy_ksp", "bcgs"), config_.get("energy_pc", "ilu")};
    equation.gradientLimiter = readGradientLimiter(config_);
    equation.secondOrderConvection = config_.boolean("energy_second_order", true);
    equation.additionalFaceFlux = compressible_ ? &additionalFaceFlux_ : nullptr;
    equation.densityField = &density;
    equation.oldDensityField = &oldDensity;
    equation.oldOldDensityField = &oldOldDensity;
    const LinearSolveInfo info = solveScalarTransport(mesh_, gradient, energySystem_, equation);
    change_ = 0.0;
    temperatureLimited_ = false;
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        if (compressible_) {
            const double recovered = (totalEnthalpy_[cell] - kineticEnergy_[cell]) / heatCapacity[cell];
            if (!std::isfinite(recovered)) {
                throw std::runtime_error("Energy equation produced a non-finite temperature");
            }
            double relaxed = previousTemperature_[cell] + relaxation_ * (recovered - previousTemperature_[cell]);
            const double minimumTemperature = std::max(1.0, 0.5 * previousTemperature_[cell]);
            const double maximumTemperature = 1.5 * previousTemperature_[cell];
            const double bounded = std::clamp(relaxed, minimumTemperature, maximumTemperature);
            temperatureLimited_ = temperatureLimited_ || bounded != relaxed;
            relaxed = bounded;
            temperature_[cell] = relaxed;
            totalEnthalpy_[cell] = heatCapacity[cell] * relaxed + kineticEnergy_[cell];
        }
        if (!(temperature_[cell] > 0.0) || !std::isfinite(temperature_[cell])) {
            throw std::runtime_error("Energy equation produced a non-positive or non-finite temperature");
        }
        change_ = std::max(change_, std::abs(temperature_[cell] - previousTemperature_[cell]));
    }
    updateBoundaries(velocity, heatCapacity);
    linearIterations_ = info.iterations;
    linearResidual_ = info.residual;
    return info;
}

void EnergyModel::updateTimeHistory() {
    oldOldTemperature_.values() = oldTemperature_.values();
    oldTemperature_.values() = temperature_.values();
    oldOldEnthalpy_.values() = oldEnthalpy_.values();
    oldEnthalpy_.values() = totalEnthalpy_.values();
}

} // namespace mp
