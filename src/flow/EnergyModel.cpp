#include "flow/EnergyModel.hpp"

#include "fvm/Operators.hpp"
#include "fvm/ScalarTransport.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace mp {

EnergyModel::EnergyModel(const Mesh& mesh, const Config& config, bool compressible)
    : mesh_(mesh), config_(config), compressible_(compressible),
      relaxation_(config.number("energy_relaxation", 0.7)),
      convectionBlending_(config.number("energy_convection_blending", config.number("convection_blending", 1.0))),
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
      diffusivity_(mesh.cells.size(), 0.0), source_(mesh.cells.size(), 0.0),
      destruction_(mesh.cells.size(), 0.0), constrainedCells_(mesh.cells.size(), 0),
      constrainedValues_(mesh.cells.size(), 0.0), energySystem_(static_cast<int>(mesh.cells.size())) {
    if (relaxation_ <= 0.0 || relaxation_ > 1.0) throw std::runtime_error("energy_relaxation must be in (0, 1]");
    if (convectionBlending_ < 0.0 || convectionBlending_ > 1.0) {
        throw std::runtime_error("energy_convection_blending must be between zero and one");
    }
    if (!(turbulentPrandtl_ > 0.0)) throw std::runtime_error("turbulent_prandtl must be positive");
    for (double value : temperature_.values()) {
        if (!(value > 0.0) || !std::isfinite(value)) throw std::runtime_error("initial_temperature must be positive");
    }
    for (const BoundaryPatch& patch : mesh.patches) boundaries_.push_back(readBoundary(patch.name));
    updateBoundaries();
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

void EnergyModel::initialize() {
    if (initialized_) return;
    previousTemperature_.values() = temperature_.values();
    oldTemperature_.values() = temperature_.values();
    oldOldTemperature_.values() = temperature_.values();
    updateBoundaries();
    initialized_ = true;
}

void EnergyModel::updateBoundaries() {
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const int owner = mesh_.faces[face].owner;
        const EnergyBoundaryCondition& condition = boundaries_[mesh_.faces[face].patch];
        temperature_.boundaryDirichlet(face) = condition.type == EnergyBoundaryType::FixedValue;
        temperature_.boundary(face) = condition.type == EnergyBoundaryType::FixedValue
            ? condition.value : temperature_[owner];
    }
}

void EnergyModel::beginIteration() {
    previousTemperature_.values() = temperature_.values();
}

LinearSolveInfo EnergyModel::solve(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                                   const GradientScheme& gradient, const std::vector<double>& density,
                                   const std::vector<double>& oldDensity, const std::vector<double>& oldOldDensity,
                                   const std::vector<double>& viscosity, const std::vector<double>& conductivity,
                                   const std::vector<double>& heatCapacity, const ScalarField* turbulentViscosity,
                                   const VectorField& velocity, const ScalarField& pressure,
                                   const ScalarField& oldPressure, const ScalarField& oldOldPressure,
                                   const Vector2& bodyForce) {
    initialize();
    updateBoundaries();
    const auto velocityGradient = gradient.vector(velocity);
    const auto pressureGradient = compressible_ ? gradient.scalar(pressure) : std::vector<ScalarGradient>{};
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOld);
    std::fill(source_.begin(), source_.end(), 0.0);
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        const double cp = heatCapacity[cell];
        if (!(cp > 0.0) || !(density[cell] > 0.0) || conductivity[cell] < 0.0) {
            throw std::runtime_error("Energy transport requires positive density and heat capacity, and non-negative conductivity");
        }
        const double eddy = turbulentViscosity ? (*turbulentViscosity)[cell] : 0.0;
        diffusivity_[cell] = conductivity[cell] / cp + eddy / turbulentPrandtl_;
        double energySource = heatSource_;
        if (compressible_) {
            double pressureRate = time.current * pressure[cell] - time.previous * oldPressure[cell]
                - time.previousPrevious * oldOldPressure[cell];
            pressureRate += dot(velocity[cell], pressureGradient[cell]);
            energySource += pressureRate + density[cell] * dot(bodyForce, velocity[cell]);
            if (viscousDissipation_) {
                const VectorGradient& derivatives = velocityGradient[cell];
                const double divergence = derivatives[0].x + derivatives[1].y;
                const double shear = 0.5 * (derivatives[0].y + derivatives[1].x);
                const double strainProduct = derivatives[0].x * derivatives[0].x
                    + derivatives[1].y * derivatives[1].y + 2.0 * shear * shear;
                energySource += viscosity[cell] * (2.0 * strainProduct - (2.0 / 3.0) * divergence * divergence);
            }
        }
        source_[cell] = energySource / cp;
    }
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const EnergyBoundaryCondition& condition = boundaries_[mesh_.faces[face].patch];
        if (condition.type != EnergyBoundaryType::FixedHeatFlux) continue;
        const int owner = mesh_.faces[face].owner;
        source_[owner] += condition.value * mesh_.faces[face].area /
            (heatCapacity[owner] * mesh_.cells[owner].volume);
    }

    ScalarTransportEquation equation{temperature_, previousTemperature_, oldTemperature_, oldOldTemperature_,
        massFlux, diffusivity_, source_, destruction_, constrainedCells_, constrainedValues_, 1.0, time,
        relaxation_, convectionBlending_, linearRelativeTolerance_, linearAbsoluteTolerance_,
        maximumLinearIterations_, config_.get("energy_ksp", "bcgs"), config_.get("energy_pc", "ilu")};
    equation.densityField = &density;
    equation.oldDensityField = &oldDensity;
    equation.oldOldDensityField = &oldOldDensity;
    const LinearSolveInfo info = solveScalarTransport(mesh_, gradient, energySystem_, equation);
    updateBoundaries();
    change_ = 0.0;
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        if (!(temperature_[cell] > 0.0) || !std::isfinite(temperature_[cell])) {
            throw std::runtime_error("Energy equation produced a non-positive or non-finite temperature");
        }
        change_ = std::max(change_, std::abs(temperature_[cell] - previousTemperature_[cell]));
    }
    linearIterations_ = info.iterations;
    linearResidual_ = info.residual;
    return info;
}

void EnergyModel::updateTimeHistory() {
    oldOldTemperature_.values() = oldTemperature_.values();
    oldTemperature_.values() = temperature_.values();
}

} // namespace mp
