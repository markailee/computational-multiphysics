#include "flow/SimpleSolver.hpp"

#include "fvm/Operators.hpp"
#include "fvm/TecplotWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace mp {
namespace {

FlowModel readFlowModel(const Config& config) {
    const std::string model = lowercase(config.get("flow_model", "incompressible"));
    if (model == "incompressible" || model == "constant_density") return FlowModel::Incompressible;
    if (model == "boussinesq" || model == "thermal" || model == "thermal_dilation") return FlowModel::Boussinesq;
    if (model == "compressible" || model == "real_density" || model == "ideal_gas") return FlowModel::Compressible;
    throw std::runtime_error("Unknown flow_model: " + model);
}

double initialFlowPressure(const Config& config) {
    if (config.has("initial_pressure")) return config.number("initial_pressure", 0.0);
    const FlowModel model = readFlowModel(config);
    return model == FlowModel::Compressible
        ? config.number("thermo_pressure", 101325.0) : 0.0;
}

} // namespace

SimpleSolver::SimpleSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory)
    : mesh_(mesh), config_(config), caseDirectory_(std::move(caseDirectory)),
      flowModel_(readFlowModel(config)), material_(makeMaterialModel(config)),
      density_(config.number("density", 1.0)), viscosity_(config.number("viscosity", 0.01)),
      velocityRelaxation_(config.number("velocity_relaxation", 0.7)),
      pressureRelaxation_(config.number("pressure_relaxation", 0.3)),
      barthJespersenLimiter_(useBarthJespersenLimiter(config)),
      outerTolerance_(config.number("tolerance", 1.0e-6)),
      continuityTolerance_(config.number("continuity_tolerance", 1.0e-8)),
      linearRelativeTolerance_(config.number("linear_tolerance", 1.0e-9)),
      linearAbsoluteTolerance_(config.number("linear_absolute_tolerance", 1.0e-12)),
      maximumLinearIterations_(config.integer("maximum_linear_iterations", 1000)),
      bodyForce_(config.vector("body_force", {})), gravity_(config.vector("gravity", {})),
      thermalExpansion_(config.number("thermal_expansion", 0.0)),
      referenceTemperature_(config.number("reference_temperature", config.number("initial_temperature", 300.0))),
      timeDerivative_(config.get("time_scheme", "steady")),
      boundaryConditions_(mesh, config),
      gradient_(makeGradientScheme(config.get("gradient", "least_squares"))),
      velocity_("U", mesh, config.vector("initial_velocity", {})),
      pressure_("p", mesh, initialFlowPressure(config)),
      pressureCorrection_("pCorrection", mesh, 0.0),
      previousVelocity_("previousU", mesh),
      oldVelocity_("oldU", mesh),
      oldOldVelocity_("oldOldU", mesh),
      oldPressure_("oldP", mesh, pressure_[0]), oldOldPressure_("oldOldP", mesh, pressure_[0]),
      massFlux_(mesh.faces.size(), 0.0), densityField_(mesh.cells.size(), density_),
      previousDensityField_(mesh.cells.size(), density_), oldDensityField_(mesh.cells.size(), density_),
      oldOldDensityField_(mesh.cells.size(), density_), viscosityField_(mesh.cells.size(), viscosity_),
      conductivityField_(mesh.cells.size(), 0.0), heatCapacityField_(mesh.cells.size(), 0.0),
      pressureDerivativeField_(mesh.cells.size(), 0.0), momentumDiagonal_(mesh.cells.size(), 1.0),
      lastImbalance_(mesh.cells.size(), 0.0),
      momentumSystem_(static_cast<int>(mesh.cells.size())), pressureSystem_(static_cast<int>(mesh.cells.size())) {
    if (density_ <= 0.0 || viscosity_ < 0.0) throw std::runtime_error("Density must be positive and viscosity non-negative");
    if (velocityRelaxation_ <= 0.0 || velocityRelaxation_ > 1.0 ||
        pressureRelaxation_ <= 0.0 || pressureRelaxation_ > 1.0) {
        throw std::runtime_error("SIMPLE relaxation factors must be in (0, 1]");
    }
    if (flowModel_ == FlowModel::Boussinesq) {
        if (material_->kind() != MaterialKind::Constant) {
            throw std::runtime_error("Boussinesq flow requires material_model = constant");
        }
        if (!(thermalExpansion_ > 0.0) || !(referenceTemperature_ > 0.0)) {
            throw std::runtime_error("Boussinesq flow requires positive thermal_expansion and reference_temperature");
        }
        energy_ = std::make_unique<EnergyModel>(mesh_, config_, false);
    }
    if (flowModel_ == FlowModel::Compressible) {
        if (material_->kind() != MaterialKind::IdealGas) {
            throw std::runtime_error("Compressible flow requires single-component material_model = ideal_gas");
        }
        if (lowercase(config.get("turbulence_model", "laminar")) != "laminar") {
            throw std::runtime_error("Compressible flow currently supports laminar calculations only");
        }
        if (config.has("species_model") && lowercase(config.get("species_model", "none")) != "none") {
            throw std::runtime_error("Compressible species coupling is not implemented");
        }
        if (timeDerivative_.isSteady() && !boundaryConditions_.hasFixedPressure()) {
            throw std::runtime_error("Steady compressible flow requires at least one fixedValue absolute-pressure boundary");
        }
        energy_ = std::make_unique<EnergyModel>(mesh_, config_, true);
        std::cout << "Compressible ideal-gas flow: pressure is absolute\n";
    }
    turbulence_ = makeTurbulenceModel(mesh_, config_, boundaryConditions_, density_, viscosity_);
    if (turbulence_) std::cout << "Turbulence model: " << turbulence_->name() << '\n';
    multiSpecies_ = makeMultiSpeciesModel(mesh_, config_, material_);
    multiphase_ = makeMultiphaseModel(mesh_, config_);
    if (material_->kind() == MaterialKind::IdealGasMixture && !multiSpecies_) {
        throw std::runtime_error("An ideal-gas mixture requires species_model = transport in the current flow solver");
    }
    if (material_->kind() == MaterialKind::VolumeFractionMixture && !multiphase_) {
        throw std::runtime_error("A volume-fraction material requires an active multiphase_model");
    }
    if (multiSpecies_) {
        std::cout << "Material model: " << material_->name() << "; passive species transport with "
                  << multiSpecies_->count() << " species\n";
    }
    if (multiphase_) std::cout << "Multiphase model: " << multiphase_->name() << " (interface only)\n";
    if (config.has("force_patch")) {
        forcePatch_ = mesh_.patchIndex(config.require("force_patch"));
        if (forcePatch_ < 0) throw std::runtime_error("Unknown force_patch: " + config.require("force_patch"));
        forceReferenceVelocity_ = config.number("force_reference_velocity", 1.0);
        forceReferenceLength_ = config.number("force_reference_length", 1.0);
        if (forceReferenceVelocity_ <= 0.0 || forceReferenceLength_ <= 0.0) {
            throw std::runtime_error("Force reference velocity and length must be positive");
        }
    }
}

void SimpleSolver::initialize() {
    if (initialized_) return;
    if (config_.has("restart")) readRestart(caseDirectory_ / config_.require("restart"));
    oldVelocity_.values() = velocity_.values();
    oldOldVelocity_.values() = velocity_.values();
    oldPressure_.values() = pressure_.values();
    oldOldPressure_.values() = pressure_.values();
    updateMaterialProperties();
    updateBoundaries();
    if (energy_) energy_->initialize(velocity_, heatCapacityField_);
    if (multiSpecies_) multiSpecies_->initialize();
    if (multiphase_) multiphase_->initialize(*material_);
    oldDensityField_ = densityField_;
    oldOldDensityField_ = densityField_;
    previousDensityField_ = densityField_;
    computeMassFlux(false);
    initialized_ = true;
    std::cout << "Mesh: " << mesh_.cells.size() << " cells, " << mesh_.faces.size() << " faces, "
              << mesh_.nodes.size() << " nodes\n";
}

void SimpleSolver::readRestart(const std::filesystem::path& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("Cannot open restart file: " + file.string());
    std::vector<std::string> names;
    std::string line;
    std::size_t cellOffset = 0;
    bool foundZone = false;
    bool foundTemperature = false;
    while (std::getline(input, line)) {
        if (line.rfind("VARIABLES", 0) == 0) {
            names.clear();
            std::size_t start = line.find('"');
            while (start != std::string::npos) {
                const std::size_t end = line.find('"', start + 1);
                if (end == std::string::npos) throw std::runtime_error("Invalid variable list in restart file");
                names.push_back(line.substr(start + 1, end - start - 1));
                start = line.find('"', end + 1);
            }
            continue;
        }
        if (line.rfind("ZONE", 0) != 0) continue;
        if (names.size() < 3) throw std::runtime_error("Restart zone has no cell-centered fields");
        const std::size_t nPosition = line.find("N=");
        const std::size_t ePosition = line.find("E=");
        if (nPosition == std::string::npos || ePosition == std::string::npos) {
            throw std::runtime_error("Cannot read restart zone dimensions");
        }
        const int nodeCount = std::stoi(line.substr(nPosition + 2));
        const int cellCount = std::stoi(line.substr(ePosition + 2));
        if (nodeCount < 0 || cellCount < 0 || cellOffset + static_cast<std::size_t>(cellCount) > mesh_.cells.size()) {
            throw std::runtime_error("Restart cell count does not match the mesh");
        }
        std::vector<double> values(nodeCount);
        for (int component = 0; component < 2; ++component) {
            for (int node = 0; node < nodeCount; ++node) input >> values[node];
        }
        for (std::size_t variable = 2; variable < names.size(); ++variable) {
            values.resize(cellCount);
            for (int cell = 0; cell < cellCount; ++cell) input >> values[cell];
            if (!input) throw std::runtime_error("Unexpected end of restart data");
            for (int cell = 0; cell < cellCount; ++cell) {
                const std::size_t target = cellOffset + static_cast<std::size_t>(cell);
                if (names[variable] == "U_x") velocity_[target].x = values[cell];
                else if (names[variable] == "U_y") velocity_[target].y = values[cell];
                else if (names[variable] == "p") pressure_[target] = values[cell];
                else if (names[variable] == "T" && energy_) {
                    energy_->temperature()[target] = values[cell];
                    foundTemperature = true;
                }
            }
        }
        cellOffset += static_cast<std::size_t>(cellCount);
        foundZone = true;
    }
    if (!foundZone || cellOffset != mesh_.cells.size()) {
        throw std::runtime_error("Restart cell count does not match the mesh");
    }
    if (energy_ && !foundTemperature) throw std::runtime_error("Restart file has no temperature field");
    std::cout << "Restart: " << file << '\n';
}

void SimpleSolver::updateBoundaries() {
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const int owner = mesh_.faces[face].owner;
        velocity_.boundary(face) = boundaryConditions_.velocity(face, velocity_[owner]);
        velocity_.boundaryDirichlet(face) = boundaryConditions_.fixedVelocity(face);
        pressure_.boundary(face) = boundaryConditions_.pressure(face, pressure_[owner]);
        pressure_.boundaryDirichlet(face) = boundaryConditions_.fixedPressure(face);
        if ((flowModel_ == FlowModel::Compressible) &&
            !(pressure_.boundary(face) > 0.0)) {
            throw std::runtime_error("Compressible boundary pressure must be absolute and positive");
        }
    }
    if (energy_) energy_->updateBoundaries(velocity_, heatCapacityField_);
    if (turbulence_) turbulence_->updateBoundaries();
}

void SimpleSolver::updateMaterialProperties() {
    MaterialState state;
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        if (material_->kind() == MaterialKind::IdealGasMixture ||
            material_->kind() == MaterialKind::VolumeFractionMixture) {
            densityField_[cell] = density_;
            viscosityField_[cell] = viscosity_;
            conductivityField_[cell] = 0.0;
            heatCapacityField_[cell] = 0.0;
            pressureDerivativeField_[cell] = 0.0;
            continue;
        }
        state.temperature = energy_ ? energy_->temperature()[cell] : config_.number("temperature", 300.0);
        state.pressure = flowModel_ == FlowModel::Compressible
            ? pressure_[cell] : config_.number("thermo_pressure", 101325.0);
        const MaterialProperties properties = material_->properties(state);
        if (!(properties.density > 0.0) || properties.viscosity < 0.0 ||
            properties.thermalConductivity < 0.0 || properties.specificHeatCapacity < 0.0) {
            throw std::runtime_error("Material returned invalid flow properties");
        }
        densityField_[cell] = flowModel_ == FlowModel::Compressible
            ? properties.density : density_;
        viscosityField_[cell] = properties.viscosity;
        conductivityField_[cell] = properties.thermalConductivity;
        heatCapacityField_[cell] = properties.specificHeatCapacity;
        pressureDerivativeField_[cell] = flowModel_ == FlowModel::Compressible
            ? densityField_[cell] / pressure_[cell] : 0.0;
    }
}

double SimpleSolver::faceDensity(std::size_t faceIndex) const {
    if (flowModel_ != FlowModel::Compressible) return density_;
    const Face& face = mesh_.faces[faceIndex];
    if (!face.isBoundary()) {
        return face.ownerWeight * densityField_[face.owner] + (1.0 - face.ownerWeight) * densityField_[face.neighbour];
    }
    MaterialState state;
    state.pressure = pressure_.boundary(faceIndex);
    state.temperature = energy_->temperature().boundary(faceIndex);
    return material_->properties(state).density;
}

void SimpleSolver::computeMassFlux(bool rhieChow) {
    std::vector<ScalarGradient> pressureGradient;
    if (rhieChow) pressureGradient = gradient_->scalar(pressure_);
    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const Vector2 interpolatedVelocity = interpolate(mesh_, velocity_, faceIndex);
        double normalVelocity = dot(interpolatedVelocity, face.normal);
        if (rhieChow) {
            const double ownerD = mesh_.cells[face.owner].volume / momentumDiagonal_[face.owner];
            const double neighbourD = mesh_.cells[face.neighbour].volume / momentumDiagonal_[face.neighbour];
            const double faceD = face.ownerWeight * ownerD + (1.0 - face.ownerWeight) * neighbourD;
            const Vector2 interpolatedGradient = face.ownerWeight * pressureGradient[face.owner] +
                (1.0 - face.ownerWeight) * pressureGradient[face.neighbour];
            const double normalPressureGradient = (pressure_[face.neighbour] - pressure_[face.owner]) / face.delta;
            normalVelocity -= faceD * (normalPressureGradient - dot(interpolatedGradient, face.normal));
        }
        massFlux_[faceIndex] = faceDensity(faceIndex) * face.area * normalVelocity;
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const Vector2 value = boundaryConditions_.velocity(faceIndex, velocity_[face.owner]);
        massFlux_[faceIndex] = faceDensity(faceIndex) * face.area * dot(value, face.normal);
    }
}

LinearSolveInfo SimpleSolver::solveMomentum(double dt) {
    const auto velocityGradient = gradient_->vector(velocity_);
    std::vector<VectorGradient> limitedVelocityGradient;
    const std::vector<VectorGradient>* convectionGradient = &velocityGradient;
    if (barthJespersenLimiter_) {
        limitedVelocityGradient = applyBarthJespersenLimiter(velocity_, velocityGradient);
        convectionGradient = &limitedVelocityGradient;
    }
    const auto pressureGradient = gradient_->scalar(pressure_);
    std::vector<ScalarGradient> kineticEnergyGradient;
    if (turbulence_) kineticEnergyGradient = gradient_->scalar(turbulence_->k());
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOldVelocity_);
    const int cellCount = static_cast<int>(mesh_.cells.size());
    SparseMatrix matrix(cellCount);
    std::vector<double> rhsX(cellCount, 0.0), rhsY(cellCount, 0.0);
    std::vector<Vector2> rhs(cellCount);

    auto turbulentStressFlux = [&](std::size_t faceIndex, const VectorGradient& faceGradient,
                                   double eddyViscosity) {
        Vector2 result;
        if (eddyViscosity <= 0.0) return result;
        const Vector2 area = faceVector(mesh_, faceIndex);
        const double divergence = faceGradient[0].x + faceGradient[1].y;
        for (int component = 0; component < 2; ++component) {
            const Vector2 transposeRow{faceGradient[0][component], faceGradient[1][component]};
            result[component] = eddyViscosity * (dot(transposeRow, area) -
                (2.0 / 3.0) * divergence * area[component]);
        }
        return result;
    };

    for (int cell = 0; cell < cellCount; ++cell) {
        const double volume = mesh_.cells[cell].volume;
        Vector2 acceleration = bodyForce_;
        if (flowModel_ == FlowModel::Compressible) acceleration += gravity_;
        rhs[cell] = acceleration * (densityField_[cell] * volume);
        if (flowModel_ == FlowModel::Boussinesq) {
            rhs[cell] -= gravity_ * (density_ * thermalExpansion_ * volume *
                (energy_->temperature()[cell] - referenceTemperature_));
        }
        matrix.add(cell, cell, densityField_[cell] * volume * time.current);
        rhs[cell] += oldVelocity_[cell] * (oldDensityField_[cell] * volume * time.previous) +
                     oldOldVelocity_[cell] * (oldOldDensityField_[cell] * volume * time.previousPrevious) -
                     pressureGradient[cell] * volume;
        if (turbulence_) {
            rhs[cell] -= kineticEnergyGradient[cell] * ((2.0 / 3.0) * densityField_[cell] * volume);
        }
    }

    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        const double eddyViscosity = turbulence_
            ? interpolate(mesh_, turbulence_->turbulentViscosity(), faceIndex) : 0.0;
        const double molecularViscosity = face.ownerWeight * viscosityField_[owner] +
            (1.0 - face.ownerWeight) * viscosityField_[neighbour];
        const double effectiveViscosity = molecularViscosity + eddyViscosity;
        const double diffusion = diffusionConductance(mesh_, faceIndex, effectiveViscosity);
        matrix.add(owner, owner, diffusion);
        matrix.add(owner, neighbour, -diffusion);
        matrix.add(neighbour, neighbour, diffusion);
        matrix.add(neighbour, owner, -diffusion);

        const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
        VectorGradient faceGradient;
        for (int component = 0; component < 2; ++component) {
            faceGradient[component] = face.ownerWeight * velocityGradient[owner][component] +
                (1.0 - face.ownerWeight) * velocityGradient[neighbour][component];
            const double correction = effectiveViscosity * dot(correctionVector, faceGradient[component]);
            rhs[owner][component] += correction;
            rhs[neighbour][component] -= correction;
        }
        const double secondaryViscosity = eddyViscosity +
            (flowModel_ == FlowModel::Compressible ? molecularViscosity : 0.0);
        const Vector2 secondaryFlux = turbulentStressFlux(faceIndex, faceGradient, secondaryViscosity);
        rhs[owner] += secondaryFlux;
        rhs[neighbour] -= secondaryFlux;

        const double flux = massFlux_[faceIndex];
        const ConvectionCoefficients convection = upwindConvection(flux);
        matrix.add(owner, owner, convection.ownerDiagonal);
        matrix.add(owner, neighbour, convection.ownerToNeighbour);
        matrix.add(neighbour, neighbour, convection.neighbourDiagonal);
        matrix.add(neighbour, owner, convection.neighbourToOwner);
        const Vector2 upwind = flux >= 0.0 ? velocity_[owner] : velocity_[neighbour];
        const Vector2 reconstructed = secondOrderUpwindFaceValue(mesh_, velocity_, *convectionGradient, faceIndex, flux);
        const Vector2 deferred = flux * (reconstructed - upwind);
        rhs[owner] -= deferred;
        rhs[neighbour] += deferred;
    }

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const Vector2 boundaryVelocity = velocity_.boundary(faceIndex);
        if (boundaryConditions_.fixedVelocity(faceIndex)) {
            if (turbulence_ && turbulence_->usesStandardWallFunction() && boundaryConditions_.noSlip(faceIndex)) {
                rhs[owner] -= turbulence_->wallShearStress(faceIndex, velocity_[owner]) * face.area;
            } else {
                const double eddyViscosity = turbulence_ ? turbulence_->faceEddyViscosity(faceIndex) : 0.0;
                const double effectiveViscosity = viscosityField_[owner] + eddyViscosity;
                const double diffusion = diffusionConductance(mesh_, faceIndex, effectiveViscosity);
                matrix.add(owner, owner, diffusion);
                const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
                for (int component = 0; component < 2; ++component) {
                    rhs[owner][component] += diffusion * boundaryVelocity[component] +
                        effectiveViscosity * dot(correctionVector, velocityGradient[owner][component]);
                }
                rhs[owner] += turbulentStressFlux(faceIndex, velocityGradient[owner], eddyViscosity +
                    (flowModel_ == FlowModel::Compressible ? viscosityField_[owner] : 0.0));
            }
        }

        const double flux = massFlux_[faceIndex];
        const Vector2 upwind = flux >= 0.0 ? velocity_[owner] : velocity_.boundary(faceIndex);
        const Vector2 reconstructed = secondOrderUpwindFaceValue(mesh_, velocity_, *convectionGradient, faceIndex, flux);
        if (flux >= 0.0) matrix.add(owner, owner, flux);
        else rhs[owner] -= flux * upwind;
        rhs[owner] -= flux * (reconstructed - upwind);
    }

    for (int cell = 0; cell < cellCount; ++cell) {
        const double diagonal = matrix.diagonal(cell);
        if (!(diagonal > 0.0)) throw std::runtime_error("Non-positive momentum diagonal; check mesh and boundary conditions");
        const double relaxation = diagonal * (1.0 / velocityRelaxation_ - 1.0);
        matrix.add(cell, cell, relaxation);
        rhs[cell] += relaxation * previousVelocity_[cell];
        momentumDiagonal_[cell] = matrix.diagonal(cell);
        rhsX[cell] = rhs[cell].x;
        rhsY[cell] = rhs[cell].y;
    }

    std::vector<double> solutionX(cellCount), solutionY(cellCount);
    for (int cell = 0; cell < cellCount; ++cell) {
        solutionX[cell] = velocity_[cell].x;
        solutionY[cell] = velocity_[cell].y;
    }
    momentumSystem_.setMatrix(matrix);
    const std::string solver = config_.get("momentum_ksp", "bcgs");
    const std::string preconditioner = config_.get("momentum_pc", "ilu");
    const LinearSolveInfo infoX = momentumSystem_.solve(rhsX, solutionX, solver, preconditioner,
        linearRelativeTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
    const LinearSolveInfo infoY = momentumSystem_.solve(rhsY, solutionY, solver, preconditioner,
        linearRelativeTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
    for (int cell = 0; cell < cellCount; ++cell) velocity_[cell] = {solutionX[cell], solutionY[cell]};
    updateBoundaries();
    computeMassFlux(true);
    return {infoX.iterations + infoY.iterations, std::max(infoX.residual, infoY.residual)};
}

LinearSolveInfo SimpleSolver::correctPressure(double dt) {
    const int cellCount = static_cast<int>(mesh_.cells.size());
    const bool compressible = flowModel_ == FlowModel::Compressible;
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOldVelocity_);
    SparseMatrix matrix(cellCount);
    std::vector<double> rhs(cellCount, 0.0);
    if (compressible) {
        for (int cell = 0; cell < cellCount; ++cell) {
            const double volume = mesh_.cells[cell].volume;
            matrix.add(cell, cell, pressureDerivativeField_[cell] * volume * time.current);
            rhs[cell] -= volume * (time.current * densityField_[cell]
                - time.previous * oldDensityField_[cell]
                - time.previousPrevious * oldOldDensityField_[cell]);
        }
    }
    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double density = faceDensity(faceIndex);
        const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density);
        matrix.add(face.owner, face.owner, coefficient);
        matrix.add(face.owner, face.neighbour, -coefficient);
        matrix.add(face.neighbour, face.neighbour, coefficient);
        matrix.add(face.neighbour, face.owner, -coefficient);
        rhs[face.owner] -= massFlux_[faceIndex];
        rhs[face.neighbour] += massFlux_[faceIndex];
        if (compressible) {
            const double volumetricFlux = massFlux_[faceIndex] / density;
            const double upwindDerivative = volumetricFlux >= 0.0
                ? pressureDerivativeField_[face.owner] : pressureDerivativeField_[face.neighbour];
            const ConvectionCoefficients convection = upwindConvection(volumetricFlux * upwindDerivative);
            matrix.add(face.owner, face.owner, convection.ownerDiagonal);
            matrix.add(face.owner, face.neighbour, convection.ownerToNeighbour);
            matrix.add(face.neighbour, face.neighbour, convection.neighbourDiagonal);
            matrix.add(face.neighbour, face.owner, convection.neighbourToOwner);
        }
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double density = faceDensity(faceIndex);
        if (boundaryConditions_.fixedPressure(faceIndex) && !boundaryConditions_.fixedVelocity(faceIndex)) {
            const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density);
            matrix.add(face.owner, face.owner, coefficient);
        } else if (compressible && !boundaryConditions_.fixedPressure(faceIndex) &&
                   !boundaryConditions_.fixedVelocity(faceIndex) && massFlux_[faceIndex] > 0.0) {
            const double coefficient = massFlux_[faceIndex] / density * pressureDerivativeField_[face.owner];
            matrix.add(face.owner, face.owner, coefficient);
        }
        rhs[face.owner] -= massFlux_[faceIndex];
    }
    if (!compressible && !boundaryConditions_.hasFixedPressure()) {
        matrix.pinZero(0);
        rhs[0] = 0.0;
    }
    pressureSystem_.setMatrix(matrix);
    const std::string defaultSolver = compressible ? "bcgs" : "cg";
    const std::string defaultPreconditioner = compressible ? "ilu" : "icc";
    const LinearSolveInfo info = pressureSystem_.solve(rhs, pressureCorrection_.values(),
        config_.get("pressure_ksp", defaultSolver), config_.get("pressure_pc", defaultPreconditioner),
        linearRelativeTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const int owner = mesh_.faces[faceIndex].owner;
        // A prescribed boundary velocity requires u_b' = 0, i.e. homogeneous Neumann for p'
        // (zero normal pressure-correction gradient), even when the static pressure is also
        // prescribed. Only a fixed-pressure boundary with an extrapolated velocity uses p_b' = 0.
        const bool correctionDirichlet = boundaryConditions_.fixedPressure(faceIndex) &&
            !boundaryConditions_.fixedVelocity(faceIndex);
        pressureCorrection_.boundary(faceIndex) = correctionDirichlet ? 0.0 : pressureCorrection_[owner];
        pressureCorrection_.boundaryDirichlet(faceIndex) = correctionDirichlet;
    }
    const auto correctionGradient = gradient_->scalar(pressureCorrection_);
    for (int cell = 0; cell < cellCount; ++cell) {
        velocity_[cell] -= correctionGradient[cell] * (mesh_.cells[cell].volume / momentumDiagonal_[cell]);
        pressure_[cell] += pressureRelaxation_ * pressureCorrection_[cell];
        if (compressible && !(pressure_[cell] > 0.0)) {
            throw std::runtime_error("Pressure correction produced a non-positive absolute pressure; reduce pressure_relaxation");
        }
    }
    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double density = faceDensity(faceIndex);
        const double predictedMassFlux = massFlux_[faceIndex];
        const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density);
        massFlux_[faceIndex] += coefficient * (pressureCorrection_[face.owner] - pressureCorrection_[face.neighbour]);
        if (compressible) {
            const double volumetricFlux = predictedMassFlux / density;
            massFlux_[faceIndex] += volumetricFlux * (volumetricFlux >= 0.0
                ? pressureCorrection_[face.owner] * pressureDerivativeField_[face.owner]
                : pressureCorrection_[face.neighbour] * pressureDerivativeField_[face.neighbour]);
        }
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const int owner = mesh_.faces[faceIndex].owner;
        const double density = faceDensity(faceIndex);
        const double predictedMassFlux = massFlux_[faceIndex];
        if (boundaryConditions_.fixedPressure(faceIndex) && !boundaryConditions_.fixedVelocity(faceIndex)) {
            const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density);
            massFlux_[faceIndex] += coefficient * pressureCorrection_[owner];
        } else if (compressible && !boundaryConditions_.fixedPressure(faceIndex) &&
                   !boundaryConditions_.fixedVelocity(faceIndex) && predictedMassFlux > 0.0) {
            massFlux_[faceIndex] += predictedMassFlux / density * pressureDerivativeField_[owner] *
                pressureCorrection_[owner];
        }
    }
    updateBoundaries();
    updateMaterialProperties();
    return info;
}

std::vector<double> SimpleSolver::cellImbalance() const {
    std::vector<double> imbalance(mesh_.cells.size(), 0.0);
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        const Cell& controlVolume = mesh_.cells[cell];
        for (std::size_t i = 0; i < controlVolume.faces.size(); ++i) {
            imbalance[cell] += controlVolume.faceSigns[i] * massFlux_[controlVolume.faces[i]];
        }
    }
    return imbalance;
}

double SimpleSolver::netBoundaryFlux() const {
    double result = 0.0;
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) result += massFlux_[face];
    return result;
}

FlowDiagnostics SimpleSolver::solveStep(double dt, int maximumIterations) {
    FlowDiagnostics diagnostics;
    const int logInterval = config_.integer("log_interval", 10);
    const double turbulenceTolerance = config_.number("turbulence_tolerance", outerTolerance_);
    const double speciesTolerance = config_.number("species_tolerance", outerTolerance_);
    const double compositionTolerance = config_.number("species_composition_tolerance", 1.0e-8);
    const double energyTolerance = config_.number("energy_tolerance", outerTolerance_);
    const double densityTolerance = config_.number("density_tolerance", outerTolerance_);
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOldVelocity_);
    for (int iteration = 1; iteration <= maximumIterations; ++iteration) {
        previousVelocity_.values() = velocity_.values();
        previousDensityField_ = densityField_;
        if (energy_) energy_->beginIteration();
        if (turbulence_) turbulence_->beginIteration();
        if (multiSpecies_) multiSpecies_->beginIteration();
        if (multiphase_) multiphase_->beginOuterIteration();
        const LinearSolveInfo momentum = solveMomentum(dt);
        const LinearSolveInfo pressure = correctPressure(dt);
        if (multiphase_) {
            MultiphaseContext context{mesh_, velocity_, pressure_, massFlux_, densityField_, viscosityField_,
                *gradient_, *material_, multiSpecies_.get(), dt, currentTime_};
            multiphase_->solve(context);
        }
        LinearSolveInfo energyInfo;
        if (energy_) {
            const Vector2 acceleration = bodyForce_ + (flowModel_ == FlowModel::Compressible ? gravity_ : Vector2{});
            energyInfo = energy_->solve(dt, hasOldOldVelocity_, massFlux_, *gradient_, densityField_,
                oldDensityField_, oldOldDensityField_, viscosityField_, conductivityField_, heatCapacityField_,
                turbulence_ ? &turbulence_->turbulentViscosity() : nullptr, velocity_, pressure_, oldPressure_,
                oldOldPressure_, acceleration);
            updateMaterialProperties();
        }
        if (turbulence_) {
            turbulence_->solve(dt, hasOldOldVelocity_, velocity_, massFlux_, *gradient_);
        }
        if (multiSpecies_) {
            multiSpecies_->solve(dt, hasOldOldVelocity_, massFlux_, *gradient_, density_,
                turbulence_ ? &turbulence_->turbulentViscosity() : nullptr);
        }
        std::vector<double> imbalance = cellImbalance();
        if (flowModel_ == FlowModel::Compressible) {
            for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
                imbalance[cell] += mesh_.cells[cell].volume * (time.current * densityField_[cell]
                    - time.previous * oldDensityField_[cell] - time.previousPrevious * oldOldDensityField_[cell]);
            }
        }
        lastImbalance_ = imbalance;
        diagnostics.iterations = iteration;
        diagnostics.momentumLinearIterations = momentum.iterations;
        diagnostics.pressureLinearIterations = pressure.iterations;
        diagnostics.turbulenceChange = turbulence_ ? turbulence_->change() : 0.0;
        diagnostics.turbulenceLinearIterations = turbulence_ ? turbulence_->linearIterations() : 0;
        diagnostics.speciesChange = multiSpecies_ ? multiSpecies_->change() : 0.0;
        diagnostics.speciesCompositionError = multiSpecies_ ? multiSpecies_->compositionError() : 0.0;
        diagnostics.minimumMassFraction = multiSpecies_ ? multiSpecies_->minimumMassFraction() : 0.0;
        diagnostics.maximumMassFraction = multiSpecies_ ? multiSpecies_->maximumMassFraction() : 1.0;
        diagnostics.speciesLinearIterations = multiSpecies_ ? multiSpecies_->linearIterations() : 0;
        diagnostics.temperatureChange = energy_ ? energy_->change() : 0.0;
        diagnostics.energyLinearIterations = energy_ ? energy_->linearIterations() : 0;
        diagnostics.velocityChange = 0.0;
        diagnostics.pressureChange = 0.0;
        diagnostics.densityChange = 0.0;
        diagnostics.maxCellImbalance = 0.0;
        for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
            diagnostics.velocityChange = std::max(diagnostics.velocityChange, norm(velocity_[cell] - previousVelocity_[cell]));
            diagnostics.pressureChange = std::max(diagnostics.pressureChange,
                std::abs(pressureRelaxation_ * pressureCorrection_[cell]));
            diagnostics.densityChange = std::max(diagnostics.densityChange,
                std::abs(densityField_[cell] - previousDensityField_[cell]));
            diagnostics.maxCellImbalance = std::max(diagnostics.maxCellImbalance, std::abs(imbalance[cell]));
        }
        diagnostics.netBoundaryFlux = netBoundaryFlux();
        diagnostics.converged = diagnostics.velocityChange < outerTolerance_ &&
            diagnostics.pressureChange < outerTolerance_ && diagnostics.maxCellImbalance < continuityTolerance_ &&
            (!turbulence_ || diagnostics.turbulenceChange < turbulenceTolerance) &&
            (!multiSpecies_ || (diagnostics.speciesChange < speciesTolerance &&
                diagnostics.speciesCompositionError < compositionTolerance &&
                multiSpecies_->bounded(compositionTolerance))) &&
            (!energy_ || (diagnostics.temperatureChange < energyTolerance && !energy_->temperatureLimited())) &&
            (flowModel_ != FlowModel::Compressible || diagnostics.densityChange < densityTolerance);
        if (iteration == 1 || iteration % logInterval == 0 || diagnostics.converged || iteration == maximumIterations) {
            std::cout << (timeDerivative_.isSteady() ? "SIMPLE " : "  outer ")
                      << std::setw(6) << iteration << "  dU=" << std::scientific << diagnostics.velocityChange
                      << "  dp=" << diagnostics.pressureChange << "  dTurb=" << diagnostics.turbulenceChange
                      << "  continuity=" << diagnostics.maxCellImbalance;
            if (energy_) {
                std::cout << "  dTemp=" << diagnostics.temperatureChange << "  drho=" << diagnostics.densityChange
                          << "  Tlimit=" << energy_->temperatureLimited();
            }
            if (multiSpecies_) {
                std::cout << "  dY=" << diagnostics.speciesChange << "  Ysum=" << diagnostics.speciesCompositionError;
            }
            std::cout << "  linear=" << diagnostics.momentumLinearIterations + diagnostics.pressureLinearIterations +
                         diagnostics.turbulenceLinearIterations + diagnostics.speciesLinearIterations +
                         diagnostics.energyLinearIterations << '\n';
        }
        if (diagnostics.converged) break;
    }
    diagnostics_ = diagnostics;
    return diagnostics;
}

void SimpleSolver::advance(double dt) {
    initialize();
    const FlowDiagnostics result = solveStep(dt, config_.integer("maximum_iterations", 100));
    currentTime_ += dt;
    updateTimeHistory();
    if (!result.converged) throw std::runtime_error("SIMPLE did not converge at time " + std::to_string(currentTime_));
}

void SimpleSolver::updateTimeHistory() {
    oldOldVelocity_.values() = oldVelocity_.values();
    oldVelocity_.values() = velocity_.values();
    oldOldPressure_.values() = oldPressure_.values();
    oldPressure_.values() = pressure_.values();
    oldOldDensityField_ = oldDensityField_;
    oldDensityField_ = densityField_;
    if (energy_) energy_->updateTimeHistory();
    if (turbulence_) turbulence_->updateTimeHistory(hasOldOldVelocity_);
    if (multiSpecies_) multiSpecies_->updateTimeHistory();
    if (multiphase_) multiphase_->updateTimeHistory();
    hasOldOldVelocity_ = true;
}

FlowForceDiagnostics SimpleSolver::forceDiagnostics() const {
    FlowForceDiagnostics result;
    if (forcePatch_ < 0) return result;
    const std::vector<VectorGradient> velocityGradient = gradient_->vector(velocity_);
    for (int faceIndex : mesh_.patches[forcePatch_].faces) {
        const Face& face = mesh_.faces[faceIndex];
        const auto& gradient = velocityGradient[face.owner];
        const double divergence = gradient[0].x + gradient[1].y;
        Vector2 viscousTraction;
        for (int component = 0; component < 2; ++component) {
            const double transpose = gradient[0][component] * face.normal.x + gradient[1][component] * face.normal.y;
            const double dilatation = flowModel_ == FlowModel::Compressible
                ? (2.0 / 3.0) * divergence * face.normal[component] : 0.0;
            viscousTraction[component] = viscosityField_[face.owner] * face.area *
                (dot(gradient[component], face.normal) + transpose - dilatation);
        }
        result.pressure += face.normal * (pressure_.boundary(faceIndex) * face.area);
        result.viscous -= viscousTraction;
    }
    result.total = result.pressure + result.viscous;
    return result;
}

std::filesystem::path SimpleSolver::forceFile() const {
    std::filesystem::path output = config_.get("force_output", "results/forces.csv");
    return output.is_absolute() ? output : caseDirectory_ / output;
}

void SimpleSolver::writeForceHistory(const std::filesystem::path& file, double time, bool append) const {
    const FlowForceDiagnostics force = forceDiagnostics();
    const double dynamicPressure = 0.5 * density_ * forceReferenceVelocity_ * forceReferenceVelocity_ * forceReferenceLength_;
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file, append ? std::ios::app : std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write flow force history: " + file.string());
    output << std::setprecision(12) << std::scientific;
    if (!append) {
        output << "time,force_x,force_y,pressure_force_x,pressure_force_y,viscous_force_x,viscous_force_y,cd,cl\n";
    }
    output << time << ',' << force.total.x << ',' << force.total.y << ','
           << force.pressure.x << ',' << force.pressure.y << ',' << force.viscous.x << ',' << force.viscous.y << ','
           << force.total.x / dynamicPressure << ',' << force.total.y / dynamicPressure << '\n';
}

std::filesystem::path SimpleSolver::outputFile() const {
    std::filesystem::path output = config_.get("output", "results/flow.plt");
    return output.is_absolute() ? output : caseDirectory_ / output;
}

void SimpleSolver::write(const std::filesystem::path& file) const {
    ScalarField imbalance("mass_imbalance", mesh_);
    imbalance.values() = lastImbalance_;
    std::map<std::string, const ScalarField*> scalars{{"p", &pressure_}, {"mass_imbalance", &imbalance}};
    ScalarField density("rho", mesh_), conductivity("lambda", mesh_), heatCapacity("cp", mesh_), mach("Mach", mesh_);
    if (flowModel_ != FlowModel::Incompressible) {
        density.values() = densityField_;
        conductivity.values() = conductivityField_;
        heatCapacity.values() = heatCapacityField_;
        scalars["rho"] = &density;
        scalars["lambda"] = &conductivity;
        scalars["cp"] = &heatCapacity;
    }
    if (flowModel_ == FlowModel::Compressible) {
        const double gasConstant = config_.number("gas_constant", 287.0);
        for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
            const double gamma = heatCapacityField_[cell] / (heatCapacityField_[cell] - gasConstant);
            mach[cell] = norm(velocity_[cell]) / std::sqrt(gamma * pressure_[cell] / densityField_[cell]);
        }
        scalars["Mach"] = &mach;
    }
    if (energy_) scalars["T"] = &energy_->temperature();
    if (flowModel_ == FlowModel::Compressible) scalars["H0"] = &energy_->totalEnthalpy();
    ScalarField vorticity("vorticity", mesh_);
    if (config_.boolean("output_vorticity", false)) {
        const std::vector<VectorGradient> gradient = gradient_->vector(velocity_);
        for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
            vorticity[cell] = gradient[cell][1].x - gradient[cell][0].y;
        }
        scalars["vorticity"] = &vorticity;
    }
    if (turbulence_) {
        scalars["k"] = &turbulence_->k();
        scalars[turbulence_->scale().name()] = &turbulence_->scale();
        scalars["mu_t"] = &turbulence_->turbulentViscosity();
        if (turbulence_->kind() == TurbulenceModelKind::SSTMenter2003) {
            scalars["F1"] = &turbulence_->blendingF1();
            scalars["F2"] = &turbulence_->blendingF2();
        }
    }
    std::map<std::string, const VectorField*> vectors{{"U", &velocity_}};
    if (multiSpecies_) multiSpecies_->appendOutputFields(scalars);
    if (multiphase_) multiphase_->appendOutputFields(scalars, vectors);
    writeTecplot(file, mesh_, scalars, vectors, name(), currentTime_);
    if (turbulence_) {
        writeWallDiagnostics(file.parent_path() / (file.stem().string() + "_wall.csv"));
    }
}

void SimpleSolver::writeSummary(const std::filesystem::path& file, double time,
                                const FlowDiagnostics& diagnostics, bool append) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file, append ? std::ios::app : std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write flow summary: " + file.string());
    output << std::setprecision(12) << std::scientific;
    if (!append) output << "time,iterations,velocity_change,pressure_change,max_cell_imbalance,net_boundary_flux,turbulence_change,turbulence_linear_iterations,species_change,species_composition_error,minimum_mass_fraction,maximum_mass_fraction,species_linear_iterations,temperature_change,density_change,energy_linear_iterations,converged\n";
    output << time << ',' << diagnostics.iterations << ',' << diagnostics.velocityChange << ','
           << diagnostics.pressureChange << ',' << diagnostics.maxCellImbalance << ','
           << diagnostics.netBoundaryFlux << ',' << diagnostics.turbulenceChange << ','
           << diagnostics.turbulenceLinearIterations << ',' << diagnostics.speciesChange << ','
           << diagnostics.speciesCompositionError << ',' << diagnostics.minimumMassFraction << ','
           << diagnostics.maximumMassFraction << ',' << diagnostics.speciesLinearIterations << ','
           << diagnostics.temperatureChange << ',' << diagnostics.densityChange << ','
           << diagnostics.energyLinearIterations << ',' << (diagnostics.converged ? 1 : 0) << '\n';
}

void SimpleSolver::writeWallDiagnostics(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot write wall diagnostics: " + file.string());
    output << std::setprecision(12) << std::scientific;
    output << "patch,x,y,distance,tau_x,tau_y,friction_velocity,y_plus\n";
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        if (!boundaryConditions_.noSlip(face)) continue;
        const TurbulenceWallDiagnostics diagnostics = turbulence_->wallDiagnostics(face, velocity_[mesh_.faces[face].owner]);
        const Face& wall = mesh_.faces[face];
        output << mesh_.patches[wall.patch].name << ',' << wall.centre.x << ',' << wall.centre.y << ','
               << diagnostics.distance << ',' << diagnostics.shear.x << ',' << diagnostics.shear.y << ','
               << diagnostics.frictionVelocity << ',' << diagnostics.yPlus << '\n';
    }
}

void SimpleSolver::run() {
    initialize();
    const std::filesystem::path output = outputFile();
    const std::filesystem::path summary = output.parent_path() / (output.stem().string() + "_summary.csv");
    if (timeDerivative_.isSteady()) {
        const FlowDiagnostics result = solveStep(0.0, config_.integer("maximum_iterations", 2000));
        write(output);
        writeSummary(summary, 0.0, result, false);
        if (forcePatch_ >= 0) writeForceHistory(forceFile(), 0.0, false);
        std::cout << "Tecplot output: " << output << '\n';
        if (!result.converged) throw std::runtime_error("Steady SIMPLE iteration did not converge; output was written for diagnosis");
        return;
    }

    const double dt = config_.number("dt", 0.01);
    const double endTime = config_.number("end_time", 1.0);
    const int maximumSteps = config_.integer("maximum_steps", 1000000);
    const int outputInterval = config_.integer("output_interval", 10);
    const bool allowIncomplete = config_.boolean("allow_incomplete_time_steps", false);
    if (dt <= 0.0 || endTime <= 0.0) throw std::runtime_error("dt and end_time must be positive");
    bool wroteSummary = false;
    if (forcePatch_ >= 0) writeForceHistory(forceFile(), currentTime_, false);
    int step = 0;
    while (currentTime_ < endTime - 1.0e-12 && step < maximumSteps) {
        if (endTime - currentTime_ <= 1.0e-8 * dt) {
            currentTime_ = endTime;
            break;
        }
        const double stepSize = std::min(dt, endTime - currentTime_);
        const bool finalStep = currentTime_ + stepSize >= endTime - 1.0e-12;
        const int iterationLimit = finalStep
            ? config_.integer("final_maximum_iterations", config_.integer("maximum_iterations", 100))
            : config_.integer("maximum_iterations", 100);
        const FlowDiagnostics result = solveStep(stepSize, iterationLimit);
        currentTime_ += stepSize;
        ++step;
        updateTimeHistory();
        writeSummary(summary, currentTime_, result, wroteSummary);
        wroteSummary = true;
        if (forcePatch_ >= 0) writeForceHistory(forceFile(), currentTime_, true);
        if (step % outputInterval == 0 || currentTime_ >= endTime - 1.0e-12) {
            const std::filesystem::path stepFile = output.parent_path() /
                (output.stem().string() + "_" + std::to_string(step) + output.extension().string());
            write(stepFile);
        }
        if (!result.converged && (!allowIncomplete || finalStep)) {
            throw std::runtime_error("SIMPLE did not converge at time " + std::to_string(currentTime_));
        }
    }
    if (currentTime_ < endTime - 1.0e-12) throw std::runtime_error("maximum_steps reached before end_time");
    std::cout << "Tecplot output: " << output.parent_path() << '\n';
}

} // namespace mp
