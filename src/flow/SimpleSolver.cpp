#include "flow/SimpleSolver.hpp"

#include "fvm/Operators.hpp"
#include "fvm/TecplotWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace mp {

SimpleSolver::SimpleSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory)
    : mesh_(mesh), config_(config), caseDirectory_(std::move(caseDirectory)),
      density_(config.number("density", 1.0)), viscosity_(config.number("viscosity", 0.01)),
      velocityRelaxation_(config.number("velocity_relaxation", 0.7)),
      pressureRelaxation_(config.number("pressure_relaxation", 0.3)),
      convectionBlending_(config.number("convection_blending", 1.0)),
      outerTolerance_(config.number("tolerance", 1.0e-6)),
      continuityTolerance_(config.number("continuity_tolerance", 1.0e-8)),
      linearRelativeTolerance_(config.number("linear_tolerance", 1.0e-9)),
      linearAbsoluteTolerance_(config.number("linear_absolute_tolerance", 1.0e-12)),
      maximumLinearIterations_(config.integer("maximum_linear_iterations", 1000)),
      bodyForce_(config.vector("body_force", {})),
      timeDerivative_(config.get("time_scheme", "steady")),
      boundaryConditions_(mesh, config),
      gradient_(makeGradientScheme(config.get("gradient", "least_squares"))),
      velocity_("U", mesh, config.vector("initial_velocity", {})),
      pressure_("p", mesh, config.number("initial_pressure", 0.0)),
      pressureCorrection_("pCorrection", mesh, 0.0),
      previousVelocity_("previousU", mesh),
      oldVelocity_("oldU", mesh),
      oldOldVelocity_("oldOldU", mesh),
      massFlux_(mesh.faces.size(), 0.0), momentumDiagonal_(mesh.cells.size(), 1.0),
      momentumSystem_(static_cast<int>(mesh.cells.size())), pressureSystem_(static_cast<int>(mesh.cells.size())) {
    if (density_ <= 0.0 || viscosity_ < 0.0) throw std::runtime_error("Density must be positive and viscosity non-negative");
    if (velocityRelaxation_ <= 0.0 || velocityRelaxation_ > 1.0 ||
        pressureRelaxation_ <= 0.0 || pressureRelaxation_ > 1.0) {
        throw std::runtime_error("SIMPLE relaxation factors must be in (0, 1]");
    }
    if (convectionBlending_ < 0.0 || convectionBlending_ > 1.0) {
        throw std::runtime_error("convection_blending must be between zero (upwind) and one (central)");
    }
    turbulence_ = makeTurbulenceModel(mesh_, config_, boundaryConditions_, density_, viscosity_);
    if (turbulence_) std::cout << "Turbulence model: " << turbulence_->name() << '\n';
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
    oldVelocity_.values() = velocity_.values();
    oldOldVelocity_.values() = velocity_.values();
    updateBoundaries();
    computeMassFlux(false);
    initialized_ = true;
    std::cout << "Mesh: " << mesh_.cells.size() << " cells, " << mesh_.faces.size() << " faces, "
              << mesh_.nodes.size() << " nodes\n";
}

void SimpleSolver::updateBoundaries() {
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const int owner = mesh_.faces[face].owner;
        velocity_.boundary(face) = boundaryConditions_.velocity(face, velocity_[owner]);
        velocity_.boundaryDirichlet(face) = boundaryConditions_.fixedVelocity(face);
        pressure_.boundary(face) = boundaryConditions_.pressure(face, pressure_[owner]);
        pressure_.boundaryDirichlet(face) = boundaryConditions_.fixedPressure(face);
    }
    if (turbulence_) turbulence_->updateBoundaries();
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
        massFlux_[faceIndex] = density_ * face.area * normalVelocity;
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const Vector2 value = boundaryConditions_.velocity(faceIndex, velocity_[face.owner]);
        massFlux_[faceIndex] = density_ * face.area * dot(value, face.normal);
    }
}

LinearSolveInfo SimpleSolver::solveMomentum(double dt) {
    const auto velocityGradient = gradient_->vector(velocity_);
    const auto pressureGradient = gradient_->scalar(pressure_);
    std::vector<ScalarGradient> kineticEnergyGradient;
    if (turbulence_) kineticEnergyGradient = gradient_->scalar(turbulence_->k());
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOldVelocity_);
    const int cellCount = static_cast<int>(mesh_.cells.size());
    SparseMatrix matrix(cellCount);
    std::vector<double> rhsX(cellCount, 0.0), rhsY(cellCount, 0.0);
    std::vector<Vector2> rhs(cellCount, bodyForce_ * density_);

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
        matrix.add(cell, cell, density_ * volume * time.current);
        rhs[cell] += oldVelocity_[cell] * (density_ * volume * time.previous) +
                     oldOldVelocity_[cell] * (density_ * volume * time.previousPrevious) -
                     pressureGradient[cell] * volume;
        if (turbulence_) {
            rhs[cell] -= kineticEnergyGradient[cell] * ((2.0 / 3.0) * density_ * volume);
        }
    }

    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        const double eddyViscosity = turbulence_
            ? interpolate(mesh_, turbulence_->turbulentViscosity(), faceIndex) : 0.0;
        const double effectiveViscosity = viscosity_ + eddyViscosity;
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
        const Vector2 turbulentFlux = turbulentStressFlux(faceIndex, faceGradient, eddyViscosity);
        rhs[owner] += turbulentFlux;
        rhs[neighbour] -= turbulentFlux;

        const double flux = massFlux_[faceIndex];
        const ConvectionCoefficients convection = upwindConvection(flux);
        matrix.add(owner, owner, convection.ownerDiagonal);
        matrix.add(owner, neighbour, convection.ownerToNeighbour);
        matrix.add(neighbour, neighbour, convection.neighbourDiagonal);
        matrix.add(neighbour, owner, convection.neighbourToOwner);
        const Vector2 upwind = flux >= 0.0 ? velocity_[owner] : velocity_[neighbour];
        const Vector2 central = centralFaceValue(mesh_, velocity_, velocityGradient, faceIndex);
        const Vector2 deferred = (convectionBlending_ * flux) * (central - upwind);
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
                const double effectiveViscosity = viscosity_ + eddyViscosity;
                const double diffusion = diffusionConductance(mesh_, faceIndex, effectiveViscosity);
                matrix.add(owner, owner, diffusion);
                const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
                for (int component = 0; component < 2; ++component) {
                    rhs[owner][component] += diffusion * boundaryVelocity[component] +
                        effectiveViscosity * dot(correctionVector, velocityGradient[owner][component]);
                }
                rhs[owner] += turbulentStressFlux(faceIndex, velocityGradient[owner], eddyViscosity);
            }
        }

        const double flux = massFlux_[faceIndex];
        const Vector2 central = centralFaceValue(mesh_, velocity_, velocityGradient, faceIndex);
        const Vector2 upwind = flux >= 0.0 ? velocity_[owner] : central;
        if (flux >= 0.0) matrix.add(owner, owner, flux);
        else rhs[owner] -= flux * central;
        const Vector2 deferred = (convectionBlending_ * flux) * (central - upwind);
        rhs[owner] -= deferred;
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

LinearSolveInfo SimpleSolver::correctPressure() {
    const int cellCount = static_cast<int>(mesh_.cells.size());
    SparseMatrix matrix(cellCount);
    std::vector<double> rhs(cellCount, 0.0);
    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density_);
        matrix.add(face.owner, face.owner, coefficient);
        matrix.add(face.owner, face.neighbour, -coefficient);
        matrix.add(face.neighbour, face.neighbour, coefficient);
        matrix.add(face.neighbour, face.owner, -coefficient);
        rhs[face.owner] -= massFlux_[faceIndex];
        rhs[face.neighbour] += massFlux_[faceIndex];
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        if (boundaryConditions_.fixedPressure(faceIndex)) {
            const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density_);
            matrix.add(face.owner, face.owner, coefficient);
        }
        rhs[face.owner] -= massFlux_[faceIndex];
    }
    if (!boundaryConditions_.hasFixedPressure()) {
        matrix.pinZero(0);
        rhs[0] = 0.0;
    }
    pressureSystem_.setMatrix(matrix);
    const LinearSolveInfo info = pressureSystem_.solve(rhs, pressureCorrection_.values(),
        config_.get("pressure_ksp", "cg"), config_.get("pressure_pc", "icc"),
        linearRelativeTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const int owner = mesh_.faces[faceIndex].owner;
        pressureCorrection_.boundary(faceIndex) = boundaryConditions_.fixedPressure(faceIndex) ? 0.0 : pressureCorrection_[owner];
        pressureCorrection_.boundaryDirichlet(faceIndex) = boundaryConditions_.fixedPressure(faceIndex);
    }
    const auto correctionGradient = gradient_->scalar(pressureCorrection_);
    for (int cell = 0; cell < cellCount; ++cell) {
        velocity_[cell] -= correctionGradient[cell] * (mesh_.cells[cell].volume / momentumDiagonal_[cell]);
        pressure_[cell] += pressureRelaxation_ * pressureCorrection_[cell];
    }
    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density_);
        massFlux_[faceIndex] += coefficient * (pressureCorrection_[face.owner] - pressureCorrection_[face.neighbour]);
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        if (boundaryConditions_.fixedPressure(faceIndex)) {
            const double coefficient = rhieChowCoefficient(mesh_, faceIndex, momentumDiagonal_, density_);
            massFlux_[faceIndex] += coefficient * pressureCorrection_[mesh_.faces[faceIndex].owner];
        }
    }
    updateBoundaries();
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
    for (int iteration = 1; iteration <= maximumIterations; ++iteration) {
        previousVelocity_.values() = velocity_.values();
        if (turbulence_) turbulence_->beginIteration();
        const LinearSolveInfo momentum = solveMomentum(dt);
        const LinearSolveInfo pressure = correctPressure();
        if (turbulence_) {
            turbulence_->solve(dt, hasOldOldVelocity_, velocity_, massFlux_, *gradient_);
        }
        const std::vector<double> imbalance = cellImbalance();
        diagnostics.iterations = iteration;
        diagnostics.momentumLinearIterations = momentum.iterations;
        diagnostics.pressureLinearIterations = pressure.iterations;
        diagnostics.turbulenceChange = turbulence_ ? turbulence_->change() : 0.0;
        diagnostics.turbulenceLinearIterations = turbulence_ ? turbulence_->linearIterations() : 0;
        diagnostics.velocityChange = 0.0;
        diagnostics.pressureChange = 0.0;
        diagnostics.maxCellImbalance = 0.0;
        for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
            diagnostics.velocityChange = std::max(diagnostics.velocityChange, norm(velocity_[cell] - previousVelocity_[cell]));
            diagnostics.pressureChange = std::max(diagnostics.pressureChange,
                std::abs(pressureRelaxation_ * pressureCorrection_[cell]));
            diagnostics.maxCellImbalance = std::max(diagnostics.maxCellImbalance, std::abs(imbalance[cell]));
        }
        diagnostics.netBoundaryFlux = netBoundaryFlux();
        diagnostics.converged = diagnostics.velocityChange < outerTolerance_ &&
            diagnostics.pressureChange < outerTolerance_ && diagnostics.maxCellImbalance < continuityTolerance_ &&
            (!turbulence_ || diagnostics.turbulenceChange < turbulenceTolerance);
        if (iteration == 1 || iteration % logInterval == 0 || diagnostics.converged || iteration == maximumIterations) {
            std::cout << (timeDerivative_.isSteady() ? "SIMPLE " : "  outer ")
                      << std::setw(6) << iteration << "  dU=" << std::scientific << diagnostics.velocityChange
                      << "  dp=" << diagnostics.pressureChange << "  dT=" << diagnostics.turbulenceChange
                      << "  continuity=" << diagnostics.maxCellImbalance
                      << "  linear=" << diagnostics.momentumLinearIterations + diagnostics.pressureLinearIterations +
                         diagnostics.turbulenceLinearIterations << '\n';
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
    if (turbulence_) turbulence_->updateTimeHistory(hasOldOldVelocity_);
    hasOldOldVelocity_ = true;
}

FlowForceDiagnostics SimpleSolver::forceDiagnostics() const {
    FlowForceDiagnostics result;
    if (forcePatch_ < 0) return result;
    const std::vector<VectorGradient> velocityGradient = gradient_->vector(velocity_);
    for (int faceIndex : mesh_.patches[forcePatch_].faces) {
        const Face& face = mesh_.faces[faceIndex];
        const auto& gradient = velocityGradient[face.owner];
        Vector2 viscousTraction;
        for (int component = 0; component < 2; ++component) {
            const double transpose = gradient[0][component] * face.normal.x + gradient[1][component] * face.normal.y;
            viscousTraction[component] = viscosity_ * face.area *
                (dot(gradient[component], face.normal) + transpose);
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
    imbalance.values() = cellImbalance();
    std::map<std::string, const ScalarField*> scalars{{"p", &pressure_}, {"mass_imbalance", &imbalance}};
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
    writeTecplot(file, mesh_, scalars, {{"U", &velocity_}}, name(), currentTime_);
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
    if (!append) output << "time,iterations,velocity_change,pressure_change,max_cell_imbalance,net_boundary_flux,turbulence_change,turbulence_linear_iterations,converged\n";
    output << time << ',' << diagnostics.iterations << ',' << diagnostics.velocityChange << ','
           << diagnostics.pressureChange << ',' << diagnostics.maxCellImbalance << ','
           << diagnostics.netBoundaryFlux << ',' << diagnostics.turbulenceChange << ','
           << diagnostics.turbulenceLinearIterations << ',' << (diagnostics.converged ? 1 : 0) << '\n';
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
        const FlowDiagnostics result = solveStep(stepSize, config_.integer("maximum_iterations", 100));
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
        if (!result.converged) throw std::runtime_error("SIMPLE did not converge at time " + std::to_string(currentTime_));
    }
    if (currentTime_ < endTime - 1.0e-12) throw std::runtime_error("maximum_steps reached before end_time");
    std::cout << "Tecplot output: " << output.parent_path() << '\n';
}

} // namespace mp
