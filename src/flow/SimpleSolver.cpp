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
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOldVelocity_);
    const int cellCount = static_cast<int>(mesh_.cells.size());
    SparseMatrix matrix(cellCount);
    std::vector<double> rhsX(cellCount, 0.0), rhsY(cellCount, 0.0);
    std::vector<Vector2> rhs(cellCount, bodyForce_ * density_);

    for (int cell = 0; cell < cellCount; ++cell) {
        const double volume = mesh_.cells[cell].volume;
        matrix.add(cell, cell, density_ * volume * time.current);
        rhs[cell] += oldVelocity_[cell] * (density_ * volume * time.previous) +
                     oldOldVelocity_[cell] * (density_ * volume * time.previousPrevious) -
                     pressureGradient[cell] * volume;
    }

    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        const double diffusion = diffusionConductance(mesh_, faceIndex, viscosity_);
        matrix.add(owner, owner, diffusion);
        matrix.add(owner, neighbour, -diffusion);
        matrix.add(neighbour, neighbour, diffusion);
        matrix.add(neighbour, owner, -diffusion);

        const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
        for (int component = 0; component < 2; ++component) {
            const Vector2 faceGradient = face.ownerWeight * velocityGradient[owner][component] +
                (1.0 - face.ownerWeight) * velocityGradient[neighbour][component];
            const double correction = viscosity_ * dot(correctionVector, faceGradient);
            rhs[owner][component] += correction;
            rhs[neighbour][component] -= correction;
        }

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
            const double diffusion = diffusionConductance(mesh_, faceIndex, viscosity_);
            matrix.add(owner, owner, diffusion);
            const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
            for (int component = 0; component < 2; ++component) {
                rhs[owner][component] += diffusion * boundaryVelocity[component] +
                    viscosity_ * dot(correctionVector, velocityGradient[owner][component]);
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
    for (int iteration = 1; iteration <= maximumIterations; ++iteration) {
        previousVelocity_.values() = velocity_.values();
        const LinearSolveInfo momentum = solveMomentum(dt);
        const LinearSolveInfo pressure = correctPressure();
        const std::vector<double> imbalance = cellImbalance();
        diagnostics.iterations = iteration;
        diagnostics.momentumLinearIterations = momentum.iterations;
        diagnostics.pressureLinearIterations = pressure.iterations;
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
            diagnostics.pressureChange < outerTolerance_ && diagnostics.maxCellImbalance < continuityTolerance_;
        if (iteration == 1 || iteration % logInterval == 0 || diagnostics.converged || iteration == maximumIterations) {
            std::cout << (timeDerivative_.isSteady() ? "SIMPLE " : "  outer ")
                      << std::setw(6) << iteration << "  dU=" << std::scientific << diagnostics.velocityChange
                      << "  dp=" << diagnostics.pressureChange << "  continuity=" << diagnostics.maxCellImbalance
                      << "  linear=" << diagnostics.momentumLinearIterations + diagnostics.pressureLinearIterations << '\n';
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
    hasOldOldVelocity_ = true;
}

std::filesystem::path SimpleSolver::outputFile() const {
    std::filesystem::path output = config_.get("output", "results/flow.plt");
    return output.is_absolute() ? output : caseDirectory_ / output;
}

void SimpleSolver::write(const std::filesystem::path& file) const {
    ScalarField imbalance("mass_imbalance", mesh_);
    imbalance.values() = cellImbalance();
    writeTecplot(file, mesh_, {{"p", &pressure_}, {"mass_imbalance", &imbalance}}, {{"U", &velocity_}}, name(), currentTime_);
}

void SimpleSolver::writeSummary(const std::filesystem::path& file, double time,
                                const FlowDiagnostics& diagnostics, bool append) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file, append ? std::ios::app : std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write flow summary: " + file.string());
    output << std::setprecision(12) << std::scientific;
    if (!append) output << "time,iterations,velocity_change,pressure_change,max_cell_imbalance,net_boundary_flux,converged\n";
    output << time << ',' << diagnostics.iterations << ',' << diagnostics.velocityChange << ','
           << diagnostics.pressureChange << ',' << diagnostics.maxCellImbalance << ','
           << diagnostics.netBoundaryFlux << ',' << (diagnostics.converged ? 1 : 0) << '\n';
}

void SimpleSolver::run() {
    initialize();
    const std::filesystem::path output = outputFile();
    const std::filesystem::path summary = output.parent_path() / (output.stem().string() + "_summary.csv");
    if (timeDerivative_.isSteady()) {
        const FlowDiagnostics result = solveStep(0.0, config_.integer("maximum_iterations", 2000));
        write(output);
        writeSummary(summary, 0.0, result, false);
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
    int step = 0;
    while (currentTime_ < endTime - 1.0e-12 && step < maximumSteps) {
        const double stepSize = std::min(dt, endTime - currentTime_);
        const FlowDiagnostics result = solveStep(stepSize, config_.integer("maximum_iterations", 100));
        currentTime_ += stepSize;
        ++step;
        updateTimeHistory();
        writeSummary(summary, currentTime_, result, wroteSummary);
        wroteSummary = true;
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
