#pragma once

#include "flow/BoundaryConditions.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/Config.hpp"
#include "fvm/Gradient.hpp"
#include "fvm/PhysicsSolver.hpp"
#include "fvm/TimeDerivative.hpp"

#include <filesystem>
#include <memory>
#include <vector>

namespace mp {

struct FlowDiagnostics {
    int iterations = 0;
    double velocityChange = 0.0;
    double pressureChange = 0.0;
    double maxCellImbalance = 0.0;
    double netBoundaryFlux = 0.0;
    int momentumLinearIterations = 0;
    int pressureLinearIterations = 0;
    bool converged = false;
};

class SimpleSolver final : public PhysicsSolver {
public:
    SimpleSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory);

    std::string name() const override { return "incompressible-simple"; }
    void initialize() override;
    void advance(double dt) override;
    void write(const std::filesystem::path& file) const override;
    void run();

    const FlowDiagnostics& diagnostics() const { return diagnostics_; }
    double currentTime() const { return currentTime_; }

private:
    FlowDiagnostics solveStep(double dt, int maximumIterations);
    void updateBoundaries();
    void computeMassFlux(bool rhieChow);
    LinearSolveInfo solveMomentum(double dt);
    LinearSolveInfo correctPressure();
    std::vector<double> cellImbalance() const;
    double netBoundaryFlux() const;
    void updateTimeHistory();
    std::filesystem::path outputFile() const;
    void writeSummary(const std::filesystem::path& file, double time,
                      const FlowDiagnostics& diagnostics, bool append) const;

    const Mesh& mesh_;
    const Config& config_;
    std::filesystem::path caseDirectory_;
    double density_;
    double viscosity_;
    double velocityRelaxation_;
    double pressureRelaxation_;
    double convectionBlending_;
    double outerTolerance_;
    double continuityTolerance_;
    double linearRelativeTolerance_;
    double linearAbsoluteTolerance_;
    int maximumLinearIterations_;
    Vector2 bodyForce_;
    TimeDerivative timeDerivative_;
    FlowBoundaryConditions boundaryConditions_;
    std::unique_ptr<GradientScheme> gradient_;
    VectorField velocity_;
    ScalarField pressure_;
    ScalarField pressureCorrection_;
    VectorField previousVelocity_;
    VectorField oldVelocity_;
    VectorField oldOldVelocity_;
    std::vector<double> massFlux_;
    std::vector<double> momentumDiagonal_;
    PetscLinearSystem momentumSystem_;
    PetscLinearSystem pressureSystem_;
    FlowDiagnostics diagnostics_;
    double currentTime_ = 0.0;
    bool initialized_ = false;
    bool hasOldOldVelocity_ = false;
};

} // namespace mp
