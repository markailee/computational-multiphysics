#pragma once

#include "flow/BoundaryConditions.hpp"
#include "flow/EnergyModel.hpp"
#include "flow/MultiSpeciesModel.hpp"
#include "flow/TurbulenceModel.hpp"
#include "flow/multiphase/MultiphaseModel.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/Config.hpp"
#include "fvm/Gradient.hpp"
#include "fvm/PhysicsSolver.hpp"
#include "fvm/TimeDerivative.hpp"
#include "material/MaterialModel.hpp"

#include <array>
#include <filesystem>
#include <memory>
#include <vector>

namespace mp {

enum class FlowModel { Incompressible, Boussinesq, Compressible };

struct FlowForceDiagnostics {
    Vector2 total;
    Vector2 pressure;
    Vector2 viscous;
};

struct FlowDiagnostics {
    int iterations = 0;
    double velocityChange = 0.0;
    double pressureChange = 0.0;
    double maxCellImbalance = 0.0;
    double netBoundaryFlux = 0.0;
    int momentumLinearIterations = 0;
    int pressureLinearIterations = 0;
    double turbulenceChange = 0.0;
    int turbulenceLinearIterations = 0;
    double speciesChange = 0.0;
    double speciesCompositionError = 0.0;
    double minimumMassFraction = 0.0;
    double maximumMassFraction = 1.0;
    int speciesLinearIterations = 0;
    double temperatureChange = 0.0;
    double densityChange = 0.0;
    int energyLinearIterations = 0;
    bool converged = false;
};

class SimpleSolver final : public PhysicsSolver {
public:
    SimpleSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory);

    std::string name() const override {
        if (flowModel_ == FlowModel::Boussinesq) return "boussinesq-simple";
        if (flowModel_ == FlowModel::Compressible) return "compressible-simple";
        return "incompressible-simple";
    }
    void initialize() override;
    void advance(double dt) override;
    void write(const std::filesystem::path& file) const override;
    void run();

    const FlowDiagnostics& diagnostics() const { return diagnostics_; }
    double currentTime() const { return currentTime_; }

private:
    FlowDiagnostics solveStep(double dt, int maximumIterations);
    void updateBoundaries();
    void readRestart(const std::filesystem::path& file);
    void updateMaterialProperties();
    double faceDensity(std::size_t faceIndex) const;
    void computeMassFlux(bool rhieChow);
    std::array<LinearSolveInfo, 2> solveMomentum(double dt);
    LinearSolveInfo correctPressure(double dt);
    std::vector<double> cellImbalance() const;
    double netBoundaryFlux() const;
    void updateTimeHistory();
    FlowForceDiagnostics forceDiagnostics() const;
    std::filesystem::path forceFile() const;
    void writeForceHistory(const std::filesystem::path& file, double time, bool append) const;
    std::filesystem::path outputFile() const;
    void writeSummary(const std::filesystem::path& file, double time,
                      const FlowDiagnostics& diagnostics, bool append) const;
    void writeWallDiagnostics(const std::filesystem::path& file) const;

    const Mesh& mesh_;
    const Config& config_;
    std::filesystem::path caseDirectory_;
    FlowModel flowModel_;
    std::shared_ptr<MaterialModel> material_;
    double density_;
    double viscosity_;
    double velocityRelaxation_;
    double pressureRelaxation_;
    bool barthJespersenLimiter_;
    double outerTolerance_;
    double continuityTolerance_;
    double linearRelativeTolerance_;
    double linearAbsoluteTolerance_;
    int maximumLinearIterations_;
    Vector2 bodyForce_;
    Vector2 gravity_;
    double thermalExpansion_;
    double referenceTemperature_;
    int forcePatch_ = -1;
    double forceReferenceVelocity_ = 1.0;
    double forceReferenceLength_ = 1.0;
    TimeDerivative timeDerivative_;
    FlowBoundaryConditions boundaryConditions_;
    std::unique_ptr<GradientScheme> gradient_;
    VectorField velocity_;
    ScalarField pressure_;
    ScalarField pressureCorrection_;
    VectorField previousVelocity_;
    VectorField oldVelocity_;
    VectorField oldOldVelocity_;
    ScalarField oldPressure_;
    ScalarField oldOldPressure_;
    std::vector<double> massFlux_;
    std::vector<double> densityField_;
    std::vector<double> previousDensityField_;
    std::vector<double> oldDensityField_;
    std::vector<double> oldOldDensityField_;
    std::vector<double> viscosityField_;
    std::vector<double> conductivityField_;
    std::vector<double> heatCapacityField_;
    std::vector<double> pressureDerivativeField_;
    std::vector<double> momentumDiagonal_;
    std::vector<double> lastImbalance_;
    PetscLinearSystem momentumSystem_;
    PetscLinearSystem pressureSystem_;
    std::unique_ptr<EnergyModel> energy_;
    std::unique_ptr<TurbulenceModel> turbulence_;
    std::unique_ptr<MultiSpeciesModel> multiSpecies_;
    std::unique_ptr<MultiphaseModel> multiphase_;
    FlowDiagnostics diagnostics_;
    double currentTime_ = 0.0;
    bool initialized_ = false;
    bool hasOldOldVelocity_ = false;
};

} // namespace mp
