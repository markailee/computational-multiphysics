#pragma once

#include "fvm/Config.hpp"
#include "fvm/Field.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/PhysicsSolver.hpp"

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace mp {

enum class NeutronModel { DiscreteOrdinates, Diffusion };
enum class NeutronProblem { FixedSource, Eigenvalue, Transient };
enum class NeutronQuadratureKind { ProductMidpoint, ProductAxis, LevelSymmetric };
enum class NeutronSnScheme { Upwind, SecondOrder };
enum class NeutronSnAcceleration { None, DSA };
enum class NeutronTimeScheme { BackwardEuler, CrankNicolson };

enum class NeutronBoundaryKind {
    Vacuum,
    Reflective,
    WhiteReflective,
    Albedo,
    IncidentFlux,
    IncidentCurrent,
    FixedScalarFlux,
    Marshak,
    ZeroCurrent,
};

struct NeutronMaterial {
    std::string name;
    std::vector<double> total;
    std::vector<double> absorption;
    std::vector<double> fission;
    std::vector<double> nuSigmaFission;
    std::vector<double> kappaSigmaFission;
    std::vector<double> chi;
    std::vector<double> diffusion;
    std::vector<double> externalSource;
    std::vector<double> scattering;   // row = destination g, col = source g'
    std::vector<double> scatteringP1; // row = destination g, col = source g'
};

struct NeutronDirection {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double weight = 0.0;
};

struct NeutronBoundaryCondition {
    NeutronBoundaryKind kind = NeutronBoundaryKind::Vacuum;
    std::vector<double> values;
    double whiteFraction = 0.0;

    double value(int group) const {
        if (values.empty()) return 0.0;
        return group < static_cast<int>(values.size()) ? values[group] : values.front();
    }
};

struct NeutronRegion {
    std::string name;
    int material = 0;
    bool isCircle = false;
    Vector2 lower;
    Vector2 upper;
    Vector2 centre;
    double radius = 0.0;
};

struct NeutronPerturbation {
    bool enabled = false;
    std::string regionName;
    std::string materialName;
    std::vector<double> deltaAbsorption;
    std::vector<double> deltaNuSigmaFission;
    std::vector<double> deltaDiffusion;
    std::vector<double> deltaScattering;
};

struct NeutronTransientEvent {
    std::string name;
    std::string regionName;
    std::string materialName;
    bool isRamp = false;
    double startTime = 0.0;
    double endTime = 0.0;
    std::vector<double> deltaAbsorption;
    std::vector<double> deltaNuSigmaFission;
    std::vector<double> deltaDiffusion;
    std::vector<double> deltaScattering;
};

struct NeutronEigenmode {
    int index = 0;
    double eigenvalue = 0.0;
    double imaginary = 0.0;
    double dominanceRatio = 1.0;
    double residual = 0.0;
    std::vector<double> scalarFlux;
};

struct NeutronKineticsSummary {
    bool adjointComputed = false;
    double adjointEigenvalue = 0.0;
    double adjointDiscrepancy = 0.0;
    bool kineticsComputed = false;
    double betaEffective = 0.0;
    double generationTime = 0.0;
    double promptLifetime = 0.0;
    std::vector<double> betaEffectiveGroups;
    bool perturbationComputed = false;
    double perturbationDeltaKOverK = 0.0;
    double perturbationDeltaRho = 0.0;
};

class NeutronSolver final : public PhysicsSolver {
public:
    NeutronSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory);

    std::string name() const override { return model_ == NeutronModel::Diffusion ? "neutron-diffusion" : "neutron-sn"; }
    void initialize() override;
    void advance(double dt) override;
    void write(const std::filesystem::path& file) const override;
    void run();

private:
    struct IterationRecord {
        int iteration = 0;
        double eigenvalue = 1.0;
        double fluxChange = 0.0;
        double eigenResidual = 0.0;
        double balance = 0.0;
        int sourceIterations = 0;
        int linearIterations = 0;
        double linearResidual = 0.0;
    };

    struct TransientRecord {
        int step = 0;
        double time = 0.0;
        double totalPower = 0.0;
        double relativePower = 1.0;
        double peakPowerDensity = 0.0;
        double fluxIntegral = 0.0;
        double precursorIntegral = 0.0;
        int linearIterations = 0;
    };

    struct Balance {
        double production = 0.0;
        double absorption = 0.0;
        double axialLeakage = 0.0;
        double leakage = 0.0;
        double incoming = 0.0;
        double externalSource = 0.0;
        double power = 0.0;
        double residual = 0.0;
        double relativeResidual = 0.0;
        double eigenvalue = 0.0;
    };

    // Setup and quadrature (NeutronQuadrature.cpp)
    void buildQuadrature();
    void buildSnSweepOrder();
    void buildReflectionMaps();

    // Core configuration and operators (NeutronSolver.cpp)
    void readConfiguration();
    void buildBoundaries();
    void buildDiffusionOperator();
    void buildDiffusionOperatorFromMaterials(const std::vector<NeutronMaterial>& materials,
                                             double invDtScale,
                                             SparseMatrix& matrix,
                                             std::vector<double>& boundaryRhs) const;
    void buildDsaOperators();

    // Solvers (NeutronSolver.cpp)
    void solveFixedSource();
    void solveEigenvalue();
    void solveSnInner(double eigenvalue, const std::vector<double>* fissionState, bool includeExternal);
    void sweepSnDirection(int group, int direction, double eigenvalue,
                          const std::vector<double>& sourceScalar,
                          const std::vector<double>& sourceCurrentX,
                          const std::vector<double>& sourceCurrentY,
                          const std::vector<double>* fissionState,
                          const std::vector<double>& previousAngular,
                          const std::vector<double>& whiteIncoming,
                          const std::vector<double>& hoCorrection,
                          bool includeExternal,
                          std::vector<double>& angularOut);
    void computeSnHighOrderCorrection(int group, int direction,
                                      const std::vector<double>& angularState,
                                      std::vector<double>& faceCorrection,
                                      std::vector<double>& cellRhsCorrection) const;
    void applyDsaCorrection(int group, const std::vector<double>& oldScalarGroup);
    void solveDiffusionStep(double eigenvalue, const std::vector<double>& fissionState, bool includeExternal);
    void assembleSnEquation(int group, int direction, double eigenvalue,
                            const std::vector<double>& sourceScalar,
                            const std::vector<double>& sourceCurrentX,
                            const std::vector<double>& sourceCurrentY,
                            const std::vector<double>* fissionState,
                            const std::vector<double>& previousAngular,
                            const std::vector<double>& whiteIncoming,
                            const std::vector<double>& hoCorrection,
                            bool includeExternal,
                            SparseMatrix& matrix, std::vector<double>& rhs) const;
    std::vector<double> diffusionRightHandSide(double eigenvalue,
                                               const std::vector<double>& fissionState,
                                               bool includeExternal) const;
    void updateScalarAndCurrentMoments();
    double fluxChange(const std::vector<double>& previous) const;
    double eigenResidual(double eigenvalue) const;
    double fissionProduction(const std::vector<double>& scalar) const;
    Balance balance(double eigenvalue, bool eigenproblem) const;
    void normalize(double factor);

    // Advanced physics (NeutronAdvanced.cpp)
    std::vector<double> applyFissionOperator(const std::vector<double>& scalar) const;
    std::vector<double> applyKOperator(const std::vector<double>& scalar);
    void solveHigherModes();
    void solveAdjoint();
    void solveDiffusionAdjoint();
    void solveSnAdjoint();
    void computeKineticsAndPerturbation();
    void solveTransient();
    std::vector<NeutronMaterial> materialsAtTime(double time) const;
    bool cellMatchesFilter(int cell, const std::string& regionName, const std::string& materialName) const;

    // Post-processing & output (NeutronSolver.cpp)
    void updateOutputFields() const;
    void writeCells(const std::filesystem::path& file) const;
    void writeSummary(const std::filesystem::path& file) const;
    void writeHistory(const std::filesystem::path& file) const;
    void writeDirections(const std::filesystem::path& file) const;
    void writeRegions(const std::filesystem::path& file) const;
    void writeModes(const std::filesystem::path& file) const;
    void writeTransient(const std::filesystem::path& file) const;
    std::filesystem::path outputPath(const std::string& suffix = {}) const;

    int scalarIndex(int group, int cell) const { return group * cellCount_ + cell; }
    int angularIndex(int group, int direction, int cell) const {
        return (group * directionCount_ + direction) * cellCount_ + cell;
    }
    int reflectionIndex(int face, int direction) const { return face * directionCount_ + direction; }

    const Mesh& mesh_;
    const Config& config_;
    std::filesystem::path caseDirectory_;
    NeutronModel model_ = NeutronModel::DiscreteOrdinates;
    NeutronProblem problem_ = NeutronProblem::FixedSource;
    NeutronQuadratureKind quadratureKind_ = NeutronQuadratureKind::ProductMidpoint;
    NeutronSnScheme snScheme_ = NeutronSnScheme::Upwind;
    NeutronSnAcceleration snAcceleration_ = NeutronSnAcceleration::None;
    NeutronTimeScheme timeScheme_ = NeutronTimeScheme::BackwardEuler;
    std::string snLimiter_ = "barth_jespersen";

    int cellCount_ = 0;
    int groupCount_ = 1;
    int directionCount_ = 0;
    int snOrder_ = 0;
    int polarCount_ = 2;
    int azimuthCount_ = 4;
    int maximumIterations_ = 100;
    int maximumSourceIterations_ = 100;
    int maximumLinearIterations_ = 1000;
    int eigenModesRequested_ = 1;
    int maximumModeIterations_ = 150;
    int logInterval_ = 1;
    int transientOutputInterval_ = 1;

    double tolerance_ = 1.0e-8;
    double eigenTolerance_ = 1.0e-10;
    double sourceTolerance_ = 1.0e-9;
    double balanceTolerance_ = 1.0e-8;
    double modeTolerance_ = 1.0e-8;
    double linearTolerance_ = 1.0e-10;
    double linearAbsoluteTolerance_ = 1.0e-12;
    double initialFlux_ = 1.0;
    double targetPower_ = 0.0;
    double axialBuckling_ = 0.0;
    double eigenvalue_ = 1.0;
    double lastEigenResidual_ = 0.0;
    double timeStep_ = 0.01;
    double endTime_ = 0.5;
    bool criticalInitialState_ = true;
    bool computeAdjoint_ = false;
    bool hasP1Scattering_ = false;

    std::string linearSolver_ = "gmres";
    std::string preconditioner_ = "ilu";
    std::string normalization_ = "fixed_source";

    std::vector<NeutronMaterial> materials_;
    std::vector<NeutronRegion> regions_;
    std::vector<int> cellMaterials_;
    std::vector<int> cellRegions_;
    std::vector<NeutronDirection> directions_;
    std::vector<NeutronBoundaryCondition> boundaries_;
    std::vector<int> reflectionDirections_;
    std::vector<double> faceIncomingWeightSum_;

    // Precomputed SN topological sweep data
    std::vector<std::vector<int>> sweepOrder_;
    std::vector<bool> directionAcyclic_;
    std::vector<double> outgoingBoundaryAngular_; // [group * dir * boundaryFace] for exact HO balance

    // Kinetics / delayed neutron data
    std::vector<double> neutronVelocity_;
    std::vector<double> delayedBeta_;
    std::vector<double> delayedLambda_;
    std::vector<double> delayedChi_; // [delayedGroup * energyGroup]
    NeutronPerturbation perturbation_;
    std::vector<NeutronTransientEvent> transientEvents_;

    // State vectors
    std::vector<double> angularFlux_;
    std::vector<double> scalarFlux_;
    std::vector<double> currentX_;
    std::vector<double> currentY_;
    std::vector<double> adjointScalarFlux_;
    std::vector<double> adjointAngularFlux_;
    std::vector<double> diffusionBoundaryRhs_;
    SparseMatrix diffusionMatrix_;
    std::unique_ptr<PetscLinearSystem> system_;
    std::vector<std::unique_ptr<PetscLinearSystem>> dsaSystems_;

    // Diagnostics & advanced results
    std::vector<IterationRecord> history_;
    std::vector<NeutronEigenmode> modes_;
    NeutronKineticsSummary kinetics_;
    std::vector<TransientRecord> transientHistory_;
    int lastSourceIterations_ = 0;
    int totalLinearIterations_ = 0;
    double lastLinearResidual_ = 0.0;
    bool initialized_ = false;
    bool solved_ = false;

    mutable std::map<std::string, std::unique_ptr<ScalarField>> outputScalars_;
    mutable std::unique_ptr<VectorField> outputCurrent_;
};

} // namespace mp
