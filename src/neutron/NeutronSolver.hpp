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
enum class NeutronProblem { FixedSource, Eigenvalue };
enum class NeutronBoundaryKind {
    Vacuum,
    Reflective,
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
    std::vector<double> scattering;
};

struct NeutronDirection {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double weight = 0.0;
};

struct NeutronBoundaryCondition {
    NeutronBoundaryKind kind = NeutronBoundaryKind::Vacuum;
    double value = 0.0;
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
        double balance = 0.0;
        int sourceIterations = 0;
        int linearIterations = 0;
        double linearResidual = 0.0;
    };

    struct Balance {
        double production = 0.0;
        double absorption = 0.0;
        double leakage = 0.0;
        double incoming = 0.0;
        double externalSource = 0.0;
        double power = 0.0;
        double residual = 0.0;
        double relativeResidual = 0.0;
        double eigenvalue = 0.0;
    };

    void readConfiguration();
    void buildQuadrature();
    void buildBoundaries();
    void buildDiffusionOperator();
    void solveFixedSource();
    void solveEigenvalue();
    void solveSnInner(double eigenvalue, const std::vector<double>* fissionState);
    void solveDiffusionStep(double eigenvalue, const std::vector<double>& fissionState);
    void assembleSnEquation(int group, int direction, double eigenvalue,
                            const std::vector<double>& sourceScalar,
                            const std::vector<double>* fissionState,
                            const std::vector<double>& previousAngular,
                            SparseMatrix& matrix, std::vector<double>& rhs) const;
    std::vector<double> diffusionRightHandSide(double eigenvalue,
                                               const std::vector<double>& fissionState) const;
    void updateScalarFlux();
    double fluxChange(const std::vector<double>& previous) const;
    double fissionProduction(const std::vector<double>& scalar) const;
    Balance balance(double eigenvalue, bool eigenproblem) const;
    void normalize(double factor);
    void updateOutputFields() const;
    void writeCells(const std::filesystem::path& file) const;
    void writeSummary(const std::filesystem::path& file) const;
    void writeHistory(const std::filesystem::path& file) const;
    void writeDirections(const std::filesystem::path& file) const;
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
    int cellCount_ = 0;
    int groupCount_ = 1;
    int directionCount_ = 0;
    int polarCount_ = 2;
    int azimuthCount_ = 4;
    int maximumIterations_ = 100;
    int maximumSourceIterations_ = 100;
    int maximumLinearIterations_ = 1000;
    int logInterval_ = 1;
    double tolerance_ = 1.0e-8;
    double eigenTolerance_ = 1.0e-10;
    double sourceTolerance_ = 1.0e-9;
    double balanceTolerance_ = 1.0e-8;
    double linearTolerance_ = 1.0e-10;
    double linearAbsoluteTolerance_ = 1.0e-12;
    double initialFlux_ = 1.0;
    double targetPower_ = 0.0;
    double eigenvalue_ = 1.0;
    std::string linearSolver_ = "gmres";
    std::string preconditioner_ = "ilu";
    std::string normalization_ = "fixed_source";
    std::vector<NeutronMaterial> materials_;
    std::vector<int> cellMaterials_;
    std::vector<NeutronDirection> directions_;
    std::vector<NeutronBoundaryCondition> boundaries_;
    std::vector<int> reflectionDirections_;
    std::vector<double> angularFlux_;
    std::vector<double> scalarFlux_;
    std::vector<double> diffusionBoundaryRhs_;
    SparseMatrix diffusionMatrix_;
    std::unique_ptr<PetscLinearSystem> system_;
    std::vector<IterationRecord> history_;
    int lastSourceIterations_ = 0;
    int totalLinearIterations_ = 0;
    double lastLinearResidual_ = 0.0;
    bool initialized_ = false;
    bool solved_ = false;
    mutable std::map<std::string, std::unique_ptr<ScalarField>> outputScalars_;
    mutable std::unique_ptr<VectorField> outputCurrent_;
};

} // namespace mp
