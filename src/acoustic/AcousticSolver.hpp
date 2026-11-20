#pragma once

#include "acoustic/AcousticBoundaryConditions.hpp"
#include "fvm/Config.hpp"
#include "fvm/PhysicsSolver.hpp"
#include "fvm/PetscLinearSystem.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace mp {

class AcousticSolver final : public PhysicsSolver {
public:
    AcousticSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory);

    std::string name() const override { return "node-centered-acoustic"; }
    void initialize() override;
    void advance(double dt) override;
    void write(const std::filesystem::path& file) const override;
    void run();

private:
    struct CellProperties {
        double density = 1.0;
        double soundSpeed = 1.0;
        Vector2 mach;
    };

    struct GeometryCheck {
        double area = 0.0;
        double subAreaSum = 0.0;
        double gradientSum = 0.0;
    };

    enum class InitialKind { Zero, CavityMode, PlaneWave, GaussianPulse };
    enum class SourceKind { None, HarmonicMonopole, FlowDipole };

    struct FlowForceSample {
        double time = 0.0;
        Vector2 force;
    };

    void assembleOperators();
    std::vector<Vector2> shapeGradients(int cellIndex, const Vector2& point) const;
    double planeWavePressure(const Vector2& point, double time) const;
    double planeWaveRate(const Vector2& point, double time) const;
    double initialPressure(const Vector2& point) const;
    double initialPressureRate(const Vector2& point) const;
    double prescribedPressure(int node, double time) const;
    double prescribedRate(int node, double time) const;
    double prescribedAcceleration(int node, double time) const;
    bool isPrescribed(int node) const { return prescribedKind_[node] >= 0; }
    void applyPrescribedPressure(SparseMatrix& matrix, std::vector<double>& rhs,
                                 const std::vector<double>& values) const;
    std::vector<double> source(double time) const;
    void loadFlowSource(const std::filesystem::path& file);
    Vector2 flowForce(double time) const;
    std::vector<double> multiply(const SparseMatrix& matrix, const std::vector<double>& values) const;
    double energy() const;
    double pressureScale() const;
    std::filesystem::path outputFile() const;
    void writeSummary(const std::filesystem::path& file, bool append, int linearIterations,
                      double linearResidual, double courant) const;
    void writeProbe(const std::filesystem::path& file, bool append) const;
    void writeProbeArray(const std::filesystem::path& file, bool append) const;
    void writeGeometry(const std::filesystem::path& file) const;

    const Mesh& mesh_;
    const Config& config_;
    std::filesystem::path caseDirectory_;
    AcousticBoundaryConditions boundaryConditions_;
    std::vector<CellProperties> properties_;
    std::vector<double> mass_;
    std::vector<double> dualArea_;
    SparseMatrix damping_;
    SparseMatrix stiffness_;
    PetscLinearSystem system_;
    std::vector<double> pressure_;
    std::vector<double> pressureRate_;
    std::vector<double> pressureAcceleration_;
    std::vector<double> previousPressure_;
    std::vector<int> prescribedKind_;
    std::vector<double> prescribedValue_;
    std::vector<double> sourceWeights_;
    std::vector<Vector2> dipoleWeights_;
    std::vector<FlowForceSample> flowForceHistory_;
    std::vector<GeometryCheck> geometryChecks_;
    InitialKind initialKind_ = InitialKind::Zero;
    SourceKind sourceKind_ = SourceKind::None;
    std::string timeScheme_;
    Vector2 meanVelocity_;
    Vector2 waveDirection_{1.0, 0.0};
    Vector2 waveVector_{6.283185307179586, 0.0};
    Vector2 cavityOrigin_;
    Vector2 cavityLengths_{1.0, 1.0};
    Vector2 modeIndices_{2.0, 1.0};
    Vector2 gaussianCenter_;
    Vector2 sourcePosition_;
    double initialAmplitude_ = 1.0;
    double initialPressureOffset_ = 0.0;
    double initialUniformRate_ = 0.0;
    double wavePhase_ = 0.0;
    double waveFrequency_ = 6.283185307179586;
    double gaussianWidth_ = 0.1;
    double newmarkBeta_ = 0.25;
    double newmarkGamma_ = 0.5;
    double soundSpeed_ = 1.0;
    double sourceAmplitude_ = 0.0;
    double sourceFrequency_ = 1.0;
    double sourcePhase_ = 0.0;
    double sourceWidth_ = 0.1;
    double sourceScale_ = 1.0;
    double sourceRampTime_ = 0.0;
    double linearRelativeTolerance_ = 1.0e-10;
    double linearAbsoluteTolerance_ = 1.0e-13;
    double currentTime_ = 0.0;
    double minimumEdge_ = 0.0;
    int maximumLinearIterations_ = 1000;
    int lastLinearIterations_ = 0;
    int probeNode_ = -1;
    std::vector<int> probeNodes_;
    double lastLinearResidual_ = 0.0;
    bool initialized_ = false;
    bool hasPreviousPressure_ = false;
};

} // namespace mp
