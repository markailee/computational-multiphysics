#pragma once

#include "flow/BoundaryConditions.hpp"
#include "fvm/Field.hpp"
#include "fvm/Gradient.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/TimeDerivative.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mp {

enum class TurbulenceModelKind { KEpsilon, SSTMenter2003 };

enum class TurbulenceBoundaryType { FixedValue, ZeroGradient, WallFunction };

struct TurbulenceBoundaryCondition {
    TurbulenceBoundaryType type = TurbulenceBoundaryType::ZeroGradient;
    double value = 0.0;
};

struct TurbulenceWallDiagnostics {
    std::size_t face = 0;
    int owner = -1;
    double distance = 0.0;
    Vector2 shear;
    double frictionVelocity = 0.0;
    double yPlus = 0.0;
};

class TurbulenceModel {
public:
    TurbulenceModel(const Mesh& mesh, const Config& config, const FlowBoundaryConditions& boundaryConditions,
                    double density, double viscosity, TurbulenceModelKind kind);

    std::string name() const;
    TurbulenceModelKind kind() const { return kind_; }
    bool usesStandardWallFunction() const { return kind_ == TurbulenceModelKind::KEpsilon; }

    void updateBoundaries();
    void beginIteration();
    void solve(double dt, bool hasOldOld, const VectorField& velocity,
               const std::vector<double>& massFlux, const GradientScheme& gradient);
    void updateTimeHistory(bool& hasOldOld);

    double kineticEnergy(std::size_t cell) const { return k_[cell]; }
    double eddyViscosity(std::size_t cell) const { return eddyViscosity_[cell]; }
    double faceEddyViscosity(std::size_t face) const { return eddyViscosity_.boundary(face); }
    double wallDistance(std::size_t cell) const { return wallDistance_[cell]; }
    Vector2 wallShearStress(std::size_t face, const Vector2& cellVelocity) const;
    TurbulenceWallDiagnostics wallDiagnostics(std::size_t face, const Vector2& cellVelocity) const;

    const ScalarField& k() const { return k_; }
    const ScalarField& scale() const { return kind_ == TurbulenceModelKind::KEpsilon ? epsilon_ : omega_; }
    const ScalarField& turbulentViscosity() const { return eddyViscosity_; }
    const ScalarField& blendingF1() const { return blendingF1_; }
    const ScalarField& blendingF2() const { return blendingF2_; }
    double change() const { return change_; }
    int linearIterations() const { return linearIterations_; }
    double linearResidual() const { return linearResidual_; }

private:
    TurbulenceBoundaryCondition readBoundary(const std::string& patch, const std::string& variable,
                                             bool noSlip) const;
    void updateFieldBoundaries(ScalarField& field, const std::vector<TurbulenceBoundaryCondition>& conditions);
    void computeProduction(const std::vector<VectorGradient>& velocityGradient);
    void updateSSTProperties(const GradientScheme& gradient,
                             const std::vector<VectorGradient>& velocityGradient);
    void solveKEpsilon(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                       const GradientScheme& gradient);
    void solveSST(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                  const GradientScheme& gradient);
    LinearSolveInfo solveField(ScalarField& field, ScalarField& previous, ScalarField& old,
                               ScalarField& oldOld, PetscLinearSystem& system,
                               const std::vector<double>& massFlux, const GradientScheme& gradient,
                               const std::vector<double>& diffusivity, const std::vector<double>& source,
                               const std::vector<double>& destruction, const std::vector<char>& constrained,
                               const std::vector<double>& constrainedValues, double dt, bool hasOldOld);
    void clip(ScalarField& field, double floor);
    void updateChange();
    double wallFunctionViscosity(std::size_t face, const Vector2& cellVelocity) const;

    const Mesh& mesh_;
    const Config& config_;
    const FlowBoundaryConditions& boundaryConditions_;
    double density_;
    double viscosity_;
    TurbulenceModelKind kind_;
    double relaxation_;
    double convectionBlending_;
    double minimumK_;
    double minimumEpsilon_;
    double minimumOmega_;
    double maximumEddyViscosityRatio_;

    ScalarField k_;
    ScalarField epsilon_;
    ScalarField omega_;
    ScalarField previousK_;
    ScalarField previousEpsilon_;
    ScalarField previousOmega_;
    ScalarField oldK_;
    ScalarField oldEpsilon_;
    ScalarField oldOmega_;
    ScalarField oldOldK_;
    ScalarField oldOldEpsilon_;
    ScalarField oldOldOmega_;
    ScalarField eddyViscosity_;
    ScalarField blendingF1_;
    ScalarField blendingF2_;
    std::vector<TurbulenceBoundaryCondition> kBoundaries_;
    std::vector<TurbulenceBoundaryCondition> epsilonBoundaries_;
    std::vector<TurbulenceBoundaryCondition> omegaBoundaries_;
    std::vector<double> wallDistance_;
    std::vector<double> production_;
    std::vector<double> strainInvariant_;
    std::vector<double> sigmaK_;
    std::vector<double> sigmaOmega_;
    std::vector<double> beta_;
    std::vector<double> gamma_;
    std::vector<char> constrainedCells_;
    std::vector<double> constrainedValues_;
    PetscLinearSystem kSystem_;
    PetscLinearSystem scaleSystem_;
    double change_ = 0.0;
    int linearIterations_ = 0;
    double linearResidual_ = 0.0;
};

std::unique_ptr<TurbulenceModel> makeTurbulenceModel(const Mesh& mesh, const Config& config,
                                                     const FlowBoundaryConditions& boundaryConditions,
                                                     double density, double viscosity);

} // namespace mp
