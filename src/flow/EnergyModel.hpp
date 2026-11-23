#pragma once

#include "fvm/Config.hpp"
#include "fvm/Field.hpp"
#include "fvm/Gradient.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/TimeDerivative.hpp"

#include <vector>

namespace mp {

enum class EnergyBoundaryType { FixedValue, ZeroGradient, FixedHeatFlux };

struct EnergyBoundaryCondition {
    EnergyBoundaryType type = EnergyBoundaryType::ZeroGradient;
    double value = 300.0;
};

class EnergyModel final {
public:
    EnergyModel(const Mesh& mesh, const Config& config, bool compressible);

    void initialize();
    void beginIteration();
    void updateBoundaries();
    LinearSolveInfo solve(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                          const GradientScheme& gradient, const std::vector<double>& density,
                          const std::vector<double>& oldDensity, const std::vector<double>& oldOldDensity,
                          const std::vector<double>& viscosity, const std::vector<double>& conductivity,
                          const std::vector<double>& heatCapacity, const ScalarField* turbulentViscosity,
                          const VectorField& velocity, const ScalarField& pressure,
                          const ScalarField& oldPressure, const ScalarField& oldOldPressure,
                          const Vector2& bodyForce);
    void updateTimeHistory();

    const ScalarField& temperature() const { return temperature_; }
    double change() const { return change_; }
    int linearIterations() const { return linearIterations_; }
    double linearResidual() const { return linearResidual_; }

private:
    EnergyBoundaryCondition readBoundary(const std::string& patch) const;

    const Mesh& mesh_;
    const Config& config_;
    bool compressible_ = false;
    double relaxation_;
    double convectionBlending_;
    double turbulentPrandtl_;
    double heatSource_;
    bool viscousDissipation_;
    double linearRelativeTolerance_;
    double linearAbsoluteTolerance_;
    int maximumLinearIterations_;
    TimeDerivative timeDerivative_;

    ScalarField temperature_;
    ScalarField previousTemperature_;
    ScalarField oldTemperature_;
    ScalarField oldOldTemperature_;
    std::vector<EnergyBoundaryCondition> boundaries_;
    std::vector<double> diffusivity_;
    std::vector<double> source_;
    std::vector<double> destruction_;
    std::vector<char> constrainedCells_;
    std::vector<double> constrainedValues_;
    PetscLinearSystem energySystem_;

    double change_ = 0.0;
    int linearIterations_ = 0;
    double linearResidual_ = 0.0;
    bool initialized_ = false;
};

} // namespace mp
