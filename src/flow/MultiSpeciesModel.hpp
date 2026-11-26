#pragma once

#include "fvm/Config.hpp"
#include "fvm/Field.hpp"
#include "fvm/Gradient.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/TimeDerivative.hpp"
#include "material/MaterialModel.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace mp {

struct SpeciesBoundaryCondition {
    bool fixedValue = false;
    std::vector<double> massFractions;
};

class MultiSpeciesModel final {
public:
    MultiSpeciesModel(const Mesh& mesh, const Config& config, std::shared_ptr<MaterialModel> material);

    void initialize();
    void updateBoundaries();
    void beginIteration();
    void solve(double dt, bool hasOldOld, const std::vector<double>& massFlux,
               const GradientScheme& gradient, double carrierDensity,
               const ScalarField* turbulentViscosity);
    void updateTimeHistory();
    void appendOutputFields(std::map<std::string, const ScalarField*>& scalars) const;

    std::size_t count() const { return fields_.size(); }
    const ScalarField& massFraction(std::size_t species) const { return fields_.at(species); }
    double change() const { return change_; }
    double compositionError() const { return compositionError_; }
    double minimumMassFraction() const { return minimumMassFraction_; }
    double maximumMassFraction() const { return maximumMassFraction_; }
    bool bounded(double tolerance) const;
    int linearIterations() const { return linearIterations_; }
    double linearResidual() const { return linearResidual_; }
    const std::vector<double>& linearResiduals() const { return linearResiduals_; }

private:
    SpeciesBoundaryCondition readBoundary(const std::string& patch) const;
    void updateFieldBoundaries(std::size_t species);
    void updateBackgroundSpecies();
    void updateProperties();
    void computeTransportDiffusivities(double carrierDensity, const ScalarField* turbulentViscosity);
    void computeCorrectionFluxes(const GradientScheme& gradient);
    void updateDiagnostics();

    const Mesh& mesh_;
    const Config& config_;
    std::shared_ptr<MaterialModel> material_;
    std::size_t speciesCount_ = 0;
    double temperature_;
    double thermodynamicPressure_;
    double relaxation_;
    double turbulentSchmidt_;
    double linearRelativeTolerance_;
    double linearAbsoluteTolerance_;
    int maximumLinearIterations_;
    TimeDerivative timeDerivative_;

    std::vector<ScalarField> fields_;
    std::vector<ScalarField> previousFields_;
    std::vector<ScalarField> oldFields_;
    std::vector<ScalarField> oldOldFields_;
    std::vector<ScalarField> molecularDiffusivityFields_;
    std::vector<SpeciesBoundaryCondition> boundaries_;
    ScalarField densityField_;
    ScalarField viscosityField_;
    ScalarField conductivityField_;
    ScalarField heatCapacityField_;

    std::vector<std::vector<double>> diffusivities_;
    std::vector<std::vector<double>> correctionFluxes_;
    std::vector<std::vector<ScalarGradient>> gradients_;
    std::vector<double> source_;
    std::vector<double> destruction_;
    std::vector<char> constrainedCells_;
    std::vector<double> constrainedValues_;
    PetscLinearSystem transportSystem_;

    double change_ = 0.0;
    double compositionError_ = 0.0;
    double minimumMassFraction_ = 0.0;
    double maximumMassFraction_ = 1.0;
    int linearIterations_ = 0;
    double linearResidual_ = 0.0;
    std::vector<double> linearResiduals_;
    bool initialized_ = false;
};

std::unique_ptr<MultiSpeciesModel> makeMultiSpeciesModel(const Mesh& mesh, const Config& config,
                                                         std::shared_ptr<MaterialModel> material);

} // namespace mp
