#pragma once

#include "fvm/Config.hpp"
#include "fvm/Field.hpp"
#include "fvm/Gradient.hpp"
#include "material/MaterialModel.hpp"

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace mp {

class MultiSpeciesModel;

struct MultiphaseContext {
    const Mesh& mesh;
    VectorField& velocity;
    ScalarField& pressure;
    const std::vector<double>& massFlux;
    const std::vector<double>& density;
    const std::vector<double>& viscosity;
    const GradientScheme& gradient;
    MaterialModel& material;
    MultiSpeciesModel* species = nullptr;
    double dt = 0.0;
    double time = 0.0;
};

class MultiphaseModel {
public:
    virtual ~MultiphaseModel() = default;
    virtual std::string name() const = 0;

    virtual void initialize(MaterialModel& material);
    virtual void beginOuterIteration();
    virtual void solve(MultiphaseContext& context);
    virtual void updateTimeHistory();
    virtual double change() const;
    virtual void appendOutputFields(std::map<std::string, const ScalarField*>& scalars,
                                    std::map<std::string, const VectorField*>& vectors) const;

protected:
    [[noreturn]] void notImplemented() const;
};

std::unique_ptr<MultiphaseModel> makeMultiphaseModel(const Mesh& mesh, const Config& config);

} // namespace mp
