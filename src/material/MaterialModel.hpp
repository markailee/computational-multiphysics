#pragma once

#include "fvm/Config.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mp {

constexpr double universalGasConstant = 8.31446261815324;

enum class MaterialKind { Constant, IdealGas, IdealGasMixture, VolumeFractionMixture };

struct SpeciesProperties {
    std::string name;
    double molarMass = 0.0;
    double viscosity = 0.0;
    double thermalConductivity = 0.0;
    double specificHeatCapacity = 0.0;
    double diffusivity = 0.0;
};

struct MaterialState {
    double temperature = 300.0;
    double pressure = 101325.0;
    std::vector<double> massFractions;
    std::vector<double> volumeFractions;
};

struct MaterialProperties {
    double density = 0.0;
    double viscosity = 0.0;
    double thermalConductivity = 0.0;
    double specificHeatCapacity = 0.0;
    std::vector<double> speciesDiffusivities;
};

class MaterialModel {
public:
    virtual ~MaterialModel() = default;
    virtual std::string name() const = 0;
    virtual MaterialKind kind() const = 0;
    virtual MaterialProperties properties(const MaterialState& state) const = 0;
    virtual const std::vector<SpeciesProperties>* speciesDefinitions() const { return nullptr; }
};

class ConstantMaterial final : public MaterialModel {
public:
    ConstantMaterial(double density, double viscosity, double thermalConductivity = 0.0,
                     double specificHeatCapacity = 0.0);

    std::string name() const override { return "constant"; }
    MaterialKind kind() const override { return MaterialKind::Constant; }
    MaterialProperties properties(const MaterialState& state) const override;

private:
    MaterialProperties properties_;
};

class IdealGas final : public MaterialModel {
public:
    explicit IdealGas(const Config& config);

    std::string name() const override { return "ideal-gas"; }
    MaterialKind kind() const override { return MaterialKind::IdealGas; }
    MaterialProperties properties(const MaterialState& state) const override;
    double gasConstant() const { return gasConstant_; }

private:
    double gasConstant_ = 0.0;
    double viscosity_ = 0.0;
    double thermalConductivity_ = 0.0;
    double specificHeatCapacity_ = 0.0;
};

std::vector<double> massToMoleFractions(const std::vector<double>& massFractions,
                                        const std::vector<SpeciesProperties>& species);
double mixtureMolarMass(const std::vector<double>& massFractions,
                        const std::vector<SpeciesProperties>& species);

std::unique_ptr<MaterialModel> makeMaterialModel(const Config& config);

} // namespace mp
