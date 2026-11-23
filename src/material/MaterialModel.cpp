#include "material/MaterialModel.hpp"

#include "material/IdealGasMixture.hpp"
#include "material/MaterialParsing.hpp"
#include "material/VolumeFractionMixture.hpp"

#include <cmath>
#include <stdexcept>

namespace mp {

ConstantMaterial::ConstantMaterial(double density, double viscosity, double thermalConductivity,
                                   double specificHeatCapacity) {
    if (!(density > 0.0) || viscosity < 0.0 || thermalConductivity < 0.0 || specificHeatCapacity < 0.0) {
        throw std::runtime_error("Constant material requires positive density and non-negative transport properties");
    }
    properties_.density = density;
    properties_.viscosity = viscosity;
    properties_.thermalConductivity = thermalConductivity;
    properties_.specificHeatCapacity = specificHeatCapacity;
}

MaterialProperties ConstantMaterial::properties(const MaterialState&) const {
    return properties_;
}

IdealGas::IdealGas(const Config& config)
    : gasConstant_(config.number("gas_constant", 287.0)),
      viscosity_(config.number("viscosity", 0.01)),
      thermalConductivity_(config.number("thermal_conductivity", 0.0)),
      specificHeatCapacity_(config.number("specific_heat", 0.0)) {
    if (!(gasConstant_ > 0.0) || viscosity_ < 0.0 || thermalConductivity_ < 0.0 ||
        !(specificHeatCapacity_ > 0.0)) {
        throw std::runtime_error("Ideal gas requires positive gas_constant and specific_heat, and non-negative transport properties");
    }
}

MaterialProperties IdealGas::properties(const MaterialState& state) const {
    if (!(state.pressure > 0.0) || !(state.temperature > 0.0)) {
        throw std::runtime_error("Ideal-gas pressure and temperature must be positive");
    }
    MaterialProperties result;
    result.density = state.pressure / (gasConstant_ * state.temperature);
    result.viscosity = viscosity_;
    result.thermalConductivity = thermalConductivity_;
    result.specificHeatCapacity = specificHeatCapacity_;
    return result;
}

std::vector<double> massToMoleFractions(const std::vector<double>& massFractions,
                                        const std::vector<SpeciesProperties>& species) {
    if (massFractions.size() != species.size() || species.empty()) {
        throw std::runtime_error("Mass-fraction and species sizes are inconsistent");
    }
    std::vector<double> result(species.size());
    double sum = 0.0;
    for (std::size_t k = 0; k < species.size(); ++k) {
        if (!(species[k].molarMass > 0.0)) throw std::runtime_error("Species molar mass must be positive");
        result[k] = massFractions[k] / species[k].molarMass;
        sum += result[k];
    }
    if (!(sum > 0.0) || !std::isfinite(sum)) throw std::runtime_error("Invalid mass fractions for mole conversion");
    for (double& value : result) value /= sum;
    return result;
}

double mixtureMolarMass(const std::vector<double>& massFractions,
                        const std::vector<SpeciesProperties>& species) {
    const std::vector<double> moleFractions = massToMoleFractions(massFractions, species);
    double result = 0.0;
    for (std::size_t k = 0; k < species.size(); ++k) result += moleFractions[k] * species[k].molarMass;
    return result;
}

std::unique_ptr<MaterialModel> makeMaterialModel(const Config& config) {
    const std::string model = lowercase(config.get("material_model", "constant"));
    if (model == "constant" || model == "pure" || model == "none") {
        return std::make_unique<ConstantMaterial>(config.number("density", 1.0),
            config.number("viscosity", 0.01), config.number("thermal_conductivity", 0.0),
            config.number("specific_heat", 0.0));
    }
    if (model == "ideal_gas" && !config.has("species")) {
        return std::make_unique<IdealGas>(config);
    }
    if (model == "ideal_gas" || model == "ideal_gas_mixture" || model == "multispecies") {
        return std::make_unique<IdealGasMixture>(config);
    }
    if (model == "volume_fraction" || model == "volume_fraction_mixture" || model == "multiphase") {
        return std::make_unique<VolumeFractionMixture>(config);
    }
    throw std::runtime_error("Unknown material_model: " + model);
}

} // namespace mp
