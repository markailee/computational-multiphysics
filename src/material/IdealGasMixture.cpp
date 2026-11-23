#include "material/IdealGasMixture.hpp"

#include "material/MaterialParsing.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace mp {

IdealGasMixture::IdealGasMixture(const Config& config) {
    const std::vector<std::string> names = material::names(config, "species");
    for (const std::string& name : names) {
        const std::string prefix = "species." + name + ".";
        SpeciesProperties properties;
        properties.name = name;
        properties.molarMass = material::number(config, prefix + "molar_mass");
        properties.viscosity = material::number(config, prefix + "viscosity");
        properties.thermalConductivity = material::number(config, prefix + "thermal_conductivity");
        properties.specificHeatCapacity = material::number(config, prefix + "specific_heat");
        properties.diffusivity = material::number(config, prefix + "diffusivity");
        species_.push_back(properties);
    }

    const std::string diffusion = lowercase(config.get("diffusion_model", "constant"));
    if (diffusion == "constant" || diffusion == "fick") {
        mixtureAveragedDiffusion_ = false;
    } else if (diffusion == "mixture_averaged" || diffusion == "wilke") {
        mixtureAveragedDiffusion_ = true;
    } else {
        throw std::runtime_error("Unknown diffusion_model: " + diffusion);
    }

    const std::size_t count = species_.size();
    binaryDiffusivities_.assign(count, std::vector<double>(count, 0.0));
    for (std::size_t k = 0; k < count; ++k) binaryDiffusivities_[k][k] = species_[k].diffusivity;
    if (!mixtureAveragedDiffusion_ || count < 2) return;
    for (std::size_t k = 0; k < count; ++k) {
        for (std::size_t l = k + 1; l < count; ++l) {
            const std::string key = "binary_diffusion." + species_[k].name + "." + species_[l].name;
            const double value = material::number(config, key);
            binaryDiffusivities_[k][l] = value;
            binaryDiffusivities_[l][k] = value;
        }
    }
}

MaterialProperties IdealGasMixture::properties(const MaterialState& state) const {
    if (!(state.temperature > 0.0) || !(state.pressure > 0.0)) {
        throw std::runtime_error("Ideal-gas temperature and thermodynamic pressure must be positive");
    }
    const std::vector<double> massFractions = material::fractions(state.massFractions, species_.size(),
                                                                  "Ideal-gas composition");
    const std::vector<double> moleFractions = massToMoleFractions(massFractions, species_);
    const std::size_t count = species_.size();
    MaterialProperties result;
    result.specificHeatCapacity = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        result.specificHeatCapacity += massFractions[k] * species_[k].specificHeatCapacity;
    }
    result.density = state.pressure * mixtureMolarMass(massFractions, species_) /
        (universalGasConstant * state.temperature);

    double arithmeticConductivity = 0.0;
    double harmonicConductivity = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
        double denominator = 0.0;
        for (std::size_t l = 0; l < count; ++l) {
            const double viscosityRatio = species_[k].viscosity / species_[l].viscosity;
            const double molarRatio = species_[l].molarMass / species_[k].molarMass;
            const double numerator = 1.0 + std::sqrt(viscosityRatio) * std::pow(molarRatio, 0.25);
            const double phi = numerator * numerator /
                std::sqrt(8.0 * (1.0 + species_[k].molarMass / species_[l].molarMass));
            denominator += moleFractions[l] * phi;
        }
        if (!(denominator > 0.0)) throw std::runtime_error("Invalid Wilke viscosity denominator");
        result.viscosity += moleFractions[k] * species_[k].viscosity / denominator;
        arithmeticConductivity += moleFractions[k] * species_[k].thermalConductivity;
        harmonicConductivity += moleFractions[k] / species_[k].thermalConductivity;
    }
    result.thermalConductivity = 0.5 * (arithmeticConductivity + 1.0 / harmonicConductivity);

    result.speciesDiffusivities.assign(count, 0.0);
    for (std::size_t k = 0; k < count; ++k) {
        double denominator = 0.0;
        if (mixtureAveragedDiffusion_ && count > 1) {
            for (std::size_t l = 0; l < count; ++l) {
                if (l != k) denominator += moleFractions[l] / binaryDiffusivities_[k][l];
            }
        }
        const double numerator = 1.0 - moleFractions[k];
        result.speciesDiffusivities[k] = mixtureAveragedDiffusion_ && denominator > 0.0 && numerator > 1.0e-14
            ? numerator / denominator : species_[k].diffusivity;
        if (!(result.speciesDiffusivities[k] > 0.0) || !std::isfinite(result.speciesDiffusivities[k])) {
            throw std::runtime_error("Invalid effective diffusivity for species " + species_[k].name);
        }
    }
    return result;
}

} // namespace mp
