#include "material/VolumeFractionMixture.hpp"

#include "material/MaterialParsing.hpp"

#include <stdexcept>

namespace mp {

VolumeFractionMixture::VolumeFractionMixture(const Config& config) {
    phaseNames_ = material::names(config, "phases");
    for (const std::string& name : phaseNames_) {
        const std::string prefix = "phase." + name + ".";
        MaterialProperties properties;
        properties.density = material::number(config, prefix + "density");
        properties.viscosity = material::number(config, prefix + "viscosity");
        properties.thermalConductivity = material::number(config, prefix + "thermal_conductivity", false);
        properties.specificHeatCapacity = material::number(config, prefix + "specific_heat", false);
        phases_.push_back(properties);
    }
}

MaterialProperties VolumeFractionMixture::properties(const MaterialState& state) const {
    const std::vector<double> alpha = material::fractions(state.volumeFractions, phases_.size(),
                                                          "Volume-fraction composition");
    MaterialProperties result;
    for (std::size_t k = 0; k < phases_.size(); ++k) {
        result.density += alpha[k] * phases_[k].density;
        result.viscosity += alpha[k] * phases_[k].viscosity;
        result.thermalConductivity += alpha[k] * phases_[k].thermalConductivity;
        result.specificHeatCapacity += alpha[k] * phases_[k].density * phases_[k].specificHeatCapacity;
    }
    if (!(result.density > 0.0)) throw std::runtime_error("Volume-fraction mixture density must be positive");
    result.specificHeatCapacity /= result.density;
    return result;
}

} // namespace mp
