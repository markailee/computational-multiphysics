#pragma once

#include "material/MaterialModel.hpp"

namespace mp {

class IdealGasMixture final : public MaterialModel {
public:
    explicit IdealGasMixture(const Config& config);

    std::string name() const override { return "ideal-gas-mixture"; }
    MaterialKind kind() const override { return MaterialKind::IdealGasMixture; }
    MaterialProperties properties(const MaterialState& state) const override;
    const std::vector<SpeciesProperties>* speciesDefinitions() const override { return &species_; }

private:
    std::vector<SpeciesProperties> species_;
    std::vector<std::vector<double>> binaryDiffusivities_;
    bool mixtureAveragedDiffusion_ = false;
};

} // namespace mp
