#pragma once

#include "material/MaterialModel.hpp"

namespace mp {

class VolumeFractionMixture final : public MaterialModel {
public:
    explicit VolumeFractionMixture(const Config& config);

    std::string name() const override { return "volume-fraction-mixture"; }
    MaterialKind kind() const override { return MaterialKind::VolumeFractionMixture; }
    MaterialProperties properties(const MaterialState& state) const override;

private:
    std::vector<std::string> phaseNames_;
    std::vector<MaterialProperties> phases_;
};

} // namespace mp
