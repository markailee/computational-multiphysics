#pragma once

#include "flow/multiphase/MultiphaseModel.hpp"

#include <cstdint>

namespace mp {

struct ParticleParcel {
    Vector2 position;
    Vector2 velocity;
    double diameter = 0.0;
    double temperature = 300.0;
    double mass = 0.0;
    std::uint64_t particleCount = 1;
    std::size_t materialIndex = 0;
};

struct DiscretePhaseSources {
    std::vector<double> mass;
    std::vector<double> energy;
    std::vector<double> speciesMass;
    std::vector<Vector2> momentum;
};

class DpmModel final : public MultiphaseModel {
public:
    std::string name() const override { return "dpm"; }
    void inject(const ParticleParcel& parcel);
    const DiscretePhaseSources& continuousPhaseSources() const;
};

} // namespace mp
