#pragma once

#include "flow/multiphase/MultiphaseModel.hpp"

namespace mp {

class VofModel final : public MultiphaseModel {
public:
    std::string name() const override { return "vof"; }
    const ScalarField& volumeFraction(std::size_t phase) const;
    const VectorField& surfaceTensionForce() const;
};

} // namespace mp
