#pragma once

#include "flow/multiphase/MultiphaseModel.hpp"

namespace mp {

class EulerEulerModel final : public MultiphaseModel {
public:
    std::string name() const override { return "euler-euler"; }
    const ScalarField& phaseFraction(std::size_t phase) const;
    const VectorField& phaseVelocity(std::size_t phase) const;
    const VectorField& interphaseMomentumTransfer(std::size_t phase) const;
};

} // namespace mp
