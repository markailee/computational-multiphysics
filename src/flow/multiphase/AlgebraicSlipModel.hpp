#pragma once

#include "flow/multiphase/MultiphaseModel.hpp"

namespace mp {

class AlgebraicSlipModel final : public MultiphaseModel {
public:
    std::string name() const override { return "algebraic-slip"; }
    const ScalarField& phaseFraction(std::size_t phase) const;
    const VectorField& slipVelocity(std::size_t phase) const;
};

} // namespace mp
