#pragma once

#include "fvm/Config.hpp"
#include "fvm/Field.hpp"

namespace mp {

enum class VelocityBoundaryType { FixedValue, NoSlip, ZeroGradient, Slip, ParabolicX };
enum class PressureBoundaryType { FixedValue, ZeroGradient };

struct PatchBoundaryCondition {
    VelocityBoundaryType velocityType = VelocityBoundaryType::NoSlip;
    PressureBoundaryType pressureType = PressureBoundaryType::ZeroGradient;
    Vector2 velocity;
    double pressure = 0.0;
    double meanVelocity = 0.0;
    double minimumY = 0.0;
    double maximumY = 0.0;
};

class FlowBoundaryConditions {
public:
    FlowBoundaryConditions(const Mesh& mesh, const Config& config);

    Vector2 velocity(std::size_t faceIndex, const Vector2& cellVelocity) const;
    double pressure(std::size_t faceIndex, double cellPressure) const;
    bool fixedVelocity(std::size_t faceIndex) const;
    bool fixedPressure(std::size_t faceIndex) const;
    bool noSlip(std::size_t faceIndex) const;
    bool slip(std::size_t faceIndex) const;
    bool hasFixedPressure() const;

private:
    const Mesh* mesh_;
    std::vector<PatchBoundaryCondition> conditions_;
    bool hasFixedPressure_ = false;
};

} // namespace mp
