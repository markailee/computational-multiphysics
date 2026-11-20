#pragma once

#include "fvm/Config.hpp"
#include "mesh/Mesh.hpp"

#include <vector>

namespace mp {

enum class AcousticBoundaryKind {
    Rigid,
    Absorbing,
    FixedValue,
    PlaneWave,
};

struct AcousticBoundaryCondition {
    AcousticBoundaryKind kind = AcousticBoundaryKind::Rigid;
    double value = 0.0;
};

class AcousticBoundaryConditions {
public:
    AcousticBoundaryConditions(const Mesh& mesh, const Config& config);

    const AcousticBoundaryCondition& face(int faceIndex) const { return conditions_[faceIndex]; }
    static bool isDirichlet(AcousticBoundaryKind kind) {
        return kind == AcousticBoundaryKind::FixedValue || kind == AcousticBoundaryKind::PlaneWave;
    }

private:
    std::vector<AcousticBoundaryCondition> conditions_;
};

} // namespace mp
