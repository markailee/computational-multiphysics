#include "acoustic/AcousticBoundaryConditions.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

namespace mp {
namespace {

AcousticBoundaryCondition parseCondition(const std::string& patch, const std::string& text) {
    std::istringstream input(text);
    std::string name;
    input >> name;
    name = lowercase(name);
    AcousticBoundaryCondition result;
    if (name == "rigid" || name == "hardwall" || name == "wall" || name == "zerogradient") {
        result.kind = AcousticBoundaryKind::Rigid;
    } else if (name == "absorbing" || name == "nonreflecting" || name == "impedance") {
        result.kind = AcousticBoundaryKind::Absorbing;
    } else if (name == "fixedvalue") {
        result.kind = AcousticBoundaryKind::FixedValue;
        if (!(input >> result.value)) {
            throw std::runtime_error("Acoustic fixedValue requires a pressure on patch: " + patch);
        }
    } else if (name == "planewave") {
        result.kind = AcousticBoundaryKind::PlaneWave;
    } else {
        throw std::runtime_error("Unknown acoustic pressure boundary '" + name + "' on patch: " + patch);
    }
    std::string extra;
    if (input >> extra) throw std::runtime_error("Unexpected acoustic boundary data on patch: " + patch);
    return result;
}

} // namespace

AcousticBoundaryConditions::AcousticBoundaryConditions(const Mesh& mesh, const Config& config)
    : conditions_(mesh.faces.size()) {
    for (std::size_t patch = 0; patch < mesh.patches.size(); ++patch) {
        const std::string& name = mesh.patches[patch].name;
        const AcousticBoundaryCondition condition = parseCondition(name, config.require("boundary." + name + ".pressure"));
        for (int face : mesh.patches[patch].faces) conditions_[face] = condition;
    }
}

} // namespace mp
