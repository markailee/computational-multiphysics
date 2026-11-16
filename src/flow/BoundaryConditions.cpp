#include "flow/BoundaryConditions.hpp"

#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace mp {

FlowBoundaryConditions::FlowBoundaryConditions(const Mesh& mesh, const Config& config) : mesh_(&mesh) {
    for (std::size_t patchIndex = 0; patchIndex < mesh.patches.size(); ++patchIndex) {
        const BoundaryPatch& patch = mesh.patches[patchIndex];
        const std::string prefix = "boundary." + patch.name + ".";
        std::istringstream velocity(config.require(prefix + "velocity"));
        std::string velocityType;
        velocity >> velocityType;
        std::istringstream pressure(config.require(prefix + "pressure"));
        std::string pressureType;
        pressure >> pressureType;

        PatchBoundaryCondition condition;
        const std::string type = lowercase(velocityType);
        if (type == "fixedvalue") {
            condition.velocityType = VelocityBoundaryType::FixedValue;
            velocity >> condition.velocity;
        } else if (type == "noslip") {
            condition.velocityType = VelocityBoundaryType::NoSlip;
        } else if (type == "zerogradient") {
            condition.velocityType = VelocityBoundaryType::ZeroGradient;
        } else if (type == "slip" || type == "symmetry") {
            condition.velocityType = VelocityBoundaryType::Slip;
        } else if (type == "parabolic" || type == "parabolicx") {
            condition.velocityType = VelocityBoundaryType::ParabolicX;
            velocity >> condition.meanVelocity;
        } else {
            throw std::runtime_error("Unknown velocity boundary condition on patch " + patch.name);
        }
        if (!velocity && type != "noslip" && type != "zerogradient" && type != "slip" && type != "symmetry") {
            throw std::runtime_error("Invalid velocity data on patch " + patch.name);
        }

        const std::string pressureName = lowercase(pressureType);
        if (pressureName == "fixedvalue") {
            condition.pressureType = PressureBoundaryType::FixedValue;
            pressure >> condition.pressure;
            hasFixedPressure_ = true;
        } else if (pressureName == "zerogradient") {
            condition.pressureType = PressureBoundaryType::ZeroGradient;
        } else {
            throw std::runtime_error("Unknown pressure boundary condition on patch " + patch.name);
        }
        if (!pressure) throw std::runtime_error("Invalid pressure data on patch " + patch.name);

        condition.minimumY = std::numeric_limits<double>::max();
        condition.maximumY = std::numeric_limits<double>::lowest();
        for (int face : patch.faces) {
            condition.minimumY = std::min({condition.minimumY, mesh.nodes[mesh.faces[face].node0].y,
                                           mesh.nodes[mesh.faces[face].node1].y});
            condition.maximumY = std::max({condition.maximumY, mesh.nodes[mesh.faces[face].node0].y,
                                           mesh.nodes[mesh.faces[face].node1].y});
        }
        if (condition.velocityType == VelocityBoundaryType::ParabolicX &&
            condition.maximumY - condition.minimumY <= 0.0) {
            throw std::runtime_error("Parabolic inlet requires a vertical patch: " + patch.name);
        }
        conditions_.push_back(condition);
    }
}

Vector2 FlowBoundaryConditions::velocity(std::size_t faceIndex, const Vector2& cellVelocity) const {
    const Face& face = mesh_->faces[faceIndex];
    const PatchBoundaryCondition& condition = conditions_[face.patch];
    switch (condition.velocityType) {
    case VelocityBoundaryType::FixedValue:
        return condition.velocity;
    case VelocityBoundaryType::NoSlip:
        return {};
    case VelocityBoundaryType::ZeroGradient:
        return cellVelocity;
    case VelocityBoundaryType::Slip:
        return cellVelocity - face.normal * dot(cellVelocity, face.normal);
    case VelocityBoundaryType::ParabolicX: {
        const double eta = (face.centre.y - condition.minimumY) / (condition.maximumY - condition.minimumY);
        return {6.0 * condition.meanVelocity * eta * (1.0 - eta), 0.0};
    }
    }
    return {};
}

double FlowBoundaryConditions::pressure(std::size_t faceIndex, double cellPressure) const {
    const Face& face = mesh_->faces[faceIndex];
    const PatchBoundaryCondition& condition = conditions_[face.patch];
    return condition.pressureType == PressureBoundaryType::FixedValue ? condition.pressure : cellPressure;
}

bool FlowBoundaryConditions::fixedVelocity(std::size_t faceIndex) const {
    const auto type = conditions_[mesh_->faces[faceIndex].patch].velocityType;
    return type == VelocityBoundaryType::FixedValue || type == VelocityBoundaryType::NoSlip ||
           type == VelocityBoundaryType::ParabolicX;
}

bool FlowBoundaryConditions::fixedPressure(std::size_t faceIndex) const {
    return conditions_[mesh_->faces[faceIndex].patch].pressureType == PressureBoundaryType::FixedValue;
}

bool FlowBoundaryConditions::slip(std::size_t faceIndex) const {
    return conditions_[mesh_->faces[faceIndex].patch].velocityType == VelocityBoundaryType::Slip;
}

bool FlowBoundaryConditions::hasFixedPressure() const {
    return hasFixedPressure_;
}

} // namespace mp
