#include "fvm/Operators.hpp"

#include <algorithm>
#include <stdexcept>

namespace mp {

Vector2 faceVector(const Mesh& mesh, std::size_t faceIndex) {
    const Face& face = mesh.faces[faceIndex];
    return face.normal * face.area;
}

Vector2 nonOrthogonalVector(const Mesh& mesh, std::size_t faceIndex) {
    const Face& face = mesh.faces[faceIndex];
    const Vector2 centreDistance = face.isBoundary()
        ? face.centre - mesh.cells[face.owner].centre
        : mesh.cells[face.neighbour].centre - mesh.cells[face.owner].centre;
    const Vector2 orthogonal = centreDistance * (face.area / face.delta);
    return faceVector(mesh, faceIndex) - orthogonal;
}

ConvectionCoefficients upwindConvection(double massFlux) {
    ConvectionCoefficients coefficients;
    coefficients.ownerDiagonal = std::max(massFlux, 0.0);
    coefficients.neighbourDiagonal = std::max(-massFlux, 0.0);
    coefficients.ownerToNeighbour = std::min(massFlux, 0.0);
    coefficients.neighbourToOwner = -std::max(massFlux, 0.0);
    return coefficients;
}

double diffusionConductance(const Mesh& mesh, std::size_t faceIndex, double diffusivity) {
    return diffusivity * mesh.faces[faceIndex].area / mesh.faces[faceIndex].delta;
}

double rhieChowCoefficient(const Mesh& mesh, std::size_t faceIndex,
                           const std::vector<double>& momentumDiagonal, double density) {
    const Face& face = mesh.faces[faceIndex];
    const double owner = mesh.cells[face.owner].volume / momentumDiagonal[face.owner];
    const double interpolated = face.isBoundary()
        ? owner
        : face.ownerWeight * owner + (1.0 - face.ownerWeight) * mesh.cells[face.neighbour].volume / momentumDiagonal[face.neighbour];
    return density * interpolated * face.area / face.delta;
}

} // namespace mp
