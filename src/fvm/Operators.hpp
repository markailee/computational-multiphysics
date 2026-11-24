#pragma once

#include "fvm/Field.hpp"

namespace mp {

template <typename Value>
Value interpolate(const Mesh& mesh, const Field<Value>& field, std::size_t faceIndex) {
    const Face& face = mesh.faces[faceIndex];
    if (face.isBoundary()) {
        return field.boundary(faceIndex);
    }
    return face.ownerWeight * field[face.owner] + (1.0 - face.ownerWeight) * field[face.neighbour];
}

template <typename Gradient>
Gradient interpolateGradient(const Mesh& mesh, const std::vector<Gradient>& gradient,
                             std::size_t faceIndex) {
    const Face& face = mesh.faces[faceIndex];
    if (face.isBoundary()) {
        return gradient[face.owner];
    }
    return face.ownerWeight * gradient[face.owner] + (1.0 - face.ownerWeight) * gradient[face.neighbour];
}

Vector2 faceVector(const Mesh& mesh, std::size_t faceIndex);
Vector2 nonOrthogonalVector(const Mesh& mesh, std::size_t faceIndex);

inline double reconstruct(double value, const Vector2& gradient, const Vector2& displacement) {
    return value + dot(gradient, displacement);
}

inline Vector2 reconstruct(const Vector2& value, const VectorGradient& gradient,
                           const Vector2& displacement) {
    return {value.x + dot(gradient[0], displacement), value.y + dot(gradient[1], displacement)};
}

template <typename Value, typename Gradient>
Value secondOrderUpwindFaceValue(const Mesh& mesh, const Field<Value>& field,
                                 const std::vector<Gradient>& gradients, std::size_t faceIndex,
                                 double massFlux) {
    const Face& face = mesh.faces[faceIndex];
    if (face.isBoundary() && field.boundaryDirichlet(faceIndex)) return field.boundary(faceIndex);
    const int upwindCell = !face.isBoundary() && massFlux < 0.0 ? face.neighbour : face.owner;
    return reconstruct(field[upwindCell], gradients[upwindCell],
                       face.centre - mesh.cells[upwindCell].centre);
}

struct ConvectionCoefficients {
    double ownerDiagonal = 0.0;
    double neighbourDiagonal = 0.0;
    double ownerToNeighbour = 0.0;
    double neighbourToOwner = 0.0;
};

ConvectionCoefficients upwindConvection(double massFlux);
double diffusionConductance(const Mesh& mesh, std::size_t faceIndex, double diffusivity);
double rhieChowCoefficient(const Mesh& mesh, std::size_t faceIndex,
                           const std::vector<double>& momentumDiagonal, double density);

} // namespace mp
