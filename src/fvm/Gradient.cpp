#include "fvm/Gradient.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace mp {
namespace {

template <typename Value>
struct GradientType;

template <>
struct GradientType<double> {
    using Type = ScalarGradient;
};

template <>
struct GradientType<Vector2> {
    using Type = VectorGradient;
};

template <typename Value>
using Gradient = typename GradientType<Value>::Type;

ScalarGradient contribution(double value, const Vector2& area) {
    return value * area;
}

VectorGradient contribution(const Vector2& value, const Vector2& area) {
    return {value.x * area, value.y * area};
}

void accumulate(ScalarGradient& target, const ScalarGradient& value) {
    target += value;
}

void accumulate(VectorGradient& target, const VectorGradient& value) {
    target[0] += value[0];
    target[1] += value[1];
}

void divide(ScalarGradient& value, double divisor) {
    value /= divisor;
}

void divide(VectorGradient& value, double divisor) {
    value[0] /= divisor;
    value[1] /= divisor;
}

double reconstruct(double value, const ScalarGradient& gradient, const Vector2& displacement) {
    return value + dot(gradient, displacement);
}

Vector2 reconstruct(const Vector2& value, const VectorGradient& gradient, const Vector2& displacement) {
    return {value.x + dot(gradient[0], displacement), value.y + dot(gradient[1], displacement)};
}

template <typename Value>
Gradient<Value> cellGauss(const Mesh& mesh, const Field<Value>& field, int cellIndex,
                          const std::vector<Gradient<Value>>* reconstruction) {
    const Cell& cell = mesh.cells[cellIndex];
    Gradient<Value> result{};
    for (std::size_t i = 0; i < cell.faces.size(); ++i) {
        const int faceIndex = cell.faces[i];
        const Face& face = mesh.faces[faceIndex];
        Value value;
        if (face.isBoundary()) {
            value = field.boundary(faceIndex);
            if (reconstruction && !field.boundaryDirichlet(faceIndex)) {
                value = reconstruct(field[cellIndex], (*reconstruction)[cellIndex], face.centre - cell.centre);
            }
        } else if (reconstruction) {
            const int neighbour = face.owner == cellIndex ? face.neighbour : face.owner;
            const Value owner = reconstruct(field[face.owner], (*reconstruction)[face.owner],
                                            face.centre - mesh.cells[face.owner].centre);
            const Value other = reconstruct(field[face.neighbour], (*reconstruction)[face.neighbour],
                                            face.centre - mesh.cells[face.neighbour].centre);
            value = face.ownerWeight * owner + (1.0 - face.ownerWeight) * other;
            (void)neighbour;
        } else {
            value = face.ownerWeight * field[face.owner] + (1.0 - face.ownerWeight) * field[face.neighbour];
        }
        const Vector2 area = face.normal * (face.area * cell.faceSigns[i]);
        accumulate(result, contribution(value, area));
    }
    divide(result, cell.volume);
    return result;
}

template <typename Value>
std::vector<Gradient<Value>> gauss(const Field<Value>& field, int correctionIterations) {
    const Mesh& mesh = field.mesh();
    std::vector<Gradient<Value>> gradient(mesh.cells.size());
    for (std::size_t cell = 0; cell < mesh.cells.size(); ++cell) {
        gradient[cell] = cellGauss(mesh, field, static_cast<int>(cell), nullptr);
    }
    for (int iteration = 0; iteration < correctionIterations; ++iteration) {
        std::vector<Gradient<Value>> corrected(mesh.cells.size());
        for (std::size_t cell = 0; cell < mesh.cells.size(); ++cell) {
            corrected[cell] = cellGauss(mesh, field, static_cast<int>(cell), &gradient);
        }
        gradient = std::move(corrected);
    }
    return gradient;
}

template <typename Value>
void addNormalEquation(double weight, const Vector2& distance, const Value& difference,
                       double (&normalMatrix)[2][2], Gradient<Value>& rightHandSide) {
    normalMatrix[0][0] += weight * distance.x * distance.x;
    normalMatrix[0][1] += weight * distance.x * distance.y;
    normalMatrix[1][0] += weight * distance.y * distance.x;
    normalMatrix[1][1] += weight * distance.y * distance.y;
    if constexpr (std::is_same_v<Value, double>) {
        rightHandSide += distance * (weight * difference);
    } else {
        rightHandSide[0] += distance * (weight * difference.x);
        rightHandSide[1] += distance * (weight * difference.y);
    }
}

template <typename Value>
std::vector<Gradient<Value>> leastSquares(const Field<Value>& field) {
    const Mesh& mesh = field.mesh();
    std::vector<Gradient<Value>> gradients(mesh.cells.size());
    for (std::size_t cellIndex = 0; cellIndex < mesh.cells.size(); ++cellIndex) {
        const Cell& cell = mesh.cells[cellIndex];
        double matrix[2][2] = {};
        Gradient<Value> rhs{};
        for (std::size_t i = 0; i < cell.faces.size(); ++i) {
            const int faceIndex = cell.faces[i];
            const Face& face = mesh.faces[faceIndex];
            if (face.isBoundary()) {
                if (!field.boundaryDirichlet(faceIndex)) continue;
                const Vector2 distance = face.centre - cell.centre;
                const Value difference = field.boundary(faceIndex) - field[cellIndex];
                addNormalEquation(1.0 / normSquared(distance), distance, difference, matrix, rhs);
            } else {
                const int neighbour = face.owner == static_cast<int>(cellIndex) ? face.neighbour : face.owner;
                const Vector2 distance = mesh.cells[neighbour].centre - cell.centre;
                const Value difference = field[neighbour] - field[cellIndex];
                addNormalEquation(1.0 / normSquared(distance), distance, difference, matrix, rhs);
            }
        }
        const double determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        const double scale = matrix[0][0] * matrix[1][1] + matrix[0][1] * matrix[1][0];
        if (std::abs(determinant) <= 1.0e-12 * scale || std::abs(determinant) < 1.0e-28) {
            gradients[cellIndex] = cellGauss(mesh, field, static_cast<int>(cellIndex), nullptr);
            continue;
        }
        if constexpr (std::is_same_v<Value, double>) {
            gradients[cellIndex] = {(matrix[1][1] * rhs.x - matrix[0][1] * rhs.y) / determinant,
                                    (matrix[0][0] * rhs.y - matrix[1][0] * rhs.x) / determinant};
        } else {
            for (int component = 0; component < 2; ++component) {
                gradients[cellIndex][component] = {
                    (matrix[1][1] * rhs[component].x - matrix[0][1] * rhs[component].y) / determinant,
                    (matrix[0][0] * rhs[component].y - matrix[1][0] * rhs[component].x) / determinant};
            }
        }
    }
    return gradients;
}

template <typename Value>
void limitReconstructedGradient(const Field<Value>& field, std::vector<Gradient<Value>>& gradients) {
    const Mesh& mesh = field.mesh();
    constexpr int componentCount = std::is_same_v<Value, double> ? 1 : 2;
    for (std::size_t cellIndex = 0; cellIndex < mesh.cells.size(); ++cellIndex) {
        const Cell& cell = mesh.cells[cellIndex];
        for (int componentIndex = 0; componentIndex < componentCount; ++componentIndex) {
            const auto component = [componentIndex](const Value& value) {
                if constexpr (std::is_same_v<Value, double>) return value;
                else return value[componentIndex];
            };
            const double centreValue = component(field[cellIndex]);
            double minimum = centreValue;
            double maximum = centreValue;
            for (int faceIndex : cell.faces) {
                const Face& face = mesh.faces[faceIndex];
                const double neighbourValue = face.isBoundary()
                    ? component(field.boundary(faceIndex))
                    : component(field[face.owner == static_cast<int>(cellIndex) ? face.neighbour : face.owner]);
                minimum = std::min(minimum, neighbourValue);
                maximum = std::max(maximum, neighbourValue);
            }

            double alpha = 1.0;
            const double tolerance = std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(centreValue), std::abs(minimum), std::abs(maximum)});
            for (int faceIndex : cell.faces) {
                const Vector2 displacement = mesh.faces[faceIndex].centre - cell.centre;
                const double delta = [&] {
                    if constexpr (std::is_same_v<Value, double>) return dot(gradients[cellIndex], displacement);
                    else return dot(gradients[cellIndex][componentIndex], displacement);
                }();
                if (delta > tolerance) {
                    alpha = std::min(alpha, std::max(0.0, (maximum - centreValue) / delta));
                } else if (delta < -tolerance) {
                    alpha = std::min(alpha, std::max(0.0, (minimum - centreValue) / delta));
                }
            }
            alpha = std::clamp(alpha, 0.0, 1.0);
            if constexpr (std::is_same_v<Value, double>) gradients[cellIndex] *= alpha;
            else gradients[cellIndex][componentIndex] *= alpha;
        }
    }
}

} // namespace

bool useBarthJespersenLimiter(const Config& config) {
    const std::string limiter = lowercase(config.get("gradient_limiter", "none"));
    if (limiter == "none" || limiter == "off" || limiter == "false") return false;
    if (limiter == "barth_jespersen" || limiter == "barth-jespersen" || limiter == "bj") return true;
    throw std::runtime_error("Unknown gradient limiter: " + limiter);
}

std::vector<ScalarGradient> applyBarthJespersenLimiter(const ScalarField& field,
                                                       std::vector<ScalarGradient> gradients) {
    limitReconstructedGradient(field, gradients);
    return gradients;
}

std::vector<VectorGradient> applyBarthJespersenLimiter(const VectorField& field,
                                                       std::vector<VectorGradient> gradients) {
    limitReconstructedGradient(field, gradients);
    return gradients;
}

std::vector<ScalarGradient> GaussGradient::scalar(const ScalarField& field) const {
    return gauss(field, correctionIterations_);
}

std::vector<VectorGradient> GaussGradient::vector(const VectorField& field) const {
    return gauss(field, correctionIterations_);
}

std::vector<ScalarGradient> LeastSquaresGradient::scalar(const ScalarField& field) const {
    return leastSquares(field);
}

std::vector<VectorGradient> LeastSquaresGradient::vector(const VectorField& field) const {
    return leastSquares(field);
}

std::unique_ptr<GradientScheme> makeGradientScheme(const std::string& name) {
    std::string value = name;
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::tolower(c); });
    if (value == "gauss" || value == "green_gauss" || value == "green-gauss") {
        return std::make_unique<GaussGradient>();
    }
    if (value == "least_squares" || value == "least-squares" || value == "lsq") {
        return std::make_unique<LeastSquaresGradient>();
    }
    throw std::runtime_error("Unknown gradient scheme: " + name);
}

} // namespace mp
