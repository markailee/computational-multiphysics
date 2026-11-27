#include "neutron/NeutronSolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>

namespace mp {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

void gaussLegendre(int count, std::vector<double>& nodes, std::vector<double>& weights) {
    nodes.assign(count, 0.0);
    weights.assign(count, 0.0);
    for (int i = 0; i < (count + 1) / 2; ++i) {
        double z = std::cos(pi * (i + 0.75) / (count + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double p0 = 1.0;
            double p1 = z;
            for (int degree = 2; degree <= count; ++degree) {
                const double p2 = ((2.0 * degree - 1.0) * z * p1 - (degree - 1.0) * p0) / degree;
                p0 = p1;
                p1 = p2;
            }
            derivative = count * (z * p1 - p0) / (z * z - 1.0);
            const double step = p1 / derivative;
            z -= step;
            if (std::abs(step) < 1.0e-15) break;
        }
        const double weight = 2.0 / ((1.0 - z * z) * derivative * derivative);
        nodes[i] = -z;
        nodes[count - 1 - i] = z;
        weights[i] = weight;
        weights[count - 1 - i] = weight;
    }
}

void addLevelSymmetricOctant(int order, std::vector<NeutronDirection>& octant) {
    octant.clear();
    struct LevelEntry {
        std::array<int, 3> idx;
        double weight;
    };
    std::vector<double> mu;
    std::vector<LevelEntry> entries;

    if (order == 2) {
        mu = {std::sqrt(1.0 / 3.0)};
        entries = {{{0, 0, 0}, 1.0}};
    } else if (order == 4) {
        const double u1 = (4.0 - std::sqrt(6.4)) / 12.0;
        const double u2 = 1.0 - 2.0 * u1;
        mu = {std::sqrt(u1), std::sqrt(u2)};
        entries = {{{0, 0, 1}, 1.0 / 3.0}};
    } else if (order == 6) {
        const double u1 = 0.2666355 * 0.2666355;
        const double d = (1.0 - 3.0 * u1) / 2.0;
        const double u2 = u1 + d;
        const double u3 = u1 + 2.0 * d;
        const double a1 = 2.0 * u1 * u1 + u3 * u3;
        const double a2 = u1 * u1 + 2.0 * u2 * u2;
        const double w1 = (0.2 - a2 / 3.0) / (a1 - a2);
        const double w2 = 1.0 / 3.0 - w1;
        mu = {std::sqrt(u1), std::sqrt(u2), std::sqrt(u3)};
        entries = {{{0, 0, 2}, w1}, {{0, 1, 1}, w2}};
    } else if (order == 8) {
        const double u1 = 1.0 / 21.0;
        const double d = (1.0 - 3.0 * u1) / 3.0;
        mu = {std::sqrt(u1), std::sqrt(u1 + d), std::sqrt(u1 + 2.0 * d), std::sqrt(u1 + 3.0 * d)};
        entries = {
            {{0, 0, 3}, 49.0 / 405.0},
            {{0, 1, 2}, 49.0 / 540.0},
            {{1, 1, 1}, 5.0 / 54.0},
        };
    } else {
        throw std::runtime_error("Supported level_symmetric sn_order values are 2, 4, 6, and 8 "
                                 "(use product_midpoint for arbitrarily high angular orders)");
    }

    for (const LevelEntry& entry : entries) {
        std::array<int, 3> p = entry.idx;
        std::sort(p.begin(), p.end());
        do {
            octant.push_back({mu[p[0]], mu[p[1]], mu[p[2]], entry.weight});
        } while (std::next_permutation(p.begin(), p.end()));
    }
}

} // namespace

void NeutronSolver::buildQuadrature() {
    directions_.clear();
    if (quadratureKind_ == NeutronQuadratureKind::LevelSymmetric) {
        std::vector<NeutronDirection> octant;
        addLevelSymmetricOctant(snOrder_, octant);
        // 2D xy symmetry: upper hemisphere (z > 0) with 4 quadrants, weights sum to 4*pi.
        constexpr std::array<std::pair<double, double>, 4> signs = {{{1.0, 1.0}, {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}}};
        for (const auto& [sx, sy] : signs) {
            for (const NeutronDirection& dir : octant) {
                directions_.push_back({sx * dir.x, sy * dir.y, dir.z, pi * dir.weight});
            }
        }
    } else {
        std::vector<double> nodes;
        std::vector<double> weights;
        gaussLegendre(polarCount_, nodes, weights);
        const bool midpoint = (quadratureKind_ == NeutronQuadratureKind::ProductMidpoint);
        // Fold symmetric +/-z into upper hemisphere (z > 0) for 2D xy transport: halves direction count exactly.
        for (int polar = 0; polar < polarCount_; ++polar) {
            if (nodes[polar] < 0.0) continue;
            const double transverse = std::sqrt(std::max(0.0, 1.0 - nodes[polar] * nodes[polar]));
            const double polarWeight = 2.0 * weights[polar];
            for (int azimuth = 0; azimuth < azimuthCount_; ++azimuth) {
                const double offset = midpoint ? 0.5 : 0.0;
                const double angle = 2.0 * pi * (azimuth + offset) / azimuthCount_;
                directions_.push_back({transverse * std::cos(angle),
                                       transverse * std::sin(angle),
                                       nodes[polar],
                                       polarWeight * 2.0 * pi / azimuthCount_});
            }
        }
    }

    directionCount_ = static_cast<int>(directions_.size());
    const double weightSum = std::accumulate(directions_.begin(), directions_.end(), 0.0,
        [](double sum, const NeutronDirection& direction) { return sum + direction.weight; });
    for (NeutronDirection& direction : directions_) {
        direction.weight *= 4.0 * pi / weightSum;
    }
}

void NeutronSolver::buildReflectionMaps() {
    if (model_ != NeutronModel::DiscreteOrdinates) return;

    reflectionDirections_.assign(mesh_.faces.size() * directionCount_, -1);
    faceIncomingWeightSum_.assign(mesh_.faces.size(), 0.0);

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const NeutronBoundaryCondition& condition = boundaries_[faceIndex];

        double incomingWeightSum = 0.0;
        for (int direction = 0; direction < directionCount_; ++direction) {
            const NeutronDirection& omega = directions_[direction];
            const double normalProjection = omega.x * face.normal.x + omega.y * face.normal.y;
            if (normalProjection < -1.0e-14) {
                incomingWeightSum += omega.weight * (-normalProjection);
            }
        }
        faceIncomingWeightSum_[faceIndex] = incomingWeightSum;

        const bool needsSpecular =
            (condition.kind == NeutronBoundaryKind::Reflective) ||
            (condition.kind == NeutronBoundaryKind::Albedo && condition.whiteFraction < 1.0 - 1.0e-12);
        if (!needsSpecular) continue;

        for (int direction = 0; direction < directionCount_; ++direction) {
            const NeutronDirection& omega = directions_[direction];
            const double normalProjection = omega.x * face.normal.x + omega.y * face.normal.y;
            if (normalProjection >= -1.0e-14) continue;
            const double reflectedX = omega.x - 2.0 * normalProjection * face.normal.x;
            const double reflectedY = omega.y - 2.0 * normalProjection * face.normal.y;
            int best = -1;
            double bestDistance = std::numeric_limits<double>::max();
            for (int candidate = 0; candidate < directionCount_; ++candidate) {
                const NeutronDirection& other = directions_[candidate];
                const double distance = (other.x - reflectedX) * (other.x - reflectedX) +
                    (other.y - reflectedY) * (other.y - reflectedY) + (other.z - omega.z) * (other.z - omega.z);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    best = candidate;
                }
            }
            if (best < 0 || bestDistance > 1.0e-18 ||
                std::abs(directions_[best].weight - omega.weight) > 1.0e-10 * omega.weight) {
                throw std::runtime_error(
                    "No exact specular reflected ordinate exists for this boundary orientation and quadrature "
                    "(use white_reflective or an x/y-aligned boundary)");
            }
            reflectionDirections_[reflectionIndex(static_cast<int>(faceIndex), direction)] = best;
        }
    }
}

void NeutronSolver::buildSnSweepOrder() {
    sweepOrder_.assign(directionCount_, {});
    directionAcyclic_.assign(directionCount_, false);

    for (int direction = 0; direction < directionCount_; ++direction) {
        const NeutronDirection& omega = directions_[direction];
        std::vector<std::vector<int>> outgoing(cellCount_);
        std::vector<int> inDegree(cellCount_, 0);

        for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
            const Face& face = mesh_.faces[faceIndex];
            const double proj = omega.x * face.normal.x + omega.y * face.normal.y;
            if (proj > 1.0e-14) {
                outgoing[face.owner].push_back(face.neighbour);
                ++inDegree[face.neighbour];
            } else if (proj < -1.0e-14) {
                outgoing[face.neighbour].push_back(face.owner);
                ++inDegree[face.owner];
            }
        }

        std::queue<int> ready;
        for (int cell = 0; cell < cellCount_; ++cell) {
            if (inDegree[cell] == 0) ready.push(cell);
        }

        std::vector<int> order;
        order.reserve(cellCount_);
        while (!ready.empty()) {
            const int current = ready.front();
            ready.pop();
            order.push_back(current);
            for (int next : outgoing[current]) {
                if (--inDegree[next] == 0) ready.push(next);
            }
        }

        if (static_cast<int>(order.size()) == cellCount_) {
            sweepOrder_[direction] = std::move(order);
            directionAcyclic_[direction] = true;
        } else {
            directionAcyclic_[direction] = false;
        }
    }
}

} // namespace mp
