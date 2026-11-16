#include "mesh/Mesh.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace mp {
namespace {

struct EdgeKey {
    int first;
    int second;

    EdgeKey(int a, int b) : first(std::min(a, b)), second(std::max(a, b)) {}
    bool operator<(const EdgeKey& other) const {
        return std::tie(first, second) < std::tie(other.first, other.second);
    }
};

void orderCounterClockwise(std::vector<int>& nodes, const std::vector<Vector2>& coordinates) {
    double twiceArea = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const Vector2& a = coordinates[nodes[i]];
        const Vector2& b = coordinates[nodes[(i + 1) % nodes.size()]];
        twiceArea += cross(a, b);
    }
    if (twiceArea < 0.0) {
        std::reverse(nodes.begin(), nodes.end());
    }
}

void cellGeometry(Cell& cell, const std::vector<Vector2>& coordinates) {
    double twiceArea = 0.0;
    Vector2 centroid;
    for (std::size_t i = 0; i < cell.nodes.size(); ++i) {
        const Vector2& a = coordinates[cell.nodes[i]];
        const Vector2& b = coordinates[cell.nodes[(i + 1) % cell.nodes.size()]];
        const double contribution = cross(a, b);
        twiceArea += contribution;
        centroid += (a + b) * contribution;
    }
    if (twiceArea <= 1.0e-14) {
        throw std::runtime_error("Degenerate or non-convex cell in mesh");
    }
    cell.volume = 0.5 * twiceArea;
    cell.centre = centroid / (3.0 * twiceArea);
}

} // namespace

void Mesh::buildGeometry(const BoundaryFaceMap& boundaryFacePatches) {
    if (nodes.empty() || cells.empty()) {
        throw std::runtime_error("Mesh contains no two-dimensional cells");
    }

    for (Cell& cell : cells) {
        orderCounterClockwise(cell.nodes, nodes);
        cellGeometry(cell, nodes);
        cell.faces.clear();
        cell.faceSigns.clear();
    }

    std::map<EdgeKey, int> edgeToFace;
    faces.clear();
    for (int cellIndex = 0; cellIndex < static_cast<int>(cells.size()); ++cellIndex) {
        Cell& cell = cells[cellIndex];
        for (std::size_t edge = 0; edge < cell.nodes.size(); ++edge) {
            const int a = cell.nodes[edge];
            const int b = cell.nodes[(edge + 1) % cell.nodes.size()];
            const EdgeKey key(a, b);
            auto [iterator, inserted] = edgeToFace.emplace(key, static_cast<int>(faces.size()));
            if (inserted) {
                Face face;
                face.node0 = a;
                face.node1 = b;
                face.owner = cellIndex;
                faces.push_back(face);
                cell.faces.push_back(iterator->second);
                cell.faceSigns.push_back(1);
            } else {
                Face& face = faces[iterator->second];
                if (face.neighbour >= 0) {
                    throw std::runtime_error("An edge is shared by more than two cells");
                }
                face.neighbour = cellIndex;
                cell.faces.push_back(iterator->second);
                cell.faceSigns.push_back(-1);
            }
        }
    }

    std::vector<int> order(faces.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [this](int a, int b) {
        return faces[a].isBoundary() < faces[b].isBoundary();
    });
    std::vector<int> oldToNew(faces.size());
    std::vector<Face> orderedFaces;
    orderedFaces.reserve(faces.size());
    for (int newIndex = 0; newIndex < static_cast<int>(order.size()); ++newIndex) {
        oldToNew[order[newIndex]] = newIndex;
        orderedFaces.push_back(faces[order[newIndex]]);
    }
    faces = std::move(orderedFaces);
    for (Cell& cell : cells) {
        for (int& face : cell.faces) {
            face = oldToNew[face];
        }
    }

    internalFaceCount_ = 0;
    for (Face& face : faces) {
        const Vector2& a = nodes[face.node0];
        const Vector2& b = nodes[face.node1];
        const Vector2 edge = b - a;
        face.area = norm(edge);
        if (face.area <= 1.0e-14) {
            throw std::runtime_error("Degenerate face in mesh");
        }
        face.normal = {edge.y / face.area, -edge.x / face.area};
        face.centre = 0.5 * (a + b);
        const Cell& owner = cells[face.owner];
        if (face.isBoundary()) {
            face.delta = dot(face.normal, face.centre - owner.centre);
            const auto patch = boundaryFacePatches.find({std::min(face.node0, face.node1),
                                                         std::max(face.node0, face.node1)});
            if (patch == boundaryFacePatches.end()) {
                throw std::runtime_error("A geometric boundary face has no Gmsh line element");
            }
            face.patch = patch->second;
        } else {
            ++internalFaceCount_;
            const Cell& neighbour = cells[face.neighbour];
            const Vector2 centreDistance = neighbour.centre - owner.centre;
            face.delta = dot(face.normal, centreDistance);
            const double ownerDistance = norm(face.centre - owner.centre);
            const double neighbourDistance = norm(face.centre - neighbour.centre);
            face.ownerWeight = neighbourDistance / (ownerDistance + neighbourDistance);
        }
        if (face.delta <= 1.0e-12 * std::sqrt(owner.volume)) {
            throw std::runtime_error("Mesh contains a strongly non-convex or inverted cell");
        }
    }

    for (BoundaryPatch& patch : patches) {
        patch.faces.clear();
    }
    for (int faceIndex = 0; faceIndex < static_cast<int>(faces.size()); ++faceIndex) {
        const Face& face = faces[faceIndex];
        if (face.isBoundary()) {
            if (face.patch < 0 || face.patch >= static_cast<int>(patches.size())) {
                throw std::runtime_error("Invalid physical boundary tag");
            }
            patches[face.patch].faces.push_back(faceIndex);
        }
    }
}

std::size_t Mesh::internalFaceCount() const {
    return internalFaceCount_;
}

int Mesh::patchIndex(const std::string& name) const {
    for (std::size_t i = 0; i < patches.size(); ++i) {
        if (patches[i].name == name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

double Mesh::minDelta() const {
    double result = faces.front().delta;
    for (const Face& face : faces) {
        result = std::min(result, face.delta);
    }
    return result;
}

} // namespace mp
