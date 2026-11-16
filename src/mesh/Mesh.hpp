#pragma once

#include "mesh/Vector2.hpp"

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace mp {

struct Cell {
    std::vector<int> nodes;
    std::vector<int> faces;
    std::vector<int> faceSigns;
    Vector2 centre;
    double volume = 0.0;
};

struct Face {
    int node0 = -1;
    int node1 = -1;
    int owner = -1;
    int neighbour = -1;
    int patch = -1;
    Vector2 centre;
    Vector2 normal;
    double area = 0.0;
    double delta = 0.0;
    double ownerWeight = 0.5;

    bool isBoundary() const { return neighbour < 0; }
};

struct BoundaryPatch {
    std::string name;
    std::vector<int> faces;
};

class Mesh {
public:
    std::vector<Vector2> nodes;
    std::vector<Cell> cells;
    std::vector<Face> faces;
    std::vector<BoundaryPatch> patches;

    using BoundaryFaceMap = std::map<std::pair<int, int>, int>;
    void buildGeometry(const BoundaryFaceMap& boundaryFacePatches = {});
    std::size_t internalFaceCount() const;
    std::size_t boundaryFaceCount() const { return faces.size() - internalFaceCount(); }
    int patchIndex(const std::string& name) const;
    double minDelta() const;

private:
    std::size_t internalFaceCount_ = 0;
};

} // namespace mp
