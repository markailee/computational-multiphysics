#include "mesh/GmshReader.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mp {
namespace {

using PhysicalNames = std::map<std::pair<int, int>, std::string>;
using EntityPhysicals = std::map<std::pair<int, int>, std::vector<int>>;

struct RawElement {
    int type = 0;
    int entityDimension = 0;
    int entity = 0;
    std::vector<int> physicals;
    std::vector<int> nodeTags;
};

int nodeCount(int type) {
    switch (type) {
    case 1: return 2;
    case 2: return 3;
    case 3: return 4;
    case 8: return 3;
    case 9: return 6;
    case 10: return 9;
    case 15: return 1;
    default: throw std::runtime_error("Unsupported Gmsh element type " + std::to_string(type));
    }
}

void requireSection(std::istream& input, const std::string& expected) {
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line == expected) return;
        break;
    }
    throw std::runtime_error("Malformed Gmsh file: expected " + expected);
}

PhysicalNames readPhysicalNames(std::istream& input) {
    int count = 0;
    input >> count;
    std::string line;
    std::getline(input, line);
    PhysicalNames names;
    for (int i = 0; i < count; ++i) {
        std::getline(input, line);
        const std::size_t firstQuote = line.find('"');
        const std::size_t lastQuote = line.rfind('"');
        if (firstQuote == std::string::npos || lastQuote == firstQuote) {
            throw std::runtime_error("Malformed $PhysicalNames section");
        }
        std::istringstream prefix(line.substr(0, firstQuote));
        int dimension = 0;
        int tag = 0;
        prefix >> dimension >> tag;
        names[{dimension, tag}] = line.substr(firstQuote + 1, lastQuote - firstQuote - 1);
    }
    requireSection(input, "$EndPhysicalNames");
    return names;
}

EntityPhysicals readEntities(std::istream& input) {
    std::size_t counts[4] = {};
    input >> counts[0] >> counts[1] >> counts[2] >> counts[3];
    EntityPhysicals physicals;
    for (int dimension = 0; dimension < 4; ++dimension) {
        for (std::size_t i = 0; i < counts[dimension]; ++i) {
            int tag = 0;
            input >> tag;
            if (dimension == 0) {
                double x = 0.0, y = 0.0, z = 0.0;
                input >> x >> y >> z;
            } else {
                double bounds = 0.0;
                for (int j = 0; j < 6; ++j) input >> bounds;
            }
            int physicalCount = 0;
            input >> physicalCount;
            std::vector<int> tags(physicalCount);
            for (int& physical : tags) input >> physical;
            physicals[{dimension, tag}] = tags;
            if (dimension > 0) {
                int boundingCount = 0;
                input >> boundingCount;
                for (int j = 0; j < boundingCount; ++j) {
                    int bounding = 0;
                    input >> bounding;
                }
            }
        }
    }
    requireSection(input, "$EndEntities");
    return physicals;
}

std::vector<int> readNodeTags(std::istream& input, int count) {
    std::vector<int> tags(count);
    for (int& tag : tags) input >> tag;
    return tags;
}

std::vector<RawElement> readElements2(std::istream& input) {
    int count = 0;
    input >> count;
    std::vector<RawElement> elements;
    elements.reserve(count);
    for (int i = 0; i < count; ++i) {
        int id = 0;
        int tagCount = 0;
        RawElement element;
        input >> id >> element.type >> tagCount;
        std::vector<int> tags(tagCount);
        for (int& tag : tags) input >> tag;
        if (!tags.empty()) element.physicals.push_back(tags[0]);
        element.nodeTags = readNodeTags(input, nodeCount(element.type));
        element.entityDimension = element.type == 1 || element.type == 8 ? 1 : 2;
        elements.push_back(std::move(element));
    }
    requireSection(input, "$EndElements");
    return elements;
}

std::vector<RawElement> readElements4(std::istream& input) {
    int blockCount = 0;
    int elementCount = 0;
    int minimumTag = 0;
    int maximumTag = 0;
    input >> blockCount >> elementCount >> minimumTag >> maximumTag;
    std::vector<RawElement> elements;
    elements.reserve(elementCount);
    for (int block = 0; block < blockCount; ++block) {
        int dimension = 0;
        int entity = 0;
        int type = 0;
        int count = 0;
        input >> dimension >> entity >> type >> count;
        for (int i = 0; i < count; ++i) {
            int id = 0;
            input >> id;
            RawElement element;
            element.type = type;
            element.entityDimension = dimension;
            element.entity = entity;
            element.nodeTags = readNodeTags(input, nodeCount(type));
            elements.push_back(std::move(element));
        }
    }
    requireSection(input, "$EndElements");
    return elements;
}

} // namespace

Mesh readGmsh(const std::filesystem::path& file) {
    std::ifstream input(file);
    if (!input) {
        throw std::runtime_error("Cannot open mesh file: " + file.string());
    }

    double version = 0.0;
    PhysicalNames physicalNames;
    EntityPhysicals entityPhysicals;
    std::vector<int> nodeTags;
    std::vector<Vector2> coordinates;
    std::vector<RawElement> elements;
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line == "$MeshFormat") {
            int binary = 0;
            int dataSize = 0;
            input >> version >> binary >> dataSize;
            if (binary != 0) {
                throw std::runtime_error("Only ASCII Gmsh files are supported; use Mesh.Binary=0");
            }
            std::getline(input, line);
            requireSection(input, "$EndMeshFormat");
        } else if (line == "$PhysicalNames") {
            physicalNames = readPhysicalNames(input);
        } else if (line == "$Entities") {
            entityPhysicals = readEntities(input);
        } else if (line == "$Nodes") {
            if (version < 4.0) {
                int count = 0;
                input >> count;
                nodeTags.resize(count);
                coordinates.resize(count);
                for (int i = 0; i < count; ++i) {
                    double z = 0.0;
                    input >> nodeTags[i] >> coordinates[i].x >> coordinates[i].y >> z;
                }
            } else {
                int blockCount = 0;
                int count = 0;
                int minimumTag = 0;
                int maximumTag = 0;
                input >> blockCount >> count >> minimumTag >> maximumTag;
                nodeTags.reserve(count);
                coordinates.reserve(count);
                for (int block = 0; block < blockCount; ++block) {
                    int dimension = 0;
                    int entity = 0;
                    int parametric = 0;
                    int blockSize = 0;
                    input >> dimension >> entity >> parametric >> blockSize;
                    std::vector<int> blockTags = readNodeTags(input, blockSize);
                    for (int tag : blockTags) {
                        Vector2 point;
                        double z = 0.0;
                        input >> point.x >> point.y >> z;
                        for (int j = 0; j < parametric * dimension; ++j) {
                            double parameter = 0.0;
                            input >> parameter;
                        }
                        nodeTags.push_back(tag);
                        coordinates.push_back(point);
                    }
                }
            }
            requireSection(input, "$EndNodes");
        } else if (line == "$Elements") {
            elements = version < 4.0 ? readElements2(input) : readElements4(input);
        }
    }

    if (version < 2.0 || version > 4.1 || nodeTags.empty() || elements.empty()) {
        throw std::runtime_error("Unsupported or incomplete Gmsh mesh: " + file.string());
    }

    std::map<int, int> tagToNode;
    for (std::size_t i = 0; i < nodeTags.size(); ++i) tagToNode[nodeTags[i]] = static_cast<int>(i);
    auto localNode = [&tagToNode](int tag) {
        const auto iterator = tagToNode.find(tag);
        if (iterator == tagToNode.end()) throw std::runtime_error("Element references an unknown node");
        return iterator->second;
    };
    auto patchName = [&physicalNames](int dimension, int tag) {
        const auto iterator = physicalNames.find({dimension, tag});
        return iterator == physicalNames.end() ? "patch_" + std::to_string(tag) : iterator->second;
    };

    Mesh mesh;
    mesh.nodes = std::move(coordinates);
    std::map<std::string, int> patchIndices;
    Mesh::BoundaryFaceMap boundaryFaces;
    for (RawElement& element : elements) {
        if (version >= 4.0) {
            const auto iterator = entityPhysicals.find({element.entityDimension, element.entity});
            if (iterator != entityPhysicals.end()) element.physicals = iterator->second;
        }
        if (element.type == 15) continue;
        if (element.type != 1 && element.type != 2 && element.type != 3) {
            throw std::runtime_error("Only linear Gmsh lines, triangles, and quadrilaterals are supported");
        }
        const int physical = element.physicals.empty() ? element.entity : element.physicals.front();
        const int dimension = element.type == 1 ? 1 : 2;
        if (dimension == 2) {
            Cell cell;
            for (int tag : element.nodeTags) cell.nodes.push_back(localNode(tag));
            mesh.cells.push_back(std::move(cell));
        } else {
            const std::string name = patchName(1, physical);
            auto [iterator, inserted] = patchIndices.emplace(name, static_cast<int>(mesh.patches.size()));
            if (inserted) mesh.patches.push_back({name, {}});
            const int a = localNode(element.nodeTags[0]);
            const int b = localNode(element.nodeTags[1]);
            boundaryFaces[{std::min(a, b), std::max(a, b)}] = iterator->second;
        }
    }
    mesh.buildGeometry(boundaryFaces);
    return mesh;
}

} // namespace mp
