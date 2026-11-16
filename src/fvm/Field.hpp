#pragma once

#include "mesh/Mesh.hpp"
#include "mesh/Vector2.hpp"

#include <array>
#include <string>
#include <vector>

namespace mp {

template <typename Value>
class Field {
public:
    Field(std::string name, const Mesh& mesh, const Value& initial = Value{})
        : name_(std::move(name)), mesh_(&mesh), values_(mesh.cells.size(), initial),
          boundary_(mesh.faces.size(), initial), boundaryDirichlet_(mesh.faces.size(), false) {}

    const std::string& name() const { return name_; }
    const Mesh& mesh() const { return *mesh_; }
    std::size_t size() const { return values_.size(); }

    Value& operator[](std::size_t cell) { return values_[cell]; }
    const Value& operator[](std::size_t cell) const { return values_[cell]; }
    Value& boundary(std::size_t face) { return boundary_[face]; }
    const Value& boundary(std::size_t face) const { return boundary_[face]; }
    char& boundaryDirichlet(std::size_t face) { return boundaryDirichlet_[face]; }
    bool boundaryDirichlet(std::size_t face) const { return boundaryDirichlet_[face]; }

    std::vector<Value>& values() { return values_; }
    const std::vector<Value>& values() const { return values_; }

private:
    std::string name_;
    const Mesh* mesh_;
    std::vector<Value> values_;
    std::vector<Value> boundary_;
    std::vector<char> boundaryDirichlet_;
};

using ScalarField = Field<double>;
using VectorField = Field<Vector2>;
using ScalarGradient = Vector2;
using VectorGradient = std::array<Vector2, 2>;

} // namespace mp
