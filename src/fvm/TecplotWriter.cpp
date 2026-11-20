#include "fvm/TecplotWriter.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace mp {

void writeTecplot(const std::filesystem::path& file, const Mesh& mesh,
                  const std::map<std::string, const ScalarField*>& scalars,
                  const std::map<std::string, const VectorField*>& vectors,
                  const std::string& zoneTitle, double solutionTime) {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create Tecplot file: " + file.string());
    output << std::setprecision(12) << std::scientific;

    std::vector<std::string> names = {"X", "Y"};
    for (const auto& [name, field] : vectors) {
        names.push_back(name + "_x");
        names.push_back(name + "_y");
    }
    for (const auto& [name, field] : scalars) names.push_back(name);
    output << "VARIABLES = ";
    for (const auto& name : names) output << '"' << name << "\" ";
    output << '\n';

    for (int cellSize = 3; cellSize <= 4; ++cellSize) {
        std::vector<int> cells;
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
            if (static_cast<int>(mesh.cells[i].nodes.size()) == cellSize) cells.push_back(static_cast<int>(i));
        }
        if (cells.empty()) continue;
        output << "ZONE T=\"" << zoneTitle << (cellSize == 3 ? "-triangles" : "-quadrilaterals")
               << "\", N=" << mesh.nodes.size() << ", E=" << cells.size()
               << ", DATAPACKING=BLOCK, ZONETYPE=" << (cellSize == 3 ? "FETRIANGLE" : "FEQUADRILATERAL")
               << ", SOLUTIONTIME=" << solutionTime;
        if (!names.empty()) output << ", VARLOCATION=([3-" << names.size() << "]=CELLCENTERED)";
        output << '\n';
        for (int component = 0; component < 2; ++component) {
            for (const Vector2& node : mesh.nodes) output << node[component] << '\n';
        }
        auto writeCellValues = [&output, &cells](const auto& values) {
            for (int cell : cells) output << values[cell] << '\n';
        };
        for (const auto& [name, field] : vectors) {
            for (int component = 0; component < 2; ++component) {
                for (int cell : cells) output << (*field)[cell][component] << '\n';
            }
        }
        for (const auto& [name, field] : scalars) writeCellValues(field->values());
        for (int cell : cells) {
            const auto& connectivity = mesh.cells[cell].nodes;
            for (int node : connectivity) output << node + 1 << ' ';
            output << '\n';
        }
    }
    if (!output) throw std::runtime_error("Failed while writing Tecplot file: " + file.string());
}

void writeNodalTecplot(const std::filesystem::path& file, const Mesh& mesh,
                       const std::map<std::string, const std::vector<double>*>& scalars,
                       const std::string& zoneTitle, double solutionTime) {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot create Tecplot file: " + file.string());
    output << std::setprecision(12) << std::scientific;

    std::vector<std::string> names = {"X", "Y"};
    for (const auto& [name, values] : scalars) {
        if (values->size() != mesh.nodes.size()) throw std::runtime_error("Nodal Tecplot field size mismatch");
        names.push_back(name);
    }
    output << "VARIABLES = ";
    for (const auto& name : names) output << '"' << name << "\" ";
    output << '\n';

    for (int cellSize = 3; cellSize <= 4; ++cellSize) {
        std::vector<int> cells;
        for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
            if (static_cast<int>(mesh.cells[i].nodes.size()) == cellSize) cells.push_back(static_cast<int>(i));
        }
        if (cells.empty()) continue;
        output << "ZONE T=\"" << zoneTitle << (cellSize == 3 ? "-triangles" : "-quadrilaterals")
               << "\", N=" << mesh.nodes.size() << ", E=" << cells.size()
               << ", DATAPACKING=BLOCK, ZONETYPE=" << (cellSize == 3 ? "FETRIANGLE" : "FEQUADRILATERAL")
               << ", SOLUTIONTIME=" << solutionTime << '\n';
        for (int component = 0; component < 2; ++component) {
            for (const Vector2& node : mesh.nodes) output << node[component] << '\n';
        }
        for (const auto& [name, values] : scalars) {
            for (double value : *values) output << value << '\n';
        }
        for (int cell : cells) {
            for (int node : mesh.cells[cell].nodes) output << node + 1 << ' ';
            output << '\n';
        }
    }
    if (!output) throw std::runtime_error("Failed while writing Tecplot file: " + file.string());
}

} // namespace mp
