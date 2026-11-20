#pragma once

#include "fvm/Field.hpp"

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace mp {

void writeTecplot(const std::filesystem::path& file, const Mesh& mesh,
                  const std::map<std::string, const ScalarField*>& scalars,
                  const std::map<std::string, const VectorField*>& vectors,
                  const std::string& zoneTitle = "flow", double solutionTime = 0.0);

void writeNodalTecplot(const std::filesystem::path& file, const Mesh& mesh,
                       const std::map<std::string, const std::vector<double>*>& scalars,
                       const std::string& zoneTitle = "acoustic", double solutionTime = 0.0);

} // namespace mp
