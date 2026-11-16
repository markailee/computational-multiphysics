#pragma once

#include "fvm/Field.hpp"

#include <filesystem>
#include <map>
#include <string>

namespace mp {

void writeTecplot(const std::filesystem::path& file, const Mesh& mesh,
                  const std::map<std::string, const ScalarField*>& scalars,
                  const std::map<std::string, const VectorField*>& vectors,
                  const std::string& zoneTitle = "flow", double solutionTime = 0.0);

} // namespace mp
