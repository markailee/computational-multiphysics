#pragma once

#include "mesh/Mesh.hpp"

#include <filesystem>

namespace mp {

Mesh readGmsh(const std::filesystem::path& file);

} // namespace mp
