#pragma once

#include <filesystem>
#include <string>

namespace mp {

class PhysicsSolver {
public:
    virtual ~PhysicsSolver() = default;
    virtual std::string name() const = 0;
    virtual void initialize() = 0;
    virtual void advance(double dt) = 0;
    virtual void write(const std::filesystem::path& file) const = 0;
};

} // namespace mp
