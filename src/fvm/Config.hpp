#pragma once

#include "mesh/Vector2.hpp"

#include <filesystem>
#include <map>
#include <string>

namespace mp {

class Config {
public:
    Config() = default;
    explicit Config(const std::filesystem::path& file);

    bool has(const std::string& key) const;
    std::string get(const std::string& key, const std::string& fallback) const;
    std::string require(const std::string& key) const;
    double number(const std::string& key, double fallback) const;
    int integer(const std::string& key, int fallback) const;
    bool boolean(const std::string& key, bool fallback) const;
    Vector2 vector(const std::string& key, const Vector2& fallback) const;

private:
    std::map<std::string, std::string> values_;
};

std::string lowercase(std::string value);

} // namespace mp
