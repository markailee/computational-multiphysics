#include "fvm/Config.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mp {
namespace {

std::string trim(const std::string& value) {
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char c) { return std::isspace(c); });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) { return std::isspace(c); }).base();
    if (first >= last) return {};
    return std::string(first, last);
}

} // namespace

std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::tolower(c); });
    return value;
}

Config::Config(const std::filesystem::path& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("Cannot open case configuration: " + file.string());
    std::string section;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::size_t comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        line = trim(line);
        if (line.empty()) continue;
        if (line.front() == '[' && line.back() == ']') {
            section = trim(line.substr(1, line.size() - 2));
            continue;
        }
        const std::size_t equals = line.find('=');
        if (equals == std::string::npos) {
            throw std::runtime_error("Invalid configuration line " + std::to_string(lineNumber) + " in " + file.string());
        }
        const std::string key = trim(line.substr(0, equals));
        std::string value = trim(line.substr(equals + 1));
        if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
            value = value.substr(1, value.size() - 2);
        }
        values_[section.empty() ? key : section + "." + key] = value;
    }
}

bool Config::has(const std::string& key) const {
    return values_.find(key) != values_.end();
}

std::string Config::get(const std::string& key, const std::string& fallback) const {
    const auto iterator = values_.find(key);
    return iterator == values_.end() ? fallback : iterator->second;
}

std::string Config::require(const std::string& key) const {
    const auto iterator = values_.find(key);
    if (iterator == values_.end()) throw std::runtime_error("Missing required case setting: " + key);
    return iterator->second;
}

double Config::number(const std::string& key, double fallback) const {
    if (!has(key)) return fallback;
    std::istringstream input(require(key));
    double value = 0.0;
    input >> value;
    if (!input || !input.eof()) throw std::runtime_error("Expected a number for setting: " + key);
    return value;
}

int Config::integer(const std::string& key, int fallback) const {
    if (!has(key)) return fallback;
    std::istringstream input(require(key));
    int value = 0;
    input >> value;
    if (!input || !input.eof()) throw std::runtime_error("Expected an integer for setting: " + key);
    return value;
}

bool Config::boolean(const std::string& key, bool fallback) const {
    const std::string value = lowercase(get(key, fallback ? "true" : "false"));
    if (value == "true" || value == "yes" || value == "1" || value == "on") return true;
    if (value == "false" || value == "no" || value == "0" || value == "off") return false;
    throw std::runtime_error("Expected a Boolean for setting: " + key);
}

Vector2 Config::vector(const std::string& key, const Vector2& fallback) const {
    if (!has(key)) return fallback;
    std::istringstream input(require(key));
    Vector2 value;
    input >> value;
    if (!input) throw std::runtime_error("Expected two components for setting: " + key);
    return value;
}

} // namespace mp
