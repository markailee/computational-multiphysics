#pragma once

#include "fvm/Config.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mp::material {

inline std::vector<std::string> names(const Config& config, const std::string& key) {
    std::istringstream input(config.require(key));
    std::vector<std::string> result;
    std::string name;
    while (input >> name) {
        if (name.empty() || !(std::isalpha(static_cast<unsigned char>(name.front())) || name.front() == '_') ||
            !std::all_of(name.begin() + 1, name.end(), [](unsigned char c) { return std::isalnum(c) || c == '_'; })) {
            throw std::runtime_error("Invalid material name in " + key + ": " + name);
        }
        if (std::find(result.begin(), result.end(), name) != result.end()) {
            throw std::runtime_error("Duplicate material name in " + key + ": " + name);
        }
        result.push_back(name);
    }
    if (result.empty()) throw std::runtime_error("At least one material name is required for " + key);
    return result;
}

inline double number(const Config& config, const std::string& key, bool positive = true) {
    std::istringstream input(config.require(key));
    double value = 0.0;
    input >> value;
    if (!input || !input.eof() || !std::isfinite(value) || (positive ? value <= 0.0 : value < 0.0)) {
        throw std::runtime_error(std::string("Expected a ") + (positive ? "positive" : "non-negative") +
                                 " number for " + key);
    }
    return value;
}

inline std::vector<double> fractions(const std::vector<double>& values, std::size_t expected,
                                     const std::string& context) {
    if (values.size() != expected) {
        throw std::runtime_error(context + " requires exactly " + std::to_string(expected) + " fractions");
    }
    double sum = 0.0;
    for (double value : values) {
        if (!std::isfinite(value) || value < -1.0e-12 || value > 1.0 + 1.0e-12) {
            throw std::runtime_error(context + " fractions must be in [0, 1]");
        }
        sum += value;
    }
    if (!(sum > 0.0) || std::abs(sum - 1.0) > 1.0e-8) {
        throw std::runtime_error(context + " fractions must sum to one within 1e-8");
    }
    std::vector<double> result = values;
    for (double& value : result) value /= sum;
    return result;
}

inline std::vector<double> readFractions(const Config& config, const std::string& key,
                                         std::size_t expected) {
    std::istringstream input(config.require(key));
    std::vector<double> values;
    double value = 0.0;
    while (input >> value) values.push_back(value);
    if (!input.eof()) throw std::runtime_error("Invalid numeric list for " + key);
    return fractions(values, expected, key);
}

} // namespace mp::material
