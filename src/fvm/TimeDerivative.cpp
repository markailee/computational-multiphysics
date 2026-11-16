#include "fvm/TimeDerivative.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace mp {

TimeDerivative::TimeDerivative(const std::string& scheme) {
    std::string value = scheme;
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return std::tolower(c); });
    if (value == "steady") scheme_ = TimeScheme::Steady;
    else if (value == "euler" || value == "backward_euler") scheme_ = TimeScheme::Euler;
    else if (value == "bdf2") scheme_ = TimeScheme::BDF2;
    else throw std::runtime_error("Unknown time scheme: " + scheme);
}

TimeCoefficients TimeDerivative::coefficients(double dt, bool hasPreviousPrevious) const {
    if (scheme_ == TimeScheme::Steady) return {};
    if (dt <= 0.0) throw std::runtime_error("Transient calculations require dt > 0");
    if (scheme_ == TimeScheme::Euler || !hasPreviousPrevious) return {1.0 / dt, 1.0 / dt, 0.0};
    return {1.5 / dt, 2.0 / dt, -0.5 / dt};
}

} // namespace mp
