#pragma once

#include <string>

namespace mp {

enum class TimeScheme { Steady, Euler, BDF2 };

struct TimeCoefficients {
    double current = 0.0;
    double previous = 0.0;
    double previousPrevious = 0.0;
};

class TimeDerivative {
public:
    explicit TimeDerivative(TimeScheme scheme = TimeScheme::Steady) : scheme_(scheme) {}
    explicit TimeDerivative(const std::string& scheme);

    TimeCoefficients coefficients(double dt, bool hasPreviousPrevious) const;
    TimeScheme scheme() const { return scheme_; }
    bool isSteady() const { return scheme_ == TimeScheme::Steady; }

private:
    TimeScheme scheme_;
};

} // namespace mp
