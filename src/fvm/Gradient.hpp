#pragma once

#include "fvm/Config.hpp"
#include "fvm/Field.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mp {

class GradientScheme {
public:
    virtual ~GradientScheme() = default;
    virtual std::vector<ScalarGradient> scalar(const ScalarField& field) const = 0;
    virtual std::vector<VectorGradient> vector(const VectorField& field) const = 0;
};

class GaussGradient final : public GradientScheme {
public:
    explicit GaussGradient(int correctionIterations = 2) : correctionIterations_(correctionIterations) {}
    std::vector<ScalarGradient> scalar(const ScalarField& field) const override;
    std::vector<VectorGradient> vector(const VectorField& field) const override;

private:
    int correctionIterations_;
};

class LeastSquaresGradient final : public GradientScheme {
public:
    std::vector<ScalarGradient> scalar(const ScalarField& field) const override;
    std::vector<VectorGradient> vector(const VectorField& field) const override;
};

enum class GradientLimiter { None, BarthJespersen, Venkatakrishnan };

GradientLimiter readGradientLimiter(const Config& config);
bool useBarthJespersenLimiter(const Config& config);
std::vector<ScalarGradient> applyGradientLimiter(const ScalarField& field,
                                                 std::vector<ScalarGradient> gradients,
                                                 GradientLimiter limiter);
std::vector<VectorGradient> applyGradientLimiter(const VectorField& field,
                                                 std::vector<VectorGradient> gradients,
                                                 GradientLimiter limiter);
std::vector<ScalarGradient> applyBarthJespersenLimiter(const ScalarField& field,
                                                       std::vector<ScalarGradient> gradients);
std::vector<VectorGradient> applyBarthJespersenLimiter(const VectorField& field,
                                                       std::vector<VectorGradient> gradients);

std::unique_ptr<GradientScheme> makeGradientScheme(const std::string& name);

} // namespace mp
