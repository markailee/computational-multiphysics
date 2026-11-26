#pragma once

#include "fvm/Field.hpp"
#include "fvm/Gradient.hpp"
#include "fvm/PetscLinearSystem.hpp"
#include "fvm/TimeDerivative.hpp"

#include <string>
#include <vector>

namespace mp {

struct ScalarTransportEquation {
    ScalarField& field;
    const ScalarField& previous;
    const ScalarField& old;
    const ScalarField& oldOld;
    const std::vector<double>& massFlux;
    const std::vector<double>& diffusivity;
    const std::vector<double>& source;
    const std::vector<double>& destruction;
    const std::vector<char>& constrainedCells;
    const std::vector<double>& constrainedValues;
    double density = 1.0;
    TimeCoefficients time;
    double relaxation = 0.7;
    double linearRelativeTolerance = 1.0e-9;
    double linearAbsoluteTolerance = 1.0e-12;
    int maximumLinearIterations = 1000;
    std::string linearSolver = "bcgs";
    std::string preconditioner = "ilu";
    const std::vector<double>* additionalFaceFlux = nullptr;
    const std::vector<double>* densityField = nullptr;
    const std::vector<double>* oldDensityField = nullptr;
    const std::vector<double>* oldOldDensityField = nullptr;
    bool barthJespersen = false;
    GradientLimiter gradientLimiter = GradientLimiter::None;
    bool secondOrderConvection = true;
};

LinearSolveInfo solveScalarTransport(const Mesh& mesh, const GradientScheme& gradient,
                                     PetscLinearSystem& system, ScalarTransportEquation& equation);

} // namespace mp
