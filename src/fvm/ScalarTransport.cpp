#include "fvm/ScalarTransport.hpp"

#include "fvm/Operators.hpp"

#include <cmath>
#include <stdexcept>

namespace mp {

LinearSolveInfo solveScalarTransport(const Mesh& mesh, const GradientScheme& gradientScheme,
                                     PetscLinearSystem& system, ScalarTransportEquation& equation) {
    const int cellCount = static_cast<int>(mesh.cells.size());
    if (equation.additionalFaceFlux && equation.additionalFaceFlux->size() != mesh.faces.size()) {
        throw std::runtime_error("Additional scalar face-flux size does not match the mesh");
    }
    if ((equation.densityField && equation.densityField->size() != mesh.cells.size()) ||
        (equation.oldDensityField && equation.oldDensityField->size() != mesh.cells.size()) ||
        (equation.oldOldDensityField && equation.oldOldDensityField->size() != mesh.cells.size())) {
        throw std::runtime_error("Scalar-transport density-field size does not match the mesh");
    }
    const auto gradient = gradientScheme.scalar(equation.field);
    std::vector<ScalarGradient> limitedGradient;
    const std::vector<ScalarGradient>* convectionGradient = &gradient;
    const GradientLimiter activeLimiter = equation.gradientLimiter != GradientLimiter::None
        ? equation.gradientLimiter
        : (equation.barthJespersen ? GradientLimiter::BarthJespersen : GradientLimiter::None);
    if (activeLimiter != GradientLimiter::None) {
        limitedGradient = applyGradientLimiter(equation.field, gradient, activeLimiter);
        convectionGradient = &limitedGradient;
    }
    SparseMatrix matrix(cellCount);
    std::vector<double> rhs(cellCount);
    for (int cell = 0; cell < cellCount; ++cell) {
        const double volume = mesh.cells[cell].volume;
        const double currentDensity = equation.densityField ? (*equation.densityField)[cell] : equation.density;
        const double oldDensity = equation.oldDensityField ? (*equation.oldDensityField)[cell] : currentDensity;
        const double oldOldDensity = equation.oldOldDensityField ? (*equation.oldOldDensityField)[cell] : oldDensity;
        matrix.add(cell, cell, currentDensity * volume * equation.time.current +
                               volume * equation.destruction[cell]);
        rhs[cell] = volume * equation.source[cell] +
                    oldDensity * volume * equation.time.previous * equation.old[cell] +
                    oldOldDensity * volume * equation.time.previousPrevious * equation.oldOld[cell];
    }

    for (std::size_t faceIndex = 0; faceIndex < mesh.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh.faces[faceIndex];
        const int owner = face.owner;
        const int neighbour = face.neighbour;
        const double faceDiffusivity = face.ownerWeight * equation.diffusivity[owner] +
            (1.0 - face.ownerWeight) * equation.diffusivity[neighbour];
        const double diffusion = diffusionConductance(mesh, faceIndex, faceDiffusivity);
        matrix.add(owner, owner, diffusion);
        matrix.add(owner, neighbour, -diffusion);
        matrix.add(neighbour, neighbour, diffusion);
        matrix.add(neighbour, owner, -diffusion);

        const Vector2 correctionVector = nonOrthogonalVector(mesh, faceIndex);
        const Vector2 faceGradient = face.ownerWeight * gradient[owner] +
            (1.0 - face.ownerWeight) * gradient[neighbour];
        const double correction = faceDiffusivity * dot(correctionVector, faceGradient);
        rhs[owner] += correction;
        rhs[neighbour] -= correction;

        const double flux = equation.massFlux[faceIndex];
        const ConvectionCoefficients convection = upwindConvection(flux);
        matrix.add(owner, owner, convection.ownerDiagonal);
        matrix.add(owner, neighbour, convection.ownerToNeighbour);
        matrix.add(neighbour, neighbour, convection.neighbourDiagonal);
        matrix.add(neighbour, owner, convection.neighbourToOwner);
        const double upwind = flux >= 0.0 ? equation.field[owner] : equation.field[neighbour];
        const double reconstructed = equation.secondOrderConvection
            ? secondOrderUpwindFaceValue(mesh, equation.field, *convectionGradient, faceIndex, flux) : upwind;
        const double deferred = flux * (reconstructed - upwind);
        rhs[owner] -= deferred;
        rhs[neighbour] += deferred;
        if (equation.additionalFaceFlux) {
            rhs[owner] -= (*equation.additionalFaceFlux)[faceIndex];
            rhs[neighbour] += (*equation.additionalFaceFlux)[faceIndex];
        }
    }

    for (std::size_t faceIndex = mesh.internalFaceCount(); faceIndex < mesh.faces.size(); ++faceIndex) {
        const Face& face = mesh.faces[faceIndex];
        const int owner = face.owner;
        if (equation.field.boundaryDirichlet(faceIndex)) {
            const double diffusion = diffusionConductance(mesh, faceIndex, equation.diffusivity[owner]);
            matrix.add(owner, owner, diffusion);
            rhs[owner] += diffusion * equation.field.boundary(faceIndex) +
                equation.diffusivity[owner] * dot(nonOrthogonalVector(mesh, faceIndex), gradient[owner]);
        }

        const double flux = equation.massFlux[faceIndex];
        const double upwind = flux >= 0.0 ? equation.field[owner] : equation.field.boundary(faceIndex);
        const double reconstructed = equation.secondOrderConvection
            ? secondOrderUpwindFaceValue(mesh, equation.field, *convectionGradient, faceIndex, flux) : upwind;
        if (flux >= 0.0) matrix.add(owner, owner, flux);
        else rhs[owner] -= flux * upwind;
        rhs[owner] -= flux * (reconstructed - upwind);
        if (equation.additionalFaceFlux) rhs[owner] -= (*equation.additionalFaceFlux)[faceIndex];
    }

    for (int cell = 0; cell < cellCount; ++cell) {
        if (equation.constrainedCells[cell]) continue;
        const Cell& controlVolume = mesh.cells[cell];
        const double currentDensity = equation.densityField ? (*equation.densityField)[cell] : equation.density;
        const double oldDensity = equation.oldDensityField ? (*equation.oldDensityField)[cell] : currentDensity;
        const double oldOldDensity = equation.oldOldDensityField ? (*equation.oldOldDensityField)[cell] : oldDensity;
        double imbalance = controlVolume.volume * (equation.time.current * currentDensity -
            equation.time.previous * oldDensity - equation.time.previousPrevious * oldOldDensity);
        for (std::size_t i = 0; i < controlVolume.faces.size(); ++i) {
            imbalance += controlVolume.faceSigns[i] * equation.massFlux[controlVolume.faces[i]];
        }
        if (imbalance < 0.0) {
            matrix.add(cell, cell, -imbalance);
        } else if (imbalance > 0.0) {
            rhs[cell] += imbalance * equation.previous[cell];
        }
    }

    for (int cell = 0; cell < cellCount; ++cell) {
        if (equation.constrainedCells[cell]) {
            for (int other = 0; other < cellCount; ++other) {
                if (other == cell) continue;
                const auto& entries = matrix.row(other);
                const auto coefficient = entries.find(cell);
                if (coefficient != entries.end()) rhs[other] -= coefficient->second * equation.constrainedValues[cell];
            }
            matrix.pinZero(cell);
            rhs[cell] = equation.constrainedValues[cell];
            continue;
        }
        const double diagonal = matrix.diagonal(cell);
        if (!(diagonal > 0.0)) throw std::runtime_error("Non-positive scalar-transport diagonal");
        const double relaxation = diagonal * (1.0 / equation.relaxation - 1.0);
        matrix.add(cell, cell, relaxation);
        rhs[cell] += relaxation * equation.previous[cell];
    }

    system.setMatrix(matrix);
    return system.solve(rhs, equation.field.values(), equation.linearSolver, equation.preconditioner,
                        equation.linearRelativeTolerance, equation.linearAbsoluteTolerance,
                        equation.maximumLinearIterations);
}

} // namespace mp
