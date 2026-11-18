#include "fvm/ScalarTransport.hpp"

#include "fvm/Operators.hpp"

#include <cmath>
#include <stdexcept>

namespace mp {

LinearSolveInfo solveScalarTransport(const Mesh& mesh, const GradientScheme& gradientScheme,
                                     PetscLinearSystem& system, ScalarTransportEquation& equation) {
    const int cellCount = static_cast<int>(mesh.cells.size());
    const auto gradient = gradientScheme.scalar(equation.field);
    SparseMatrix matrix(cellCount);
    std::vector<double> rhs(cellCount);
    for (int cell = 0; cell < cellCount; ++cell) {
        const double volume = mesh.cells[cell].volume;
        matrix.add(cell, cell, equation.density * volume * equation.time.current +
                               volume * equation.destruction[cell]);
        rhs[cell] = volume * equation.source[cell] +
                    equation.density * volume * equation.time.previous * equation.old[cell] +
                    equation.density * volume * equation.time.previousPrevious * equation.oldOld[cell];
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
        const double central = centralFaceValue(mesh, equation.field, gradient, faceIndex);
        const double deferred = equation.convectionBlending * flux * (central - upwind);
        rhs[owner] -= deferred;
        rhs[neighbour] += deferred;
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
        const double central = centralFaceValue(mesh, equation.field, gradient, faceIndex);
        const double upwind = flux >= 0.0 ? equation.field[owner] : central;
        if (flux >= 0.0) matrix.add(owner, owner, flux);
        else rhs[owner] -= flux * central;
        rhs[owner] -= equation.convectionBlending * flux * (central - upwind);
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
