#include "fvm/PetscLinearSystem.hpp"

#include <stdexcept>

namespace mp {
namespace {

void check(PetscErrorCode status) {
    if (status != PETSC_SUCCESS) throw std::runtime_error("PETSc operation failed with error code " + std::to_string(status));
}

} // namespace

double SparseMatrix::diagonal(int row) const {
    const auto iterator = rows_[row].find(row);
    return iterator == rows_[row].end() ? 0.0 : iterator->second;
}

void SparseMatrix::pinZero(int pinnedRow) {
    for (auto& row : rows_) row.erase(pinnedRow);
    rows_[pinnedRow].clear();
    rows_[pinnedRow][pinnedRow] = 1.0;
}

PetscLinearSystem::PetscLinearSystem(int size) : size_(size) {
    check(MatCreate(PETSC_COMM_SELF, &matrix_));
    check(MatSetType(matrix_, MATSEQAIJ));
    check(MatSetSizes(matrix_, size, size, size, size));
    check(MatSeqAIJSetPreallocation(matrix_, 10, PETSC_NULLPTR));
    check(MatSetOption(matrix_, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    check(VecCreate(PETSC_COMM_SELF, &rhs_));
    check(VecSetSizes(rhs_, size, size));
    check(VecSetType(rhs_, VECSEQ));
    check(VecDuplicate(rhs_, &solution_));
    check(KSPCreate(PETSC_COMM_SELF, &ksp_));
    check(KSPSetOperators(ksp_, matrix_, matrix_));
    check(KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE));
}

PetscLinearSystem::~PetscLinearSystem() {
    KSPDestroy(&ksp_);
    VecDestroy(&solution_);
    VecDestroy(&rhs_);
    MatDestroy(&matrix_);
}

void PetscLinearSystem::setMatrix(const SparseMatrix& matrix) {
    if (matrix.size() != size_) throw std::runtime_error("Linear system size mismatch");
    check(MatZeroEntries(matrix_));
    for (int row = 0; row < size_; ++row) {
        const auto& entries = matrix.row(row);
        std::vector<PetscInt> columns;
        std::vector<PetscScalar> values;
        columns.reserve(entries.size());
        values.reserve(entries.size());
        for (const auto& [column, value] : entries) {
            columns.push_back(column);
            values.push_back(value);
        }
        check(MatSetValues(matrix_, 1, &row, static_cast<PetscInt>(columns.size()), columns.data(), values.data(), INSERT_VALUES));
    }
    check(MatAssemblyBegin(matrix_, MAT_FINAL_ASSEMBLY));
    check(MatAssemblyEnd(matrix_, MAT_FINAL_ASSEMBLY));
}

LinearSolveInfo PetscLinearSystem::solve(const std::vector<double>& rightHandSide, std::vector<double>& solution,
                                         const std::string& solver, const std::string& preconditioner,
                                         double relativeTolerance, double absoluteTolerance, int maximumIterations) {
    if (static_cast<int>(rightHandSide.size()) != size_ || static_cast<int>(solution.size()) != size_) {
        throw std::runtime_error("Linear vector size mismatch");
    }
    for (int i = 0; i < size_; ++i) {
        check(VecSetValue(rhs_, i, rightHandSide[i], INSERT_VALUES));
        check(VecSetValue(solution_, i, solution[i], INSERT_VALUES));
    }
    check(VecAssemblyBegin(rhs_));
    check(VecAssemblyEnd(rhs_));
    check(VecAssemblyBegin(solution_));
    check(VecAssemblyEnd(solution_));
    check(KSPSetType(ksp_, solver.c_str()));
    PC preconditionerObject = nullptr;
    check(KSPGetPC(ksp_, &preconditionerObject));
    check(PCSetType(preconditionerObject, preconditioner.c_str()));
    check(KSPSetTolerances(ksp_, relativeTolerance, absoluteTolerance, PETSC_DEFAULT, maximumIterations));
    check(KSPSetFromOptions(ksp_));
    check(KSPSetUp(ksp_));
    check(KSPSolve(ksp_, rhs_, solution_));
    KSPConvergedReason reason;
    check(KSPGetConvergedReason(ksp_, &reason));
    if (reason < 0) {
        throw std::runtime_error("PETSc " + solver + " solve failed with convergence reason " + std::to_string(reason));
    }
    LinearSolveInfo info;
    check(KSPGetIterationNumber(ksp_, &info.iterations));
    check(KSPGetResidualNorm(ksp_, &info.residual));
    const PetscScalar* values = nullptr;
    check(VecGetArrayRead(solution_, &values));
    for (int i = 0; i < size_; ++i) solution[i] = static_cast<double>(values[i]);
    check(VecRestoreArrayRead(solution_, &values));
    return info;
}

} // namespace mp
