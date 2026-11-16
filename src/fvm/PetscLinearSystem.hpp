#pragma once

#include <petscksp.h>

#include <map>
#include <string>
#include <vector>

namespace mp {

class SparseMatrix {
public:
    explicit SparseMatrix(int size = 0) : rows_(size) {}
    int size() const { return static_cast<int>(rows_.size()); }
    void add(int row, int column, double value) { rows_[row][column] += value; }
    double diagonal(int row) const;
    void pinZero(int row);
    const std::map<int, double>& row(int rowIndex) const { return rows_[rowIndex]; }

private:
    std::vector<std::map<int, double>> rows_;
};

struct LinearSolveInfo {
    int iterations = 0;
    double residual = 0.0;
};

class PetscLinearSystem {
public:
    explicit PetscLinearSystem(int size);
    ~PetscLinearSystem();

    PetscLinearSystem(const PetscLinearSystem&) = delete;
    PetscLinearSystem& operator=(const PetscLinearSystem&) = delete;

    void setMatrix(const SparseMatrix& matrix);
    LinearSolveInfo solve(const std::vector<double>& rightHandSide, std::vector<double>& solution,
                          const std::string& solver, const std::string& preconditioner,
                          double relativeTolerance, double absoluteTolerance, int maximumIterations);

private:
    Mat matrix_ = nullptr;
    Vec rhs_ = nullptr;
    Vec solution_ = nullptr;
    KSP ksp_ = nullptr;
    int size_ = 0;
};

} // namespace mp
