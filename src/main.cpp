#include "flow/SimpleSolver.hpp"
#include "mesh/GmshReader.hpp"

#include <petsc.h>

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    PetscErrorCode status = PetscInitialize(&argc, &argv, nullptr,
        "2D unstructured finite-volume multiphysics solver\n"
        "Usage: multiphysics <case-directory-or-case.in> [PETSc options]\n");
    if (status != PETSC_SUCCESS) return static_cast<int>(status);
    int result = 0;
    try {
        int ranks = 1;
        MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
        if (ranks != 1) throw std::runtime_error("This build currently supports serial execution; run without mpirun");
        if (argc < 2) throw std::runtime_error("Usage: multiphysics <case-directory-or-case.in> [PETSc options]");
        std::filesystem::path input = std::filesystem::absolute(argv[1]);
        const std::filesystem::path caseFile = std::filesystem::is_directory(input) ? input / "case.in" : input;
        const std::filesystem::path caseDirectory = caseFile.parent_path();
        mp::Config config(caseFile);
        std::filesystem::path meshFile = config.require("mesh");
        if (meshFile.is_relative()) meshFile = caseDirectory / meshFile;
        mp::Mesh mesh = mp::readGmsh(meshFile);
        mp::SimpleSolver solver(mesh, config, caseDirectory);
        solver.run();
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        result = 1;
    }
    status = PetscFinalize();
    return result != 0 ? result : static_cast<int>(status);
}
