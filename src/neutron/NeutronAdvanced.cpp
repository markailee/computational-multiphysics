#include "neutron/NeutronSolver.hpp"

#include <petscblaslapack.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace mp {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

double dotVec(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) sum += a[i] * b[i];
    return sum;
}

double normVec(const std::vector<double>& a) {
    return std::sqrt(dotVec(a, a));
}

void orthonormalizeColumns(std::vector<std::vector<double>>& basis) {
    // Two-pass Modified Gram-Schmidt for machine-precision orthogonality
    std::vector<std::vector<double>> kept;
    kept.reserve(basis.size());
    for (std::vector<double>& vec : basis) {
        for (int pass = 0; pass < 2; ++pass) {
            for (const std::vector<double>& q : kept) {
                const double h = dotVec(q, vec);
                for (std::size_t i = 0; i < vec.size(); ++i) vec[i] -= h * q[i];
            }
        }
        const double nrm = normVec(vec);
        if (nrm > 1.0e-12) {
            for (double& x : vec) x /= nrm;
            kept.push_back(std::move(vec));
        }
    }
    basis = std::move(kept);
}

struct SmallEigen {
    std::vector<double> wr;
    std::vector<double> wi;
    std::vector<std::vector<double>> vr; // column-major eigenvectors (real part for real eigenvalues)
};

SmallEigen solveSmallEigen(int m, const std::vector<double>& colMajorMatrix) {
    PetscBLASInt n = static_cast<PetscBLASInt>(m);
    PetscBLASInt lda = n;
    PetscBLASInt ldvl = 1;
    PetscBLASInt ldvr = n;
    PetscBLASInt info = 0;
    std::vector<PetscScalar> a(colMajorMatrix.begin(), colMajorMatrix.end());
    std::vector<PetscReal> wr(m, 0.0);
    std::vector<PetscReal> wi(m, 0.0);
    std::vector<PetscScalar> vl(1, 0.0);
    std::vector<PetscScalar> vr(m * m, 0.0);

    PetscBLASInt lwork = -1;
    PetscScalar workQuery = 0.0;
    LAPACKgeev_("N", "V", &n, a.data(), &lda, wr.data(), wi.data(),
                vl.data(), &ldvl, vr.data(), &ldvr, &workQuery, &lwork, &info);
    if (info != 0) throw std::runtime_error("LAPACKgeev workspace query failed");
    lwork = std::max<PetscBLASInt>(1, static_cast<PetscBLASInt>(PetscRealPart(workQuery)));
    std::vector<PetscScalar> work(lwork, 0.0);
    LAPACKgeev_("N", "V", &n, a.data(), &lda, wr.data(), wi.data(),
                vl.data(), &ldvl, vr.data(), &ldvr, work.data(), &lwork, &info);
    if (info != 0) throw std::runtime_error("LAPACKgeev failed with info = " + std::to_string(info));

    SmallEigen result;
    result.wr.assign(wr.begin(), wr.end());
    result.wi.assign(wi.begin(), wi.end());
    result.vr.assign(m, std::vector<double>(m, 0.0));
    for (int col = 0; col < m; ++col) {
        for (int row = 0; row < m; ++row) {
            result.vr[col][row] = static_cast<double>(PetscRealPart(vr[col * m + row]));
        }
    }
    return result;
}

} // namespace

bool NeutronSolver::cellMatchesFilter(int cell, const std::string& regionName, const std::string& materialName) const {
    if (!regionName.empty() && regionName != "all") {
        if (cellRegions_[cell] < 0 || regions_[cellRegions_[cell]].name != regionName) return false;
    }
    if (!materialName.empty() && materialName != "all") {
        if (materials_[cellMaterials_[cell]].name != materialName) return false;
    }
    return true;
}

std::vector<double> NeutronSolver::applyFissionOperator(const std::vector<double>& scalar) const {
    std::vector<double> result(groupCount_ * cellCount_, 0.0);
    for (int cell = 0; cell < cellCount_; ++cell) {
        const NeutronMaterial& mat = materials_[cellMaterials_[cell]];
        double fission = 0.0;
        for (int g = 0; g < groupCount_; ++g) {
            fission += mat.nuSigmaFission[g] * scalar[scalarIndex(g, cell)];
        }
        for (int g = 0; g < groupCount_; ++g) {
            result[scalarIndex(g, cell)] = mesh_.cells[cell].volume * mat.chi[g] * fission;
        }
    }
    return result;
}

std::vector<double> NeutronSolver::applyKOperator(const std::vector<double>& scalar) {
    // Applies K(x) = P * A^{-1} * F * x on scalar flux vector x (without modifying solver's converged state)
    const int n = groupCount_ * cellCount_;
    if (model_ == NeutronModel::Diffusion) {
        std::vector<double> rhs = applyFissionOperator(scalar);
        std::vector<double> sol(n, 0.0);
        const LinearSolveInfo info = system_->solve(rhs, sol, linearSolver_, preconditioner_,
            linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
        totalLinearIterations_ += info.iterations;
        lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);
        return sol;
    }

    // For SN: solve A * psi = F * scalar using zero initial guess for signed higher modes
    const std::vector<double> savedScalar = scalarFlux_;
    const std::vector<double> savedJx = currentX_;
    const std::vector<double> savedJy = currentY_;
    const std::vector<double> savedAngular = angularFlux_;
    const std::vector<double> savedOut = outgoingBoundaryAngular_;

    scalarFlux_.assign(n, 0.0);
    currentX_.assign(n, 0.0);
    currentY_.assign(n, 0.0);
    angularFlux_.assign(groupCount_ * directionCount_ * cellCount_, 0.0);
    std::fill(outgoingBoundaryAngular_.begin(), outgoingBoundaryAngular_.end(), 0.0);

    // Note: higher modes can have both positive and negative values; use linear upwind sweep for exact linearity
    const NeutronSnScheme savedScheme = snScheme_;
    snScheme_ = NeutronSnScheme::Upwind;
    solveSnInner(1.0, &scalar, false);
    snScheme_ = savedScheme;

    std::vector<double> result = scalarFlux_;
    scalarFlux_ = savedScalar;
    currentX_ = savedJx;
    currentY_ = savedJy;
    angularFlux_ = savedAngular;
    outgoingBoundaryAngular_ = savedOut;
    return result;
}

void NeutronSolver::solveHigherModes() {
    const int n = groupCount_ * cellCount_;
    const int requested = std::min(eigenModesRequested_, n);
    if (requested <= 1) return;

    const int blockSize = std::min(n, std::max(requested + 6, 12));

    // Domain bounding box for structured symmetry-breaking initial basis
    double xmin = mesh_.cells[0].centre.x, xmax = xmin;
    double ymin = mesh_.cells[0].centre.y, ymax = ymin;
    for (int c = 1; c < cellCount_; ++c) {
        xmin = std::min(xmin, mesh_.cells[c].centre.x);
        xmax = std::max(xmax, mesh_.cells[c].centre.x);
        ymin = std::min(ymin, mesh_.cells[c].centre.y);
        ymax = std::max(ymax, mesh_.cells[c].centre.y);
    }
    const double lx = std::max(xmax - xmin, 1.0e-12);
    const double ly = std::max(ymax - ymin, 1.0e-12);

    std::vector<std::vector<double>> basis;
    basis.reserve(blockSize);
    basis.push_back(scalarFlux_); // column 0 = converged fundamental mode

    // Deterministic trigonometric + pseudo-random basis spanning all symmetries (captures degenerate modes!)
    const std::vector<std::pair<int, int>> harmonics = {
        {1, 0}, {0, 1}, {1, 1}, {2, 0}, {0, 2}, {2, 1}, {1, 2}, {2, 2}, {3, 0}, {0, 3}, {3, 1}, {1, 3}
    };
    unsigned long long seed = 13579246801ULL;
    auto nextRand = [&seed]() {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((seed >> 33) & 0x7fffffffULL) / 2147483647.0 - 0.5;
    };

    for (int col = 1; col < blockSize; ++col) {
        std::vector<double> v(n, 0.0);
        const auto [mx, my] = harmonics[(col - 1) % harmonics.size()];
        for (int c = 0; c < cellCount_; ++c) {
            const double tx = (mesh_.cells[c].centre.x - xmin) / lx;
            const double ty = (mesh_.cells[c].centre.y - ymin) / ly;
            const double spatial = std::cos(pi * mx * tx) * std::cos(pi * my * ty);
            for (int g = 0; g < groupCount_; ++g) {
                const double gWeight = 1.0 + 0.25 * ((g + col) % 3);
                v[scalarIndex(g, c)] = spatial * gWeight + 0.05 * nextRand();
            }
        }
        basis.push_back(std::move(v));
    }
    orthonormalizeColumns(basis);

    std::vector<NeutronEigenmode> latestModes;
    for (int iter = 1; iter <= maximumModeIterations_; ++iter) {
        const int b = static_cast<int>(basis.size());
        std::vector<std::vector<double>> z(b);
        for (int j = 0; j < b; ++j) {
            z[j] = applyKOperator(basis[j]);
        }

        // Rayleigh-Ritz projection B = Q^T * Z (column-major for LAPACK)
        std::vector<double> matB(b * b, 0.0);
        for (int col = 0; col < b; ++col) {
            for (int row = 0; row < b; ++row) {
                matB[col * b + row] = dotVec(basis[row], z[col]);
            }
        }
        const SmallEigen eig = solveSmallEigen(b, matB);

        // Sort eigenvalues by descending real part (for real eigenvalues; fallback to magnitude)
        std::vector<int> order(b);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int c) {
            const double ra = (std::abs(eig.wi[a]) < 1.0e-8) ? eig.wr[a] : -1.0e30 + std::hypot(eig.wr[a], eig.wi[a]);
            const double rc = (std::abs(eig.wi[c]) < 1.0e-8) ? eig.wr[c] : -1.0e30 + std::hypot(eig.wr[c], eig.wi[c]);
            return ra > rc;
        });

        latestModes.clear();
        double maxResidual = 0.0;
        const double k0 = std::max(std::abs(eig.wr[order[0]]), 1.0e-30);

        for (int idx = 0; idx < requested && idx < b; ++idx) {
            const int col = order[idx];
            const double lam = eig.wr[col];
            const double imag = eig.wi[col];
            std::vector<double> ritzVec(n, 0.0);
            std::vector<double> kzVec(n, 0.0);
            for (int j = 0; j < b; ++j) {
                const double yj = eig.vr[col][j];
                for (int i = 0; i < n; ++i) {
                    ritzVec[i] += yj * basis[j][i];
                    kzVec[i] += yj * z[j][i];
                }
            }
            const double vNorm = std::max(normVec(ritzVec), 1.0e-30);
            double resSq = 0.0;
            for (int i = 0; i < n; ++i) {
                const double d = kzVec[i] - lam * ritzVec[i];
                resSq += d * d;
            }
            const double residual = std::sqrt(resSq) / std::max(std::abs(lam) * vNorm, 1.0e-30);
            maxResidual = std::max(maxResidual, residual);

            // Normalize mode cleanly
            if (idx == 0) {
                ritzVec = scalarFlux_;
            } else {
                double maxAbs = 0.0;
                double signRef = 1.0;
                for (int c = 0; c < cellCount_; ++c) {
                    double tot = 0.0;
                    for (int g = 0; g < groupCount_; ++g) tot += ritzVec[scalarIndex(g, c)];
                    if (std::abs(tot) > maxAbs) {
                        maxAbs = std::abs(tot);
                        signRef = (tot >= 0.0) ? 1.0 : -1.0;
                    }
                }
                if (maxAbs > 1.0e-30) {
                    const double scale = signRef / maxAbs;
                    for (double& x : ritzVec) x *= scale;
                }
            }

            NeutronEigenmode mode;
            mode.index = idx;
            mode.eigenvalue = lam;
            mode.imaginary = imag;
            mode.dominanceRatio = std::hypot(lam, imag) / k0;
            mode.residual = residual;
            mode.scalarFlux = std::move(ritzVec);
            latestModes.push_back(std::move(mode));
        }

        if (iter == 1 || iter % std::max(logInterval_, 5) == 0 || maxResidual < modeTolerance_) {
            std::cout << "Neutron subspace iteration " << iter << ": max mode residual " << maxResidual;
            for (int m = 0; m < std::min<int>(requested, 4); ++m) {
                std::cout << ", k" << m << "=" << std::setprecision(8) << latestModes[m].eigenvalue;
            }
            std::cout << '\n';
        }

        if (maxResidual < modeTolerance_) {
            modes_ = std::move(latestModes);
            return;
        }

        // Next block subspace = orthonormalized Z (with column 0 anchored to converged fundamental mode)
        z[0] = scalarFlux_;
        orthonormalizeColumns(z);
        basis = std::move(z);
    }
    modes_ = std::move(latestModes);
    std::cout << "Warning: Neutron higher-mode subspace iteration reached maximum_mode_iterations ("
              << maximumModeIterations_ << ")\n";
}

void NeutronSolver::solveAdjoint() {
    if (model_ == NeutronModel::Diffusion) {
        solveDiffusionAdjoint();
    } else {
        solveSnAdjoint();
    }
}

void NeutronSolver::solveDiffusionAdjoint() {
    const int n = groupCount_ * cellCount_;
    SparseMatrix transposeMatrix(n);
    for (int row = 0; row < n; ++row) {
        for (const auto& [col, val] : diffusionMatrix_.row(row)) {
            transposeMatrix.add(col, row, val);
        }
    }
    PetscLinearSystem adjSystem(n);
    adjSystem.setMatrix(transposeMatrix);

    adjointScalarFlux_.assign(n, 1.0);
    double kAdj = eigenvalue_;

    auto adjointFission = [&](const std::vector<double>& adj) {
        std::vector<double> rhs(n, 0.0);
        for (int cell = 0; cell < cellCount_; ++cell) {
            const NeutronMaterial& mat = materials_[cellMaterials_[cell]];
            double chiDotAdj = 0.0;
            for (int g = 0; g < groupCount_; ++g) {
                chiDotAdj += mat.chi[g] * adj[scalarIndex(g, cell)];
            }
            for (int g = 0; g < groupCount_; ++g) {
                rhs[scalarIndex(g, cell)] = mesh_.cells[cell].volume * mat.nuSigmaFission[g] * chiDotAdj;
            }
        }
        return rhs;
    };

    auto adjointScore = [&](const std::vector<double>& adj) {
        double score = 0.0;
        for (int cell = 0; cell < cellCount_; ++cell) {
            const NeutronMaterial& mat = materials_[cellMaterials_[cell]];
            for (int g = 0; g < groupCount_; ++g) {
                score += mesh_.cells[cell].volume * mat.chi[g] * adj[scalarIndex(g, cell)];
            }
        }
        return score;
    };

    for (int iter = 1; iter <= maximumIterations_; ++iter) {
        const double oldK = kAdj;
        const std::vector<double> prev = adjointScalarFlux_;
        const double prevScore = adjointScore(prev);
        std::vector<double> rhs = adjointFission(prev);
        for (double& v : rhs) v /= oldK;
        const LinearSolveInfo info = adjSystem.solve(rhs, adjointScalarFlux_, linearSolver_, preconditioner_,
            linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
        totalLinearIterations_ += info.iterations;
        const double newScore = adjointScore(adjointScalarFlux_);
        if (!(newScore > 0.0)) throw std::runtime_error("Adjoint diffusion iteration collapsed to zero");
        kAdj = oldK * newScore / prevScore;
        for (double& v : adjointScalarFlux_) v /= newScore;

        double diffSq = 0.0, normSq = 0.0;
        for (int i = 0; i < n; ++i) {
            const double d = adjointScalarFlux_[i] - prev[i] / prevScore;
            diffSq += d * d;
            normSq += adjointScalarFlux_[i] * adjointScalarFlux_[i];
        }
        const double relDiff = std::sqrt(diffSq / std::max(normSq, 1.0e-30));
        const double kDiff = std::abs(kAdj - oldK) / std::abs(kAdj);
        if (relDiff < tolerance_ && kDiff < eigenTolerance_) break;
    }

    // Normalize adjoint scalar flux to volume-averaged total = 1
    double integral = 0.0;
    double volume = 0.0;
    for (int c = 0; c < cellCount_; ++c) {
        volume += mesh_.cells[c].volume;
        for (int g = 0; g < groupCount_; ++g) integral += mesh_.cells[c].volume * adjointScalarFlux_[scalarIndex(g, c)];
    }
    if (integral > 0.0) {
        const double scale = volume / integral;
        for (double& v : adjointScalarFlux_) v *= scale;
    }

    kinetics_.adjointComputed = true;
    kinetics_.adjointEigenvalue = kAdj;
    kinetics_.adjointDiscrepancy = std::abs(kAdj - eigenvalue_) / std::max(std::abs(eigenvalue_), 1.0e-30);
}

void NeutronSolver::solveSnAdjoint() {
    // Exact algebraic transpose of the 1st-order upwind discrete SN operator:
    // Forward row (g, m, c):
    //   (T_{g,m} u_{g,m})_c - refl(u) - (V_c/4pi) sum_{g',n} w_n [S0_{g'->g} + 3 S1_{g'->g} (Om_m . Om_n)] u_{g',n,c}
    //   = (1/k) * (V_c * chi_g / 4pi) * sum_{g',n} w_n nuSf_{g'} u_{g',n,c}
    // Adjoint variable v_{g,m,c} satisfies A^T v = (1/k_adj) F^T v.
    // T_{g,m}^T is solved directly in EXACT REVERSE topological order!
    const std::size_t nAng = static_cast<std::size_t>(groupCount_) * directionCount_ * cellCount_;
    const std::size_t nInt = mesh_.internalFaceCount();
    const std::size_t bCount = mesh_.boundaryFaceCount();

    adjointAngularFlux_.assign(nAng, 0.0);
    for (int g = 0; g < groupCount_; ++g) {
        for (int m = 0; m < directionCount_; ++m) {
            for (int c = 0; c < cellCount_; ++c) {
                adjointAngularFlux_[angularIndex(g, m, c)] = mesh_.cells[c].volume * directions_[m].weight;
            }
        }
    }

    // Build reverse reflection lookup: for boundary face f and outgoing direction m*, which incoming m reflects to m*?
    std::vector<std::vector<int>> reverseReflect(bCount * directionCount_);
    for (std::size_t fIdx = nInt; fIdx < mesh_.faces.size(); ++fIdx) {
        const NeutronBoundaryCondition& cond = boundaries_[fIdx];
        const bool usesSpecular = (cond.kind == NeutronBoundaryKind::Reflective) ||
                                  (cond.kind == NeutronBoundaryKind::Albedo && cond.whiteFraction < 1.0 - 1.0e-12);
        if (!usesSpecular) continue;
        for (int m = 0; m < directionCount_; ++m) {
            const int refl = reflectionDirections_[reflectionIndex(static_cast<int>(fIdx), m)];
            if (refl >= 0) {
                reverseReflect[(fIdx - nInt) * directionCount_ + refl].push_back(m);
            }
        }
    }

    auto computeAdjointMoments = [&](const std::vector<double>& v,
                                     std::vector<double>& v0,
                                     std::vector<double>& vX,
                                     std::vector<double>& vY) {
        v0.assign(groupCount_ * cellCount_, 0.0);
        vX.assign(groupCount_ * cellCount_, 0.0);
        vY.assign(groupCount_ * cellCount_, 0.0);
        for (int g = 0; g < groupCount_; ++g) {
            for (int m = 0; m < directionCount_; ++m) {
                const NeutronDirection& omega = directions_[m];
                for (int c = 0; c < cellCount_; ++c) {
                    const double val = v[angularIndex(g, m, c)];
                    const int sIdx = scalarIndex(g, c);
                    v0[sIdx] += val;
                    vX[sIdx] += val * omega.x;
                    vY[sIdx] += val * omega.y;
                }
            }
        }
    };

    auto adjointFissionNorm = [&](const std::vector<double>& v0) {
        double score = 0.0;
        for (int c = 0; c < cellCount_; ++c) {
            const NeutronMaterial& mat = materials_[cellMaterials_[c]];
            for (int g = 0; g < groupCount_; ++g) {
                score += mat.chi[g] * v0[scalarIndex(g, c)];
            }
        }
        return score;
    };

    double kAdj = eigenvalue_;
    std::vector<double> v0, vX, vY;
    computeAdjointMoments(adjointAngularFlux_, v0, vX, vY);

    for (int outer = 1; outer <= maximumIterations_; ++outer) {
        const double oldK = kAdj;
        const std::vector<double> prevV = adjointAngularFlux_;
        const std::vector<double> fissionV0 = v0;
        const double prevScore = adjointFissionNorm(fissionV0);

        // Inner adjoint source iteration
        for (int inner = 1; inner <= maximumSourceIterations_; ++inner) {
            const std::vector<double> prevInnerV = adjointAngularFlux_;
            std::vector<double> curV0, curVX, curVY;
            computeAdjointMoments(prevInnerV, curV0, curVX, curVY);

            for (int g = groupCount_ - 1; g >= 0; --g) {
                // Precompute white-reflection adjoint incoming sum per boundary face for group g
                std::vector<double> whiteAdjSum(bCount, 0.0);
                for (std::size_t fIdx = nInt; fIdx < mesh_.faces.size(); ++fIdx) {
                    const NeutronBoundaryCondition& cond = boundaries_[fIdx];
                    const bool usesWhite = (cond.kind == NeutronBoundaryKind::WhiteReflective) ||
                                           (cond.kind == NeutronBoundaryKind::Albedo && cond.whiteFraction > 1.0e-12);
                    if (!usesWhite || faceIncomingWeightSum_[fIdx] <= 1.0e-15) continue;
                    const Face& face = mesh_.faces[fIdx];
                    const int owner = face.owner;
                    double acc = 0.0;
                    for (int m = 0; m < directionCount_; ++m) {
                        const NeutronDirection& omega = directions_[m];
                        const double flux = face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
                        if (flux < -1.0e-14) {
                            acc += (-flux) * prevInnerV[angularIndex(g, m, owner)];
                        }
                    }
                    whiteAdjSum[fIdx - nInt] = cond.value(g) * cond.whiteFraction * acc / faceIncomingWeightSum_[fIdx];
                }

                for (int m = 0; m < directionCount_; ++m) {
                    const NeutronDirection& omega = directions_[m];
                    if (!directionAcyclic_[m]) {
                        throw std::runtime_error("Cyclic mesh SN adjoint sweep not supported");
                    }
                    const int baseOffset = (g * directionCount_ + m) * cellCount_;
                    // Reverse topological order!
                    for (int idx = cellCount_ - 1; idx >= 0; --idx) {
                        const int cell = sweepOrder_[m][idx];
                        const Cell& cellGeom = mesh_.cells[cell];
                        const NeutronMaterial& mat = materials_[cellMaterials_[cell]];

                        double scatAdj = 0.0;
                        double fisAdj = 0.0;
                        for (int dst = 0; dst < groupCount_; ++dst) {
                            const int dIdx = scalarIndex(dst, cell);
                            scatAdj += mat.scattering[dst * groupCount_ + g] * curV0[dIdx];
                            if (hasP1Scattering_) {
                                const double s1 = mat.scatteringP1[dst * groupCount_ + g];
                                if (s1 != 0.0) {
                                    scatAdj += 3.0 * s1 * (omega.x * curVX[dIdx] + omega.y * curVY[dIdx]);
                                }
                            }
                            fisAdj += mat.chi[dst] * fissionV0[dIdx];
                        }
                        double rhs = (omega.weight / (4.0 * pi)) * cellGeom.volume *
                            (scatAdj + mat.nuSigmaFission[g] * fisAdj / oldK);

                        double diag = cellGeom.volume * (mat.total[g] + mat.diffusion[g] * axialBuckling_);
                        for (std::size_t local = 0; local < cellGeom.faces.size(); ++local) {
                            const int fIdx = cellGeom.faces[local];
                            const int sign = cellGeom.faceSigns[local];
                            const Face& face = mesh_.faces[fIdx];
                            const double flux = sign * face.area * (omega.x * face.normal.x + omega.y * face.normal.y);
                            if (flux >= 0.0) {
                                diag += flux;
                                if (!face.isBoundary()) {
                                    // Downwind neighbour N in forward sweep is upstream in adjoint sweep!
                                    const int nbr = (face.owner == cell) ? face.neighbour : face.owner;
                                    rhs += flux * adjointAngularFlux_[baseOffset + nbr];
                                } else {
                                    // Outgoing boundary face: receives transposed reflection from incoming directions
                                    const NeutronBoundaryCondition& cond = boundaries_[fIdx];
                                    const std::size_t bIdx = fIdx - nInt;
                                    if (cond.kind == NeutronBoundaryKind::WhiteReflective) {
                                        rhs += omega.weight * (omega.x * face.normal.x + omega.y * face.normal.y) *
                                               whiteAdjSum[bIdx];
                                    } else if (cond.kind == NeutronBoundaryKind::Reflective ||
                                               cond.kind == NeutronBoundaryKind::Albedo) {
                                        const double beta = cond.value(g);
                                        const double aw = cond.whiteFraction;
                                        if (aw < 1.0 - 1.0e-12) {
                                            for (int incDir : reverseReflect[bIdx * directionCount_ + m]) {
                                                const NeutronDirection& incOm = directions_[incDir];
                                                const double incFlux = face.area * (incOm.x * face.normal.x + incOm.y * face.normal.y);
                                                rhs += (-incFlux) * beta * (1.0 - aw) * prevInnerV[angularIndex(g, incDir, cell)];
                                            }
                                        }
                                        if (aw > 1.0e-12) {
                                            rhs += omega.weight * (omega.x * face.normal.x + omega.y * face.normal.y) *
                                                   whiteAdjSum[bIdx];
                                        }
                                    }
                                }
                            }
                        }
                        adjointAngularFlux_[baseOffset + cell] = rhs / diag;
                    }
                }
            }

            double maxDiff = 0.0, maxVal = 0.0;
            for (std::size_t i = 0; i < nAng; ++i) {
                maxDiff = std::max(maxDiff, std::abs(adjointAngularFlux_[i] - prevInnerV[i]));
                maxVal = std::max(maxVal, std::abs(adjointAngularFlux_[i]));
            }
            if (maxDiff / std::max(maxVal, 1.0e-30) < sourceTolerance_) break;
        }

        computeAdjointMoments(adjointAngularFlux_, v0, vX, vY);
        const double newScore = adjointFissionNorm(v0);
        if (!(newScore > 0.0)) throw std::runtime_error("Adjoint SN iteration collapsed to zero");
        kAdj = oldK * newScore / prevScore;
        for (double& x : adjointAngularFlux_) x /= newScore;
        for (double& x : v0) x /= newScore;

        double diffSq = 0.0, normSq = 0.0;
        for (std::size_t i = 0; i < nAng; ++i) {
            const double d = adjointAngularFlux_[i] - prevV[i] / prevScore;
            diffSq += d * d;
            normSq += adjointAngularFlux_[i] * adjointAngularFlux_[i];
        }
        const double relDiff = std::sqrt(diffSq / std::max(normSq, 1.0e-30));
        const double kDiff = std::abs(kAdj - oldK) / std::abs(kAdj);
        if (relDiff < tolerance_ && kDiff < eigenTolerance_) break;
    }

    // Physical adjoint scalar flux phi^dagger_{g,c} = (1 / V_c) * sum_m v_{g,m,c}, normalized to mean = 1
    adjointScalarFlux_.assign(groupCount_ * cellCount_, 0.0);
    double integral = 0.0;
    double volume = 0.0;
    for (int c = 0; c < cellCount_; ++c) {
        volume += mesh_.cells[c].volume;
        for (int g = 0; g < groupCount_; ++g) {
            double sum = 0.0;
            for (int m = 0; m < directionCount_; ++m) sum += adjointAngularFlux_[angularIndex(g, m, c)];
            const double physAdj = sum / mesh_.cells[c].volume;
            adjointScalarFlux_[scalarIndex(g, c)] = physAdj;
            integral += mesh_.cells[c].volume * physAdj;
        }
    }
    if (integral > 0.0) {
        const double scale = volume / integral;
        for (double& x : adjointScalarFlux_) x *= scale;
        for (double& x : adjointAngularFlux_) x *= scale;
    }

    kinetics_.adjointComputed = true;
    kinetics_.adjointEigenvalue = kAdj;
    kinetics_.adjointDiscrepancy = std::abs(kAdj - eigenvalue_) / std::max(std::abs(eigenvalue_), 1.0e-30);
}

void NeutronSolver::computeKineticsAndPerturbation() {
    if (!kinetics_.adjointComputed) return;

    // Adjoint-weighted denominator <Phi^dagger, F Phi>
    double denom = 0.0;
    if (model_ == NeutronModel::Diffusion) {
        const std::vector<double> fPhi = applyFissionOperator(scalarFlux_);
        denom = dotVec(adjointScalarFlux_, fPhi);
    } else {
        for (int c = 0; c < cellCount_; ++c) {
            const NeutronMaterial& mat = materials_[cellMaterials_[c]];
            double fis = 0.0;
            for (int g = 0; g < groupCount_; ++g) fis += mat.nuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
            for (int g = 0; g < groupCount_; ++g) {
                double sumV = 0.0;
                for (int m = 0; m < directionCount_; ++m) sumV += adjointAngularFlux_[angularIndex(g, m, c)];
                denom += (mesh_.cells[c].volume * mat.chi[g] / (4.0 * pi)) * fis * sumV;
            }
        }
    }
    if (std::abs(denom) <= 1.0e-30) return;

    // Kinetics parameters (Eq. 5.312): beta_eff and Lambda_eff
    if (!delayedBeta_.empty() || !neutronVelocity_.empty()) {
        kinetics_.kineticsComputed = true;
        if (!delayedBeta_.empty()) {
            const int nDel = static_cast<int>(delayedBeta_.size());
            kinetics_.betaEffectiveGroups.assign(nDel, 0.0);
            double betaTot = 0.0;
            for (int i = 0; i < nDel; ++i) {
                double num = 0.0;
                for (int c = 0; c < cellCount_; ++c) {
                    const NeutronMaterial& mat = materials_[cellMaterials_[c]];
                    double fis = 0.0;
                    for (int g = 0; g < groupCount_; ++g) fis += mat.nuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
                    for (int g = 0; g < groupCount_; ++g) {
                        const double chiD = delayedChi_.empty() ? mat.chi[g] : delayedChi_[i * groupCount_ + g];
                        if (model_ == NeutronModel::Diffusion) {
                            num += adjointScalarFlux_[scalarIndex(g, c)] * mesh_.cells[c].volume * chiD * delayedBeta_[i] * fis;
                        } else {
                            double sumV = 0.0;
                            for (int m = 0; m < directionCount_; ++m) sumV += adjointAngularFlux_[angularIndex(g, m, c)];
                            num += (mesh_.cells[c].volume * chiD / (4.0 * pi)) * delayedBeta_[i] * fis * sumV;
                        }
                    }
                }
                kinetics_.betaEffectiveGroups[i] = num / denom;
                betaTot += kinetics_.betaEffectiveGroups[i];
            }
            kinetics_.betaEffective = betaTot;
        }
        if (!neutronVelocity_.empty()) {
            double mvNum = 0.0;
            for (int c = 0; c < cellCount_; ++c) {
                const double vol = mesh_.cells[c].volume;
                for (int g = 0; g < groupCount_; ++g) {
                    if (model_ == NeutronModel::Diffusion) {
                        mvNum += adjointScalarFlux_[scalarIndex(g, c)] * (vol / neutronVelocity_[g]) * scalarFlux_[scalarIndex(g, c)];
                    } else {
                        for (int m = 0; m < directionCount_; ++m) {
                            mvNum += adjointAngularFlux_[angularIndex(g, m, c)] *
                                     (vol / neutronVelocity_[g]) * angularFlux_[angularIndex(g, m, c)];
                        }
                    }
                }
            }
            kinetics_.generationTime = mvNum / denom;
            kinetics_.promptLifetime = eigenvalue_ * kinetics_.generationTime;
        }
    }

    // First-order perturbation theory (Eq. 5.313): dk/k = <Phi^+, (dF - k dA) Phi> / <Phi^+, F Phi>
    if (perturbation_.enabled) {
        double num = 0.0;
        if (model_ == NeutronModel::Diffusion) {
            // Build perturbed cell-wise materials and exact discrete deltaA = A_pert - A_base
            std::vector<NeutronMaterial> cellPertMats;
            cellPertMats.reserve(cellCount_);
            std::vector<int> savedCellMats = cellMaterials_;
            for (int c = 0; c < cellCount_; ++c) {
                NeutronMaterial m = materials_[cellMaterials_[c]];
                if (cellMatchesFilter(c, perturbation_.regionName, perturbation_.materialName)) {
                    for (int g = 0; g < groupCount_; ++g) {
                        m.absorption[g] += perturbation_.deltaAbsorption[g];
                        m.nuSigmaFission[g] += perturbation_.deltaNuSigmaFission[g];
                        m.diffusion[g] += perturbation_.deltaDiffusion[g];
                        double outAdd = 0.0;
                        for (int dst = 0; dst < groupCount_; ++dst) {
                            const double ds = perturbation_.deltaScattering[dst * groupCount_ + g];
                            m.scattering[dst * groupCount_ + g] += ds;
                            outAdd += ds;
                        }
                        m.total[g] += perturbation_.deltaAbsorption[g] + outAdd;
                    }
                }
                cellPertMats.push_back(std::move(m));
            }
            // Temporarily assign 1-to-1 cell materials to assemble exact perturbed matrix A_pert
            std::vector<int> identity(cellCount_);
            std::iota(identity.begin(), identity.end(), 0);
            cellMaterials_ = identity;
            SparseMatrix pertMatrix;
            std::vector<double> dummyRhs;
            buildDiffusionOperatorFromMaterials(cellPertMats, 0.0, pertMatrix, dummyRhs);
            cellMaterials_ = savedCellMats;

            const int n = groupCount_ * cellCount_;
            for (int row = 0; row < n; ++row) {
                double aPertX = 0.0;
                for (const auto& [col, val] : pertMatrix.row(row)) aPertX += val * scalarFlux_[col];
                double aBaseX = 0.0;
                for (const auto& [col, val] : diffusionMatrix_.row(row)) aBaseX += val * scalarFlux_[col];
                num -= eigenvalue_ * adjointScalarFlux_[row] * (aPertX - aBaseX);
            }
            for (int c = 0; c < cellCount_; ++c) {
                if (!cellMatchesFilter(c, perturbation_.regionName, perturbation_.materialName)) continue;
                const NeutronMaterial& mat = materials_[cellMaterials_[c]];
                double dFis = 0.0;
                for (int g = 0; g < groupCount_; ++g) {
                    dFis += perturbation_.deltaNuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
                }
                for (int g = 0; g < groupCount_; ++g) {
                    num += adjointScalarFlux_[scalarIndex(g, c)] * mesh_.cells[c].volume * mat.chi[g] * dFis;
                }
            }
        } else {
            for (int c = 0; c < cellCount_; ++c) {
                if (!cellMatchesFilter(c, perturbation_.regionName, perturbation_.materialName)) continue;
                const NeutronMaterial& mat = materials_[cellMaterials_[c]];
                const double vol = mesh_.cells[c].volume;
                double dFis = 0.0;
                for (int g = 0; g < groupCount_; ++g) {
                    dFis += perturbation_.deltaNuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
                }
                for (int g = 0; g < groupCount_; ++g) {
                    double outAdd = 0.0;
                    double inAdd = 0.0;
                    for (int src = 0; src < groupCount_; ++src) {
                        outAdd += perturbation_.deltaScattering[src * groupCount_ + g];
                        inAdd += perturbation_.deltaScattering[g * groupCount_ + src] * scalarFlux_[scalarIndex(src, c)];
                    }
                    const double dTotal = perturbation_.deltaAbsorption[g] + outAdd;
                    double sumV = 0.0;
                    double vDotU = 0.0;
                    for (int m = 0; m < directionCount_; ++m) {
                        const double v = adjointAngularFlux_[angularIndex(g, m, c)];
                        sumV += v;
                        vDotU += v * angularFlux_[angularIndex(g, m, c)];
                    }
                    num += (vol * mat.chi[g] / (4.0 * pi)) * dFis * sumV;
                    num -= eigenvalue_ * (vol * dTotal * vDotU - (vol / (4.0 * pi)) * inAdd * sumV);
                }
            }
        }
        kinetics_.perturbationComputed = true;
        kinetics_.perturbationDeltaKOverK = num / denom;
        kinetics_.perturbationDeltaRho = kinetics_.perturbationDeltaKOverK / eigenvalue_;
    }
}

std::vector<NeutronMaterial> NeutronSolver::materialsAtTime(double time) const {
    std::vector<NeutronMaterial> cellMats;
    cellMats.reserve(cellCount_);
    for (int c = 0; c < cellCount_; ++c) {
        NeutronMaterial mat = materials_[cellMaterials_[c]];
        for (const NeutronTransientEvent& ev : transientEvents_) {
            if (!cellMatchesFilter(c, ev.regionName, ev.materialName)) continue;
            double frac = 0.0;
            if (!ev.isRamp) {
                frac = (time >= ev.startTime - 1.0e-12) ? 1.0 : 0.0;
            } else {
                if (time <= ev.startTime) frac = 0.0;
                else if (time >= ev.endTime) frac = 1.0;
                else frac = (time - ev.startTime) / (ev.endTime - ev.startTime);
            }
            if (frac == 0.0) continue;
            for (int g = 0; g < groupCount_; ++g) {
                mat.absorption[g] += frac * ev.deltaAbsorption[g];
                mat.nuSigmaFission[g] += frac * ev.deltaNuSigmaFission[g];
                mat.diffusion[g] += frac * ev.deltaDiffusion[g];
                double outAdd = 0.0;
                for (int dst = 0; dst < groupCount_; ++dst) {
                    const double ds = frac * ev.deltaScattering[dst * groupCount_ + g];
                    mat.scattering[dst * groupCount_ + g] += ds;
                    outAdd += ds;
                }
                mat.total[g] += frac * ev.deltaAbsorption[g] + outAdd;
            }
        }
        cellMats.push_back(std::move(mat));
    }
    return cellMats;
}

void NeutronSolver::solveTransient() {
    if (neutronVelocity_.empty()) {
        throw std::runtime_error("Transient neutron problem requires neutron_velocity");
    }
    if (model_ != NeutronModel::Diffusion) {
        throw std::runtime_error("Transient neutron kinetics currently uses the multigroup diffusion model");
    }

    // 1. Solve initial steady-state eigenvalue problem
    solveEigenvalue();
    const double initialPower = (targetPower_ > 0.0) ? targetPower_ : 1.0;
    const double currentPower = balance(eigenvalue_, true).power;
    if (currentPower <= 0.0) throw std::runtime_error("Cannot start transient from zero-power steady state");
    normalize(initialPower / currentPower);
    normalization_ = (targetPower_ > 0.0) ? "target_power" : "initial_power_1";

    // Scale nu*Sigma_f by 1/k_0 if critical_initial_state is enabled (standard TWIGL / LRA convention)
    if (criticalInitialState_) {
        for (NeutronMaterial& mat : materials_) {
            for (double& v : mat.nuSigmaFission) v /= eigenvalue_;
        }
    }

    const int nDel = static_cast<int>(delayedBeta_.size());
    const double betaTotal = std::accumulate(delayedBeta_.begin(), delayedBeta_.end(), 0.0);
    std::vector<double> precursors(nDel * cellCount_, 0.0);

    // Initialize steady-state precursors: C_i(0) = (beta_i / lambda_i) * sum_g nu*Sigma_f,g * phi_g(0)
    for (int c = 0; c < cellCount_; ++c) {
        const NeutronMaterial& mat = materials_[cellMaterials_[c]];
        double fisRate = 0.0;
        for (int g = 0; g < groupCount_; ++g) {
            fisRate += mat.nuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
        }
        for (int i = 0; i < nDel; ++i) {
            precursors[i * cellCount_ + c] = delayedBeta_[i] * fisRate / delayedLambda_[i];
        }
    }

    auto recordStep = [&](int step, double time, const std::vector<NeutronMaterial>& cellMats, int linIters) {
        double power = 0.0;
        double peakPower = 0.0;
        double fluxInt = 0.0;
        double precInt = 0.0;
        for (int c = 0; c < cellCount_; ++c) {
            const double vol = mesh_.cells[c].volume;
            double cellP = 0.0;
            for (int g = 0; g < groupCount_; ++g) {
                const double phi = scalarFlux_[scalarIndex(g, c)];
                fluxInt += vol * phi;
                cellP += cellMats[c].kappaSigmaFission[g] * phi;
            }
            power += vol * cellP;
            peakPower = std::max(peakPower, cellP);
            for (int i = 0; i < nDel; ++i) precInt += vol * precursors[i * cellCount_ + c];
        }
        transientHistory_.push_back({step, time, power, power / initialPower, peakPower, fluxInt, precInt, linIters});
    };

    transientHistory_.clear();
    const std::vector<int> savedCellMats = cellMaterials_;
    std::vector<int> identity(cellCount_);
    std::iota(identity.begin(), identity.end(), 0);

    // Evaluate t = 0 base state (before any t = 0 step discontinuity takes effect on n -> n+1)
    std::vector<NeutronMaterial> baseMats(cellCount_);
    for (int c = 0; c < cellCount_; ++c) baseMats[c] = materials_[savedCellMats[c]];
    recordStep(0, 0.0, baseMats, 0);

    const double theta = (timeScheme_ == NeutronTimeScheme::CrankNicolson) ? 0.5 : 1.0;
    const int nSteps = static_cast<int>(std::llround(endTime_ / timeStep_));
    const double dt = endTime_ / nSteps;
    const int n = groupCount_ * cellCount_;

    auto assembleOperatorA = [&](const std::vector<NeutronMaterial>& cellMats,
                                 SparseMatrix& matA,
                                 std::vector<double>& bRhs) {
        cellMaterials_ = identity;
        buildDiffusionOperatorFromMaterials(cellMats, 0.0, matA, bRhs);
        cellMaterials_ = savedCellMats;
    };

    SparseMatrix matA_n;
    std::vector<double> bRhs_n;
    // For step events at t = 0, use t = 0+ for the operator on [0, dt] in CN or BE
    std::vector<NeutronMaterial> mats_n = materialsAtTime(0.0 + 1.0e-15);
    assembleOperatorA(mats_n, matA_n, bRhs_n);

    for (int step = 1; step <= nSteps; ++step) {
        const double tNext = step * dt;
        std::vector<NeutronMaterial> mats_next = materialsAtTime(tNext);
        SparseMatrix matA_next;
        std::vector<double> bRhs_next;
        assembleOperatorA(mats_next, matA_next, bRhs_next);

        // Coupled implicit system:
        // [ M_v/dt + theta*(A^{n+1} - F_eff^{n+1}) ] phi^{n+1} = rhs
        SparseMatrix lhs(n);
        std::vector<double> rhs(n, 0.0);

        // Add theta * A^{n+1} to LHS
        for (int row = 0; row < n; ++row) {
            for (const auto& [col, val] : matA_next.row(row)) {
                lhs.add(row, col, theta * val);
            }
            rhs[row] += theta * bRhs_next[row];
        }

        // If Crank-Nicolson (theta < 1), subtract (1 - theta) * A^n * phi^n on RHS
        if (theta < 1.0) {
            for (int row = 0; row < n; ++row) {
                double ax = 0.0;
                for (const auto& [col, val] : matA_n.row(row)) ax += val * scalarFlux_[col];
                rhs[row] += (1.0 - theta) * (bRhs_n[row] - ax);
            }
        }

        // Add M_v/dt, prompt fission, and exact eliminated delayed neutron precursor coupling
        for (int c = 0; c < cellCount_; ++c) {
            const double vol = mesh_.cells[c].volume;
            const NeutronMaterial& mNext = mats_next[c];
            const NeutronMaterial& mPrev = mats_n[c];

            double fisPrev = 0.0;
            for (int g = 0; g < groupCount_; ++g) {
                fisPrev += mPrev.nuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
            }

            for (int g = 0; g < groupCount_; ++g) {
                const int row = scalarIndex(g, c);
                const double mv = vol / (neutronVelocity_[g] * dt);
                lhs.add(row, row, mv);
                rhs[row] += mv * scalarFlux_[row];

                // External source
                rhs[row] += vol * (theta * mNext.externalSource[g] + (1.0 - theta) * mPrev.externalSource[g]);

                // Explicit (1 - theta) prompt fission at step n
                if (theta < 1.0) {
                    rhs[row] += (1.0 - theta) * vol * mPrev.chi[g] * (1.0 - betaTotal) * fisPrev;
                }

                // Implicit prompt + delayed effective fission coupling at step n+1
                double delayedFactor = 0.0;
                double delayedHistory = 0.0;
                for (int i = 0; i < nDel; ++i) {
                    const double lam = delayedLambda_[i];
                    const double beta = delayedBeta_[i];
                    const double chiD = delayedChi_.empty() ? mNext.chi[g] : delayedChi_[i * groupCount_ + g];
                    const double chiDPrev = delayedChi_.empty() ? mPrev.chi[g] : delayedChi_[i * groupCount_ + g];
                    const double denom = 1.0 + theta * lam * dt;
                    const double cPrev = precursors[i * cellCount_ + c];

                    // C_i^{n+1} = a_i * C_i^n + b_i * P_{d,i}^n + c_i * P_{d,i}^{n+1}
                    const double a_i = (1.0 - (1.0 - theta) * lam * dt) / denom;
                    const double b_i = ((1.0 - theta) * dt) / denom;
                    const double c_i = (theta * dt) / denom;

                    delayedFactor += chiD * lam * c_i * beta;
                    delayedHistory += theta * chiD * lam * (a_i * cPrev + b_i * beta * fisPrev)
                                    + (1.0 - theta) * chiDPrev * lam * cPrev;
                }
                rhs[row] += vol * delayedHistory;

                for (int src = 0; src < groupCount_; ++src) {
                    const double effCoeff = theta * vol *
                        (mNext.chi[g] * (1.0 - betaTotal) + delayedFactor) * mNext.nuSigmaFission[src];
                    if (effCoeff != 0.0) lhs.add(row, scalarIndex(src, c), -effCoeff);
                }
            }
        }

        const std::vector<double> prevFlux = scalarFlux_;
        system_->setMatrix(lhs);
        const LinearSolveInfo info = system_->solve(rhs, scalarFlux_, linearSolver_, preconditioner_,
            linearTolerance_, linearAbsoluteTolerance_, maximumLinearIterations_);
        totalLinearIterations_ += info.iterations;
        lastLinearResidual_ = std::max(lastLinearResidual_, info.residual);

        // Update precursors exactly using converged phi^{n+1}
        for (int c = 0; c < cellCount_; ++c) {
            double fisPrev = 0.0;
            double fisNext = 0.0;
            for (int g = 0; g < groupCount_; ++g) {
                fisPrev += mats_n[c].nuSigmaFission[g] * prevFlux[scalarIndex(g, c)];
                fisNext += mats_next[c].nuSigmaFission[g] * scalarFlux_[scalarIndex(g, c)];
            }
            for (int i = 0; i < nDel; ++i) {
                const double lam = delayedLambda_[i];
                const double beta = delayedBeta_[i];
                const double denom = 1.0 + theta * lam * dt;
                const double cOld = precursors[i * cellCount_ + c];
                precursors[i * cellCount_ + c] =
                    ((1.0 - (1.0 - theta) * lam * dt) * cOld +
                     (1.0 - theta) * dt * beta * fisPrev +
                     theta * dt * beta * fisNext) / denom;
            }
        }

        if (step % transientOutputInterval_ == 0 || step == nSteps) {
            recordStep(step, tNext, mats_next, info.iterations);
        }
        if (step == 1 || step % logInterval_ == 0 || step == nSteps) {
            double stepPower = 0.0;
            for (int c = 0; c < cellCount_; ++c) {
                for (int g = 0; g < groupCount_; ++g) {
                    stepPower += mesh_.cells[c].volume * mats_next[c].kappaSigmaFission[g] *
                                 scalarFlux_[scalarIndex(g, c)];
                }
            }
            std::cout << "Neutron transient step " << step << " (t=" << std::setprecision(5) << tNext
                      << " s): relative power " << std::setprecision(8) << (stepPower / initialPower) << '\n';
        }

        mats_n = std::move(mats_next);
        matA_n = std::move(matA_next);
        bRhs_n = std::move(bRhs_next);
    }
}

} // namespace mp
