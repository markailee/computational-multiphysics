#include "acoustic/AcousticSolver.hpp"

#include "fvm/TecplotWriter.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace mp {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

double polygonArea(const std::vector<Vector2>& points) {
    double twiceArea = 0.0;
    for (std::size_t i = 0; i < points.size(); ++i) {
        twiceArea += cross(points[i], points[(i + 1) % points.size()]);
    }
    return 0.5 * twiceArea;
}

bool diagonalMatrix(const SparseMatrix& matrix) {
    for (int row = 0; row < matrix.size(); ++row) {
        for (const auto& [column, value] : matrix.row(row)) {
            if (column != row && value != 0.0) return false;
        }
    }
    return true;
}

} // namespace

AcousticSolver::AcousticSolver(const Mesh& mesh, const Config& config, std::filesystem::path caseDirectory)
    : mesh_(mesh), config_(config), caseDirectory_(std::move(caseDirectory)),
      boundaryConditions_(mesh, config), properties_(mesh.cells.size()),
      mass_(mesh.nodes.size(), 0.0), dualArea_(mesh.nodes.size(), 0.0),
      damping_(static_cast<int>(mesh.nodes.size())), stiffness_(static_cast<int>(mesh.nodes.size())),
      system_(static_cast<int>(mesh.nodes.size())), pressure_(mesh.nodes.size(), 0.0),
      pressureRate_(mesh.nodes.size(), 0.0), pressureAcceleration_(mesh.nodes.size(), 0.0),
      previousPressure_(mesh.nodes.size(), 0.0), prescribedKind_(mesh.nodes.size(), -1),
      prescribedValue_(mesh.nodes.size(), 0.0), sourceWeights_(mesh.nodes.size(), 0.0),
      dipoleWeights_(mesh.nodes.size()), geometryChecks_(mesh.cells.size()) {
    soundSpeed_ = config.number("sound_speed", 1.0);
    const double density = config.number("density", 1.0);
    meanVelocity_ = config.vector("mean_velocity", {});
    if (soundSpeed_ <= 0.0 || density <= 0.0) {
        throw std::runtime_error("Acoustic density and sound_speed must be positive");
    }
    for (CellProperties& property : properties_) {
        property.density = density;
        property.soundSpeed = soundSpeed_;
        property.mach = meanVelocity_ / soundSpeed_;
        if (normSquared(property.mach) >= 1.0) {
            throw std::runtime_error("The acoustic formulation requires a subsonic background: |M| < 1");
        }
    }

    cavityOrigin_ = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vector2 upper{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for (const Vector2& node : mesh.nodes) {
        cavityOrigin_.x = std::min(cavityOrigin_.x, node.x);
        cavityOrigin_.y = std::min(cavityOrigin_.y, node.y);
        upper.x = std::max(upper.x, node.x);
        upper.y = std::max(upper.y, node.y);
    }
    cavityLengths_ = upper - cavityOrigin_;
    if (cavityLengths_.x <= 0.0 || cavityLengths_.y <= 0.0) {
        throw std::runtime_error("Acoustic cases require a non-degenerate two-dimensional domain");
    }

    const std::string initial = lowercase(config.get("initial_condition", "zero"));
    if (initial == "zero" || initial == "constant") initialKind_ = InitialKind::Zero;
    else if (initial == "cavity_mode" || initial == "cavitymode") initialKind_ = InitialKind::CavityMode;
    else if (initial == "plane_wave" || initial == "planewave") initialKind_ = InitialKind::PlaneWave;
    else if (initial == "gaussian_pulse" || initial == "gaussianpulse") initialKind_ = InitialKind::GaussianPulse;
    else throw std::runtime_error("Unknown acoustic initial_condition: " + initial);

    initialAmplitude_ = config.number("initial_amplitude", 1.0);
    initialPressureOffset_ = config.number("initial_pressure", 0.0);
    initialUniformRate_ = config.number("initial_pressure_rate", 0.0);
    modeIndices_ = config.vector("mode_indices", {2.0, 1.0});
    if (modeIndices_.x < 0.0 || modeIndices_.y < 0.0 ||
        std::abs(modeIndices_.x - std::round(modeIndices_.x)) > 1.0e-10 ||
        std::abs(modeIndices_.y - std::round(modeIndices_.y)) > 1.0e-10 ||
        (modeIndices_.x == 0.0 && modeIndices_.y == 0.0)) {
        throw std::runtime_error("mode_indices must be non-negative integers and cannot both be zero");
    }

    waveDirection_ = config.vector("wave_direction", {1.0, 0.0});
    if (normSquared(waveDirection_) <= 0.0) throw std::runtime_error("wave_direction cannot be zero");
    waveDirection_ = normalized(waveDirection_);
    const double waveNumber = config.number("wave_number", 2.0 * pi);
    if (waveNumber <= 0.0) throw std::runtime_error("wave_number must be positive");
    waveVector_ = waveDirection_ * waveNumber;
    wavePhase_ = config.number("wave_phase", 0.0);
    waveFrequency_ = soundSpeed_ * waveNumber + dot(meanVelocity_, waveVector_);

    gaussianCenter_ = config.vector("gaussian_center", cavityOrigin_ + 0.5 * cavityLengths_);
    gaussianWidth_ = config.number("gaussian_width", 0.1);
    if (gaussianWidth_ <= 0.0) throw std::runtime_error("gaussian_width must be positive");
    if (initialKind_ == InitialKind::GaussianPulse && soundSpeed_ + dot(meanVelocity_, waveDirection_) <= 0.0) {
        throw std::runtime_error("The Gaussian pulse propagation speed must be positive");
    }

    timeScheme_ = lowercase(config.get("time_scheme", "newmark"));
    if (timeScheme_ == "central" || timeScheme_ == "explicit" || timeScheme_ == "centraldifference") {
        timeScheme_ = "central_difference";
    }
    if (timeScheme_ != "newmark" && timeScheme_ != "central_difference") {
        throw std::runtime_error("Acoustic time_scheme must be newmark or central_difference");
    }
    newmarkBeta_ = config.number("newmark_beta", 0.25);
    newmarkGamma_ = config.number("newmark_gamma", 0.5);
    if ((timeScheme_ == "newmark" && newmarkBeta_ <= 0.0) || newmarkGamma_ <= 0.0) {
        throw std::runtime_error("Newmark beta and gamma must be positive");
    }

    sourceAmplitude_ = config.number("source_amplitude", 0.0);
    sourceFrequency_ = config.number("source_frequency", 1.0);
    sourcePhase_ = config.number("source_phase", 0.0);
    sourcePosition_ = config.vector("source_position", cavityOrigin_ + 0.5 * cavityLengths_);
    sourceWidth_ = config.number("source_width", 0.1);
    sourceScale_ = config.number("source_scale", 1.0);
    sourceRampTime_ = config.number("source_ramp_time", 0.0);
    const std::string sourceType = lowercase(config.get("source_type", sourceAmplitude_ == 0.0 ? "none" : "harmonic"));
    if (sourceType == "none") sourceKind_ = SourceKind::None;
    else if (sourceType == "harmonic" || sourceType == "monopole") sourceKind_ = SourceKind::HarmonicMonopole;
    else if (sourceType == "flow_dipole" || sourceType == "flowdipole") sourceKind_ = SourceKind::FlowDipole;
    else throw std::runtime_error("Unknown acoustic source_type: " + sourceType);
    if (sourceRampTime_ < 0.0) throw std::runtime_error("source_ramp_time cannot be negative");
    if (sourceKind_ == SourceKind::HarmonicMonopole && (sourceFrequency_ <= 0.0 || sourceWidth_ <= 0.0)) {
        throw std::runtime_error("A harmonic source requires positive source_frequency and source_width");
    }
    if (sourceKind_ == SourceKind::FlowDipole) {
        if (sourceWidth_ <= 0.0) throw std::runtime_error("A flow dipole requires positive source_width");
        std::filesystem::path sourceFile = config.require("source_file");
        if (sourceFile.is_relative()) sourceFile = caseDirectory_ / sourceFile;
        loadFlowSource(sourceFile);
    }

    linearRelativeTolerance_ = config.number("linear_tolerance", 1.0e-10);
    linearAbsoluteTolerance_ = config.number("linear_absolute_tolerance", 1.0e-13);
    maximumLinearIterations_ = config.integer("maximum_linear_iterations", 1000);
    if (linearRelativeTolerance_ <= 0.0 || linearAbsoluteTolerance_ <= 0.0 || maximumLinearIterations_ <= 0) {
        throw std::runtime_error("Invalid acoustic linear-solver settings");
    }
    if (config.has("probe_position")) {
        const Vector2 probe = config.vector("probe_position", {});
        probeNode_ = static_cast<int>(std::min_element(mesh.nodes.begin(), mesh.nodes.end(),
            [&probe](const Vector2& a, const Vector2& b) {
                return normSquared(a - probe) < normSquared(b - probe);
            }) - mesh.nodes.begin());
    }
    if (config.has("probe_positions")) {
        std::istringstream input(config.require("probe_positions"));
        std::vector<double> coordinates;
        double value = 0.0;
        while (input >> value) coordinates.push_back(value);
        if (!input.eof() || coordinates.empty() || coordinates.size() % 2 != 0) {
            throw std::runtime_error("probe_positions requires x y pairs");
        }
        for (std::size_t index = 0; index < coordinates.size(); index += 2) {
            const Vector2 probe{coordinates[index], coordinates[index + 1]};
            probeNodes_.push_back(static_cast<int>(std::min_element(mesh.nodes.begin(), mesh.nodes.end(),
                [&probe](const Vector2& a, const Vector2& b) {
                    return normSquared(a - probe) < normSquared(b - probe);
                }) - mesh.nodes.begin()));
        }
    }
}

std::vector<Vector2> AcousticSolver::shapeGradients(int cellIndex, const Vector2& point) const {
    const Cell& cell = mesh_.cells[cellIndex];
    const std::size_t count = cell.nodes.size();
    std::vector<Vector2> gradients(count);
    if (count == 3) {
        const Vector2& a = mesh_.nodes[cell.nodes[0]];
        const Vector2& b = mesh_.nodes[cell.nodes[1]];
        const Vector2& c = mesh_.nodes[cell.nodes[2]];
        const double twiceArea = cross(b - a, c - a);
        gradients[0] = {b.y - c.y, c.x - b.x};
        gradients[1] = {c.y - a.y, a.x - c.x};
        gradients[2] = {a.y - b.y, b.x - a.x};
        for (Vector2& gradient : gradients) gradient /= twiceArea;
        return gradients;
    }
    if (count != 4) throw std::runtime_error("Acoustic elements must be linear triangles or quadrilaterals");

    constexpr std::array<double, 4> signs{-1.0, 1.0, 1.0, -1.0};
    double xi = 0.0;
    double eta = 0.0;
    double determinant = 0.0;
    for (int iteration = 0; iteration < 20; ++iteration) {
        Vector2 mapped;
        Vector2 xXi;
        Vector2 xEta;
        for (std::size_t i = 0; i < 4; ++i) {
            const double shape = 0.25 * (1.0 + signs[i] * xi) * (1.0 + signs[(i + 3) % 4] * eta);
            const double dXi = 0.25 * signs[i] * (1.0 + signs[(i + 3) % 4] * eta);
            const double dEta = 0.25 * signs[(i + 3) % 4] * (1.0 + signs[i] * xi);
            const Vector2& coordinate = mesh_.nodes[cell.nodes[i]];
            mapped += coordinate * shape;
            xXi += coordinate * dXi;
            xEta += coordinate * dEta;
        }
        const Vector2 residual = mapped - point;
        determinant = xXi.x * xEta.y - xEta.x * xXi.y;
        if (std::abs(determinant) <= 1.0e-15) {
            throw std::runtime_error("Singular quadrilateral Jacobian in acoustic inverse mapping");
        }
        const double deltaXi = (xEta.y * residual.x - xEta.x * residual.y) / determinant;
        const double deltaEta = (-xXi.y * residual.x + xXi.x * residual.y) / determinant;
        xi -= deltaXi;
        eta -= deltaEta;
        if (std::abs(deltaXi) + std::abs(deltaEta) < 1.0e-13) break;
        if (iteration == 19) throw std::runtime_error("Quadrilateral inverse mapping did not converge");
    }

    Vector2 xXi;
    Vector2 xEta;
    std::array<double, 4> dXi{};
    std::array<double, 4> dEta{};
    for (std::size_t i = 0; i < 4; ++i) {
        dXi[i] = 0.25 * signs[i] * (1.0 + signs[(i + 3) % 4] * eta);
        dEta[i] = 0.25 * signs[(i + 3) % 4] * (1.0 + signs[i] * xi);
        const Vector2& coordinate = mesh_.nodes[cell.nodes[i]];
        xXi += coordinate * dXi[i];
        xEta += coordinate * dEta[i];
    }
    determinant = xXi.x * xEta.y - xEta.x * xXi.y;
    if (std::abs(determinant) <= 1.0e-15) throw std::runtime_error("Singular quadrilateral Jacobian");
    for (std::size_t i = 0; i < 4; ++i) {
        gradients[i] = {(xEta.y * dXi[i] - xXi.y * dEta[i]) / determinant,
                        (-xEta.x * dXi[i] + xXi.x * dEta[i]) / determinant};
    }
    return gradients;
}

void AcousticSolver::assembleOperators() {
    const int nodeCount = static_cast<int>(mesh_.nodes.size());
    std::fill(mass_.begin(), mass_.end(), 0.0);
    std::fill(dualArea_.begin(), dualArea_.end(), 0.0);
    damping_ = SparseMatrix(nodeCount);
    stiffness_ = SparseMatrix(nodeCount);
    minimumEdge_ = std::numeric_limits<double>::max();

    for (std::size_t cellIndex = 0; cellIndex < mesh_.cells.size(); ++cellIndex) {
        const Cell& cell = mesh_.cells[cellIndex];
        const CellProperties& property = properties_[cellIndex];
        const std::size_t count = cell.nodes.size();
        const std::vector<Vector2> centreGradients = shapeGradients(static_cast<int>(cellIndex), cell.centre);
        Vector2 centreGradientSum;
        for (const Vector2& gradient : centreGradients) centreGradientSum += gradient;
        double gradientSum = norm(centreGradientSum);
        double subAreaSum = 0.0;
        for (std::size_t edge = 0; edge < count; ++edge) {
            minimumEdge_ = std::min(minimumEdge_, norm(mesh_.nodes[cell.nodes[(edge + 1) % count]] -
                                                       mesh_.nodes[cell.nodes[edge]]));
        }

        for (std::size_t a = 0; a < count; ++a) {
            const int row = cell.nodes[a];
            const Vector2& vertex = mesh_.nodes[row];
            const Vector2 previousMiddle = 0.5 * (vertex + mesh_.nodes[cell.nodes[(a + count - 1) % count]]);
            const Vector2 nextMiddle = 0.5 * (vertex + mesh_.nodes[cell.nodes[(a + 1) % count]]);
            const double subArea = polygonArea({vertex, nextMiddle, cell.centre, previousMiddle});
            if (subArea <= 0.0) throw std::runtime_error("Non-positive acoustic median-dual sub-area");
            subAreaSum += subArea;
            dualArea_[row] += subArea;
            mass_[row] += subArea / (property.density * property.soundSpeed * property.soundSpeed);
            for (std::size_t b = 0; b < count; ++b) {
                damping_.add(row, cell.nodes[b], 2.0 * subArea / (property.density * property.soundSpeed) *
                    dot(property.mach, centreGradients[b]));
            }

            const std::array<std::pair<Vector2, Vector2>, 2> segments{{
                {nextMiddle, cell.centre}, {cell.centre, previousMiddle}
            }};
            for (const auto& [start, end] : segments) {
                const Vector2 areaVector{end.y - start.y, start.x - end.x};
                const Vector2 middle = 0.5 * (start + end);
                const std::vector<Vector2> gradients = shapeGradients(static_cast<int>(cellIndex), middle);
                Vector2 sum;
                for (const Vector2& gradient : gradients) sum += gradient;
                gradientSum = std::max(gradientSum, norm(sum));
                for (std::size_t b = 0; b < count; ++b) {
                    const Vector2 tensorGradient = gradients[b] - property.mach * dot(property.mach, gradients[b]);
                    stiffness_.add(row, cell.nodes[b], -dot(tensorGradient, areaVector) / property.density);
                }
            }
        }
        geometryChecks_[cellIndex] = {cell.volume, subAreaSum, gradientSum};
        if (std::abs(subAreaSum - cell.volume) > 1.0e-10 * cell.volume) {
            throw std::runtime_error("Acoustic dual areas do not reproduce element area");
        }
        if (gradientSum > 1.0e-9 / std::sqrt(cell.volume)) {
            throw std::runtime_error("Acoustic shape functions violate the zero-gradient partition check");
        }
    }
    for (double value : mass_) {
        if (!(value > 0.0)) throw std::runtime_error("Acoustic lumped mass is not positive at every node");
    }

    std::fill(prescribedKind_.begin(), prescribedKind_.end(), -1);
    std::fill(prescribedValue_.begin(), prescribedValue_.end(), 0.0);
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const AcousticBoundaryCondition& condition = boundaryConditions_.face(static_cast<int>(faceIndex));
        const CellProperties& property = properties_[face.owner];
        const double normalMach = dot(property.mach, face.normal);
        if (AcousticBoundaryConditions::isDirichlet(condition.kind)) {
            const int kind = condition.kind == AcousticBoundaryKind::FixedValue ? 0 : 1;
            for (int node : {face.node0, face.node1}) {
                if (prescribedKind_[node] >= 0 && prescribedKind_[node] != kind) {
                    throw std::runtime_error("Conflicting prescribed acoustic pressures at a boundary node");
                }
                if (kind == 0 && prescribedKind_[node] == 0 &&
                    std::abs(prescribedValue_[node] - condition.value) > 1.0e-12) {
                    throw std::runtime_error("Conflicting fixed acoustic pressures at a boundary node");
                }
                prescribedKind_[node] = kind;
                prescribedValue_[node] = condition.value;
            }
            continue;
        }
        if (condition.kind == AcousticBoundaryKind::Rigid) continue;

        const double outgoingSpeed = property.soundSpeed * (1.0 + normalMach);
        if (outgoingSpeed <= 0.0 || 1.0 - normalMach * normalMach <= 0.0) {
            throw std::runtime_error("Invalid outgoing characteristic speed on an absorbing boundary");
        }
        const Vector2 tensorNormal = face.normal - property.mach * normalMach;
        const std::vector<Vector2> gradients = shapeGradients(face.owner, face.centre);
        for (int node : {face.node0, face.node1}) {
            const double halfLength = 0.5 * face.area;
            damping_.add(node, node, halfLength * (1.0 - normalMach * normalMach) /
                (property.density * outgoingSpeed));
            const auto& connectivity = mesh_.cells[face.owner].nodes;
            for (std::size_t b = 0; b < connectivity.size(); ++b) {
                const Vector2 tangentialGradient = gradients[b] - face.normal * dot(face.normal, gradients[b]);
                stiffness_.add(node, connectivity[b], -halfLength * dot(tensorNormal, tangentialGradient) /
                    property.density);
            }
        }
    }

    std::fill(sourceWeights_.begin(), sourceWeights_.end(), 0.0);
    std::fill(dipoleWeights_.begin(), dipoleWeights_.end(), Vector2{});
    if (sourceKind_ == SourceKind::HarmonicMonopole) {
        double total = 0.0;
        for (std::size_t node = 0; node < mesh_.nodes.size(); ++node) {
            const double distance2 = normSquared(mesh_.nodes[node] - sourcePosition_);
            sourceWeights_[node] = dualArea_[node] * std::exp(-distance2 / (2.0 * sourceWidth_ * sourceWidth_));
            total += sourceWeights_[node];
        }
        if (!(total > 0.0)) throw std::runtime_error("Acoustic source normalization failed");
        for (double& weight : sourceWeights_) weight /= total;
    } else if (sourceKind_ == SourceKind::FlowDipole) {
        double total = 0.0;
        std::vector<double> gaussian(mesh_.nodes.size(), 0.0);
        for (std::size_t node = 0; node < mesh_.nodes.size(); ++node) {
            gaussian[node] = std::exp(-normSquared(mesh_.nodes[node] - sourcePosition_) /
                (2.0 * sourceWidth_ * sourceWidth_));
            total += dualArea_[node] * gaussian[node];
        }
        if (!(total > 0.0)) throw std::runtime_error("Flow-dipole normalization failed");
        for (std::size_t node = 0; node < mesh_.nodes.size(); ++node) {
            const Vector2 gradient = (mesh_.nodes[node] - sourcePosition_) *
                (-gaussian[node] / (total * sourceWidth_ * sourceWidth_));
            dipoleWeights_[node] = gradient * (dualArea_[node] / properties_[0].density);
        }
    }
}

double AcousticSolver::planeWavePressure(const Vector2& point, double time) const {
    return initialAmplitude_ * std::sin(dot(waveVector_, point) - waveFrequency_ * time + wavePhase_);
}

double AcousticSolver::planeWaveRate(const Vector2& point, double time) const {
    return -initialAmplitude_ * waveFrequency_ * std::cos(dot(waveVector_, point) - waveFrequency_ * time + wavePhase_);
}

double AcousticSolver::initialPressure(const Vector2& point) const {
    switch (initialKind_) {
        case InitialKind::Zero:
            return initialPressureOffset_;
        case InitialKind::CavityMode:
            return initialAmplitude_ *
                std::cos(modeIndices_.x * pi * (point.x - cavityOrigin_.x) / cavityLengths_.x) *
                std::cos(modeIndices_.y * pi * (point.y - cavityOrigin_.y) / cavityLengths_.y);
        case InitialKind::PlaneWave:
            return planeWavePressure(point, 0.0);
        case InitialKind::GaussianPulse: {
            const double distance = dot(point - gaussianCenter_, waveDirection_);
            return initialAmplitude_ * std::exp(-distance * distance / (2.0 * gaussianWidth_ * gaussianWidth_));
        }
    }
    return 0.0;
}

double AcousticSolver::initialPressureRate(const Vector2& point) const {
    switch (initialKind_) {
        case InitialKind::Zero:
        case InitialKind::CavityMode:
            return initialUniformRate_;
        case InitialKind::PlaneWave:
            return planeWaveRate(point, 0.0);
        case InitialKind::GaussianPulse: {
            const double distance = dot(point - gaussianCenter_, waveDirection_);
            const double pressure = initialPressure(point);
            const double propagationSpeed = soundSpeed_ + dot(meanVelocity_, waveDirection_);
            return propagationSpeed * distance * pressure / (gaussianWidth_ * gaussianWidth_);
        }
    }
    return 0.0;
}

double AcousticSolver::prescribedPressure(int node, double time) const {
    return prescribedKind_[node] == 0 ? prescribedValue_[node] : planeWavePressure(mesh_.nodes[node], time);
}

double AcousticSolver::prescribedRate(int node, double time) const {
    return prescribedKind_[node] == 0 ? 0.0 : planeWaveRate(mesh_.nodes[node], time);
}

double AcousticSolver::prescribedAcceleration(int node, double time) const {
    if (prescribedKind_[node] == 0) return 0.0;
    return -waveFrequency_ * waveFrequency_ * planeWavePressure(mesh_.nodes[node], time);
}

void AcousticSolver::initialize() {
    if (initialized_) return;
    assembleOperators();
    for (std::size_t node = 0; node < mesh_.nodes.size(); ++node) {
        pressure_[node] = initialPressure(mesh_.nodes[node]);
        pressureRate_[node] = initialPressureRate(mesh_.nodes[node]);
        if (isPrescribed(static_cast<int>(node))) {
            pressure_[node] = prescribedPressure(static_cast<int>(node), 0.0);
            pressureRate_[node] = prescribedRate(static_cast<int>(node), 0.0);
        }
    }
    const std::vector<double> dampingRate = multiply(damping_, pressureRate_);
    const std::vector<double> stiffnessPressure = multiply(stiffness_, pressure_);
    const std::vector<double> forcing = source(0.0);
    for (std::size_t node = 0; node < mesh_.nodes.size(); ++node) {
        pressureAcceleration_[node] = isPrescribed(static_cast<int>(node))
            ? prescribedAcceleration(static_cast<int>(node), 0.0)
            : (forcing[node] - dampingRate[node] - stiffnessPressure[node]) / mass_[node];
    }
    previousPressure_ = pressure_;
    initialized_ = true;

    int triangles = 0;
    int quadrilaterals = 0;
    double areaError = 0.0;
    double gradientSum = 0.0;
    for (const GeometryCheck& check : geometryChecks_) {
        areaError = std::max(areaError, std::abs(check.subAreaSum - check.area));
        gradientSum = std::max(gradientSum, check.gradientSum);
    }
    for (const Cell& cell : mesh_.cells) {
        if (cell.nodes.size() == 3) ++triangles;
        else if (cell.nodes.size() == 4) ++quadrilaterals;
    }
    std::cout << "Acoustic mesh: " << mesh_.nodes.size() << " nodes, " << triangles << " triangles, "
              << quadrilaterals << " quadrilaterals\n"
              << "Geometry checks: max area error " << areaError << ", max |sum(grad N)| " << gradientSum << '\n';
}

std::vector<double> AcousticSolver::multiply(const SparseMatrix& matrix,
                                             const std::vector<double>& values) const {
    std::vector<double> result(values.size(), 0.0);
    for (int row = 0; row < matrix.size(); ++row) {
        for (const auto& [column, coefficient] : matrix.row(row)) result[row] += coefficient * values[column];
    }
    return result;
}

void AcousticSolver::loadFlowSource(const std::filesystem::path& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("Cannot open flow-dipole source file: " + file.string());
    const auto split = [](const std::string& line) {
        std::vector<std::string> fields;
        std::istringstream stream(line);
        std::string field;
        while (std::getline(stream, field, ',')) fields.push_back(field);
        return fields;
    };
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Empty flow-dipole source file");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::vector<std::string> header = split(line);
    int timeColumn = -1, xColumn = -1, yColumn = -1;
    for (std::size_t column = 0; column < header.size(); ++column) {
        if (header[column] == "time") timeColumn = static_cast<int>(column);
        else if (header[column] == "force_x") xColumn = static_cast<int>(column);
        else if (header[column] == "force_y") yColumn = static_cast<int>(column);
    }
    if (timeColumn < 0 || xColumn < 0 || yColumn < 0) {
        throw std::runtime_error("Flow-dipole CSV requires time, force_x, and force_y columns");
    }
    const int requiredColumns = std::max(timeColumn, std::max(xColumn, yColumn));
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const std::vector<std::string> fields = split(line);
        if (static_cast<int>(fields.size()) <= requiredColumns) {
            throw std::runtime_error("Incomplete row in flow-dipole source file");
        }
        FlowForceSample sample;
        try {
            sample.time = std::stod(fields[timeColumn]);
            sample.force = {std::stod(fields[xColumn]), std::stod(fields[yColumn])};
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid number in flow-dipole source file");
        }
        if (!flowForceHistory_.empty() && sample.time <= flowForceHistory_.back().time) {
            throw std::runtime_error("Flow-dipole source times must be strictly increasing");
        }
        flowForceHistory_.push_back(sample);
    }
    if (flowForceHistory_.size() < 2) throw std::runtime_error("Flow dipole requires at least two force samples");
}

Vector2 AcousticSolver::flowForce(double time) const {
    const double tolerance = 1.0e-10;
    if (time < flowForceHistory_.front().time - tolerance || time > flowForceHistory_.back().time + tolerance) {
        throw std::runtime_error("Acoustic time is outside the flow-force history");
    }
    if (time <= flowForceHistory_.front().time) return flowForceHistory_.front().force;
    if (time >= flowForceHistory_.back().time) return flowForceHistory_.back().force;
    const auto upper = std::lower_bound(flowForceHistory_.begin(), flowForceHistory_.end(), time,
        [](const FlowForceSample& sample, double value) { return sample.time < value; });
    const auto lower = upper - 1;
    const double fraction = (time - lower->time) / (upper->time - lower->time);
    return lower->force * (1.0 - fraction) + upper->force * fraction;
}

std::vector<double> AcousticSolver::source(double time) const {
    std::vector<double> result(mesh_.nodes.size(), 0.0);
    if (sourceKind_ == SourceKind::HarmonicMonopole) {
        const double temporal = sourceAmplitude_ * std::sin(2.0 * pi * sourceFrequency_ * time + sourcePhase_);
        for (std::size_t node = 0; node < result.size(); ++node) result[node] = temporal * sourceWeights_[node];
    } else if (sourceKind_ == SourceKind::FlowDipole) {
        const Vector2 force = flowForce(time) * sourceScale_;
        const double ramp = sourceRampTime_ > 0.0 && time < sourceRampTime_
            ? 0.5 * (1.0 - std::cos(pi * time / sourceRampTime_)) : 1.0;
        for (std::size_t node = 0; node < result.size(); ++node) {
            result[node] = ramp * dot(dipoleWeights_[node], force);
        }
    }
    return result;
}

void AcousticSolver::applyPrescribedPressure(SparseMatrix& matrix, std::vector<double>& rhs,
                                             const std::vector<double>& values) const {
    for (int node = 0; node < matrix.size(); ++node) {
        if (!isPrescribed(node)) continue;
        const double value = values[node];
        for (int row = 0; row < matrix.size(); ++row) {
            if (row == node || isPrescribed(row)) continue;
            const auto iterator = matrix.row(row).find(node);
            if (iterator != matrix.row(row).end()) rhs[row] -= iterator->second * value;
        }
        matrix.pinZero(node);
        rhs[node] = value;
    }
}

void AcousticSolver::advance(double dt) {
    if (!initialized_) initialize();
    if (!(dt > 0.0)) throw std::runtime_error("Acoustic time step must be positive");
    const int count = static_cast<int>(mesh_.nodes.size());
    const double targetTime = currentTime_ + dt;
    std::vector<double> nextPressure = pressure_;
    SparseMatrix effective(count);
    std::vector<double> rhs;
    double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0;

    if (timeScheme_ == "newmark") {
        a0 = 1.0 / (newmarkBeta_ * dt * dt);
        a1 = newmarkGamma_ / (newmarkBeta_ * dt);
        a2 = 1.0 / (newmarkBeta_ * dt);
        a3 = 1.0 / (2.0 * newmarkBeta_) - 1.0;
        a4 = newmarkGamma_ / newmarkBeta_ - 1.0;
        a5 = 0.5 * dt * (newmarkGamma_ / newmarkBeta_ - 2.0);
        effective = stiffness_;
        for (int row = 0; row < count; ++row) {
            effective.add(row, row, a0 * mass_[row]);
            for (const auto& [column, value] : damping_.row(row)) effective.add(row, column, a1 * value);
        }
        rhs = source(targetTime);
        std::vector<double> dampingHistory(count, 0.0);
        for (int node = 0; node < count; ++node) {
            rhs[node] += mass_[node] * (a0 * pressure_[node] + a2 * pressureRate_[node] +
                                        a3 * pressureAcceleration_[node]);
            dampingHistory[node] = a1 * pressure_[node] + a4 * pressureRate_[node] +
                                   a5 * pressureAcceleration_[node];
        }
        const std::vector<double> dampingContribution = multiply(damping_, dampingHistory);
        for (int node = 0; node < count; ++node) rhs[node] += dampingContribution[node];
    } else {
        if (!hasPreviousPressure_) {
            for (int node = 0; node < count; ++node) {
                previousPressure_[node] = pressure_[node] - dt * pressureRate_[node] +
                    0.5 * dt * dt * pressureAcceleration_[node];
                if (isPrescribed(node)) previousPressure_[node] = prescribedPressure(node, currentTime_ - dt);
            }
        }
        for (int row = 0; row < count; ++row) {
            effective.add(row, row, mass_[row] / (dt * dt));
            for (const auto& [column, value] : damping_.row(row)) effective.add(row, column, value / (2.0 * dt));
        }
        rhs = source(currentTime_);
        const std::vector<double> stiffnessPressure = multiply(stiffness_, pressure_);
        const std::vector<double> dampingPrevious = multiply(damping_, previousPressure_);
        for (int node = 0; node < count; ++node) {
            rhs[node] += -stiffnessPressure[node] + mass_[node] * (2.0 * pressure_[node] - previousPressure_[node]) /
                (dt * dt) + dampingPrevious[node] / (2.0 * dt);
        }
    }

    for (int node = 0; node < count; ++node) {
        if (isPrescribed(node)) nextPressure[node] = prescribedPressure(node, targetTime);
    }
    lastLinearIterations_ = 0;
    lastLinearResidual_ = 0.0;
    if (diagonalMatrix(effective)) {
        for (int node = 0; node < count; ++node) {
            if (isPrescribed(node)) continue;
            const double diagonal = effective.diagonal(node);
            if (std::abs(diagonal) <= 1.0e-30) throw std::runtime_error("Singular diagonal acoustic update");
            nextPressure[node] = rhs[node] / diagonal;
        }
    } else {
        applyPrescribedPressure(effective, rhs, nextPressure);
        system_.setMatrix(effective);
        const LinearSolveInfo info = system_.solve(rhs, nextPressure, config_.get("acoustic_ksp", "bcgs"),
            config_.get("acoustic_pc", "ilu"), linearRelativeTolerance_, linearAbsoluteTolerance_,
            maximumLinearIterations_);
        lastLinearIterations_ = info.iterations;
        lastLinearResidual_ = info.residual;
    }

    if (timeScheme_ == "newmark") {
        std::vector<double> nextRate(count, 0.0);
        std::vector<double> nextAcceleration(count, 0.0);
        for (int node = 0; node < count; ++node) {
            nextAcceleration[node] = a0 * (nextPressure[node] - pressure_[node]) - a2 * pressureRate_[node] -
                a3 * pressureAcceleration_[node];
            nextRate[node] = a1 * (nextPressure[node] - pressure_[node]) - a4 * pressureRate_[node] -
                a5 * pressureAcceleration_[node];
        }
        pressure_ = std::move(nextPressure);
        pressureRate_ = std::move(nextRate);
        pressureAcceleration_ = std::move(nextAcceleration);
    } else {
        std::vector<double> nextRate(count, 0.0);
        std::vector<double> nextAcceleration(count, 0.0);
        for (int node = 0; node < count; ++node) {
            nextAcceleration[node] = (nextPressure[node] - 2.0 * pressure_[node] + previousPressure_[node]) / (dt * dt);
            nextRate[node] = (nextPressure[node] - previousPressure_[node]) / (2.0 * dt);
        }
        previousPressure_ = pressure_;
        pressure_ = std::move(nextPressure);
        pressureRate_ = std::move(nextRate);
        pressureAcceleration_ = std::move(nextAcceleration);
        hasPreviousPressure_ = true;
    }
    currentTime_ = targetTime;
    for (int node = 0; node < count; ++node) {
        if (isPrescribed(node)) {
            pressureRate_[node] = prescribedRate(node, currentTime_);
            pressureAcceleration_[node] = prescribedAcceleration(node, currentTime_);
        }
        if (!std::isfinite(pressure_[node]) || !std::isfinite(pressureRate_[node]) ||
            !std::isfinite(pressureAcceleration_[node])) {
            throw std::runtime_error("Non-finite acoustic solution at time " + std::to_string(currentTime_));
        }
    }
}

double AcousticSolver::energy() const {
    const std::vector<double> stiffnessPressure = multiply(stiffness_, pressure_);
    double result = 0.0;
    for (std::size_t node = 0; node < pressure_.size(); ++node) {
        result += 0.5 * (mass_[node] * pressureRate_[node] * pressureRate_[node] +
                         pressure_[node] * stiffnessPressure[node]);
    }
    return result;
}

double AcousticSolver::pressureScale() const {
    double result = 0.0;
    for (double value : pressure_) result = std::max(result, std::abs(value));
    return result;
}

std::filesystem::path AcousticSolver::outputFile() const {
    std::filesystem::path output = config_.get("output", "results/acoustic.plt");
    return output.is_absolute() ? output : caseDirectory_ / output;
}

void AcousticSolver::write(const std::filesystem::path& file) const {
    writeNodalTecplot(file, mesh_, {{"p", &pressure_}, {"p_t", &pressureRate_},
        {"p_tt", &pressureAcceleration_}, {"dual_area", &dualArea_}}, name(), currentTime_);
}

void AcousticSolver::writeSummary(const std::filesystem::path& file, bool append, int linearIterations,
                                  double linearResidual, double courant) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file, append ? std::ios::app : std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write acoustic summary: " + file.string());
    output << std::setprecision(12) << std::scientific;
    if (!append) output << "time,min_pressure,max_pressure,energy,linear_iterations,linear_residual,courant\n";
    const auto [minimum, maximum] = std::minmax_element(pressure_.begin(), pressure_.end());
    output << currentTime_ << ',' << *minimum << ',' << *maximum << ',' << energy() << ','
           << linearIterations << ',' << linearResidual << ',' << courant << '\n';
}

void AcousticSolver::writeProbe(const std::filesystem::path& file, bool append) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file, append ? std::ios::app : std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write acoustic probe: " + file.string());
    output << std::setprecision(12) << std::scientific;
    if (!append) output << "time,p,p_t,p_tt\n";
    output << currentTime_ << ',' << pressure_[probeNode_] << ',' << pressureRate_[probeNode_] << ','
           << pressureAcceleration_[probeNode_] << '\n';
}

void AcousticSolver::writeProbeArray(const std::filesystem::path& file, bool append) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file, append ? std::ios::app : std::ios::trunc);
    if (!output) throw std::runtime_error("Cannot write acoustic probe array: " + file.string());
    output << std::setprecision(12) << std::scientific;
    if (!append) output << "time,probe,x,y,p,p_t,p_tt\n";
    for (std::size_t probe = 0; probe < probeNodes_.size(); ++probe) {
        const int node = probeNodes_[probe];
        output << currentTime_ << ',' << probe << ',' << mesh_.nodes[node].x << ',' << mesh_.nodes[node].y << ','
               << pressure_[node] << ',' << pressureRate_[node] << ',' << pressureAcceleration_[node] << '\n';
    }
}

void AcousticSolver::writeGeometry(const std::filesystem::path& file) const {
    std::filesystem::create_directories(file.parent_path());
    std::ofstream output(file);
    if (!output) throw std::runtime_error("Cannot write acoustic geometry diagnostics: " + file.string());
    output << std::setprecision(12) << std::scientific << "cell,area,subarea_sum,area_error,gradient_sum\n";
    for (std::size_t cell = 0; cell < geometryChecks_.size(); ++cell) {
        const GeometryCheck& check = geometryChecks_[cell];
        output << cell << ',' << check.area << ',' << check.subAreaSum << ','
               << std::abs(check.area - check.subAreaSum) << ',' << check.gradientSum << '\n';
    }
}

void AcousticSolver::run() {
    initialize();
    const double endTime = config_.number("end_time", sourceKind_ == SourceKind::FlowDipole
        ? flowForceHistory_.back().time : 1.0);
    if (endTime <= 0.0) throw std::runtime_error("Acoustic end_time must be positive");
    double dt = config_.number("dt", 0.0);
    const double propagationSpeed = soundSpeed_ + norm(meanVelocity_);
    if (dt <= 0.0) {
        const double cfl = config_.number("cfl", 0.4);
        if (cfl <= 0.0) throw std::runtime_error("Acoustic cfl must be positive");
        dt = cfl * minimumEdge_ / propagationSpeed;
    }
    const int maximumSteps = config_.integer("maximum_steps", 1000000);
    const int steps = static_cast<int>(std::ceil(endTime / dt - 1.0e-12));
    if (steps <= 0 || steps > maximumSteps) throw std::runtime_error("Acoustic step count exceeds maximum_steps");
    dt = endTime / static_cast<double>(steps);
    const double courant = propagationSpeed * dt / minimumEdge_;
    const int outputInterval = config_.integer("output_interval", steps);
    const int logInterval = config_.integer("log_interval", 100);
    const std::filesystem::path output = outputFile();
    const std::filesystem::path summary = output.parent_path() / (output.stem().string() + "_summary.csv");
    const std::filesystem::path probe = output.parent_path() / (output.stem().string() + "_probe.csv");
    const std::filesystem::path probeArray = output.parent_path() / (output.stem().string() + "_probes.csv");
    const std::filesystem::path geometry = output.parent_path() / (output.stem().string() + "_geometry.csv");
    writeGeometry(geometry);
    writeSummary(summary, false, 0, 0.0, courant);
    if (probeNode_ >= 0) writeProbe(probe, false);
    if (!probeNodes_.empty()) writeProbeArray(probeArray, false);
    std::cout << "Acoustic time scheme: " << timeScheme_ << ", steps " << steps << ", dt " << dt
              << ", Courant " << courant << '\n';

    for (int step = 1; step <= steps; ++step) {
        advance(dt);
        writeSummary(summary, true, lastLinearIterations_, lastLinearResidual_, courant);
        if (probeNode_ >= 0) writeProbe(probe, true);
        if (!probeNodes_.empty()) writeProbeArray(probeArray, true);
        if ((outputInterval > 0 && step % outputInterval == 0) || step == steps) {
            const std::filesystem::path stepFile = output.parent_path() /
                (output.stem().string() + "_" + std::to_string(step) + output.extension().string());
            write(stepFile);
        }
        if (step == 1 || step == steps || (logInterval > 0 && step % logInterval == 0)) {
            std::cout << "step " << step << ", time " << currentTime_ << ", max|p| " << pressureScale()
                      << ", energy " << energy() << ", linear " << lastLinearIterations_ << '\n';
        }
    }
    std::cout << "Tecplot output: " << output.parent_path() << '\n';
}

} // namespace mp
