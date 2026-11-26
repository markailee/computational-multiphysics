#include "flow/MultiSpeciesModel.hpp"

#include "fvm/Operators.hpp"
#include "fvm/ScalarTransport.hpp"
#include "material/MaterialParsing.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace mp {

MultiSpeciesModel::MultiSpeciesModel(const Mesh& mesh, const Config& config,
                                     std::shared_ptr<MaterialModel> material)
    : mesh_(mesh), config_(config), material_(std::move(material)),
      temperature_(config.number("temperature", 300.0)),
      thermodynamicPressure_(config.number("thermo_pressure", 101325.0)),
      relaxation_(config.number("species_relaxation", 0.7)),
      turbulentSchmidt_(config.number("turbulent_schmidt", 0.7)),
      linearRelativeTolerance_(config.number("linear_tolerance", 1.0e-9)),
      linearAbsoluteTolerance_(config.number("linear_absolute_tolerance", 1.0e-12)),
      maximumLinearIterations_(config.integer("maximum_linear_iterations", 1000)),
      timeDerivative_(config.get("time_scheme", "steady")),
      densityField_("rho_mix", mesh), viscosityField_("mu_mix", mesh),
      conductivityField_("lambda_mix", mesh), heatCapacityField_("cp_mix", mesh),
      source_(mesh.cells.size(), 0.0), destruction_(mesh.cells.size(), 0.0),
      constrainedCells_(mesh.cells.size(), 0), constrainedValues_(mesh.cells.size(), 0.0),
      transportSystem_(static_cast<int>(mesh.cells.size())) {
    const auto* definitions = material_->speciesDefinitions();
    if (!definitions || definitions->empty()) {
        throw std::runtime_error("species_model = transport requires an ideal-gas material with species definitions");
    }
    speciesCount_ = definitions->size();
    if (speciesCount_ < 2) throw std::runtime_error("Multispecies transport requires at least two species");
    if (!(temperature_ > 0.0) || !(thermodynamicPressure_ > 0.0)) {
        throw std::runtime_error("temperature and thermo_pressure must be positive");
    }
    if (relaxation_ <= 0.0 || relaxation_ > 1.0) throw std::runtime_error("species_relaxation must be in (0, 1]");
    if (!(turbulentSchmidt_ > 0.0)) throw std::runtime_error("turbulent_schmidt must be positive");

    const std::vector<double> initial = material::readFractions(config, "initial_mass_fractions", speciesCount_);
    fields_.reserve(speciesCount_);
    previousFields_.reserve(speciesCount_);
    oldFields_.reserve(speciesCount_);
    oldOldFields_.reserve(speciesCount_);
    molecularDiffusivityFields_.reserve(speciesCount_);
    for (std::size_t k = 0; k < speciesCount_; ++k) {
        fields_.emplace_back("Y_" + (*definitions)[k].name, mesh, initial[k]);
        previousFields_.emplace_back("previousY_" + (*definitions)[k].name, mesh, initial[k]);
        oldFields_.emplace_back("oldY_" + (*definitions)[k].name, mesh, initial[k]);
        oldOldFields_.emplace_back("oldOldY_" + (*definitions)[k].name, mesh, initial[k]);
        molecularDiffusivityFields_.emplace_back("D_" + (*definitions)[k].name, mesh,
                                                 (*definitions)[k].diffusivity);
    }
    diffusivities_.assign(speciesCount_, std::vector<double>(mesh.cells.size(), 0.0));
    correctionFluxes_.assign(speciesCount_, std::vector<double>(mesh.faces.size(), 0.0));
    gradients_.resize(speciesCount_);
    for (const BoundaryPatch& patch : mesh.patches) boundaries_.push_back(readBoundary(patch.name));
    updateBoundaries();
    updateProperties();
    updateDiagnostics();
}

SpeciesBoundaryCondition MultiSpeciesModel::readBoundary(const std::string& patch) const {
    const std::string key = "boundary." + patch + ".mass_fractions";
    std::istringstream input(config_.require(key));
    std::string typeName;
    input >> typeName;
    const std::string type = lowercase(typeName);
    SpeciesBoundaryCondition condition;
    if (type == "fixedvalue") {
        condition.fixedValue = true;
        double value = 0.0;
        while (input >> value) condition.massFractions.push_back(value);
        if (!input.eof()) throw std::runtime_error("Invalid mass-fraction boundary data for " + key);
        condition.massFractions = material::fractions(condition.massFractions, speciesCount_, key);
    } else if (type == "zerogradient") {
        condition.massFractions.assign(speciesCount_, 0.0);
    } else {
        throw std::runtime_error("Unknown mass_fractions boundary condition for " + key);
    }
    return condition;
}

void MultiSpeciesModel::initialize() {
    if (initialized_) return;
    for (std::size_t k = 0; k < speciesCount_; ++k) {
        previousFields_[k].values() = fields_[k].values();
        oldFields_[k].values() = fields_[k].values();
        oldOldFields_[k].values() = fields_[k].values();
    }
    updateBoundaries();
    updateProperties();
    initialized_ = true;
}

void MultiSpeciesModel::updateFieldBoundaries(std::size_t k) {
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const SpeciesBoundaryCondition& condition = boundaries_[face.patch];
        fields_[k].boundaryDirichlet(faceIndex) = condition.fixedValue;
        fields_[k].boundary(faceIndex) = condition.fixedValue
            ? condition.massFractions[k] : fields_[k][face.owner];
    }
}

void MultiSpeciesModel::updateBackgroundSpecies() {
    const std::size_t background = speciesCount_ - 1;
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        double sum = 0.0;
        for (std::size_t k = 0; k < background; ++k) sum += fields_[k][cell];
        fields_[background][cell] = 1.0 - sum;
    }
    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const bool fixedValue = boundaries_[face.patch].fixedValue;
        double boundarySum = 0.0;
        if (fixedValue) {
            for (std::size_t k = 0; k < background; ++k) boundarySum += fields_[k].boundary(faceIndex);
        }
        fields_[background].boundaryDirichlet(faceIndex) = fixedValue;
        fields_[background].boundary(faceIndex) = fixedValue
            ? 1.0 - boundarySum : fields_[background][face.owner];
    }
}

void MultiSpeciesModel::updateBoundaries() {
    for (std::size_t k = 0; k + 1 < speciesCount_; ++k) updateFieldBoundaries(k);
    updateBackgroundSpecies();
}

void MultiSpeciesModel::beginIteration() {
    for (std::size_t k = 0; k < speciesCount_; ++k) previousFields_[k].values() = fields_[k].values();
}

void MultiSpeciesModel::updateProperties() {
    MaterialState state;
    state.temperature = temperature_;
    state.pressure = thermodynamicPressure_;
    state.massFractions.resize(speciesCount_);
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        for (std::size_t k = 0; k < speciesCount_; ++k) state.massFractions[k] = fields_[k][cell];
        const MaterialProperties properties = material_->properties(state);
        if (properties.speciesDiffusivities.size() != speciesCount_) {
            throw std::runtime_error("Material returned an inconsistent species-diffusivity count");
        }
        densityField_[cell] = properties.density;
        viscosityField_[cell] = properties.viscosity;
        conductivityField_[cell] = properties.thermalConductivity;
        heatCapacityField_[cell] = properties.specificHeatCapacity;
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            molecularDiffusivityFields_[k][cell] = properties.speciesDiffusivities[k];
        }
    }
    for (std::size_t face = mesh_.internalFaceCount(); face < mesh_.faces.size(); ++face) {
        const int owner = mesh_.faces[face].owner;
        densityField_.boundary(face) = densityField_[owner];
        viscosityField_.boundary(face) = viscosityField_[owner];
        conductivityField_.boundary(face) = conductivityField_[owner];
        heatCapacityField_.boundary(face) = heatCapacityField_[owner];
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            molecularDiffusivityFields_[k].boundary(face) = molecularDiffusivityFields_[k][owner];
        }
    }
}

void MultiSpeciesModel::computeTransportDiffusivities(double carrierDensity,
                                                       const ScalarField* turbulentViscosity) {
    if (!(carrierDensity > 0.0)) throw std::runtime_error("Species carrier density must be positive");
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        const double turbulent = turbulentViscosity
            ? (*turbulentViscosity)[cell] / turbulentSchmidt_ : 0.0;
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            diffusivities_[k][cell] = carrierDensity * molecularDiffusivityFields_[k][cell] + turbulent;
        }
    }
}

void MultiSpeciesModel::computeCorrectionFluxes(const GradientScheme& gradient) {
    for (std::size_t k = 0; k < speciesCount_; ++k) gradients_[k] = gradient.scalar(fields_[k]);
    for (auto& fluxes : correctionFluxes_) std::fill(fluxes.begin(), fluxes.end(), 0.0);

    for (std::size_t faceIndex = 0; faceIndex < mesh_.internalFaceCount(); ++faceIndex) {
        const Face& face = mesh_.faces[faceIndex];
        const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
        double uncorrectedSum = 0.0;
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            const double diffusivity = face.ownerWeight * diffusivities_[k][face.owner] +
                (1.0 - face.ownerWeight) * diffusivities_[k][face.neighbour];
            const ScalarGradient faceGradient = face.ownerWeight * gradients_[k][face.owner] +
                (1.0 - face.ownerWeight) * gradients_[k][face.neighbour];
            uncorrectedSum += diffusivity * (face.area / face.delta *
                (fields_[k][face.owner] - fields_[k][face.neighbour]) - dot(correctionVector, faceGradient));
        }
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            const double faceFraction = face.ownerWeight * fields_[k][face.owner] +
                (1.0 - face.ownerWeight) * fields_[k][face.neighbour];
            correctionFluxes_[k][faceIndex] = -faceFraction * uncorrectedSum;
        }
    }

    for (std::size_t faceIndex = mesh_.internalFaceCount(); faceIndex < mesh_.faces.size(); ++faceIndex) {
        if (!fields_[0].boundaryDirichlet(faceIndex)) continue;
        const Face& face = mesh_.faces[faceIndex];
        const Vector2 correctionVector = nonOrthogonalVector(mesh_, faceIndex);
        double uncorrectedSum = 0.0;
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            uncorrectedSum += diffusivities_[k][face.owner] * (face.area / face.delta *
                (fields_[k][face.owner] - fields_[k].boundary(faceIndex)) -
                dot(correctionVector, gradients_[k][face.owner]));
        }
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            correctionFluxes_[k][faceIndex] = -fields_[k].boundary(faceIndex) * uncorrectedSum;
        }
    }
}

void MultiSpeciesModel::solve(double dt, bool hasOldOld, const std::vector<double>& massFlux,
                              const GradientScheme& gradient, double carrierDensity,
                              const ScalarField* turbulentViscosity) {
    initialize();
    if (compositionError_ < 1.0e-7 && bounded(1.0e-8)) updateProperties();
    computeTransportDiffusivities(carrierDensity, turbulentViscosity);
    updateBoundaries();
    computeCorrectionFluxes(gradient);
    const TimeCoefficients time = timeDerivative_.coefficients(dt, hasOldOld);
    linearIterations_ = 0;
    linearResidual_ = 0.0;
    linearResiduals_.clear();
    for (std::size_t k = 0; k + 1 < speciesCount_; ++k) {
        ScalarTransportEquation equation{fields_[k], previousFields_[k], oldFields_[k], oldOldFields_[k],
            massFlux, diffusivities_[k], source_, destruction_, constrainedCells_, constrainedValues_,
            carrierDensity, time, relaxation_, linearRelativeTolerance_,
            linearAbsoluteTolerance_, maximumLinearIterations_, config_.get("species_ksp", "bcgs"),
            config_.get("species_pc", "ilu")};
        equation.barthJespersen = useBarthJespersenLimiter(config_);
        equation.additionalFaceFlux = &correctionFluxes_[k];
        const LinearSolveInfo info = solveScalarTransport(mesh_, gradient, transportSystem_, equation);
        linearIterations_ += info.iterations;
        linearResidual_ = std::max(linearResidual_, info.residual);
        linearResiduals_.push_back(info.residual);
        for (double value : fields_[k].values()) {
            if (!std::isfinite(value)) throw std::runtime_error("Non-finite value in species field " + fields_[k].name());
        }
        updateFieldBoundaries(k);
    }
    updateBackgroundSpecies();
    updateDiagnostics();
    if (compositionError_ < 1.0e-8 && bounded(1.0e-8)) updateProperties();
}

void MultiSpeciesModel::updateDiagnostics() {
    change_ = 0.0;
    compositionError_ = 0.0;
    minimumMassFraction_ = std::numeric_limits<double>::max();
    maximumMassFraction_ = std::numeric_limits<double>::lowest();
    for (std::size_t cell = 0; cell < mesh_.cells.size(); ++cell) {
        double sum = 0.0;
        for (std::size_t k = 0; k < speciesCount_; ++k) {
            sum += fields_[k][cell];
            change_ = std::max(change_, std::abs(fields_[k][cell] - previousFields_[k][cell]));
            minimumMassFraction_ = std::min(minimumMassFraction_, fields_[k][cell]);
            maximumMassFraction_ = std::max(maximumMassFraction_, fields_[k][cell]);
        }
        compositionError_ = std::max(compositionError_, std::abs(sum - 1.0));
    }
}

bool MultiSpeciesModel::bounded(double tolerance) const {
    return minimumMassFraction_ >= -tolerance && maximumMassFraction_ <= 1.0 + tolerance;
}

void MultiSpeciesModel::updateTimeHistory() {
    for (std::size_t k = 0; k < speciesCount_; ++k) {
        oldOldFields_[k].values() = oldFields_[k].values();
        oldFields_[k].values() = fields_[k].values();
    }
}

void MultiSpeciesModel::appendOutputFields(std::map<std::string, const ScalarField*>& scalars) const {
    scalars["rho_mix"] = &densityField_;
    scalars["mu_mix"] = &viscosityField_;
    scalars["lambda_mix"] = &conductivityField_;
    scalars["cp_mix"] = &heatCapacityField_;
    for (std::size_t k = 0; k < speciesCount_; ++k) {
        scalars[fields_[k].name()] = &fields_[k];
        scalars[molecularDiffusivityFields_[k].name()] = &molecularDiffusivityFields_[k];
    }
}

std::unique_ptr<MultiSpeciesModel> makeMultiSpeciesModel(const Mesh& mesh, const Config& config,
                                                         std::shared_ptr<MaterialModel> material) {
    const std::string model = lowercase(config.get("species_model", "none"));
    if (model == "none" || model == "off") return nullptr;
    if (model != "transport" && model != "multispecies" && model != "multi_species") {
        throw std::runtime_error("Unknown species_model: " + model);
    }
    const std::string coupling = lowercase(config.get("flow_coupling", "passive"));
    if (coupling != "passive" && coupling != "constant_density") {
        throw std::runtime_error("Only passive constant-density species coupling is implemented; low-Mach variable-density coupling is future work");
    }
    return std::make_unique<MultiSpeciesModel>(mesh, config, std::move(material));
}

} // namespace mp
