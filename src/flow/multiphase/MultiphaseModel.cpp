#include "flow/multiphase/MultiphaseModel.hpp"

#include "flow/multiphase/AlgebraicSlipModel.hpp"
#include "flow/multiphase/DpmModel.hpp"
#include "flow/multiphase/EulerEulerModel.hpp"
#include "flow/multiphase/VofModel.hpp"

#include <stdexcept>

namespace mp {

void MultiphaseModel::notImplemented() const {
    throw std::logic_error("Multiphase model '" + name() + "' is an interface only; its solve routines are not implemented");
}

void MultiphaseModel::initialize(MaterialModel&) { notImplemented(); }
void MultiphaseModel::beginOuterIteration() { notImplemented(); }
void MultiphaseModel::solve(MultiphaseContext&) { notImplemented(); }
void MultiphaseModel::updateTimeHistory() { notImplemented(); }
double MultiphaseModel::change() const { notImplemented(); }
void MultiphaseModel::appendOutputFields(std::map<std::string, const ScalarField*>&,
                                         std::map<std::string, const VectorField*>&) const {}

const ScalarField& VofModel::volumeFraction(std::size_t) const { notImplemented(); }
const VectorField& VofModel::surfaceTensionForce() const { notImplemented(); }
void DpmModel::inject(const ParticleParcel&) { notImplemented(); }
const DiscretePhaseSources& DpmModel::continuousPhaseSources() const { notImplemented(); }
const ScalarField& EulerEulerModel::phaseFraction(std::size_t) const { notImplemented(); }
const VectorField& EulerEulerModel::phaseVelocity(std::size_t) const { notImplemented(); }
const VectorField& EulerEulerModel::interphaseMomentumTransfer(std::size_t) const { notImplemented(); }
const ScalarField& AlgebraicSlipModel::phaseFraction(std::size_t) const { notImplemented(); }
const VectorField& AlgebraicSlipModel::slipVelocity(std::size_t) const { notImplemented(); }

std::unique_ptr<MultiphaseModel> makeMultiphaseModel(const Mesh&, const Config& config) {
    const std::string model = lowercase(config.get("multiphase_model", "none"));
    if (model == "none" || model == "off") return nullptr;
    if (model == "vof") return std::make_unique<VofModel>();
    if (model == "dpm" || model == "euler_lagrange" || model == "euler-lagrange") {
        return std::make_unique<DpmModel>();
    }
    if (model == "euler_euler" || model == "euler-euler" || model == "multi_fluid") {
        return std::make_unique<EulerEulerModel>();
    }
    if (model == "asm" || model == "algebraic_slip" || model == "algebraic-slip") {
        return std::make_unique<AlgebraicSlipModel>();
    }
    throw std::runtime_error("Unknown multiphase_model: " + model);
}

} // namespace mp
