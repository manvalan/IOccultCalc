/**
 * @file astdyn_propagation_helper.cpp
 * @brief Implementazione AstDynPropagationHelper
 */

#include "ioccultcalc/astdyn_propagation_helper.h"
#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/time_utils.h"
#include <cmath>

namespace ioccultcalc {

constexpr double MJD_TO_JD = 2400000.5;
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;

AstDynPropagationHelper::AstDynPropagationHelper()
    : propagator_(std::make_unique<AstDynPropagator>(1e-12))
    , tolerance_(1e-12)
{
}

AstDynPropagationHelper& AstDynPropagationHelper::getInstance() {
    static AstDynPropagationHelper instance;
    return instance;
}

AstDySElements AstDynPropagationHelper::propagate(
    const AstDySElements& elements,
    double target_mjd_tdb) {
    return propagator_->propagate(elements, target_mjd_tdb);
}

std::pair<double, double> AstDynPropagationHelper::getRADec(
    const AstDySElements& elements,
    double target_mjd_tdb) {
    return propagator_->getRADec(elements, target_mjd_tdb);
}

std::pair<double, double> AstDynPropagationHelper::propagateToRADec(
    const EquinoctialElements& elements,
    const JulianDate& target_jd) {
    AstDySElements astdys = convertFromEquinoctial(elements);
    double target_mjd_tdb = target_jd.jd - MJD_TO_JD;
    return getRADec(astdys, target_mjd_tdb);
}

std::vector<std::tuple<double, double, double>> 
AstDynPropagationHelper::propagateRange(
    const AstDySElements& elements,
    double start_mjd_tdb,
    double end_mjd_tdb,
    double step_days) {
    std::vector<std::tuple<double, double, double>> results;
    for (double mjd = start_mjd_tdb; mjd <= end_mjd_tdb; mjd += step_days) {
        auto radec = getRADec(elements, mjd);
        results.emplace_back(mjd, radec.first, radec.second);
    }
    return results;
}

AstDySElements AstDynPropagationHelper::convertFromEquinoctial(
    const EquinoctialElements& eq) {
    AstDySElements astdys;
    astdys.name = eq.designation;
    astdys.number = eq.number;
    astdys.a = eq.a;
    double e = sqrt(eq.h * eq.h + eq.k * eq.k);
    astdys.e = e;
    double i = 2.0 * atan(sqrt(eq.p * eq.p + eq.q * eq.q));
    astdys.i = i * RAD_TO_DEG;
    double Omega = atan2(eq.p, eq.q);
    if (Omega < 0) Omega += 2.0 * M_PI;
    astdys.Omega = Omega * RAD_TO_DEG;
    double omega_plus_Omega = atan2(eq.h, eq.k);
    double omega = omega_plus_Omega - Omega;
    if (omega < 0) omega += 2.0 * M_PI;
    astdys.omega = omega * RAD_TO_DEG;
    double lambda = eq.lambda;
    double M = lambda - omega_plus_Omega;
    while (M < 0) M += 2.0 * M_PI;
    while (M >= 2.0 * M_PI) M -= 2.0 * M_PI;
    astdys.M = M * RAD_TO_DEG;
    astdys.epoch_mjd = eq.epoch.jd - MJD_TO_JD;
    astdys.H = 0.0;
    astdys.G = 0.15;
    astdys.has_covariance = false;
    return astdys;
}

AstDySElements AstDynPropagationHelper::convertFromOrbital(
    const OrbitalElements& orb) {
    AstDySElements astdys;
    astdys.name = orb.designation;
    astdys.number = orb.number;
    astdys.a = orb.a;
    astdys.e = orb.e;
    astdys.i = orb.i * RAD_TO_DEG;
    astdys.Omega = orb.Omega * RAD_TO_DEG;
    astdys.omega = orb.omega * RAD_TO_DEG;
    astdys.M = orb.M * RAD_TO_DEG;
    astdys.epoch_mjd = orb.epoch.jd - MJD_TO_JD;
    astdys.H = orb.H;
    astdys.G = orb.G;
    astdys.has_covariance = false;
    return astdys;
}

void AstDynPropagationHelper::setTolerance(double tolerance) {
    tolerance_ = tolerance;
    propagator_ = std::make_unique<AstDynPropagator>(tolerance);
}

void AstDynPropagationHelper::setPerturbations(bool planets, bool asteroids, bool relativity) {
    // Nota: AstDynPropagator potrebbe non supportare configurazione runtime
}

namespace astdyn_utils {

AstDySElements quickPropagate(const AstDySElements& elements, double target_mjd_tdb) {
    AstDynPropagator prop(1e-12);
    return prop.propagate(elements, target_mjd_tdb);
}

std::pair<double, double> quickRADec(const AstDySElements& elements, double target_mjd_tdb) {
    AstDynPropagator prop(1e-12);
    return prop.getRADec(elements, target_mjd_tdb);
}

AstDySElements toAstDyS(const EquinoctialElements& eq) {
    return AstDynPropagationHelper::convertFromEquinoctial(eq);
}

AstDySElements toAstDyS(const OrbitalElements& orb) {
    return AstDynPropagationHelper::convertFromOrbital(orb);
}

} // namespace astdyn_utils

} // namespace ioccultcalc
