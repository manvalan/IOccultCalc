/**
 * @file italoccult_integration.cpp
 * @brief Implementazione integrazione ITALOccultLibrary in IOccultCalc
 * @author IOccultCalc Team
 * @date 1 Dicembre 2025
 */

#include "italoccult_integration.h"
#include <cmath>
#include <stdexcept>

namespace ioccultcalc {

ITALOccultIntegration::ITALOccultIntegration(bool use_high_accuracy)
    : high_accuracy_(use_high_accuracy)
    , current_epoch_(0.0)
    , initialized_(false)
{
    // Crea wrapper con configurazione appropriata
    auto settings = high_accuracy_ 
        ? PropagationSettings::highAccuracy()
        : PropagationSettings::fast();
    
    wrapper_ = std::make_unique<AstDynWrapper>(settings);
}

ITALOccultIntegration::~ITALOccultIntegration() = default;

bool ITALOccultIntegration::loadAsteroidFromEQ1(const std::string& filepath) {
    bool success = wrapper_->loadFromEQ1File(filepath);
    
    if (success) {
        asteroid_name_ = wrapper_->getObjectName();
        current_epoch_ = wrapper_->getCurrentEpoch();
        initialized_ = true;
    } else {
        initialized_ = false;
    }
    
    return success;
}

void ITALOccultIntegration::setOrbitalElements(
    const std::string& name,
    double a, double e, double i,
    double Omega, double omega, double M,
    double epoch_mjd_tdb)
{
    // Converti angoli da gradi a radianti
    constexpr double DEG_TO_RAD = M_PI / 180.0;
    
    wrapper_->setKeplerianElements(
        a, e, i * DEG_TO_RAD,
        Omega * DEG_TO_RAD, omega * DEG_TO_RAD, M * DEG_TO_RAD,
        epoch_mjd_tdb, name
    );
    
    asteroid_name_ = name;
    current_epoch_ = epoch_mjd_tdb;
    initialized_ = true;
}

AsteroidState ITALOccultIntegration::propagateToEpoch(double target_mjd_tdb) {
    if (!initialized_) {
        throw std::runtime_error("ITALOccultIntegration: nessun asteroide caricato. "
                                "Chiamare loadAsteroidFromEQ1() o setOrbitalElements() prima.");
    }
    
    // Propaga con wrapper AstDyn
    auto icrf_state = wrapper_->propagateToEpoch(target_mjd_tdb);
    
    // Converti in formato IOccultCalc
    return convertToAsteroidState(icrf_state, asteroid_name_);
}

std::vector<AsteroidState> ITALOccultIntegration::propagateToEpochs(
    const std::vector<double>& target_mjd_tdbs)
{
    std::vector<AsteroidState> results;
    results.reserve(target_mjd_tdbs.size());
    
    for (double mjd : target_mjd_tdbs) {
        results.push_back(propagateToEpoch(mjd));
    }
    
    return results;
}

double ITALOccultIntegration::getCurrentEpoch() const {
    return current_epoch_;
}

std::string ITALOccultIntegration::getAsteroidName() const {
    return asteroid_name_;
}

bool ITALOccultIntegration::isInitialized() const {
    return initialized_;
}

std::string ITALOccultIntegration::getLastPropagationInfo() const {
    return wrapper_->getLastPropagationStats();
}

void ITALOccultIntegration::setHighAccuracy(bool high_accuracy) {
    if (high_accuracy != high_accuracy_) {
        high_accuracy_ = high_accuracy;
        
        // Ricrea wrapper con nuova configurazione
        auto settings = high_accuracy_ 
            ? PropagationSettings::highAccuracy()
            : PropagationSettings::fast();
        
        wrapper_ = std::make_unique<AstDynWrapper>(settings);
        
        // Nota: elementi orbitali vanno ricaricati dopo questo cambio
        initialized_ = false;
    }
}

AsteroidState ITALOccultIntegration::convertToAsteroidState(
    const CartesianStateICRF& icrf_state,
    const std::string& name) const
{
    AsteroidState state;
    state.name = name;
    state.epoch_mjd_tdb = icrf_state.epoch_mjd_tdb;
    
    // Copia posizione e velocità
    state.x_au = icrf_state.position.x();
    state.y_au = icrf_state.position.y();
    state.z_au = icrf_state.position.z();
    state.vx_au_day = icrf_state.velocity.x();
    state.vy_au_day = icrf_state.velocity.y();
    state.vz_au_day = icrf_state.velocity.z();
    
    // Calcola parametri orbitali approssimativi (per reference)
    // Nota: per valori precisi, usare OrbitalConversions
    double r = std::sqrt(state.x_au * state.x_au + 
                        state.y_au * state.y_au + 
                        state.z_au * state.z_au);
    double v = std::sqrt(state.vx_au_day * state.vx_au_day + 
                        state.vy_au_day * state.vy_au_day + 
                        state.vz_au_day * state.vz_au_day);
    
    // Energia specifica e semiasse (approssimato)
    constexpr double GM_SUN = 0.0002959122082855911;  // AU³/day²
    double energy = 0.5 * v * v - GM_SUN / r;
    state.a_au = -GM_SUN / (2.0 * energy);
    
    // Eccentricità (approssimata da momento angolare)
    double hx = state.y_au * state.vz_au_day - state.z_au * state.vy_au_day;
    double hy = state.z_au * state.vx_au_day - state.x_au * state.vz_au_day;
    double hz = state.x_au * state.vy_au_day - state.y_au * state.vx_au_day;
    double h = std::sqrt(hx * hx + hy * hy + hz * hz);
    
    state.e = std::sqrt(1.0 + 2.0 * energy * h * h / (GM_SUN * GM_SUN));
    
    // Inclinazione (da componente z del momento angolare)
    state.i_deg = std::acos(hz / h) * 180.0 / M_PI;
    
    return state;
}

// Funzione helper
AsteroidState quickPropagateFromEQ1(const std::string& eq1_file,
                                    double target_mjd_tdb,
                                    bool high_accuracy)
{
    ITALOccultIntegration integrator(high_accuracy);
    
    if (!integrator.loadAsteroidFromEQ1(eq1_file)) {
        throw std::runtime_error("Impossibile caricare file .eq1: " + eq1_file);
    }
    
    return integrator.propagateToEpoch(target_mjd_tdb);
}

} // namespace ioccultcalc
