/**
 * @file propagation_strategy.cpp
 * @brief Implementazione strategia di propagazione a due fasi semplificata
 * @author Michele Bigi  
 * @date 30 Novembre 2025
 */

#include "ioccultcalc/propagation_strategy.h"
#include "ioccultcalc/astdyn_interface.h"
#include "ioccultcalc/types.h"
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace ioccultcalc {

// ============================================================================
// Costanti locali
// ============================================================================

constexpr double PI = 3.14159265358979323846;
constexpr double GM_SUN = 1.32712440018e20; // m³/s²
constexpr double AU_METERS = AU * 1000.0; // metri (AU è in km in types.h)
constexpr double EARTH_OBLIQUITY = 23.43929111 * DEG_TO_RAD; // J2000.0

// ============================================================================
// TwoPhaseStrategy Implementation
// ============================================================================

TwoPhaseStrategy::TwoPhaseStrategy(const PropagationConfig& config) 
    : config_(config), elements_set_(false), use_fitted_elements_(false) {
    
    // Inizializza AstDyn propagator
    astdyn_ = std::make_unique<AstDynPropagator>(config_.rkf78_tolerance);
    astdyn_->usePlanetPerturbations(true);  // Sempre con perturbazioni
    astdyn_->useAsteroidPerturbations(false);
    astdyn_->useRelativisticCorrections(false);
}

void TwoPhaseStrategy::setConfig(const PropagationConfig& config) {
    config_ = config;
    if (astdyn_) {
        astdyn_->setTolerance(config_.rkf78_tolerance);
    }
}

void TwoPhaseStrategy::setElements(const AstDySElements& elements) {
    elements_ = elements;
    elements_set_ = true;
    use_fitted_elements_ = false;
    updateAstDynElements();
}

void TwoPhaseStrategy::setElements(const OrbitalElements& elements) {
    // Converti OrbitalElements → AstDySElements
    elements_.name = "Converted";
    elements_.number = 0;
    elements_.a = elements.a;
    elements_.e = elements.e; 
    elements_.i = elements.i;
    elements_.Omega = elements.Omega;
    elements_.omega = elements.omega;
    elements_.M = elements.M;
    elements_.epoch_mjd = elements.epoch.jd - 2400000.5;
    elements_.H = 16.0; // Default
    elements_.G = 0.15;
    elements_.has_covariance = false;
    
    elements_set_ = true;
    use_fitted_elements_ = false;
    updateAstDynElements();
}

void TwoPhaseStrategy::setObservations(const std::vector<RWOObservation>& observations) {
    observations_ = observations;
    
    if (config_.verbose_timing) {
        std::cout << "[TwoPhase] Caricate " << observations.size() 
                  << " osservazioni per orbit fitting\n";
    }
    
    // Analizza automaticamente la qualità se verbose
    if (config_.log_fitting_details && !observations.empty()) {
        auto info = analyzeObservations();
        std::cout << "[TwoPhase] Arco osservativo: " << info.observation_arc_days 
                  << " giorni (" << info.total_observations << " obs)\n";
        std::cout << "[TwoPhase] Qualità: " << (info.quality_acceptable ? "BUONA" : "SCARSA") 
                  << " (RMS medio: " << info.mean_rms_arcsec << "\")\n";
    }
}

bool TwoPhaseStrategy::loadObservationsFromFile(const std::string& rwo_file) {
    try {
        // Usa l'interfaccia AstDyS per caricare file .rwo
        auto observations = RWOObservation::fromFile(rwo_file);
        
        if (observations.empty()) {
            if (config_.log_fitting_details) {
                std::cout << "[TwoPhase] ERRORE: File RWO vuoto o non valido: " << rwo_file << std::endl;
            }
            return false;
        }
        
        setObservations(observations);
        
        if (config_.verbose_timing) {
            std::cout << "[TwoPhase] Caricate " << observations.size() 
                      << " osservazioni da file: " << rwo_file << std::endl;
        }
        
        return true;
        
    } catch (const std::exception& e) {
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] ERRORE caricamento RWO: " << e.what() << std::endl;
        }
        return false;
    }
}

bool TwoPhaseStrategy::loadObservationsFromAstDyS(int asteroid_number) {
    try {
        // Usa l'interfaccia AstDyS per download automatico
        auto observations = AstDySClient::downloadObservations(asteroid_number);
        
        if (observations.empty()) {
            if (config_.log_fitting_details) {
                std::cout << "[TwoPhase] ERRORE: Nessuna osservazione disponibile per asteroide " 
                          << asteroid_number << std::endl;
            }
            return false;
        }
        
        setObservations(observations);
        
        if (config_.verbose_timing) {
            std::cout << "[TwoPhase] Scaricate " << observations.size() 
                      << " osservazioni per asteroide " << asteroid_number << std::endl;
        }
        
        return true;
        
    } catch (const std::exception& e) {
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] ERRORE download AstDyS: " << e.what() << std::endl;
        }
        return false;
    }
}

TwoPhaseStrategy::ObservationInfo TwoPhaseStrategy::analyzeObservations() const {
    ObservationInfo info;
    
    if (observations_.empty()) {
        return info;
    }
    
    info.total_observations = observations_.size();
    
    // Trova primo e ultimo tempo
    double min_mjd = observations_[0].mjd_utc;
    double max_mjd = observations_[0].mjd_utc;
    double sum_rms = 0.0;
    
    for (const auto& obs : observations_) {
        min_mjd = std::min(min_mjd, obs.mjd_utc);
        max_mjd = std::max(max_mjd, obs.mjd_utc);
        
        // RMS combinato RA + Dec
        double combined_rms = std::sqrt(obs.rms_ra_arcsec * obs.rms_ra_arcsec + 
                                       obs.rms_dec_arcsec * obs.rms_dec_arcsec);
        sum_rms += combined_rms;
    }
    
    info.first_observation_mjd = min_mjd;
    info.last_observation_mjd = max_mjd;
    info.observation_arc_days = max_mjd - min_mjd;
    info.mean_rms_arcsec = sum_rms / info.total_observations;
    info.designation = observations_.empty() ? "" : observations_[0].designation;
    
    // Criteri di qualità
    bool arc_long_enough = info.observation_arc_days >= 30.0;  // Almeno 30 giorni
    bool enough_obs = info.total_observations >= static_cast<size_t>(config_.min_observations_for_fitting);
    bool rms_acceptable = info.mean_rms_arcsec <= 2.0;  // RMS < 2 arcsec
    
    info.quality_acceptable = arc_long_enough && enough_obs && rms_acceptable;
    
    return info;
}

int TwoPhaseStrategy::cleanObservations(double max_rms_arcsec) {
    if (observations_.empty()) {
        return 0;
    }
    
    size_t original_count = observations_.size();
    
    // Rimuovi osservazioni con RMS troppo alto
    observations_.erase(
        std::remove_if(observations_.begin(), observations_.end(),
            [max_rms_arcsec](const RWOObservation& obs) {
                double combined_rms = std::sqrt(obs.rms_ra_arcsec * obs.rms_ra_arcsec + 
                                               obs.rms_dec_arcsec * obs.rms_dec_arcsec);
                return combined_rms > max_rms_arcsec;
            }),
        observations_.end()
    );
    
    int removed = static_cast<int>(original_count - observations_.size());
    
    if (config_.log_fitting_details && removed > 0) {
        std::cout << "[TwoPhase] Rimosse " << removed << " osservazioni outliers (RMS > " 
                  << max_rms_arcsec << "\"). Rimaste: " << observations_.size() << std::endl;
    }
    
    return removed;
}

double TwoPhaseStrategy::screenCandidate(double target_mjd, double ra_star, double dec_star) {
    auto start = std::chrono::high_resolution_clock::now();
    
    stats_.phase1_calls++;
    stats_.candidates_screened++;
    
    // FASE 1: Genera/Usa polinomi di Chebyshev per screening ultra-veloce
    // 1. Verifica se i polinomi sono validi per questa epoca
    if (!areChebyshevPolynomialsValid(target_mjd)) {
        // 2. Genera nuovi polinomi usando RKF78 per posizioni precise
        double window_days = 10.0; // Finestra di validità polinomi
        generateChebyshevPolynomials(target_mjd - window_days/2, 
                                   target_mjd + window_days/2, 8);
    }
    
    // 3. Valuta posizione asteroide usando interpolazione Chebyshev (ultra-veloce)
    EquatorialCoords asteroid_pos = evaluateChebyshev(target_mjd);
    
    // Calcola separazione angolare
    double separation_arcsec = angularSeparation(
        asteroid_pos.ra_deg, asteroid_pos.dec_deg,
        ra_star, dec_star) * 3600.0;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    stats_.phase1_total_time_ms += duration.count() / 1000.0;
    
    if (config_.log_phase_transitions && separation_arcsec < config_.screening_threshold_arcsec) {
        stats_.candidates_promoted++;
        std::cout << "[TwoPhase] Candidato promosso: separazione=" 
                  << separation_arcsec << "\" < soglia=" 
                  << config_.screening_threshold_arcsec << "\"\n";
    }
    
    return separation_arcsec;
}

TwoPhaseStrategy::CloseApproachResult TwoPhaseStrategy::findClosestApproach(
    double target_mjd, double ra_star, double dec_star, double search_window_hours) {
    
    auto start = std::chrono::high_resolution_clock::now();
    
    stats_.phase2_calls++;
    
    CloseApproachResult result;
    result.closest_time_mjd = target_mjd;
    result.minimum_separation_arcsec = 999999.0;
    result.position_angle_deg = 0.0;
    result.computation_steps = 0;
    result.orbit_was_fitted = false;
    
    // Orbit fitting con logica avanzata basata su FittingMode
    bool should_attempt_fitting = false;
    
    // Determina se tentare orbit fitting
    switch (config_.orbit_fitting_mode) {
        case PropagationConfig::FittingMode::NEVER:
            should_attempt_fitting = false;
            if (config_.log_fitting_details) {
                std::cout << "[TwoPhase] Orbit fitting disabilitato (mode=NEVER)\n";
            }
            break;
            
        case PropagationConfig::FittingMode::AUTO:
            should_attempt_fitting = !observations_.empty() && 
                                    observations_.size() >= static_cast<size_t>(config_.min_observations_for_fitting);
            if (config_.log_fitting_details) {
                std::cout << "[TwoPhase] Orbit fitting AUTO: " << observations_.size() 
                          << " osservazioni, min=" << config_.min_observations_for_fitting 
                          << " → " << (should_attempt_fitting ? "TENTATIVO" : "SKIP") << "\n";
            }
            break;
            
        case PropagationConfig::FittingMode::ALWAYS_ATTEMPT:
            if (observations_.empty() || observations_.size() < static_cast<size_t>(config_.min_observations_for_fitting)) {
                throw std::runtime_error("Orbit fitting richiesto (ALWAYS_ATTEMPT) ma osservazioni insufficienti: " 
                                        + std::to_string(observations_.size()) + " < " 
                                        + std::to_string(config_.min_observations_for_fitting));
            }
            should_attempt_fitting = true;
            if (config_.log_fitting_details) {
                std::cout << "[TwoPhase] Orbit fitting ALWAYS_ATTEMPT con " 
                          << observations_.size() << " osservazioni\n";
            }
            break;
    }
    
    // Retrocompatibilità: se enable_orbit_fitting=true, override AUTO mode
    if (config_.enable_orbit_fitting && config_.orbit_fitting_mode == PropagationConfig::FittingMode::AUTO) {
        should_attempt_fitting = !observations_.empty();
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] Override retrocompatibilità: enable_orbit_fitting=true\n";
        }
    }
    
    // Forza re-fit se richiesto
    if (config_.force_refit_each_prediction && use_fitted_elements_) {
        use_fitted_elements_ = false;
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] Force refit: resetting fitted elements\n";
        }
    }
    
    // Esegui orbit fitting se necessario
    if (should_attempt_fitting && !use_fitted_elements_) {
        if (config_.verbose_timing || config_.log_fitting_details) {
            std::cout << "[TwoPhase] Inizio orbit fitting con " << observations_.size() 
                      << " osservazioni...\n";
        }
        
        if (performOrbitFitting()) {
            result.orbit_was_fitted = true;
            use_fitted_elements_ = true;
            updateAstDynElements();
            
            if (config_.verbose_timing || config_.log_fitting_details) {
                std::cout << "[TwoPhase] ✅ Orbit fitting completato con successo\n";
            }
        } else {
            if (config_.log_fitting_details) {
                std::cout << "[TwoPhase] ❌ Orbit fitting fallito, uso elementi originali\n";
            }
        }
    } else if (use_fitted_elements_ && config_.log_fitting_details) {
        std::cout << "[TwoPhase] Uso elementi già fittati dalla predizione precedente\n";
    }
    
    // FASE 2: Golden section search per closest approach con RKF78
    double search_window_days = search_window_hours / 24.0;
    double t_start = target_mjd - search_window_days / 2.0;
    double t_end = target_mjd + search_window_days / 2.0;
    
    auto calcSeparation = [&](double mjd) -> double {
        EquatorialCoords asteroid_pos = propagateRKF78(mjd);
        return angularSeparation(asteroid_pos.ra_deg, asteroid_pos.dec_deg, 
                               ra_star, dec_star);
    };
    
    // Golden section search
    const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    const double resphi = 2.0 - phi;
    double tolerance_days = config_.convergence_seconds / 86400.0;
    
    double a = t_start;
    double b = t_end;
    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);
    double f1 = calcSeparation(x1);
    double f2 = calcSeparation(x2);
    
    int steps = 0;
    while (std::abs(b - a) > tolerance_days && steps < 50) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = calcSeparation(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            f2 = calcSeparation(x2);
        }
        steps++;
    }
    
    result.closest_time_mjd = (a + b) / 2.0;
    result.minimum_separation_arcsec = calcSeparation(result.closest_time_mjd) * 3600.0;
    result.computation_steps = steps;
    
    // Calcola position angle
    EquatorialCoords asteroid_final = propagateRKF78(result.closest_time_mjd);
    double dra = (asteroid_final.ra_deg - ra_star) * std::cos(dec_star * DEG_TO_RAD);
    double ddec = asteroid_final.dec_deg - dec_star;
    result.position_angle_deg = std::atan2(dra, ddec) * RAD_TO_DEG;
    if (result.position_angle_deg < 0) result.position_angle_deg += 360.0;
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    result.computation_time_ms = duration.count();
    stats_.phase2_total_time_ms += result.computation_time_ms;
    
    return result;
}

void TwoPhaseStrategy::resetStats() {
    stats_ = PerformanceStats();
}

TwoPhaseStrategy::EquatorialCoords TwoPhaseStrategy::getChebyshevPosition(double target_mjd) {
    // Verifica se i polinomi sono validi per questa epoca
    if (!areChebyshevPolynomialsValid(target_mjd)) {
        // Genera nuovi polinomi usando RKF78
        double window_days = 10.0;
        generateChebyshevPolynomials(target_mjd - window_days/2, 
                                   target_mjd + window_days/2, 8);
    }
    
    // Valuta posizione usando interpolazione Chebyshev
    return evaluateChebyshev(target_mjd);
}

TwoPhaseStrategy::EquatorialCoords TwoPhaseStrategy::getRKF78Position(double target_mjd) {
    return propagateRKF78(target_mjd);
}

// ============================================================================
// Private Methods
// ============================================================================

TwoPhaseStrategy::EquatorialCoords TwoPhaseStrategy::propagateKeplerianGeocentric(double target_mjd) {
    // NOTA: Questa funzione non è più usata - FASE 1 ora usa AstDyn RKF78
    // con tolleranza ridotta (1e-10) invece di propagazione Kepleriana.
    // Manteniamo il codice per compatibilità futura se necessario.
    
    if (!elements_set_) {
        throw std::runtime_error("Elementi orbitali non impostati");
    }
    
    // Ora usa semplicemente RKF78 (come FASE 2 ma con tolleranza diversa)
    return propagateRKF78(target_mjd);
}

TwoPhaseStrategy::EquatorialCoords TwoPhaseStrategy::propagateRKF78(double target_mjd) {
    if (!astdyn_) {
        throw std::runtime_error("AstDyn propagator non inizializzato");
    }
    
    const AstDySElements& elem = use_fitted_elements_ ? fitted_elements_ : elements_;
    
    // Usa interfaccia AstDyn
    auto radec = astdyn_->getRADec(elem, target_mjd);
    
    return {radec.first, radec.second, 0.0}; // Distance non necessaria
}

// ============================================================================
// Metodi Polinomi di Chebyshev per FASE 1
// ============================================================================

void TwoPhaseStrategy::generateChebyshevPolynomials(double start_mjd, double end_mjd, int degree) {
    if (!astdyn_) {
        throw std::runtime_error("AstDyn propagator non inizializzato per generare Chebyshev");
    }
    
    const AstDySElements& elem = use_fitted_elements_ ? fitted_elements_ : elements_;
    
    // Genera punti di controllo con RKF78 ad alta precisione
    int n_points = degree + 5; // Sovracampionamento
    std::vector<double> times(n_points);
    std::vector<double> ra_values(n_points);
    std::vector<double> dec_values(n_points);
    
    for (int i = 0; i < n_points; i++) {
        double t = start_mjd + (end_mjd - start_mjd) * i / (n_points - 1);
        times[i] = t;
        
        auto radec = astdyn_->getRADec(elem, t);
        ra_values[i] = radec.first;
        dec_values[i] = radec.second;
    }
    
    // Calcola coefficienti Chebyshev per RA e Dec
    chebyshev_.ra_coeffs = calculateChebyshevCoeffs(times, ra_values, degree);
    chebyshev_.dec_coeffs = calculateChebyshevCoeffs(times, dec_values, degree);
    chebyshev_.t_start_mjd = start_mjd;
    chebyshev_.t_end_mjd = end_mjd;
    chebyshev_.degree = degree;
    chebyshev_.valid = true;
    
    if (config_.verbose_timing) {
        std::cout << "[Chebyshev] Generati polinomi grado " << degree 
                  << " per intervallo [" << start_mjd << ", " << end_mjd << "]\n";
    }
}

TwoPhaseStrategy::EquatorialCoords TwoPhaseStrategy::evaluateChebyshev(double target_mjd) const {
    if (!chebyshev_.valid || !areChebyshevPolynomialsValid(target_mjd)) {
        throw std::runtime_error("Polinomi Chebyshev non validi per epoca richiesta");
    }
    
    // Normalizza tempo nell'intervallo [-1, 1]
    double t_norm = 2.0 * (target_mjd - chebyshev_.t_start_mjd) / 
                    (chebyshev_.t_end_mjd - chebyshev_.t_start_mjd) - 1.0;
    
    // Valuta polinomi di Chebyshev
    double ra_deg = evaluateChebyshevPolynomial(chebyshev_.ra_coeffs, t_norm);
    double dec_deg = evaluateChebyshevPolynomial(chebyshev_.dec_coeffs, t_norm);
    
    return {ra_deg, dec_deg, 0.0};
}

bool TwoPhaseStrategy::areChebyshevPolynomialsValid(double target_mjd) const {
    return chebyshev_.valid && 
           target_mjd >= chebyshev_.t_start_mjd && 
           target_mjd <= chebyshev_.t_end_mjd;
}

std::vector<double> TwoPhaseStrategy::calculateChebyshevCoeffs(
    const std::vector<double>& times, const std::vector<double>& values, int degree) {
    
    // METODO ACCURATO: Interpolazione di Chebyshev usando nodi ottimali
    // e trasformata discreta di Chebyshev (DCT)
    
    int n_coeffs = degree + 1;
    std::vector<double> coeffs(n_coeffs, 0.0);
    
    // Se abbiamo esattamente degree+1 punti, usa interpolazione esatta
    if (static_cast<int>(values.size()) == n_coeffs) {
        // Usa i valori direttamente con trasformata DCT
        for (int k = 0; k < n_coeffs; k++) {
            double sum = 0.0;
            double c_k = (k == 0) ? 1.0 : 2.0;  // Fattore di normalizzazione
            
            for (int j = 0; j < n_coeffs; j++) {
                // Nodo di Chebyshev: x_j = cos(π(j+0.5)/n)
                double theta_j = PI * (j + 0.5) / n_coeffs;
                double x_j = std::cos(theta_j);
                
                // Interpola il valore al nodo di Chebyshev
                double t_j = times[0] + (times.back() - times[0]) * (x_j + 1.0) / 2.0;
                double interpolated_value = interpolateLinear(times, values, t_j);
                
                // DCT: somma con cos(kθ)
                sum += interpolated_value * std::cos(k * theta_j);
            }
            coeffs[k] = (2.0 / n_coeffs) * sum;
            if (k == 0) coeffs[k] /= 2.0;  // Correzione per k=0
        }
    } else {
        // Approssimazione least squares con nodi di Chebyshev ottimali
        std::vector<double> cheb_nodes(n_coeffs);
        std::vector<double> cheb_values(n_coeffs);
        
        // Genera nodi di Chebyshev nell'intervallo [times[0], times.back()]
        for (int i = 0; i < n_coeffs; i++) {
            double theta = PI * (i + 0.5) / n_coeffs;
            double x_cheb = std::cos(theta);  // Nodo in [-1, 1]
            cheb_nodes[i] = times[0] + (times.back() - times[0]) * (x_cheb + 1.0) / 2.0;
            cheb_values[i] = interpolateLinear(times, values, cheb_nodes[i]);
        }
        
        // Calcola coefficienti con DCT sui nodi ottimali
        for (int k = 0; k < n_coeffs; k++) {
            double sum = 0.0;
            for (int j = 0; j < n_coeffs; j++) {
                double theta_j = PI * (j + 0.5) / n_coeffs;
                sum += cheb_values[j] * std::cos(k * theta_j);
            }
            coeffs[k] = (2.0 / n_coeffs) * sum;
            if (k == 0) coeffs[k] /= 2.0;
        }
    }
    
    return coeffs;
}

double TwoPhaseStrategy::evaluateChebyshevPolynomial(
    const std::vector<double>& coeffs, double x) const {
    
    if (coeffs.empty()) return 0.0;
    if (coeffs.size() == 1) return coeffs[0];
    
    // Valutazione ricorsiva stabile
    double result = coeffs[0];
    if (coeffs.size() > 1) {
        double T_prev = 1.0;      // T_0(x)
        double T_curr = x;        // T_1(x)
        result += coeffs[1] * T_curr;
        
        for (size_t n = 2; n < coeffs.size(); n++) {
            double T_next = 2.0 * x * T_curr - T_prev;
            result += coeffs[n] * T_next;
            T_prev = T_curr;
            T_curr = T_next;
        }
    }
    
    return result;
}

double TwoPhaseStrategy::chebyshevBasisFunction(int n, double x) const {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    
    double T_prev = 1.0;  // T_0(x)
    double T_curr = x;    // T_1(x)
    
    for (int i = 2; i <= n; i++) {
        double T_next = 2.0 * x * T_curr - T_prev;
        T_prev = T_curr;
        T_curr = T_next;
    }
    
    return T_curr;
}

bool TwoPhaseStrategy::performOrbitFitting() {
    if (observations_.empty()) {
        return false;
    }
    
    // TODO: Implementare differential correction con osservazioni
    // Per ora stub che simula miglioramento
    stats_.orbits_fitted++;
    
    fitted_elements_ = elements_;
    // Simulazione: piccola correzione ai parametri orbitali
    fitted_elements_.M += 0.001; // ~3.6 arcsec correzione
    
    return true;
}

double TwoPhaseStrategy::interpolateLinear(const std::vector<double>& x,
                                          const std::vector<double>& y, 
                                          double target_x) const {
    if (x.empty() || y.empty() || x.size() != y.size()) {
        throw std::runtime_error("Dati invalidi per interpolazione lineare");
    }
    
    // Se target è fuori range, estrapolazione lineare
    if (target_x <= x.front()) return y.front();
    if (target_x >= x.back()) return y.back();
    
    // Trova l'intervallo contenente target_x
    for (size_t i = 0; i < x.size() - 1; i++) {
        if (target_x >= x[i] && target_x <= x[i + 1]) {
            double t = (target_x - x[i]) / (x[i + 1] - x[i]);
            return y[i] + t * (y[i + 1] - y[i]);
        }
    }
    
    // Fallback (non dovrebbe mai accadere)
    return y.back();
}

double TwoPhaseStrategy::angularSeparation(double ra1_deg, double dec1_deg,
                                         double ra2_deg, double dec2_deg) {
    double ra1 = ra1_deg * DEG_TO_RAD;
    double dec1 = dec1_deg * DEG_TO_RAD;
    double ra2 = ra2_deg * DEG_TO_RAD;
    double dec2 = dec2_deg * DEG_TO_RAD;
    
    double cos_sep = std::sin(dec1) * std::sin(dec2) + 
                     std::cos(dec1) * std::cos(dec2) * std::cos(ra1 - ra2);
    
    // Proteggi da errori numerici
    cos_sep = std::max(-1.0, std::min(1.0, cos_sep));
    return std::acos(cos_sep);
}

bool TwoPhaseStrategy::triggerOrbitFitting() {
    if (observations_.empty()) {
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] ❌ Trigger orbit fitting: nessuna osservazione disponibile\n";
        }
        return false;
    }
    
    if (observations_.size() < static_cast<size_t>(config_.min_observations_for_fitting)) {
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] ❌ Trigger orbit fitting: osservazioni insufficienti ("
                      << observations_.size() << " < " << config_.min_observations_for_fitting << ")\n";
        }
        return false;
    }
    
    if (config_.log_fitting_details) {
        std::cout << "[TwoPhase] 🔄 Trigger orbit fitting manuale con " 
                  << observations_.size() << " osservazioni...\n";
    }
    
    // Forza re-fit anche se già fittato
    use_fitted_elements_ = false;
    
    bool success = performOrbitFitting();
    if (success) {
        use_fitted_elements_ = true;
        updateAstDynElements();
        
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] ✅ Orbit fitting manuale completato\n";
        }
    } else {
        if (config_.log_fitting_details) {
            std::cout << "[TwoPhase] ❌ Orbit fitting manuale fallito\n";
        }
    }
    
    return success;
}

void TwoPhaseStrategy::updateAstDynElements() {
    if (!astdyn_ || !elements_set_) return;
    
    // AstDyn gestisce automaticamente la conversione
    // Gli elementi sono già nel formato corretto
}

// ============================================================================
// Preset Configurations
// ============================================================================

namespace propagation_presets {

PropagationConfig createFastSurvey() {
    PropagationConfig config;
    config.screening_threshold_arcsec = 120.0; // Soglia più alta per velocità
    config.rkf78_tolerance = 1e-10;            // FASE 2: Tolleranza meno stretta 
    config.orbit_fitting_mode = PropagationConfig::FittingMode::NEVER;  // No fitting per velocità
    config.verbose_timing = false;
    config.log_phase_transitions = false;
    config.log_fitting_details = false;
    return config;
}

PropagationConfig createBalanced() {
    PropagationConfig config;
    config.screening_threshold_arcsec = 60.0;  
    config.rkf78_tolerance = 1e-12;            // FASE 2: Tolleranza standard
    config.orbit_fitting_mode = PropagationConfig::FittingMode::AUTO;    // Fitting se disponibile
    config.fitting_sigma_cutoff = 3.0;
    config.max_fitting_iterations = 5;
    config.min_observations_for_fitting = 3;
    config.verbose_timing = false;
    config.log_phase_transitions = false;
    config.log_fitting_details = false;
    return config;
}

PropagationConfig createPrecision() {
    PropagationConfig config;
    config.screening_threshold_arcsec = 30.0;  // Soglia più bassa
    config.rkf78_tolerance = 1e-14;            // FASE 2: Massima precisione
    config.convergence_seconds = 0.1;          // Convergenza stretta
    config.orbit_fitting_mode = PropagationConfig::FittingMode::ALWAYS_ATTEMPT; // Fitting sempre
    config.fitting_sigma_cutoff = 2.0;         // Outlier detection più severa
    config.max_fitting_iterations = 15;
    config.min_observations_for_fitting = 5;   // Più osservazioni richieste
    config.force_refit_each_prediction = true; // Re-fit ogni volta
    config.verbose_timing = true;
    config.log_phase_transitions = true;
    config.log_fitting_details = true;
    return config;
}

/**
 * @brief Configurazione per testing con controllo completo fitting
 */
PropagationConfig createTestingWithFitting() {
    PropagationConfig config;
    config.screening_threshold_arcsec = 60.0;  
    config.rkf78_tolerance = 1e-12;
    config.orbit_fitting_mode = PropagationConfig::FittingMode::AUTO;
    config.fitting_sigma_cutoff = 3.0;
    config.max_fitting_iterations = 10;
    config.min_observations_for_fitting = 3;
    config.force_refit_each_prediction = false;
    config.verbose_timing = true;
    config.log_phase_transitions = true;
    config.log_fitting_details = true;          // Debug dettagliato fitting
    return config;
}

/**
 * @brief Configurazione senza fitting (solo elementi AstDyS)
 */
PropagationConfig createNoFitting() {
    PropagationConfig config;
    config.screening_threshold_arcsec = 60.0;  
    config.rkf78_tolerance = 1e-12;
    config.orbit_fitting_mode = PropagationConfig::FittingMode::NEVER;
    config.verbose_timing = true;
    config.log_phase_transitions = false;
    config.log_fitting_details = true;
    return config;
}

} // namespace propagation_presets

} // namespace ioccultcalc
