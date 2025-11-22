#include "ioccultcalc/orbit_fitter.h"
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/force_model.h"
#include "ioccultcalc/orbit_propagator.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>

namespace ioccultcalc {

// Implementazione interna
class OrbitFitter::Impl {
public:
    Ephemeris ephemeris;  // Deprecated - usa solo 2-body
    ForceModel forceModel;
    OrbitPropagator propagator;  // Usa questo per propagazione con perturbazioni
    EquinoctialElements currentEqElements;  // Elementi correnti per propagazione
    
    Impl() {
        // Configura OrbitPropagator con RA15 e perturbazioni planetarie
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RA15;
        opts.usePlanetaryPerturbations = true;
        opts.useNonGravitational = false;
        opts.stepSize = 10.0;  // 10 giorni
        opts.tolerance = 1e-12;
        
        propagator = OrbitPropagator(opts);
    }
};

OrbitFitter::OrbitFitter() : pImpl(new Impl()) {}
OrbitFitter::~OrbitFitter() { delete pImpl; }

OrbitFitResult OrbitFitter::fit(const OrbitalElements& initialElements,
                                const ObservationSet& observations,
                                const OrbitFitOptions& options) {
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║          Orbit Fitting - Differential Correction          ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    OrbitFitResult result;
    result.fittedElements = initialElements;
    result.observations = observations;
    
    // Filtra osservazioni valide
    ObservationSet validObs;
    validObs.objectDesignation = observations.objectDesignation;
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) {
            validObs.observations.push_back(obs);
        }
    }
    validObs.computeStatistics();
    
    std::cout << "→ Osservazioni disponibili: " << validObs.numberOfObservations << "\n";
    std::cout << "→ Arc temporale: " << validObs.arcLength << " giorni\n";
    std::cout << "→ Epoca: " << observations.firstObservation.jd << " - " 
              << observations.lastObservation.jd << "\n";
    
    std::cout << "→ Epoca elementi nominali: JD " << initialElements.epoch.jd << "\n\n";
    
    // IMPORTANTE: Ephemeris.compute() propaga automaticamente dall'epoca degli elementi
    // all'epoca richiesta per ogni osservazione. Tuttavia usa solo propagazione Keplerian
    // 2-body (senza perturbazioni planetarie). Per lunghi archi temporali (>1 anno)
    // questo causa errori di propagazione.
    // TODO: Sostituire Ephemeris con OrbitPropagator che usa RA15 + perturbazioni complete
    
    OrbitalElements currentElements = initialElements;
    
    // Setup force model - currently ForceModel doesn't have runtime configuration
    // Planetary perturbations are fixed in the model
    // TODO: Add API to ForceModel for enabling/disabling specific bodies
    
    double previousRMS = 1e10;
    
    // Iterazione least-squares
    for (int iter = 0; iter < options.maxIterations; iter++) {
        std::cout << "→ Iterazione " << (iter + 1) << ":\n";
        
        // 1. Calcola effemeridi per tutte le osservazioni
        std::vector<EquatorialCoordinates> computed;
        computeEphemerides(currentElements, validObs, computed);
        
        // 2. Aggiorna residui nelle osservazioni
        for (size_t i = 0; i < validObs.observations.size(); i++) {
            auto& obs = validObs.observations[i];
            obs.computed = computed[i];
            
            // Calcola residui in arcsec
            double dRA = (obs.obs.ra - obs.computed.ra) * cos(obs.obs.dec) * RAD_TO_DEG * 3600.0;
            double dDec = (obs.obs.dec - obs.computed.dec) * RAD_TO_DEG * 3600.0;
            
            obs.raResidual = dRA;
            obs.decResidual = dDec;
            obs.totalResidual = sqrt(dRA * dRA + dDec * dDec);
        }
        
        // 3. Calcola RMS dei residui
        double sumRA = 0, sumDec = 0;
        int nValid = 0;
        for (const auto& obs : validObs.observations) {
            if (!obs.outlier) {
                sumRA += obs.raResidual * obs.raResidual;
                sumDec += obs.decResidual * obs.decResidual;
                nValid++;
            }
        }
        
        double rmsRA = sqrt(sumRA / nValid);
        double rmsDec = sqrt(sumDec / nValid);
        double rmsTotal = sqrt((sumRA + sumDec) / (2 * nValid));
        
        std::cout << "   RMS: RA=" << std::fixed << std::setprecision(3) 
                  << rmsRA << "\" Dec=" << rmsDec << "\" Total=" << rmsTotal << "\"\n";
        
        result.rmsResidual = rmsTotal;
        result.nObservations = nValid;
        
        // 4. Test convergenza
        if (fabs(rmsTotal - previousRMS) < options.convergenceThreshold) {
            std::cout << "\n✓ Convergenza raggiunta! ΔRM S < " 
                      << options.convergenceThreshold << "\"\n";
            result.converged = true;
            result.convergenceMessage = "Converged: RMS change below threshold";
            break;
        }
        
        // 5. Rigetta outlier se richiesto
        if (options.rejectOutliers && iter > 0) {
            int nOutliersBefore = 0;
            for (auto& obs : validObs.observations) {
                if (obs.outlier) nOutliersBefore++;
            }
            
            detectOutliers(validObs, options.outlierSigma);
            
            int nOutliersAfter = 0;
            for (auto& obs : validObs.observations) {
                if (obs.outlier) nOutliersAfter++;
            }
            
            if (nOutliersAfter > nOutliersBefore) {
                std::cout << "   Rigettati " << (nOutliersAfter - nOutliersBefore) 
                          << " outlier (>" << options.outlierSigma << "σ)\n";
                result.nOutliers = nOutliersAfter;
            }
        }
        
        // 6. Calcola matrice Jacobiana (derivate parziali)
        std::vector<std::vector<double>> jacobian;
        computeJacobian(currentElements, validObs, jacobian);
        
        // 7. Costruisci sistema normale A^T A x = A^T b
        const int nParams = 6;  // a, e, i, Omega, omega, M
        std::vector<std::vector<double>> ATA(nParams, std::vector<double>(nParams, 0.0));
        std::vector<double> ATb(nParams, 0.0);
        
        for (size_t i = 0; i < validObs.observations.size(); i++) {
            if (validObs.observations[i].outlier) continue;
            
            double weight = options.useWeights ? validObs.observations[i].getWeight() : 1.0;
            double resRA = validObs.observations[i].raResidual;
            double resDec = validObs.observations[i].decResidual;
            
            for (int j = 0; j < nParams; j++) {
                ATb[j] += weight * (jacobian[2*i][j] * resRA + jacobian[2*i+1][j] * resDec);
                
                for (int k = 0; k < nParams; k++) {
                    ATA[j][k] += weight * (jacobian[2*i][j] * jacobian[2*i][k] + 
                                          jacobian[2*i+1][j] * jacobian[2*i+1][k]);
                }
            }
        }
        
        // 8. Risolvi sistema normale
        std::vector<double> corrections(nParams);
        if (!solveNormalEquations(ATA, ATb, corrections)) {
            std::cout << "\n⚠  Sistema normale singolare - stop iterazioni\n";
            result.converged = false;
            result.convergenceMessage = "Failed: Singular normal equations";
            break;
        }
        
        std::cout << "   Correzioni: ";
        std::cout << "Δa=" << std::scientific << std::setprecision(2) << corrections[0] << " AU ";
        std::cout << "Δe=" << corrections[1] << " ";
        std::cout << "Δi=" << corrections[2] * RAD_TO_DEG << "°\n";
        
        // 9. Aggiorna elementi
        currentElements = updateElements(currentElements, corrections);
        
        // 10. Callback progresso
        result.iterations = iter + 1;
        result.fittedElements = currentElements;
        if (progressCallback_) {
            progressCallback_(iter + 1, result);
        }
        
        previousRMS = rmsTotal;
    }
    
    // Calcola matrice di covarianza finale
    if (result.converged) {
        result.covariance = computeCovariance(currentElements, validObs);
        
        // Estrai incertezze formali
        if (!result.covariance.empty()) {
            result.sigma_a = sqrt(result.covariance[0][0]);
            result.sigma_e = sqrt(result.covariance[1][1]);
            result.sigma_i = sqrt(result.covariance[2][2]);
            result.sigma_Omega = sqrt(result.covariance[3][3]);
            result.sigma_omega = sqrt(result.covariance[4][4]);
            result.sigma_M = sqrt(result.covariance[5][5]);
        }
    }
    
    result.observations = validObs;
    
    std::cout << "\n═══════════════════════════════════════════════════════════\n";
    std::cout << "RISULTATI FINALI:\n";
    std::cout << "  Iterazioni: " << result.iterations << "\n";
    std::cout << "  RMS finale: " << std::fixed << std::setprecision(3) 
              << result.rmsResidual << "\"\n";
    std::cout << "  Osservazioni: " << result.nObservations << " (outlier: " 
              << result.nOutliers << ")\n";
    std::cout << "  Convergenza: " << (result.converged ? "SÌ" : "NO") << "\n";
    std::cout << "═══════════════════════════════════════════════════════════\n\n";
    
    return result;
}

void OrbitFitter::computeEphemerides(const OrbitalElements& elements,
                                    const ObservationSet& observations,
                                    std::vector<EquatorialCoordinates>& computed) {
    computed.clear();
    computed.reserve(observations.observations.size());
    
    // Converti elementi in EquinoctialElements e salva per propagazione
    EquinoctialElements eqElements = elements.toEquinoctial();
    pImpl->currentEqElements = eqElements;
    
    // DEBUG: conta osservatori unici (prima iterazione)
    static bool debug_printed = false;
    if (!debug_printed) {
        std::map<std::string, int> obsCounts;
        for (const auto& obs : observations.observations) {
            obsCounts[obs.observatoryCode]++;
        }
        std::cout << "\nDEBUG: Codici osservatorio nelle osservazioni:\n";
        for (const auto& [code, count] : obsCounts) {
            Observatory o = Observatory::fromMPCCode(code);
            std::cout << "  " << code << ": " << count << " obs";
            if (o.rho_cos_phi != 0.0 || o.rho_sin_phi != 0.0) {
                std::cout << " (ρcos=" << o.rho_cos_phi << ", ρsin=" << o.rho_sin_phi << ")";
            } else {
                std::cout << " [NO COORDS]";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        debug_printed = true;
    }
    
    // PERFORMANCE: Per evitare di propagare da epoca nominale a ogni osservazione
    // (molto lento con 2215 obs e RA15), usiamo un approccio più efficiente:
    // 1. Prima iterazione: usa Ephemeris (veloce ma impreciso)
    // 2. Iterazioni successive: potremmo usare OrbitPropagator
    // Per ora usiamo sempre Ephemeris per velocità
    
    pImpl->ephemeris.setElements(eqElements);
    
    for (const auto& obs : observations.observations) {
        // Usa Ephemeris (Keplerian 2-body) per velocità
        // TODO: Implementare propagazione efficiente con restart per usare RA15
        EphemerisData ephData = pImpl->ephemeris.compute(obs.epoch);
        
        // Posizione eliocentrica asteroide
        Vector3D asteroidHelio = ephData.heliocentricPos;
        
        // Posizione Terra
        Vector3D earthPos = Ephemeris::getEarthPosition(obs.epoch);
        
        // Posizione geocentrica asteroide
        Vector3D asteroidGeo;
        asteroidGeo.x = asteroidHelio.x - earthPos.x;
        asteroidGeo.y = asteroidHelio.y - earthPos.y;
        asteroidGeo.z = asteroidHelio.z - earthPos.z;
        
        // Distanza geocentrica
        double rho_geo = std::sqrt(asteroidGeo.x*asteroidGeo.x + 
                                   asteroidGeo.y*asteroidGeo.y + 
                                   asteroidGeo.z*asteroidGeo.z);
        
        // RA/Dec geocentriche
        double ra_geo = std::atan2(asteroidGeo.y, asteroidGeo.x);
        if (ra_geo < 0) ra_geo += 2.0 * M_PI;
        double dec_geo = std::asin(asteroidGeo.z / rho_geo);
        
        EquatorialCoordinates geocentricPos;
        geocentricPos.ra = ra_geo;
        geocentricPos.dec = dec_geo;
        
        // Correzione topocentrica se osservatorio non è geocentro
        if (obs.observatoryCode != "500") {  // 500 = geocentro
            try {
                // Ottieni coordinate osservatorio dal database MPC
                Observatory observatory = Observatory::fromMPCCode(obs.observatoryCode);
                
                if (!observatory.isSpacecraft && 
                    (observatory.rho_cos_phi != 0.0 || observatory.rho_sin_phi != 0.0)) {
                    
                    // MPC fornisce coordinate geocentriche già in frame equatoriale:
                    // ρ·cos(φ') = componente nel piano equatoriale
                    // ρ·sin(φ') = componente lungo asse polare
                    // Unità = raggio equatoriale terrestre (a_earth = 6378.137 km)
                    
                    constexpr double a_earth = 6378.137;  // km
                    constexpr double KM_TO_AU = 1.0 / 149597870.7;
                    
                    // Calcola GMST (Greenwich Mean Sidereal Time) per rotazione Terra
                    double T = (obs.epoch.jd - 2451545.0) / 36525.0;  // secoli da J2000
                    double GMST_deg = 280.46061837 + 360.98564736629 * (obs.epoch.jd - 2451545.0)
                                    + 0.000387933 * T * T - T * T * T / 38710000.0;
                    double GMST_rad = std::fmod(GMST_deg, 360.0) * M_PI / 180.0;
                    
                    // Longitudine osservatorio
                    double lon_rad = observatory.location.longitude * M_PI / 180.0;
                    
                    // Vettore posizione osservatorio in ICRF (J2000)
                    double cos_gmst = std::cos(GMST_rad + lon_rad);
                    double sin_gmst = std::sin(GMST_rad + lon_rad);
                    
                    double x_obs_earth_radii = observatory.rho_cos_phi * cos_gmst;
                    double y_obs_earth_radii = observatory.rho_cos_phi * sin_gmst;
                    double z_obs_earth_radii = observatory.rho_sin_phi;
                    
                    // Converti in AU
                    double x_obs_au = x_obs_earth_radii * a_earth * KM_TO_AU;
                    double y_obs_au = y_obs_earth_radii * a_earth * KM_TO_AU;
                    double z_obs_au = z_obs_earth_radii * a_earth * KM_TO_AU;
                    
                    // Posizione asteroide topocentrica (già abbiamo asteroidGeo calcolato sopra)
                    Vector3D asteroidTopo;
                    asteroidTopo.x = asteroidGeo.x - x_obs_au;
                    asteroidTopo.y = asteroidGeo.y - y_obs_au;
                    asteroidTopo.z = asteroidGeo.z - z_obs_au;
                    
                    // Calcola distanza topocentrica
                    double rho_topo = std::sqrt(asteroidTopo.x*asteroidTopo.x + 
                                               asteroidTopo.y*asteroidTopo.y + 
                                               asteroidTopo.z*asteroidTopo.z);
                    
                    // Converti in coordinate sferiche (RA/Dec topocentriche)
                    // Salva valori geocentrici per debug
                    double ra_geo_before = geocentricPos.ra;
                    double dec_geo_before = geocentricPos.dec;
                    
                    geocentricPos.ra = std::atan2(asteroidTopo.y, asteroidTopo.x);
                    geocentricPos.dec = std::asin(asteroidTopo.z / rho_topo);
                    
                    // Normalizza RA in [0, 2π]
                    if (geocentricPos.ra < 0) geocentricPos.ra += 2.0 * M_PI;
                }
            } catch (...) {
                // Se osservatorio non trovato, usa posizione geocentrica
            }
        }
        
        // TODO: Correzione per aberrazione stellare (~20 arcsec)
        // TODO: Light-time correction (iterativo, converge in 2-3 step)
        
        computed.push_back(geocentricPos);
    }
}

void OrbitFitter::computeJacobian(const OrbitalElements& elements,
                                 const ObservationSet& observations,
                                 std::vector<std::vector<double>>& jacobian) {
    // Matrice Jacobiana: derivate parziali di (RA, Dec) rispetto agli elementi
    // Dimensioni: (2*nObs) x 6
    
    const int nObs = observations.observations.size();
    const int nParams = 6;
    jacobian.clear();
    jacobian.resize(2 * nObs, std::vector<double>(nParams, 0.0));
    
    // Calcola derivate numeriche con differenze finite
    const double delta_a = elements.a * 1e-8;
    const double delta_e = 1e-8;
    const double delta_i = 1e-8;
    const double delta_Omega = 1e-8;
    const double delta_omega = 1e-8;
    const double delta_M = 1e-8;
    
    std::vector<double> deltas = {delta_a, delta_e, delta_i, delta_Omega, delta_omega, delta_M};
    
    // Effemeridi nominali
    std::vector<EquatorialCoordinates> nominal;
    computeEphemerides(elements, observations, nominal);
    
    // Per ogni parametro
    for (int param = 0; param < nParams; param++) {
        OrbitalElements perturbed = elements;
        
        // Perturba il parametro
        switch (param) {
            case 0: perturbed.a += deltas[param]; break;
            case 1: perturbed.e += deltas[param]; break;
            case 2: perturbed.i += deltas[param]; break;
            case 3: perturbed.Omega += deltas[param]; break;
            case 4: perturbed.omega += deltas[param]; break;
            case 5: perturbed.M += deltas[param]; break;
        }
        
        // Calcola effemeridi perturbate
        std::vector<EquatorialCoordinates> perturbed_eph;
        computeEphemerides(perturbed, observations, perturbed_eph);
        
        // Calcola derivate
        for (int i = 0; i < nObs; i++) {
            double dRA = (perturbed_eph[i].ra - nominal[i].ra) / deltas[param];
            double dDec = (perturbed_eph[i].dec - nominal[i].dec) / deltas[param];
            
            // Converti in arcsec
            dRA *= cos(nominal[i].dec) * RAD_TO_DEG * 3600.0;
            dDec *= RAD_TO_DEG * 3600.0;
            
            jacobian[2*i][param] = dRA;
            jacobian[2*i+1][param] = dDec;
        }
    }
}

bool OrbitFitter::solveNormalEquations(const std::vector<std::vector<double>>& A,
                                      const std::vector<double>& b,
                                      std::vector<double>& x) {
    // Risolve A x = b usando eliminazione di Gauss con pivoting parziale
    
    int n = b.size();
    x.resize(n);
    
    // Copia matrice aumentata [A|b]
    std::vector<std::vector<double>> aug(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][n] = b[i];
    }
    
    // Eliminazione forward con pivoting
    for (int k = 0; k < n; k++) {
        // Trova pivot
        int maxRow = k;
        double maxVal = fabs(aug[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(aug[i][k]) > maxVal) {
                maxVal = fabs(aug[i][k]);
                maxRow = i;
            }
        }
        
        if (maxVal < 1e-15) {
            return false;  // Matrice singolare
        }
        
        // Swap righe
        if (maxRow != k) {
            std::swap(aug[k], aug[maxRow]);
        }
        
        // Eliminazione
        for (int i = k + 1; i < n; i++) {
            double factor = aug[i][k] / aug[k][k];
            for (int j = k; j <= n; j++) {
                aug[i][j] -= factor * aug[k][j];
            }
        }
    }
    
    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        x[i] = aug[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    
    return true;
}

OrbitalElements OrbitFitter::updateElements(const OrbitalElements& elements,
                                           const std::vector<double>& corrections) {
    OrbitalElements updated = elements;
    
    updated.a += corrections[0];
    updated.e += corrections[1];
    updated.i += corrections[2];
    updated.Omega += corrections[3];
    updated.omega += corrections[4];
    updated.M += corrections[5];
    
    // Normalizza angoli
    while (updated.Omega < 0) updated.Omega += 2.0 * M_PI;
    while (updated.Omega >= 2.0 * M_PI) updated.Omega -= 2.0 * M_PI;
    
    while (updated.omega < 0) updated.omega += 2.0 * M_PI;
    while (updated.omega >= 2.0 * M_PI) updated.omega -= 2.0 * M_PI;
    
    while (updated.M < 0) updated.M += 2.0 * M_PI;
    while (updated.M >= 2.0 * M_PI) updated.M -= 2.0 * M_PI;
    
    return updated;
}

void OrbitFitter::detectOutliers(ObservationSet& observations, double sigmaCutoff) {
    // Calcola RMS
    double sumRA = 0, sumDec = 0;
    int nValid = 0;
    
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) {
            sumRA += obs.raResidual * obs.raResidual;
            sumDec += obs.decResidual * obs.decResidual;
            nValid++;
        }
    }
    
    double rmsRA = sqrt(sumRA / nValid);
    double rmsDec = sqrt(sumDec / nValid);
    
    // Marca outlier
    for (auto& obs : observations.observations) {
        if (obs.outlier) continue;
        
        double sigmaRA = fabs(obs.raResidual) / rmsRA;
        double sigmaDec = fabs(obs.decResidual) / rmsDec;
        
        if (sigmaRA > sigmaCutoff || sigmaDec > sigmaCutoff) {
            obs.outlier = true;
        }
    }
}

std::vector<std::vector<double>> OrbitFitter::computeCovariance(
    const OrbitalElements& elements,
    const ObservationSet& observations) {
    
    // Calcola matrice di covarianza C = (A^T W A)^-1
    // dove W è la matrice dei pesi
    
    std::vector<std::vector<double>> jacobian;
    computeJacobian(elements, observations, jacobian);
    
    const int nParams = 6;
    std::vector<std::vector<double>> ATA(nParams, std::vector<double>(nParams, 0.0));
    
    for (size_t i = 0; i < observations.observations.size(); i++) {
        if (observations.observations[i].outlier) continue;
        
        double weight = observations.observations[i].getWeight();
        
        for (int j = 0; j < nParams; j++) {
            for (int k = 0; k < nParams; k++) {
                ATA[j][k] += weight * (jacobian[2*i][j] * jacobian[2*i][k] + 
                                      jacobian[2*i+1][j] * jacobian[2*i+1][k]);
            }
        }
    }
    
    // Inverti ATA per ottenere covarianza
    // TODO: Implementare inversione matrice
    // Per ora ritorna matrice vuota
    
    return std::vector<std::vector<double>>();
}

void OrbitFitter::computeResiduals(const OrbitalElements& elements,
                                  ObservationSet& observations) {
    std::vector<EquatorialCoordinates> computed;
    computeEphemerides(elements, observations, computed);
    
    for (size_t i = 0; i < observations.observations.size(); i++) {
        auto& obs = observations.observations[i];
        obs.computed = computed[i];
        
        double dRA = (obs.obs.ra - obs.computed.ra) * cos(obs.obs.dec) * RAD_TO_DEG * 3600.0;
        double dDec = (obs.obs.dec - obs.computed.dec) * RAD_TO_DEG * 3600.0;
        
        obs.raResidual = dRA;
        obs.decResidual = dDec;
        obs.totalResidual = sqrt(dRA * dRA + dDec * dDec);
    }
}

OrbitalElements OrbitFitter::propagate(const OrbitalElements& elements,
                                      const JulianDate& targetEpoch,
                                      bool includePerturbations) {
    // TODO: Implementare propagazione completa con perturbazioni
    // Per ora usa propagazione Kepleriana semplice
    
    // Converti GM_SUN da km³/s² a AU³/day²
    const double AU_TO_KM = 149597870.7;
    const double DAY_TO_SEC = 86400.0;
    const double k_gauss = PhysicalConstants::GM_SUN / (AU_TO_KM * AU_TO_KM * AU_TO_KM) * (DAY_TO_SEC * DAY_TO_SEC);
    
    double dt = targetEpoch.jd - elements.epoch.jd;
    double n = sqrt(k_gauss / (elements.a * elements.a * elements.a));  // mean motion (rad/day)
    
    OrbitalElements propagated = elements;
    propagated.M += n * dt;  // dt in giorni, n in rad/day
    propagated.epoch = targetEpoch;
    
    // Normalizza
    while (propagated.M >= 2.0 * M_PI) propagated.M -= 2.0 * M_PI;
    while (propagated.M < 0) propagated.M += 2.0 * M_PI;
    
    return propagated;
}

// ResidualAnalysis implementation
void ResidualAnalysis::computeStatistics(const ObservationSet& observations,
                                        double& meanRA, double& meanDec,
                                        double& rmsRA, double& rmsDec,
                                        double& rmsTotal) {
    meanRA = meanDec = 0;
    rmsRA = rmsDec = rmsTotal = 0;
    int nValid = 0;
    
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) {
            meanRA += obs.raResidual;
            meanDec += obs.decResidual;
            nValid++;
        }
    }
    
    if (nValid == 0) return;
    
    meanRA /= nValid;
    meanDec /= nValid;
    
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) {
            rmsRA += obs.raResidual * obs.raResidual;
            rmsDec += obs.decResidual * obs.decResidual;
            rmsTotal += obs.totalResidual * obs.totalResidual;
        }
    }
    
    rmsRA = sqrt(rmsRA / nValid);
    rmsDec = sqrt(rmsDec / nValid);
    rmsTotal = sqrt(rmsTotal / nValid);
}

bool ResidualAnalysis::detectSystematicBias(const ObservationSet& observations,
                                           double& biasRA, double& biasDec) {
    double meanRA, meanDec, rmsRA, rmsDec, rmsTotal;
    computeStatistics(observations, meanRA, meanDec, rmsRA, rmsDec, rmsTotal);
    
    biasRA = meanRA;
    biasDec = meanDec;
    
    // Test t: bias significativo se |mean| > 2*rms/sqrt(n)
    int n = 0;
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) n++;
    }
    
    double thresholdRA = 2.0 * rmsRA / sqrt(n);
    double thresholdDec = 2.0 * rmsDec / sqrt(n);
    
    return (fabs(meanRA) > thresholdRA || fabs(meanDec) > thresholdDec);
}

std::vector<int> ResidualAnalysis::residualHistogram(const ObservationSet& observations,
                                                     int nBins, double& binWidth) {
    std::vector<int> histogram(nBins, 0);
    
    // Trova min/max residuo
    double minRes = 1e10, maxRes = -1e10;
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) {
            minRes = std::min(minRes, obs.totalResidual);
            maxRes = std::max(maxRes, obs.totalResidual);
        }
    }
    
    binWidth = (maxRes - minRes) / nBins;
    
    for (const auto& obs : observations.observations) {
        if (!obs.outlier) {
            int bin = (int)((obs.totalResidual - minRes) / binWidth);
            if (bin >= nBins) bin = nBins - 1;
            if (bin < 0) bin = 0;
            histogram[bin]++;
        }
    }
    
    return histogram;
}

} // namespace ioccultcalc
