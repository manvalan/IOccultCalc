/**
 * @file test_orbfit_propagation.cpp
 * @brief Test propagazione con modello forze OrbFit
 * 
 * Usa direttamente OrbFitForceModel con integratore RK4 semplice
 * per verificare che il modello di forze sia corretto.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "ioccultcalc/orbfit_force_model.h"
#include "ioccultcalc/orbital_elements.h"

using namespace ioccultcalc;

const double RAD2DEG = 180.0 / M_PI;
const double DEG2RAD = M_PI / 180.0;
const double AU_KM = 149597870.7;

// Integratore RK4 semplice
struct IntegratorState {
    Vector3D pos;  // Posizione eclittica [AU]
    Vector3D vel;  // Velocità eclittica [AU/day]
};

IntegratorState rk4Step(const OrbFitForceModel& model, 
                        const IntegratorState& state,
                        double jd,
                        double dt) {
    // k1
    Vector3D a1 = model.computeAcceleration(state.pos, state.vel, jd);
    Vector3D v1 = state.vel;
    
    // k2
    Vector3D pos2 = state.pos + v1 * (dt/2);
    Vector3D vel2 = state.vel + a1 * (dt/2);
    Vector3D a2 = model.computeAcceleration(pos2, vel2, jd + dt/2);
    
    // k3
    Vector3D pos3 = state.pos + vel2 * (dt/2);
    Vector3D vel3 = state.vel + a2 * (dt/2);
    Vector3D a3 = model.computeAcceleration(pos3, vel3, jd + dt/2);
    
    // k4
    Vector3D pos4 = state.pos + vel3 * dt;
    Vector3D vel4 = state.vel + a3 * dt;
    Vector3D a4 = model.computeAcceleration(pos4, vel4, jd + dt);
    
    // Combina
    IntegratorState newState;
    newState.pos = state.pos + (v1 + vel2 * 2.0 + vel3 * 2.0 + vel4) * (dt / 6.0);
    newState.vel = state.vel + (a1 + a2 * 2.0 + a3 * 2.0 + a4) * (dt / 6.0);
    
    return newState;
}

// Propaga per un intervallo temporale
IntegratorState propagate(const OrbFitForceModel& model,
                          const IntegratorState& initial,
                          double jdStart,
                          double jdEnd,
                          double stepSize = 0.1) {
    IntegratorState state = initial;
    double jd = jdStart;
    double dt = (jdEnd > jdStart) ? stepSize : -stepSize;
    int nSteps = 0;
    
    while ((dt > 0 && jd < jdEnd) || (dt < 0 && jd > jdEnd)) {
        // Ultimo step esatto
        if ((dt > 0 && jd + dt > jdEnd) || (dt < 0 && jd + dt < jdEnd)) {
            dt = jdEnd - jd;
        }
        
        state = rk4Step(model, state, jd, dt);
        jd += dt;
        nSteps++;
        
        if (nSteps % 1000 == 0) {
            std::cout << "  Step " << nSteps << ", JD = " << jd << "\r" << std::flush;
        }
    }
    
    std::cout << "  Completati " << nSteps << " steps\n";
    return state;
}

// Converti elementi Kepleriani a stato cartesiano (eclittico)
IntegratorState elementsToState(const OrbitalElements& elem) {
    using namespace OrbFitConstants;
    
    // Anomalia eccentrica da anomalia media
    double M = elem.M;
    double E = M;
    for (int i = 0; i < 20; ++i) {
        E = M + elem.e * sin(E);
    }
    
    // Posizione e velocità nel piano orbitale
    double cosE = cos(E);
    double sinE = sin(E);
    double r = elem.a * (1.0 - elem.e * cosE);
    double x_orb = elem.a * (cosE - elem.e);
    double y_orb = elem.a * sqrt(1.0 - elem.e * elem.e) * sinE;
    
    // Velocità nel piano orbitale
    double n = sqrt(GM0 / (elem.a * elem.a * elem.a));  // Moto medio
    double factor = n * elem.a / r;
    double vx_orb = -factor * elem.a * sinE;
    double vy_orb = factor * elem.a * sqrt(1.0 - elem.e * elem.e) * cosE;
    
    // Matrici di rotazione
    double cosW = cos(elem.omega);
    double sinW = sin(elem.omega);
    double cosO = cos(elem.Omega);
    double sinO = sin(elem.Omega);
    double cosI = cos(elem.i);
    double sinI = sin(elem.i);
    
    // Rotazione a eclittico
    IntegratorState state;
    
    // Posizione
    state.pos.x = (cosO * cosW - sinO * sinW * cosI) * x_orb 
                + (-cosO * sinW - sinO * cosW * cosI) * y_orb;
    state.pos.y = (sinO * cosW + cosO * sinW * cosI) * x_orb 
                + (-sinO * sinW + cosO * cosW * cosI) * y_orb;
    state.pos.z = sinW * sinI * x_orb + cosW * sinI * y_orb;
    
    // Velocità
    state.vel.x = (cosO * cosW - sinO * sinW * cosI) * vx_orb 
                + (-cosO * sinW - sinO * cosW * cosI) * vy_orb;
    state.vel.y = (sinO * cosW + cosO * sinW * cosI) * vx_orb 
                + (-sinO * sinW + cosO * cosW * cosI) * vy_orb;
    state.vel.z = sinW * sinI * vx_orb + cosW * sinI * vy_orb;
    
    return state;
}

int main() {
    std::cout << "=============================================\n";
    std::cout << " TEST PROPAGAZIONE DIRETTA ORBFIT\n";
    std::cout << "=============================================\n\n";
    
    // Inizializza modello forze
    OrbFitForceOptions opts = OrbFitForceOptions::highPrecisionConfig();
    opts.relativityLevel = 1;  // Solo Schwarzschild per ora (più stabile)
    opts.includeJ2Sun = false;
    
    OrbFitForceModel forceModel(opts);
    if (!forceModel.initializeEphemeris()) {
        std::cerr << "ERRORE: Impossibile inizializzare effemeridi\n";
        return 1;
    }
    
    std::cout << "Modello forze:\n";
    std::cout << "  Relatività: Schwarzschild (solo Sole)\n";
    std::cout << "  Pianeti: tutti tranne Plutone\n\n";
    
    // Elementi Ceres MPC (epoca Nov 2026)
    OrbitalElements ceres;
    ceres.epoch = JulianDate(2461000.5);  // 18 Nov 2026
    ceres.a = 2.7656157;
    ceres.e = 0.0795763;
    ceres.i = 10.58789 * DEG2RAD;
    ceres.Omega = 80.24963 * DEG2RAD;
    ceres.omega = 73.29974 * DEG2RAD;
    ceres.M = 231.53975 * DEG2RAD;
    
    std::cout << "Elementi Ceres (epoca Nov 2026):\n";
    std::cout << "  a = " << ceres.a << " AU\n";
    std::cout << "  e = " << ceres.e << "\n";
    std::cout << "  i = " << (ceres.i * RAD2DEG) << " deg\n\n";
    
    // Stato iniziale (eclittico)
    IntegratorState stateEpoch = elementsToState(ceres);
    
    std::cout << "Stato all'epoca (eclittico):\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Pos: [" << stateEpoch.pos.x << ", " << stateEpoch.pos.y 
              << ", " << stateEpoch.pos.z << "] AU\n";
    std::cout << "  Vel: [" << stateEpoch.vel.x << ", " << stateEpoch.vel.y 
              << ", " << stateEpoch.vel.z << "] AU/day\n";
    std::cout << "  r = " << stateEpoch.pos.magnitude() << " AU\n\n";
    
    // === TEST 1: Round-trip +30 / -30 giorni ===
    std::cout << "=== TEST 1: Round-trip +30 / -30 giorni ===\n";
    
    IntegratorState state30 = propagate(forceModel, stateEpoch, 
                                        ceres.epoch.jd, ceres.epoch.jd + 30, 0.1);
    IntegratorState stateBack = propagate(forceModel, state30, 
                                          ceres.epoch.jd + 30, ceres.epoch.jd, 0.1);
    
    double errRT = (stateBack.pos - stateEpoch.pos).magnitude();
    std::cout << "Errore round-trip: " << (errRT * AU_KM) << " km\n";
    
    if (errRT * AU_KM < 100) {
        std::cout << "[OK] Round-trip corretto\n\n";
    } else {
        std::cout << "[WARN] Errore round-trip significativo\n\n";
    }
    
    // === TEST 2: Propagazione -354 giorni (verso Nov 2025) ===
    std::cout << "=== TEST 2: Propagazione -354 giorni ===\n";
    
    double jdTarget = 2460646.5;  // 29 Nov 2025
    double deltaT = jdTarget - ceres.epoch.jd;
    std::cout << "Delta tempo: " << deltaT << " giorni\n";
    
    IntegratorState stateNov2025 = propagate(forceModel, stateEpoch,
                                             ceres.epoch.jd, jdTarget, 0.1);
    
    std::cout << "\nPosizione calcolata (29 Nov 2025):\n";
    std::cout << "  X = " << stateNov2025.pos.x << " AU\n";
    std::cout << "  Y = " << stateNov2025.pos.y << " AU\n";
    std::cout << "  Z = " << stateNov2025.pos.z << " AU\n";
    std::cout << "  r = " << stateNov2025.pos.magnitude() << " AU\n";
    
    // Riferimento JPL Horizons (29 Nov 2025, eclittico J2000)
    Vector3D jplPos(1.8629918, -2.1633047, -0.4890291);
    
    std::cout << "\nRiferimento JPL Horizons:\n";
    std::cout << "  X = " << jplPos.x << " AU\n";
    std::cout << "  Y = " << jplPos.y << " AU\n";
    std::cout << "  Z = " << jplPos.z << " AU\n";
    
    Vector3D errVec = stateNov2025.pos - jplPos;
    double errJPL = errVec.magnitude();
    
    std::cout << "\nErrore vs JPL:\n";
    std::cout << "  dX = " << errVec.x << " AU = " << (errVec.x * AU_KM) << " km\n";
    std::cout << "  dY = " << errVec.y << " AU = " << (errVec.y * AU_KM) << " km\n";
    std::cout << "  dZ = " << errVec.z << " AU = " << (errVec.z * AU_KM) << " km\n";
    std::cout << "  |dR| = " << errJPL << " AU = " << (errJPL * AU_KM) << " km\n";
    
    // === RIEPILOGO ===
    std::cout << "\n=============================================\n";
    std::cout << " RIEPILOGO\n";
    std::cout << "=============================================\n";
    std::cout << "Round-trip error: " << (errRT * AU_KM) << " km\n";
    std::cout << "Errore vs JPL: " << (errJPL * AU_KM) << " km\n";
    
    if (errJPL * AU_KM < 50000) {
        std::cout << "\n[OK] Propagazione migliorata!\n";
    } else if (errJPL * AU_KM < 100000) {
        std::cout << "\n[WARN] Errore ridotto ma ancora significativo\n";
    } else {
        std::cout << "\n[INFO] L'errore persiste - il problema è negli elementi MPC epoca futura\n";
        std::cout << "       Per validare, serve testare con elementi epoca passata.\n";
    }
    
    return 0;
}
