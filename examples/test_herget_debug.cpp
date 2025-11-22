/**
 * @file test_herget_debug.cpp
 * @brief Debug dettagliato del metodo Herget
 * 
 * Traccia ogni passaggio per trovare dove si introduce l'errore di inclinazione
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/mpc_client.h"
#include "ioccultcalc/ephemeris.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

constexpr double PI = 3.14159265358979323846;
constexpr double DEG = PI / 180.0;

void printVector(const std::string& label, const Vector3D& v) {
    std::cout << std::setw(30) << std::left << label << ": [" 
              << std::fixed << std::setprecision(8)
              << std::setw(12) << v.x << ", " 
              << std::setw(12) << v.y << ", " 
              << std::setw(12) << v.z << "]"
              << "  |r|=" << std::setprecision(6) << v.magnitude() << std::endl;
}

void printElements(const std::string& label, const double* state, double mu = 0.01720209895 * 0.01720209895) {
    // Calcola elementi manualmente
    Vector3D r(state[0], state[1], state[2]);
    Vector3D v(state[3], state[4], state[5]);
    
    double rmag = r.magnitude();
    double vmag = v.magnitude();
    
    Vector3D h = r.cross(v);
    double hmag = h.magnitude();
    
    double energy = vmag * vmag / 2.0 - mu / rmag;
    double a = -mu / (2.0 * energy);
    
    Vector3D evec = (v.cross(h)) / mu - r / rmag;
    double e = evec.magnitude();
    
    double i = acos(h.z / hmag);
    
    std::cout << "\n" << label << ":\n";
    std::cout << "  a = " << std::fixed << std::setprecision(6) << a << " AU\n";
    std::cout << "  e = " << e << "\n";
    std::cout << "  i = " << (i * 180.0 / PI) << "°\n";
    std::cout << "  |h| = " << hmag << " AU²/day\n";
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║     DEBUG DETTAGLIATO METODO HERGET                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";
    
    // Download osservazioni
    MPCClient mpc;
    ObservationSet obs_set = mpc.getObservations("433");
    
    if (obs_set.observations.empty()) {
        std::cerr << "✗ Errore download osservazioni\n";
        return 1;
    }
    
    auto& observations = obs_set.observations;
    std::cout << "✓ Scaricate " << observations.size() << " osservazioni\n\n";
    
    // Seleziona 2 osservazioni recenti (span ~50 giorni)
    size_t n = observations.size();
    size_t idx1 = n - 60;
    size_t idx2 = n - 10;
    
    auto& obs1 = observations[idx1];
    auto& obs2 = observations[idx2];
    
    double span = obs2.epoch.jd - obs1.epoch.jd;
    
    std::cout << "============================================================\n";
    std::cout << "OSSERVAZIONI SELEZIONATE\n";
    std::cout << "============================================================\n\n";
    std::cout << "Osservazione 1:\n";
    std::cout << "  JD:  " << std::fixed << std::setprecision(6) << obs1.epoch.jd << "\n";
    std::cout << "  RA:  " << (obs1.obs.ra * 180.0 / PI) << "°\n";
    std::cout << "  Dec: " << (obs1.obs.dec * 180.0 / PI) << "°\n";
    std::cout << "  Obs: " << obs1.observatoryCode << "\n\n";
    
    std::cout << "Osservazione 2:\n";
    std::cout << "  JD:  " << obs2.epoch.jd << "\n";
    std::cout << "  RA:  " << (obs2.obs.ra * 180.0 / PI) << "°\n";
    std::cout << "  Dec: " << (obs2.obs.dec * 180.0 / PI) << "°\n";
    std::cout << "  Obs: " << obs2.observatoryCode << "\n\n";
    
    std::cout << "Span: " << span << " giorni\n\n";
    
    std::cout << "============================================================\n";
    std::cout << "STEP 1: POSIZIONI OSSERVATORI (EQUATORIALI)\n";
    std::cout << "============================================================\n\n";
    
    Vector3D obs1_pos_eq = Coordinates::observerPosition(obs1.observatoryCode, obs1.epoch);
    Vector3D obs2_pos_eq = Coordinates::observerPosition(obs2.observatoryCode, obs2.epoch);
    
    printVector("Osservatore 1 (equatoriali)", obs1_pos_eq);
    printVector("Osservatore 2 (equatoriali)", obs2_pos_eq);
    
    std::cout << "\n============================================================\n";
    std::cout << "STEP 2: DIREZIONI VERSO OGGETTO (EQUATORIALI)\n";
    std::cout << "============================================================\n\n";
    
    Vector3D dir1_eq = Coordinates::raDecToUnitVector(obs1.obs.ra, obs1.obs.dec);
    Vector3D dir2_eq = Coordinates::raDecToUnitVector(obs2.obs.ra, obs2.obs.dec);
    
    printVector("Direzione 1 (equatoriali)", dir1_eq);
    printVector("Direzione 2 (equatoriali)", dir2_eq);
    
    std::cout << "\n============================================================\n";
    std::cout << "STEP 3: CONVERSIONE A COORDINATE ECLITTICHE\n";
    std::cout << "============================================================\n\n";
    
    Vector3D obs1_pos_ecl = Coordinates::equatorialToEcliptic(obs1_pos_eq);
    Vector3D obs2_pos_ecl = Coordinates::equatorialToEcliptic(obs2_pos_eq);
    Vector3D dir1_ecl = Coordinates::equatorialToEcliptic(dir1_eq);
    Vector3D dir2_ecl = Coordinates::equatorialToEcliptic(dir2_eq);
    
    printVector("Osservatore 1 (eclittiche)", obs1_pos_ecl);
    printVector("Osservatore 2 (eclittiche)", obs2_pos_ecl);
    std::cout << "\n";
    printVector("Direzione 1 (eclittiche)", dir1_ecl);
    printVector("Direzione 2 (eclittiche)", dir2_ecl);
    
    std::cout << "\nComponente Z (dovrebbe essere ~0 per eclittiche):\n";
    std::cout << "  Osservatore 1 z: " << obs1_pos_ecl.z << " AU\n";
    std::cout << "  Osservatore 2 z: " << obs2_pos_ecl.z << " AU\n";
    
    std::cout << "\n============================================================\n";
    std::cout << "STEP 4: POSIZIONI HELIOCENTRICHE CON DISTANZE GUESS\n";
    std::cout << "============================================================\n\n";
    
    // Test con distanze tipiche per Eros (main belt interno)
    double r1_guess = 1.4; // AU
    double r2_guess = 1.5; // AU
    
    std::cout << "Distanze guess: R1=" << r1_guess << " AU, R2=" << r2_guess << " AU\n\n";
    
    Vector3D pos1_ecl = obs1_pos_ecl + dir1_ecl * r1_guess;
    Vector3D pos2_ecl = obs2_pos_ecl + dir2_ecl * r2_guess;
    
    printVector("Posizione 1 heliocentrica (ecl)", pos1_ecl);
    printVector("Posizione 2 heliocentrica (ecl)", pos2_ecl);
    
    std::cout << "\n============================================================\n";
    std::cout << "STEP 5: LAMBERT SOLVER\n";
    std::cout << "============================================================\n\n";
    
    double mu = 0.01720209895 * 0.01720209895;
    Vector3D v1_ecl, v2_ecl;
    
    int result = InitialOrbit::solveLambert(pos1_ecl, pos2_ecl, span, mu, true, v1_ecl, v2_ecl);
    
    if (result == 0) {
        std::cout << "✓ Lambert solver convergenza OK\n\n";
        printVector("Velocità 1 (eclittiche)", v1_ecl);
        printVector("Velocità 2 (eclittiche)", v2_ecl);
        
        // Elementi orbitali in eclittiche
        double state_ecl[6] = {pos1_ecl.x, pos1_ecl.y, pos1_ecl.z, v1_ecl.x, v1_ecl.y, v1_ecl.z};
        printElements("Elementi orbitali (ECLITTICHE)", state_ecl, mu);
        
        std::cout << "\n============================================================\n";
        std::cout << "STEP 6: VERIFICA - CONVERSIONE IN EQUATORIALI\n";
        std::cout << "============================================================\n\n";
        
        // Converti in equatoriali per vedere cosa succede
        Vector3D pos1_eq_back = Coordinates::eclipticToEquatorial(pos1_ecl);
        Vector3D v1_eq_back = Coordinates::eclipticToEquatorial(v1_ecl);
        
        printVector("Posizione 1 (eq da ecl)", pos1_eq_back);
        printVector("Velocità 1 (eq da ecl)", v1_eq_back);
        
        double state_eq[6] = {pos1_eq_back.x, pos1_eq_back.y, pos1_eq_back.z, 
                               v1_eq_back.x, v1_eq_back.y, v1_eq_back.z};
        printElements("Elementi orbitali (EQUATORIALI)", state_eq, mu);
        
        std::cout << "\n============================================================\n";
        std::cout << "CONFRONTO CON EROS REALE\n";
        std::cout << "============================================================\n\n";
        
        std::cout << "Eros reale (eclittiche):\n";
        std::cout << "  a = 1.458 AU\n";
        std::cout << "  e = 0.223\n";
        std::cout << "  i = 10.83°\n\n";
        
        // Calcola errori
        Vector3D r_test(state_ecl[0], state_ecl[1], state_ecl[2]);
        Vector3D v_test(state_ecl[3], state_ecl[4], state_ecl[5]);
        Vector3D h_test = r_test.cross(v_test);
        double i_test = acos(h_test.z / h_test.magnitude()) * 180.0 / PI;
        
        std::cout << "Inclinazione calcolata (eclittiche): " << i_test << "°\n";
        std::cout << "Errore: " << std::abs(i_test - 10.83) << "°\n\n";
        
        if (std::abs(i_test - 10.83) < 5.0) {
            std::cout << "✓ ECCELLENTE: Inclinazione corretta!\n";
        } else if (std::abs(i_test - 10.83) < 15.0) {
            std::cout << "✓ BUONO: Inclinazione ragionevole\n";
        } else {
            std::cout << "✗ ERRORE: Inclinazione troppo diversa\n";
            std::cout << "\nPossibili cause:\n";
            std::cout << "  1. Direzioni non convertite correttamente\n";
            std::cout << "  2. Posizioni osservatore sbagliate\n";
            std::cout << "  3. Distanze guess non appropriate\n";
        }
        
    } else {
        std::cout << "✗ Lambert solver FALLITO\n";
        std::cout << "   Prova distanze diverse o span più lungo\n";
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "============================================================\n";
    
    return 0;
}
