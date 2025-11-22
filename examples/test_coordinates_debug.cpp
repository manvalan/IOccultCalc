/**
 * @file test_coordinates_debug.cpp
 * @brief Debug delle coordinate e sistema di riferimento
 * 
 * Verifica che le osservazioni e le posizioni dell'osservatore siano
 * nel sistema di coordinate corretto (eclittica J2000 per orbite heliocentriche)
 */

#include "ioccultcalc/initial_orbit.h"
#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/mpc_client.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

constexpr double PI = 3.14159265358979323846;
constexpr double DEG = PI / 180.0;

void printVector(const std::string& label, const Vector3D& v) {
    std::cout << std::setw(25) << std::left << label << ": [" 
              << std::fixed << std::setprecision(6)
              << std::setw(10) << v.x << ", " 
              << std::setw(10) << v.y << ", " 
              << std::setw(10) << v.z << "]"
              << "  |r|=" << v.magnitude() << std::endl;
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║     DEBUG SISTEMA COORDINATE                           ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n\n";
    
    // Download osservazioni Eros
    std::cout << "============================================================\n";
    std::cout << "1. DOWNLOAD OSSERVAZIONI EROS\n";
    std::cout << "============================================================\n\n";
    
    MPCClient mpc;
    ObservationSet obs_set = mpc.getObservations("433");
    
    if (obs_set.observations.empty()) {
        std::cerr << "✗ Errore download osservazioni\n";
        return 1;
    }
    
    auto& observations = obs_set.observations;
    std::cout << "✓ Scaricate " << observations.size() << " osservazioni\n\n";
    
    // Prendi 3 osservazioni recenti ben distribuite
    if (observations.size() < 100) {
        std::cerr << "✗ Troppe poche osservazioni\n";
        return 1;
    }
    
    // Ultime 100 osservazioni, prendi prima/metà/ultima
    size_t n = observations.size();
    size_t idx1 = n - 100;
    size_t idx2 = n - 50;
    size_t idx3 = n - 1;
    
    std::cout << "============================================================\n";
    std::cout << "2. ANALISI OSSERVAZIONI SELEZIONATE\n";
    std::cout << "============================================================\n\n";
    
    for (size_t idx : {idx1, idx2, idx3}) {
        auto& obs = observations[idx];
        
        std::cout << "--- Osservazione " << (idx == idx1 ? "1" : idx == idx2 ? "2" : "3") << " ---\n";
        std::cout << "  Epoca:       JD " << std::fixed << std::setprecision(6) << obs.epoch.jd << "\n";
        std::cout << "  Observatory: " << obs.observatoryCode << "\n";
        std::cout << "  RA:          " << (obs.obs.ra * 180.0 / PI) << "° = " 
                  << (obs.obs.ra * 12.0 / PI) << " h\n";
        std::cout << "  Dec:         " << (obs.obs.dec * 180.0 / PI) << "°\n";
        
        // Verifica che RA sia in range [0, 2π)
        if (obs.obs.ra < 0 || obs.obs.ra >= 2.0 * PI) {
            std::cout << "  ⚠️  WARNING: RA fuori range [0, 2π)\n";
        }
        
        // Verifica che Dec sia in range [-π/2, π/2]
        if (obs.obs.dec < -PI/2 || obs.obs.dec > PI/2) {
            std::cout << "  ⚠️  WARNING: Dec fuori range [-π/2, π/2]\n";
        }
        
        // Posizione osservatore
        Vector3D obs_pos = Coordinates::observerPosition(obs.observatoryCode, obs.epoch);
        printVector("  Pos osservatore", obs_pos);
        
        // La posizione osservatore dovrebbe essere ~1 AU dal Sole
        double r_obs = obs_pos.magnitude();
        if (r_obs < 0.98 || r_obs > 1.02) {
            std::cout << "  ⚠️  WARNING: Distanza osservatore dal Sole anomala: " << r_obs << " AU\n";
        }
        
        // Direzione verso oggetto (vettore unitario)
        Vector3D unit_vec = Coordinates::raDecToUnitVector(obs.obs.ra, obs.obs.dec);
        printVector("  Direzione (unit vec)", unit_vec);
        
        // Verifica che sia unitario
        double mag = unit_vec.magnitude();
        if (std::abs(mag - 1.0) > 1e-6) {
            std::cout << "  ⚠️  WARNING: Vettore non unitario: |v| = " << mag << "\n";
        }
        
        std::cout << "\n";
    }
    
    std::cout << "============================================================\n";
    std::cout << "3. VERIFICA COERENZA GEOMETRICA\n";
    std::cout << "============================================================\n\n";
    
    // Per Eros dovremmo avere:
    // - Distanza dal Sole: 1.13 - 1.78 AU (perielio - afelio)
    // - Inclinazione orbitale: ~10.8°
    
    std::cout << "Eros - parametri noti:\n";
    std::cout << "  a = 1.458 AU\n";
    std::cout << "  e = 0.223\n";
    std::cout << "  q = a(1-e) = 1.133 AU (perielio)\n";
    std::cout << "  Q = a(1+e) = 1.783 AU (afelio)\n";
    std::cout << "  i = 10.83°\n\n";
    
    // Stima distanze geocentriche assumendo Eros a ~1.5 AU
    auto& obs1 = observations[idx1];
    auto& obs2 = observations[idx2];
    auto& obs3 = observations[idx3];
    
    Vector3D obs1_pos = Coordinates::observerPosition(obs1.observatoryCode, obs1.epoch);
    Vector3D obs2_pos = Coordinates::observerPosition(obs2.observatoryCode, obs2.epoch);
    Vector3D obs3_pos = Coordinates::observerPosition(obs3.observatoryCode, obs3.epoch);
    
    Vector3D dir1 = Coordinates::raDecToUnitVector(obs1.obs.ra, obs1.obs.dec);
    Vector3D dir2 = Coordinates::raDecToUnitVector(obs2.obs.ra, obs2.obs.dec);
    Vector3D dir3 = Coordinates::raDecToUnitVector(obs3.obs.ra, obs3.obs.dec);
    
    // Assumiamo Eros a distanza heliocentriche ~1.5 AU
    // Distanza geocentrica può essere 0.5 - 2.5 AU a seconda della configurazione
    
    std::cout << "Angoli tra direzioni:\n";
    double angle12 = acos(dir1.dot(dir2)) * 180.0 / PI;
    double angle23 = acos(dir2.dot(dir3)) * 180.0 / PI;
    double angle13 = acos(dir1.dot(dir3)) * 180.0 / PI;
    
    std::cout << "  θ(obs1, obs2) = " << angle12 << "°\n";
    std::cout << "  θ(obs2, obs3) = " << angle23 << "°\n";
    std::cout << "  θ(obs1, obs3) = " << angle13 << "°\n\n";
    
    if (angle13 < 1.0) {
        std::cout << "  ⚠️  WARNING: Arco angolare molto piccolo (<1°)\n";
        std::cout << "     Difficile determinare orbita con precisione\n\n";
    }
    
    std::cout << "============================================================\n";
    std::cout << "4. TEST CONVERSIONE COORDINATE\n";
    std::cout << "============================================================\n\n";
    
    // Test: converti RA/Dec in vettore e riconverti
    double test_ra = 45.0 * DEG;
    double test_dec = 30.0 * DEG;
    
    std::cout << "Test conversione RA/Dec <-> Vector:\n";
    std::cout << "  Input:  RA = " << (test_ra * 180.0 / PI) << "°, Dec = " << (test_dec * 180.0 / PI) << "°\n";
    
    Vector3D test_vec = Coordinates::raDecToUnitVector(test_ra, test_dec);
    printVector("  Vector", test_vec);
    
    // Riconversione (se esiste la funzione)
    // Per ora verifichiamo manualmente
    double ra_back = atan2(test_vec.y, test_vec.x);
    if (ra_back < 0) ra_back += 2.0 * PI;
    double dec_back = asin(test_vec.z);
    
    std::cout << "  Output: RA = " << (ra_back * 180.0 / PI) << "°, Dec = " << (dec_back * 180.0 / PI) << "°\n";
    std::cout << "  Δ RA  = " << std::abs(ra_back - test_ra) * 180.0 / PI << "°\n";
    std::cout << "  Δ Dec = " << std::abs(dec_back - test_dec) * 180.0 / PI << "°\n";
    
    if (std::abs(ra_back - test_ra) > 1e-10 || std::abs(dec_back - test_dec) > 1e-10) {
        std::cout << "  ✗ ERRORE nella conversione!\n";
    } else {
        std::cout << "  ✓ Conversione OK\n";
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "5. VERIFICA PIANO ECLITTICA vs EQUATORIALE\n";
    std::cout << "============================================================\n\n";
    
    // Obliquità dell'eclittica ε ≈ 23.44°
    double epsilon = 23.44 * DEG;
    
    std::cout << "Sistema coordinate:\n";
    std::cout << "  - RA/Dec sono in coordinate EQUATORIALI (riferimento: equatore terrestre)\n";
    std::cout << "  - Orbite heliocentriche sono in coordinate ECLITTICHE (riferimento: piano orbitale Terra)\n";
    std::cout << "  - Obliquità eclittica ε = " << (epsilon * 180.0 / PI) << "°\n\n";
    
    std::cout << "Per Eros con i=10.83° rispetto all'eclittica:\n";
    std::cout << "  - Se osservazioni in coord. equatoriali e orbita in eclittiche\n";
    std::cout << "  - Serve rotazione di ε=23.44° per passare da equatoriali a eclittiche\n";
    std::cout << "  - Inclinazione apparente può sembrare diversa\n\n";
    
    // Verifica posizione osservatore
    std::cout << "Posizione osservatore (Terra):\n";
    std::cout << "  - Se in coord. eclittiche: z dovrebbe essere ~0 (Terra su piano eclittica)\n";
    std::cout << "  - Se in coord. equatoriali: z può essere != 0\n\n";
    
    std::cout << "Componente Z posizione osservatore:\n";
    std::cout << "  Obs1: z = " << obs1_pos.z << " AU\n";
    std::cout << "  Obs2: z = " << obs2_pos.z << " AU\n";
    std::cout << "  Obs3: z = " << obs3_pos.z << " AU\n\n";
    
    double avg_z = (std::abs(obs1_pos.z) + std::abs(obs2_pos.z) + std::abs(obs3_pos.z)) / 3.0;
    
    if (avg_z < 0.01) {
        std::cout << "  ✓ Posizioni osservatore nel piano eclittica (z~0)\n";
        std::cout << "    → Sistema ECLITTICO\n";
    } else {
        std::cout << "  ⚠️  Posizioni osservatore fuori dal piano eclittica\n";
        std::cout << "    → Possibile sistema EQUATORIALE o MISTO\n";
        std::cout << "    → Serve conversione coordinate!\n";
    }
    
    std::cout << "\n============================================================\n";
    std::cout << "6. DIAGNOSI PROBLEMA INCLINAZIONE\n";
    std::cout << "============================================================\n\n";
    
    std::cout << "Problema osservato:\n";
    std::cout << "  - Herget calcola i ≈ 63° invece di i ≈ 11°\n";
    std::cout << "  - Differenza Δi ≈ 52°\n\n";
    
    std::cout << "Possibili cause:\n";
    std::cout << "  1. Osservazioni in coord. equatoriali, orbita calcolata in equatoriali\n";
    std::cout << "     ma poi confrontata con riferimento eclittico\n";
    std::cout << "  2. Posizione osservatore in coord. equatoriali ma direzione RA/Dec\n";
    std::cout << "     usata come se fosse eclittica\n";
    std::cout << "  3. Missing rotazione equatoriale → eclittica dopo calcolo orbita\n\n";
    
    std::cout << "Nota: 63° - 23.44° = 39.56° (ancora lontano da 11°)\n";
    std::cout << "      Oppure: combinazione di rotazioni\n\n";
    
    std::cout << "============================================================\n";
    std::cout << "RACCOMANDAZIONI\n";
    std::cout << "============================================================\n\n";
    
    std::cout << "1. Verificare che observerPosition() ritorni coordinate ECLITTICHE J2000\n";
    std::cout << "2. Verificare che raDecToUnitVector() produca vettori in coord. ECLITTICHE\n";
    std::cout << "3. Se osservazioni MPC sono in equatoriali, serve conversione:\n";
    std::cout << "   - Ruotare direzioni RA/Dec da equatoriali a eclittiche\n";
    std::cout << "   - Rotazione di -ε attorno all'asse X\n";
    std::cout << "4. Dopo calcolo orbita, convertire elementi da eclittica a display\n\n";
    
    std::cout << "============================================================\n";
    std::cout << "TEST COMPLETATO\n";
    std::cout << "============================================================\n";
    
    return 0;
}
