// Test per calcolare la differenza tra frame HELIOCENTRIC e BARYCENTRIC
// Questa potrebbe essere la causa dell'errore di 807 km!

#include "spk_reader.h"
#include "jd_time.h"
#include "vector3d.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    try {
        SPKReader spk;
        
        // MJD 55918.807 (Dec 23, 2011)
        double mjd_initial = 55918.807;
        double jd_initial = mjd_initial + 2400000.5;
        
        // MJD 60300.0 (Jan 01, 2024) - Epoca finale test
        double mjd_final = 60300.0;
        double jd_final = mjd_final + 2400000.5;
        
        std::cout << std::setprecision(15);
        std::cout << "=== DIFFERENZA HELIOCENTRIC vs BARYCENTRIC ===\n\n";
        
        // Posizione OrbFit fitted @ T=0 (Equatoriale J2000)
        double eros_eq_x = -0.223090281851654;
        double eros_eq_y =  1.12387595809716;
        double eros_eq_z =  0.0860668959807618;
        
        // Conversione a Eclittica
        const double eps = 23.439291111 * M_PI / 180.0;
        double eros_ecl_x = eros_eq_x;
        double eros_ecl_y = eros_eq_y * std::cos(eps) + eros_eq_z * std::sin(eps);
        double eros_ecl_z = -eros_eq_y * std::sin(eps) + eros_eq_z * std::cos(eps);
        
        std::cout << "1. POSIZIONE EROS @ MJD " << mjd_initial << " (Eclittica J2000):\n";
        std::cout << "   r = (" << eros_ecl_x << ", " << eros_ecl_y << ", " << eros_ecl_z << ") AU\n\n";
        
        // Posizione Sole rispetto al baricentro @ T=0
        Vector3D sun_bary_t0 = spk.getPosition(10, jd_initial, 0);  // 0 = SSB
        std::cout << "2. SOLE rispetto BARICENTRO @ MJD " << mjd_initial << ":\n";
        std::cout << "   Sun-SSB = (" << sun_bary_t0.x << ", " << sun_bary_t0.y << ", " << sun_bary_t0.z << ") AU\n";
        std::cout << "   |Sun-SSB| = " << sun_bary_t0.magnitude() << " AU = " 
                  << (sun_bary_t0.magnitude() * 1.496e8) << " km\n\n";
        
        // Se OrbFit fornisce coordinate BARYCENTRICHE, allora:
        // r_helio = r_bary - r_sun_bary
        Vector3D eros_bary_t0(eros_ecl_x, eros_ecl_y, eros_ecl_z);
        Vector3D eros_helio_t0 = eros_bary_t0 - sun_bary_t0;
        
        std::cout << "3. SE OrbFit usa BARYCENTRIC → Correzione HELIOCENTRIC:\n";
        std::cout << "   r_helio = r_bary - r_sun_bary\n";
        std::cout << "   r_helio = (" << eros_helio_t0.x << ", " << eros_helio_t0.y << ", " << eros_helio_t0.z << ") AU\n";
        std::cout << "   Differenza: " << (eros_bary_t0 - eros_helio_t0).magnitude() * 1.496e8 << " km\n\n";
        
        // Posizione Sole @ T=finale
        Vector3D sun_bary_tf = spk.getPosition(10, jd_final, 0);
        std::cout << "4. SOLE rispetto BARICENTRO @ MJD " << mjd_final << ":\n";
        std::cout << "   Sun-SSB = (" << sun_bary_tf.x << ", " << sun_bary_tf.y << ", " << sun_bary_tf.z << ") AU\n";
        std::cout << "   |Sun-SSB| = " << sun_bary_tf.magnitude() << " AU = " 
                  << (sun_bary_tf.magnitude() * 1.496e8) << " km\n\n";
        
        // Differenza tra offset baricentrico a T=0 e T=finale
        Vector3D delta_sun = sun_bary_tf - sun_bary_t0;
        std::cout << "5. VARIAZIONE offset baricentrico in 12 anni:\n";
        std::cout << "   Δ(Sun-SSB) = (" << delta_sun.x << ", " << delta_sun.y << ", " << delta_sun.z << ") AU\n";
        std::cout << "   |Δ(Sun-SSB)| = " << delta_sun.magnitude() << " AU = " 
                  << (delta_sun.magnitude() * 1.496e8) << " km\n\n";
        
        std::cout << "6. CONFRONTO CON ERRORE MISURATO:\n";
        std::cout << "   Errore test: 807,346 km\n";
        std::cout << "   Offset baricentrico @ T=0: " << (sun_bary_t0.magnitude() * 1.496e8) << " km\n";
        std::cout << "   Variazione offset in 12 anni: " << (delta_sun.magnitude() * 1.496e8) << " km\n\n";
        
        std::cout << "CONCLUSIONE:\n";
        if (sun_bary_t0.magnitude() * 1.496e8 > 500000.0) {
            std::cout << "   ✅ Offset baricentrico SIGNIFICATIVO!\n";
            std::cout << "   ⚠️  Se OrbFit usa frame BARYCENTRIC e IOccultCalc usa HELIOCENTRIC,\n";
            std::cout << "      questo spiega l'errore di ~800 km!\n";
        } else {
            std::cout << "   ❌ Offset baricentrico troppo piccolo per spiegare 807 km\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
