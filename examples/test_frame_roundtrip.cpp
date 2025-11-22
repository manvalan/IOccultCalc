// Test roundtrip conversione Equatoriale <-> Eclittica
// Per verificare che non ci siano errori sistematici nelle conversioni di frame

#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    // Obliquity J2000
    const double eps = 23.439291111 * M_PI / 180.0;
    const double cos_eps = std::cos(eps);
    const double sin_eps = std::sin(eps);
    
    // Stato OrbFit @ MJD 55918.807 (Equatoriale J2000)
    double r_eq_orig_x = -2.23090281851654E-01;
    double r_eq_orig_y =  1.12387595809716E+00;
    double r_eq_orig_z =  8.60668959807618E-02;
    
    double v_eq_orig_x = -1.68323419811707E-02;
    double v_eq_orig_y = -4.35994327067165E-03;
    double v_eq_orig_z = -3.12865629583693E-03;
    
    std::cout << std::setprecision(15);
    std::cout << "=== TEST ROUNDTRIP CONVERSIONE FRAME ===\n\n";
    
    std::cout << "1. STATO ORIGINALE (Equatoriale J2000):\n";
    std::cout << "   r = (" << r_eq_orig_x << ", " << r_eq_orig_y << ", " << r_eq_orig_z << ") AU\n";
    std::cout << "   v = (" << v_eq_orig_x << ", " << v_eq_orig_y << ", " << v_eq_orig_z << ") AU/d\n\n";
    
    // Conversione Equatoriale -> Eclittica
    double r_ecl_x = r_eq_orig_x;
    double r_ecl_y = r_eq_orig_y * cos_eps + r_eq_orig_z * sin_eps;
    double r_ecl_z = -r_eq_orig_y * sin_eps + r_eq_orig_z * cos_eps;
    
    double v_ecl_x = v_eq_orig_x;
    double v_ecl_y = v_eq_orig_y * cos_eps + v_eq_orig_z * sin_eps;
    double v_ecl_z = -v_eq_orig_y * sin_eps + v_eq_orig_z * cos_eps;
    
    std::cout << "2. CONVERSIONE -> Eclittica J2000:\n";
    std::cout << "   r = (" << r_ecl_x << ", " << r_ecl_y << ", " << r_ecl_z << ") AU\n";
    std::cout << "   v = (" << v_ecl_x << ", " << v_ecl_y << ", " << v_ecl_z << ") AU/d\n\n";
    
    // Conversione Eclittica -> Equatoriale (roundtrip)
    double r_eq_back_x = r_ecl_x;
    double r_eq_back_y = r_ecl_y * cos_eps - r_ecl_z * sin_eps;
    double r_eq_back_z = r_ecl_y * sin_eps + r_ecl_z * cos_eps;
    
    double v_eq_back_x = v_ecl_x;
    double v_eq_back_y = v_ecl_y * cos_eps - v_ecl_z * sin_eps;
    double v_eq_back_z = v_ecl_y * sin_eps + v_ecl_z * cos_eps;
    
    std::cout << "3. CONVERSIONE <- Equatoriale (ROUNDTRIP):\n";
    std::cout << "   r = (" << r_eq_back_x << ", " << r_eq_back_y << ", " << r_eq_back_z << ") AU\n";
    std::cout << "   v = (" << v_eq_back_x << ", " << v_eq_back_y << ", " << v_eq_back_z << ") AU/d\n\n";
    
    // Calcolo differenze
    double dr_x = (r_eq_back_x - r_eq_orig_x) * 1.496e8;  // AU -> km
    double dr_y = (r_eq_back_y - r_eq_orig_y) * 1.496e8;
    double dr_z = (r_eq_back_z - r_eq_orig_z) * 1.496e8;
    double dr_mag = std::sqrt(dr_x*dr_x + dr_y*dr_y + dr_z*dr_z);
    
    double dv_x = (v_eq_back_x - v_eq_orig_x) * 1.496e8 / 86400.0;  // AU/d -> m/s
    double dv_y = (v_eq_back_y - v_eq_orig_y) * 1.496e8 / 86400.0;
    double dv_z = (v_eq_back_z - v_eq_orig_z) * 1.496e8 / 86400.0;
    double dv_mag = std::sqrt(dv_x*dv_x + dv_y*dv_y + dv_z*dv_z);
    
    std::cout << "4. ERRORE ROUNDTRIP:\n";
    std::cout << "   Δr = (" << dr_x << ", " << dr_y << ", " << dr_z << ") km\n";
    std::cout << "   |Δr| = " << dr_mag << " km\n";
    std::cout << "   Δv = (" << dv_x << ", " << dv_y << ", " << dv_z << ") m/s\n";
    std::cout << "   |Δv| = " << dv_mag << " m/s\n\n";
    
    // Verifica precisione macchina
    double epsilon_pos = dr_mag / (std::sqrt(r_eq_orig_x*r_eq_orig_x + 
                                              r_eq_orig_y*r_eq_orig_y + 
                                              r_eq_orig_z*r_eq_orig_z) * 1.496e8);
    
    std::cout << "5. VALUTAZIONE:\n";
    std::cout << "   Errore relativo: " << epsilon_pos << "\n";
    
    if (dr_mag < 1e-6 && dv_mag < 1e-9) {
        std::cout << "   ✅ CONVERSIONE CORRETTA (errore < precisione macchina)\n";
    } else if (dr_mag < 1.0) {
        std::cout << "   ⚠️  Errore trascurabile ma misurabile\n";
    } else {
        std::cout << "   ❌ ERRORE SIGNIFICATIVO nella conversione!\n";
    }
    
    return 0;
}
