/**
 * @file test_perturbations_sierks.cpp
 * @brief Confronta precisione con/senza perturbazioni per evento (17030) Sierks
 * 
 * Evento di riferimento IOTA:
 * - Asteroide: (17030) Sierks
 * - Stella: UCAC4 552-011427 / Gaia DR3 3411546266140512128
 * - RA: 73.416106°, Dec: +20.331662°
 * - Data IOTA: 28 Nov 2025 00:35 UT
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "ioccultcalc/ephemeris.h"
#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/coordinates.h"

using namespace ioccultcalc;

// Costanti per la stella di riferimento
constexpr double STAR_RA_DEG = 73.416106;   // gradi
constexpr double STAR_DEC_DEG = 20.331662;  // gradi
constexpr double IOTA_TIME_JD = 2460645.524305556; // 28 Nov 2025 00:35 UT

void printTimeUTC(double jd) {
    // Converti JD in data/ora UTC
    double z = std::floor(jd + 0.5);
    double f = jd + 0.5 - z;
    
    double alpha = std::floor((z - 1867216.25) / 36524.25);
    double A = z + 1 + alpha - std::floor(alpha / 4.0);
    double B = A + 1524;
    double C = std::floor((B - 122.1) / 365.25);
    double D = std::floor(365.25 * C);
    double E = std::floor((B - D) / 30.6001);
    
    int day = static_cast<int>(B - D - std::floor(30.6001 * E));
    int month = static_cast<int>(E < 14 ? E - 1 : E - 13);
    int year = static_cast<int>(month > 2 ? C - 4716 : C - 4715);
    
    double hours_total = f * 24.0;
    int hours = static_cast<int>(hours_total);
    double mins_total = (hours_total - hours) * 60.0;
    int mins = static_cast<int>(mins_total);
    double secs = (mins_total - mins) * 60.0;
    
    std::cout << year << "-" << std::setw(2) << std::setfill('0') << month 
              << "-" << std::setw(2) << day << " "
              << std::setw(2) << hours << ":" 
              << std::setw(2) << mins << ":"
              << std::setw(5) << std::setfill('0') << std::fixed << std::setprecision(2) << secs
              << " UT" << std::setfill(' ');
}

double computeSeparation(double ra1, double dec1, double ra2, double dec2) {
    // Calcola separazione angolare in arcsec
    double ra1_rad = ra1 * M_PI / 180.0;
    double dec1_rad = dec1 * M_PI / 180.0;
    double ra2_rad = ra2 * M_PI / 180.0;
    double dec2_rad = dec2 * M_PI / 180.0;
    
    double cos_sep = std::sin(dec1_rad) * std::sin(dec2_rad) +
                     std::cos(dec1_rad) * std::cos(dec2_rad) * std::cos(ra1_rad - ra2_rad);
    
    if (cos_sep > 1.0) cos_sep = 1.0;
    if (cos_sep < -1.0) cos_sep = -1.0;
    
    double sep_rad = std::acos(cos_sep);
    return sep_rad * 180.0 / M_PI * 3600.0; // arcsec
}

struct SearchResult {
    double jd_min;
    double sep_min_arcsec;
    bool found;
};

SearchResult findMinimumSeparation(Ephemeris& eph, 
                                    double jd_start, double jd_end, 
                                    double step_days,
                                    double star_ra, double star_dec) {
    SearchResult result;
    result.found = false;
    result.sep_min_arcsec = 1e9;
    result.jd_min = 0;
    
    double jd_min_coarse = jd_start;
    double sep_min_coarse = 1e9;
    
    // Ricerca grossolana
    for (double jd = jd_start; jd <= jd_end; jd += step_days) {
        EphemerisData data = eph.compute(JulianDate(jd));
        double sep = computeSeparation(data.geocentricPos.ra, data.geocentricPos.dec,
                                       star_ra, star_dec);
        if (sep < sep_min_coarse) {
            sep_min_coarse = sep;
            jd_min_coarse = jd;
        }
    }
    
    // Ricerca fine attorno al minimo
    double jd_fine_start = jd_min_coarse - step_days * 2;
    double jd_fine_end = jd_min_coarse + step_days * 2;
    double fine_step = step_days / 100.0;
    
    for (double jd = jd_fine_start; jd <= jd_fine_end; jd += fine_step) {
        EphemerisData data = eph.compute(JulianDate(jd));
        double sep = computeSeparation(data.geocentricPos.ra, data.geocentricPos.dec,
                                       star_ra, star_dec);
        if (sep < result.sep_min_arcsec) {
            result.sep_min_arcsec = sep;
            result.jd_min = jd;
            result.found = true;
        }
    }
    
    return result;
}

int main() {
    std::cout << "======================================================\n";
    std::cout << "  TEST PERTURBAZIONI: (17030) Sierks\n";
    std::cout << "======================================================\n\n";
    
    std::cout << "Evento di riferimento IOTA:\n";
    std::cout << "  - Data prevista: 28 Nov 2025 00:35 UT\n";
    std::cout << "  - Stella: RA " << STAR_RA_DEG << "°, Dec +" << STAR_DEC_DEG << "°\n\n";
    
    // Scarica elementi orbitali
    std::cout << "Scaricamento elementi orbitali da AstDyS...\n";
    AstDysClient client;
    EquinoctialElements elements;
    
    try {
        elements = client.getElements("17030");
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: impossibile scaricare elementi per (17030) Sierks: " << e.what() << "\n";
        return 1;
    }
    
    std::cout << "  a = " << elements.a << " AU\n";
    std::cout << "  h = " << elements.h << " (e·sin(ω))\n";
    std::cout << "  k = " << elements.k << " (e·cos(ω))\n";
    std::cout << "  epoca = JD " << elements.epoch.jd << "\n\n";
    
    // Finestra di ricerca: ±2 ore da IOTA prediction
    double jd_start = IOTA_TIME_JD - 2.0/24.0;  // 2 ore prima
    double jd_end = IOTA_TIME_JD + 2.0/24.0;    // 2 ore dopo
    
    // =========================================
    // TEST 1: Senza perturbazioni (Keplero puro)
    // =========================================
    std::cout << "--- TEST 1: Senza perturbazioni (Keplero puro) ---\n";
    Ephemeris ephKeplero(elements);
    ephKeplero.setOptions(EphemerisOptions::fast());  // Keplero puro
    
    SearchResult resKeplero = findMinimumSeparation(ephKeplero, jd_start, jd_end, 
                                                     1.0/1440.0, // 1 minuto
                                                     STAR_RA_DEG, STAR_DEC_DEG);
    
    std::cout << "  Separazione minima: " << std::fixed << std::setprecision(2) 
              << resKeplero.sep_min_arcsec << " arcsec\n";
    std::cout << "  Istante minimo:     ";
    printTimeUTC(resKeplero.jd_min);
    std::cout << "\n";
    
    double delta_t_keplero = (resKeplero.jd_min - IOTA_TIME_JD) * 24.0 * 60.0; // minuti
    std::cout << "  Differenza da IOTA: " << std::showpos << std::fixed << std::setprecision(1)
              << delta_t_keplero << std::noshowpos << " minuti\n\n";
    
    // =========================================
    // TEST 2: Con perturbazioni planetarie
    // =========================================
    std::cout << "--- TEST 2: Con perturbazioni planetarie ---\n";
    Ephemeris ephPerturbato(elements);
    ephPerturbato.setOptions(EphemerisOptions::standard());  // Con perturbazioni
    
    SearchResult resPerturbato = findMinimumSeparation(ephPerturbato, jd_start, jd_end,
                                                        1.0/1440.0, // 1 minuto
                                                        STAR_RA_DEG, STAR_DEC_DEG);
    
    std::cout << "  Separazione minima: " << std::fixed << std::setprecision(2)
              << resPerturbato.sep_min_arcsec << " arcsec\n";
    std::cout << "  Istante minimo:     ";
    printTimeUTC(resPerturbato.jd_min);
    std::cout << "\n";
    
    double delta_t_perturbato = (resPerturbato.jd_min - IOTA_TIME_JD) * 24.0 * 60.0;
    std::cout << "  Differenza da IOTA: " << std::showpos << std::fixed << std::setprecision(1)
              << delta_t_perturbato << std::noshowpos << " minuti\n\n";
    
    // =========================================
    // TEST 3: Alta precisione (perturbazioni + relatività)
    // =========================================
    std::cout << "--- TEST 3: Alta precisione (perturbazioni + relatività) ---\n";
    Ephemeris ephHP(elements);
    ephHP.setOptions(EphemerisOptions::highPrecision());  // Tutto abilitato
    
    SearchResult resHP = findMinimumSeparation(ephHP, jd_start, jd_end,
                                                1.0/1440.0,
                                                STAR_RA_DEG, STAR_DEC_DEG);
    
    std::cout << "  Separazione minima: " << std::fixed << std::setprecision(2)
              << resHP.sep_min_arcsec << " arcsec\n";
    std::cout << "  Istante minimo:     ";
    printTimeUTC(resHP.jd_min);
    std::cout << "\n";
    
    double delta_t_hp = (resHP.jd_min - IOTA_TIME_JD) * 24.0 * 60.0;
    std::cout << "  Differenza da IOTA: " << std::showpos << std::fixed << std::setprecision(1)
              << delta_t_hp << std::noshowpos << " minuti\n\n";
    
    // =========================================
    // RIEPILOGO
    // =========================================
    std::cout << "======================================================\n";
    std::cout << "  RIEPILOGO\n";
    std::cout << "======================================================\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "                     Separazione    ΔT da IOTA\n";
    std::cout << "  Keplero:           " << std::setw(8) << resKeplero.sep_min_arcsec 
              << " \"     " << std::showpos << std::setw(7) << std::setprecision(1) 
              << delta_t_keplero << std::noshowpos << " min\n";
    std::cout << "  Perturbato:        " << std::setw(8) << std::setprecision(2) << resPerturbato.sep_min_arcsec
              << " \"     " << std::showpos << std::setw(7) << std::setprecision(1) 
              << delta_t_perturbato << std::noshowpos << " min\n";
    std::cout << "  Alta precisione:   " << std::setw(8) << std::setprecision(2) << resHP.sep_min_arcsec
              << " \"     " << std::showpos << std::setw(7) << std::setprecision(1)
              << delta_t_hp << std::noshowpos << " min\n\n";
    
    std::cout << "Miglioramento separazione: "
              << std::fixed << std::setprecision(1)
              << (resKeplero.sep_min_arcsec - resHP.sep_min_arcsec) << " arcsec\n";
    std::cout << "Miglioramento tempo:       "
              << std::abs(delta_t_keplero) - std::abs(delta_t_hp) << " minuti\n";
    
    return 0;
}
