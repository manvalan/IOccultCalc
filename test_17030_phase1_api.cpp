/**
 * @file test_17030_phase1_api.cpp
 * @brief Test Phase1CandidateScreening usando l'API screenCandidates()
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <algorithm>

#include "phase1_candidate_screening.h"
#include "ioc_gaialib/unified_gaia_catalog.h"

using namespace std;
using namespace ioccultcalc;

constexpr double MJD_28_NOV_2025 = 61007.0;

int main() {
    cout << "\n";
    cout << "=== TEST PHASE1 API - 16 PUNTI OTTIMIZZATI ===\n";
    cout << "Asteroide: (17030) 1999 CQ3\n";
    cout << "Data: 28 Novembre 2025\n\n";
    
    auto t_total_start = chrono::high_resolution_clock::now();
    
    // Init Gaia
    const char* home = getenv("HOME");
    if (!home) return 1;
    
    string json_config = R"({"catalog_type": "multifile_v2", "multifile_directory": ")" + 
                         string(home) + R"(/.catalog/gaia_mag18_v2_multifile"})";
    
    if (!ioc::gaia::UnifiedGaiaCatalog::initialize(json_config)) return 1;
    cout << "Gaia OK\n";
    
    // Load elements
    Phase1CandidateScreening screener;
    string eq1_path = "17030_astdys.eq1";
    { ifstream test(eq1_path); if (!test.good()) eq1_path = "../17030_astdys.eq1"; }
    
    if (!screener.loadAsteroidFromEQ1(eq1_path)) return 1;
    cout << "Elementi OK\n\n";
    
    // Config
    Phase1Config config;
    config.start_mjd_tdb = MJD_28_NOV_2025;
    config.end_mjd_tdb = MJD_28_NOV_2025 + 1.0;
    config.corridor_width_deg = 0.05;
    config.max_magnitude = 18.0;
    config.closest_approach_threshold_arcsec = 10.0;
    
    // Run
    cout << "Esecuzione screening...\n";
    auto t_start = chrono::high_resolution_clock::now();
    Phase1Results results = screener.screenCandidates(config);
    auto t_end = chrono::high_resolution_clock::now();
    
    double total_ms = chrono::duration<double, milli>(t_end - t_start).count();
    
    // Results
    cout << "\n=== TIMING ===\n";
    cout << "Path (16 punti):      " << fixed << setprecision(1) << results.propagation_time_ms << " ms\n";
    cout << "Query Gaia corridor:  " << results.corridor_query_time_ms << " ms\n";
    cout << "Closest approach:     " << results.closest_approach_calc_time_ms << " ms\n";
    cout << "TOTALE:               " << total_ms << " ms\n\n";
    
    cout << "=== STATISTICHE ===\n";
    cout << "Punti path:           " << results.num_path_points << "\n";
    cout << "Stelle corridor:      " << results.num_stars_in_corridor << "\n";
    cout << "Candidate (<10 as):   " << results.num_candidates_filtered << "\n\n";
    
    if (!results.candidates.empty()) {
        cout << "=== TOP CANDIDATE ===\n";
        auto candidates = results.candidates;
        sort(candidates.begin(), candidates.end(), 
             [](const CandidateStar& a, const CandidateStar& b) {
                 return a.closest_approach_arcsec < b.closest_approach_arcsec;
             });
        
        for (int i = 0; i < min(5, (int)candidates.size()); i++) {
            const auto& s = candidates[i];
            cout << i+1 << ". Source " << s.source_id 
                 << " Mag=" << setprecision(1) << s.phot_g_mean_mag
                 << " CA=" << setprecision(3) << s.closest_approach_arcsec << " arcsec\n";
        }
    }
    
    return 0;
}
