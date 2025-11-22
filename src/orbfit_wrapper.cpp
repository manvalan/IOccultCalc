/**
 * @file orbfit_wrapper.cpp
 * @brief Implementazione wrapper C++ per OrbFit
 */

#include "ioccultcalc/orbfit_wrapper.h"
#include "ioccultcalc/coordinates.h"
#include <cmath>
#include <iostream>
#include <cstring>

namespace ioccultcalc {

// Static members
bool OrbFitWrapper::s_initialized = false;
std::string OrbFitWrapper::s_lib_path;

// ============================================================================
// C INTERFACE DECLARATIONS (from Fortran wrapper with bind(C))
// ============================================================================

extern "C" {
    
    /**
     * @brief Gauss method from OrbFit (via Fortran wrapper)
     * bind(C, name="orbfit_gauss_c")
     */
    void orbfit_gauss_c(
        double* tobs,       // [3] Time of observations (MJD, TDT)
        double* alpha,      // [3] Right ascension (rad)
        double* delta,      // [3] Declination (rad)
        int* obscod,        // [3] Observatory codes
        double* elem,       // [6,3] Orbital elements output (column-major!)
        char* eletyp,       // [4,3] Element types (column-major!)
        double* t0,         // [3] Epoch of elements
        int* nroots,        // Number of polynomial roots
        int* nsol,          // Number of valid solutions
        double* rtop,       // [3] Topocentric distances
        int* fail,          // Failure flag
        char* errmsg,       // [200] Error message
        int* debug          // Debug flag
    );

    /**
     * @brief Vaisala method from OrbFit (via Fortran wrapper)
     * bind(C, name="orbfit_vaisala_c")
     */
    void orbfit_vaisala_c(
        double* tobs,       // [2] Time of observations (MJD, TDT)
        double* alpha,      // [2] Right ascension (rad)
        double* delta,      // [2] Declination (rad)
        int* obscod,        // [2] Observatory codes
        double* solar_r,    // Initial heliocentric distance guess
        double* elem,       // [6] Orbital elements output
        char* eletyp,       // [4] Element type
        double* t0,         // Epoch of elements
        double* rtop,       // Topocentric distance
        int* fail,          // Failure flag
        char* errmsg,       // [200] Error message
        int* debug          // Debug flag
    );

    /**
     * @brief OrbFit initialization (loads data files)
     * bind(C, name="orbfit_init_c")
     */
    void orbfit_init_c(
        const char* lib_path,  // Path to OrbFit lib directory
        int lib_path_len       // Length of path string
    );

    /**
     * @brief RA15 propagator initialization
     * bind(C, name="orbfit_ra15_init_c")
     */
    void orbfit_ra15_init_c();

    /**
     * @brief RA15 propagator - DISABLED
     * NOTE: propag() requires orbit_elem type which is too complex to wrap
     * bind(C, name="orbfit_ra15_propag_c")
     */
    // void orbfit_ra15_propag_c(...); // DISABLED
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

namespace {

/**
 * @brief Converte Modified Julian Date (MJD) in Julian Date (JD)
 */
inline double mjd_to_jd(double mjd) {
    return mjd + 2400000.5;
}

/**
 * @brief Converte Julian Date (JD) in Modified Julian Date (MJD)
 */
inline double jd_to_mjd(double jd) {
    return jd - 2400000.5;
}

/**
 * @brief Converte elementi orbitali da formato OrbFit a OrbitalElements
 * 
 * OrbFit usa:
 * - KEP: a, e, i, Omega, omega, M (mean ecliptic J2000)
 * - COM: q, e, i, Omega, omega, T (cometary)
 * - EQU: equinoctial elements
 * - CAR: cartesian state
 */
OrbitalElements convertFromOrbFit(const double* elem, const char* eltype, double epoch_mjd) {
    OrbitalElements result;
    result.epoch = JulianDate(mjd_to_jd(epoch_mjd));
    
    // OrbFit può ritornare diversi tipi di elementi
    std::string type(eltype, 3);
    
    if (type == "KEP" || type == "kep") {
        // Keplerian: a, e, i, Omega, omega, M
        result.a = elem[0];
        result.e = elem[1];
        result.i = elem[2] * 180.0 / M_PI;  // rad → deg
        result.Omega = elem[3] * 180.0 / M_PI;
        result.omega = elem[4] * 180.0 / M_PI;
        result.M = elem[5] * 180.0 / M_PI;
    }
    else if (type == "COM" || type == "com") {
        // Cometary: q, e, i, Omega, omega, T_peri
        double q = elem[0];
        result.e = elem[1];
        result.a = q / (1.0 - result.e);  // a = q / (1-e)
        result.i = elem[2] * 180.0 / M_PI;
        result.Omega = elem[3] * 180.0 / M_PI;
        result.omega = elem[4] * 180.0 / M_PI;
        double T_peri = elem[5];  // MJD
        // Calcola M dalla T_peri (semplificato)
        double n = std::sqrt(0.01720209895 * 0.01720209895 / (result.a * result.a * result.a));
        result.M = n * (epoch_mjd - T_peri) * 180.0 / M_PI;
        while (result.M < 0) result.M += 360.0;
        while (result.M >= 360.0) result.M -= 360.0;
    }
    else {
        std::cerr << "Warning: Unknown OrbFit element type '" << type << "', assuming KEP" << std::endl;
        result.a = elem[0];
        result.e = elem[1];
        result.i = elem[2] * 180.0 / M_PI;
        result.Omega = elem[3] * 180.0 / M_PI;
        result.omega = elem[4] * 180.0 / M_PI;
        result.M = elem[5] * 180.0 / M_PI;
    }
    
    return result;
}

/**
 * @brief Converte OrbitalElements in formato OrbFit
 */
void convertToOrbFit(const OrbitalElements& elements, double* elem, double& epoch_mjd) {
    // OrbFit KEP format: a, e, i, Omega, omega, M (radianti, eclittica J2000)
    elem[0] = elements.a;
    elem[1] = elements.e;
    elem[2] = elements.i * M_PI / 180.0;
    elem[3] = elements.Omega * M_PI / 180.0;
    elem[4] = elements.omega * M_PI / 180.0;
    elem[5] = elements.M * M_PI / 180.0;
    
    epoch_mjd = jd_to_mjd(elements.epoch.jd);
}

} // anonymous namespace

// ============================================================================
// PUBLIC METHODS
// ============================================================================

bool OrbFitWrapper::initialize(const std::string& orbfit_lib_path) {
    if (s_initialized) {
        return true;  // Already initialized
    }
    
    s_lib_path = orbfit_lib_path;
    
    // Set environment variable for OrbFit to find data files
    setenv("ORBFIT_DATA", orbfit_lib_path.c_str(), 1);
    
    // Call Fortran wrapper initialization
    orbfit_init_c(orbfit_lib_path.c_str(), orbfit_lib_path.length());
    
    s_initialized = true;
    std::cout << "OrbFit wrapper initialized (lib path: " << orbfit_lib_path << ")" << std::endl;
    
    return true;
}

void OrbFitWrapper::cleanup() {
    s_initialized = false;
}

std::vector<OrbFitWrapper::Solution> OrbFitWrapper::gauss(
    const AstrometricObservation& obs1,
    const AstrometricObservation& obs2,
    const AstrometricObservation& obs3,
    bool debug)
{
    std::vector<Solution> solutions;
    
    // Prepare input arrays
    double tobs[3] = {
        jd_to_mjd(obs1.epoch.jd),
        jd_to_mjd(obs2.epoch.jd),
        jd_to_mjd(obs3.epoch.jd)
    };
    
    double alpha[3] = {
        obs1.obs.ra * M_PI / 180.0,  // deg → rad
        obs2.obs.ra * M_PI / 180.0,
        obs3.obs.ra * M_PI / 180.0
    };
    
    double delta[3] = {
        obs1.obs.dec * M_PI / 180.0,
        obs2.obs.dec * M_PI / 180.0,
        obs3.obs.dec * M_PI / 180.0
    };
    
    // Observatory codes - OrbFit usa codici numerici
    // Per ora usiamo 500 (geocentrico) come default
    // TODO: mappare observatory codes MPC → OrbFit numeric
    int obscod[3] = {500, 500, 500};
    
    // Output arrays
    double elem[18];  // [6,3] - 3 soluzioni possibili (column-major!)
    char eletyp[12];  // [4,3] - tipo elementi per ogni soluzione (column-major!)
    double t0[3];     // Epoch per ogni soluzione
    double rtop[3];   // Distanze topocentric
    int nroots, nsol;
    int fail = 0;
    char errmsg[200];
    int debug_flag = debug ? 1 : 0;
    
    // Initialize arrays
    std::memset(elem, 0, sizeof(elem));
    std::memset(eletyp, 0, sizeof(eletyp));
    std::memset(errmsg, 0, sizeof(errmsg));
    
    // Call Fortran wrapper (C interface)
    orbfit_gauss_c(tobs, alpha, delta, obscod, 
                   elem, eletyp, t0, &nroots, &nsol, rtop,
                   &fail, errmsg, &debug_flag);
    
    if (fail != 0) {
        std::string err_str(errmsg);
        std::cerr << "OrbFit gaussn failed: " << err_str << std::endl;
        return solutions;
    }
    
    // Convert solutions
    // NOTA: Fortran usa column-major, quindi elem è [6,3] ma in memoria:
    // elem[0..5] = soluzione 1, elem[6..11] = soluzione 2, elem[12..17] = soluzione 3
    for (int i = 0; i < nsol; i++) {
        Solution sol;
        
        // Extract element type (4 chars) - column-major: eletyp[0..3]=sol1, [4..7]=sol2, [8..11]=sol3
        char eltype[5];
        std::strncpy(eltype, eletyp + i*4, 4);
        eltype[4] = '\0';
        
        // Convert elements (6 per solution)
        sol.elements = convertFromOrbFit(&elem[i*6], eltype, t0[i]);
        sol.topocentric_distance = rtop[i];
        sol.root_index = i;
        sol.is_valid = true;
        sol.rms_residual = 0.0;  // TODO: compute from residuals
        
        solutions.push_back(sol);
    }
    
    return solutions;
}

OrbFitWrapper::Solution OrbFitWrapper::vaisala(
    const AstrometricObservation& obs1,
    const AstrometricObservation& obs2,
    double heliocentric_distance,
    bool debug)
{
    Solution solution;
    solution.is_valid = false;
    
    // Prepare input arrays
    double tobs[2] = {
        jd_to_mjd(obs1.epoch.jd),
        jd_to_mjd(obs2.epoch.jd)
    };
    
    double alpha[2] = {
        obs1.obs.ra * M_PI / 180.0,
        obs2.obs.ra * M_PI / 180.0
    };
    
    double delta[2] = {
        obs1.obs.dec * M_PI / 180.0,
        obs2.obs.dec * M_PI / 180.0
    };
    
    int obscod[2] = {500, 500};  // Geocentric
    
    // Output arrays
    double elem[6];
    char eletyp[4];
    double t0;
    double rtop;
    int fail = 0;
    char errmsg[200];
    int debug_flag = debug ? 1 : 0;
    
    double solar_r = heliocentric_distance;
    
    // Initialize
    std::memset(elem, 0, sizeof(elem));
    std::memset(eletyp, 0, sizeof(eletyp));
    std::memset(errmsg, 0, sizeof(errmsg));
    
    // Call Fortran wrapper (C interface)
    orbfit_vaisala_c(tobs, alpha, delta, obscod, &solar_r,
                     elem, eletyp, &t0, &rtop, &fail, errmsg, &debug_flag);
    
    if (fail != 0) {
        std::string err_str(errmsg);
        std::cerr << "OrbFit vaisala failed: " << err_str << std::endl;
        return solution;
    }
    
    // Convert solution
    char eltype_cstr[5];
    std::strncpy(eltype_cstr, eletyp, 4);
    eltype_cstr[4] = '\0';
    solution.elements = convertFromOrbFit(elem, eltype_cstr, t0);
    solution.topocentric_distance = rtop;
    solution.root_index = 0;
    solution.is_valid = true;
    solution.rms_residual = 0.0;
    
    return solution;
}

// DISABLED: propagate() requires orbit_elem type which is too complex to wrap
// For orbit propagation, use our native RA15Integrator or TwoBodyPropagator
// OrbFit wrapper only provides preliminary orbit determination (Gauss, Vaisala)
/*
bool OrbFitWrapper::propagate(
    const OrbitalElements& elements,
    double epoch_start,
    double epoch_end,
    double state_out[6],
    bool use_nongrav)
{
    // DISABLED - use native propagators instead
    std::cerr << "OrbFit propagate() is disabled. Use native RA15Integrator or TwoBodyPropagator." << std::endl;
    return false;
}
*/

// fitObservations() - NOT IMPLEMENTED
// Requires complex file I/O and configuration
/*
bool OrbFitWrapper::fitObservations(...)
{
    std::cerr << "OrbFit fitObservations not yet implemented" << std::endl;
    return false;
}
*/

} // namespace ioccultcalc
