/**
 * @file iers_data.h
 * @brief Earth Orientation Parameters from IERS
 * 
 * Implements EOP corrections for high-precision astrometry:
 * - UT1-UTC (Earth rotation)
 * - Polar motion (x, y coordinates)
 * - Nutation corrections (dψ, dε)
 * 
 * Data sources:
 * - IERS Bulletin A (predictions): https://datacenter.iers.org/data/9/finals2000A.all
 * - IERS Bulletin B (historical): https://datacenter.iers.org/data/224/eopc04.62-now
 */

#ifndef IOCCULTCALC_IERS_DATA_H
#define IOCCULTCALC_IERS_DATA_H

#include "ioccultcalc/time_utils.h"
#include <string>
#include <vector>
#include <map>

namespace ioccultcalc {

/**
 * @brief Earth Orientation Parameters for a specific date
 */
struct EarthOrientationParams {
    double mjd;          ///< Modified Julian Date
    double ut1_utc;      ///< UT1-UTC in seconds
    double x_pole;       ///< Polar motion X in arcseconds
    double y_pole;       ///< Polar motion Y in arcseconds
    double dPsi;         ///< Nutation correction in longitude (arcsec)
    double dEps;         ///< Nutation correction in obliquity (arcsec)
    
    // Uncertainties (if available)
    double sigma_ut1;    ///< Uncertainty in UT1-UTC (seconds)
    double sigma_pole;   ///< Uncertainty in polar motion (arcsec)
};

/**
 * @brief Manager for IERS Earth Orientation Parameters
 */
class IERSData {
public:
    /**
     * @brief Initialize IERS data manager
     * @param cache_dir Directory for caching IERS files
     */
    explicit IERSData(const std::string& cache_dir = "");
    
    /**
     * @brief Get EOP for a specific Julian Date
     * @param jd Julian Date (UTC)
     * @return EOP structure with interpolated values
     * 
     * If data not available, downloads from IERS.
     * Uses linear interpolation between tabulated values.
     */
    EarthOrientationParams getEOP(double jd);
    
    /**
     * @brief Download latest IERS data
     * @param force Force re-download even if cache exists
     * @return true if successful
     * 
     * Downloads finals2000A.all (Bulletin A) covering 1992-present.
     */
    bool downloadLatestData(bool force = false);
    
    /**
     * @brief Load EOP from cached file
     * @param filename Path to IERS data file
     * @return Number of records loaded
     */
    int loadFromFile(const std::string& filename);
    
    /**
     * @brief Check if EOP data available for date
     * @param jd Julian Date
     * @return true if data available (no download needed)
     */
    bool hasDataFor(double jd) const;
    
    /**
     * @brief Get date range covered by loaded data
     * @param mjd_start Output: start MJD
     * @param mjd_end Output: end MJD
     * @return true if data loaded
     */
    bool getDataRange(double& mjd_start, double& mjd_end) const;

private:
    std::string cache_dir_;
    std::map<int, EarthOrientationParams> eop_data_;  // Keyed by MJD (integer part)
    
    // Linear interpolation between two EOPs
    EarthOrientationParams interpolate(
        const EarthOrientationParams& eop1,
        const EarthOrientationParams& eop2,
        double mjd) const;
    
    // Parse IERS finals2000A format
    bool parseFinals2000A(const std::string& line, EarthOrientationParams& eop);
    
    // Default cache directory
    std::string getDefaultCacheDir() const;
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_IERS_DATA_H
