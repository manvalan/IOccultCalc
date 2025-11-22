/**
 * @file iers_data.cpp
 * @brief Implementation of IERS Earth Orientation Parameters manager
 */

#include "ioccultcalc/iers_data.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sys/stat.h>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace ioccultcalc {

IERSData::IERSData(const std::string& cache_dir)
    : cache_dir_(cache_dir.empty() ? getDefaultCacheDir() : cache_dir)
{
    // Try to load existing data
    std::string finals_file = cache_dir_ + "/finals2000A.all";
    struct stat buffer;
    if (stat(finals_file.c_str(), &buffer) == 0) {
        loadFromFile(finals_file);
    }
}

std::string IERSData::getDefaultCacheDir() const {
    // Use ~/.ioccultcalc/iers/
    const char* home = getenv("HOME");
    if (!home) home = "/tmp";
    
    std::string dir = std::string(home) + "/.ioccultcalc/iers";
    
    // Create directory if doesn't exist
#ifdef _WIN32
    _mkdir(dir.c_str());
#else
    mkdir(dir.c_str(), 0755);
#endif
    
    return dir;
}

EarthOrientationParams IERSData::getEOP(double jd) {
    // Convert JD to MJD
    double mjd = jd - 2400000.5;
    int mjd_int = static_cast<int>(std::floor(mjd));
    
    // Check if we have data
    if (eop_data_.empty() || !hasDataFor(jd)) {
        // Try to download
        if (!downloadLatestData(false)) {
            // Return zero EOP if download fails
            EarthOrientationParams zero = {mjd, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            return zero;
        }
    }
    
    // Find surrounding data points for interpolation
    auto it_after = eop_data_.lower_bound(mjd_int);
    
    if (it_after == eop_data_.end()) {
        // Beyond data range - use last point
        return eop_data_.rbegin()->second;
    }
    
    if (it_after == eop_data_.begin()) {
        // Before data range - use first point
        return it_after->second;
    }
    
    // Interpolate between two points
    auto it_before = std::prev(it_after);
    return interpolate(it_before->second, it_after->second, mjd);
}

EarthOrientationParams IERSData::interpolate(
    const EarthOrientationParams& eop1,
    const EarthOrientationParams& eop2,
    double mjd) const
{
    // Linear interpolation
    double f = (mjd - eop1.mjd) / (eop2.mjd - eop1.mjd);
    
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    
    EarthOrientationParams result;
    result.mjd = mjd;
    result.ut1_utc = eop1.ut1_utc + f * (eop2.ut1_utc - eop1.ut1_utc);
    result.x_pole = eop1.x_pole + f * (eop2.x_pole - eop1.x_pole);
    result.y_pole = eop1.y_pole + f * (eop2.y_pole - eop1.y_pole);
    result.dPsi = eop1.dPsi + f * (eop2.dPsi - eop1.dPsi);
    result.dEps = eop1.dEps + f * (eop2.dEps - eop1.dEps);
    result.sigma_ut1 = std::max(eop1.sigma_ut1, eop2.sigma_ut1);
    result.sigma_pole = std::max(eop1.sigma_pole, eop2.sigma_pole);
    
    return result;
}

bool IERSData::downloadLatestData(bool force) {
    std::string finals_file = cache_dir_ + "/finals2000A.all";
    
    // Check if file exists and is recent (< 7 days old)
    struct stat buffer;
    if (!force && stat(finals_file.c_str(), &buffer) == 0) {
        time_t now = time(nullptr);
        double age_days = difftime(now, buffer.st_mtime) / 86400.0;
        if (age_days < 7.0) {
            // File is recent enough
            return loadFromFile(finals_file) > 0;
        }
    }
    
    // Download using curl
    std::string url = "https://datacenter.iers.org/data/9/finals2000A.all";
    std::string cmd = "curl -s -o \"" + finals_file + "\" \"" + url + "\"";
    
    int result = system(cmd.c_str());
    if (result != 0) {
        return false;
    }
    
    // Load the downloaded file
    return loadFromFile(finals_file) > 0;
}

int IERSData::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return 0;
    }
    
    eop_data_.clear();
    std::string line;
    int count = 0;
    
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        EarthOrientationParams eop;
        if (parseFinals2000A(line, eop)) {
            int mjd_key = static_cast<int>(std::floor(eop.mjd));
            eop_data_[mjd_key] = eop;
            count++;
        }
    }
    
    return count;
}

bool IERSData::parseFinals2000A(const std::string& line, EarthOrientationParams& eop) {
    // Format: finals2000A.all (fixed-width columns)
    // Col 1-2:   Year (last 2 digits)
    // Col 3-4:   Month
    // Col 5-6:   Day
    // Col 8-15:  MJD
    // Col 19-27: x_pole (arcsec)
    // Col 38-46: y_pole (arcsec)
    // Col 59-68: UT1-UTC (seconds)
    // Col 98-106: dPsi (milliarcsec)
    // Col 117-125: dEps (milliarcsec)
    
    if (line.length() < 125) return false;
    
    try {
        // MJD
        std::string mjd_str = line.substr(7, 8);
        eop.mjd = std::stod(mjd_str);
        
        // x_pole (arcsec)
        std::string x_str = line.substr(18, 9);
        eop.x_pole = std::stod(x_str);
        
        // y_pole (arcsec)
        std::string y_str = line.substr(37, 9);
        eop.y_pole = std::stod(y_str);
        
        // UT1-UTC (seconds)
        std::string ut1_str = line.substr(58, 10);
        eop.ut1_utc = std::stod(ut1_str);
        
        // dPsi (convert milliarcsec to arcsec)
        if (line.length() >= 106) {
            std::string dpsi_str = line.substr(97, 9);
            eop.dPsi = std::stod(dpsi_str) / 1000.0;
        } else {
            eop.dPsi = 0.0;
        }
        
        // dEps (convert milliarcsec to arcsec)
        if (line.length() >= 125) {
            std::string deps_str = line.substr(116, 9);
            eop.dEps = std::stod(deps_str) / 1000.0;
        } else {
            eop.dEps = 0.0;
        }
        
        // Default uncertainties
        eop.sigma_ut1 = 0.0001;  // 0.1 ms
        eop.sigma_pole = 0.001;  // 1 mas
        
        return true;
    } catch (...) {
        return false;
    }
}

bool IERSData::hasDataFor(double jd) const {
    if (eop_data_.empty()) return false;
    
    double mjd = jd - 2400000.5;
    int mjd_int = static_cast<int>(std::floor(mjd));
    
    int mjd_min = eop_data_.begin()->first;
    int mjd_max = eop_data_.rbegin()->first;
    
    return (mjd_int >= mjd_min && mjd_int <= mjd_max);
}

bool IERSData::getDataRange(double& mjd_start, double& mjd_end) const {
    if (eop_data_.empty()) return false;
    
    mjd_start = eop_data_.begin()->second.mjd;
    mjd_end = eop_data_.rbegin()->second.mjd;
    
    return true;
}

} // namespace ioccultcalc
