/**
 * @file test_frame_detection.cpp
 * @brief Testa quale frame danno risultati accurati: equatoriale o eclittico
 */

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const double M_PI = 3.14159265358979323846;
const double AU_TO_KM = 149597870.7;
const double DEG_TO_RAD = M_PI / 180.0;
const double RAD_TO_DEG = 180.0 / M_PI;
const double EPSILON_J2000 = 23.4392911 * DEG_TO_RAD;

struct CartesianState {
    double x, y, z;
    double distance() const {
        return std::sqrt(x*x + y*y + z*z);
    }
};

void getRaDec(double x, double y, double z, double& ra, double& dec) {
    ra = std::atan2(y, x);
    if (ra < 0) ra += 2 * M_PI;
    dec = std::atan2(z, std::sqrt(x*x + y*y));
}

CartesianState eclipticToEquatorial(double x, double y, double z) {
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    
    return {
        x,
        y * cos_eps - z * sin_eps,
        y * sin_eps + z * cos_eps
    };
}

CartesianState equatorialToEcliptic(double x, double y, double z) {
    double cos_eps = std::cos(EPSILON_J2000);
    double sin_eps = std::sin(EPSILON_J2000);
    
    return {
        x,
        y * cos_eps + z * sin_eps,
        -y * sin_eps + z * cos_eps
    };
}

int main() {
    const double JPL_X = 1.08361;
    const double JPL_Y = 3.08629;
    const double JPL_Z = -0.09162;
    
    std::cout << "\n=======================================================\n";
    std::cout << "Frame Detection Test - Asteroid 17030\n";
    std::cout << "=======================================================\n";
    
    std::cout << "\nJPL HORIZONS REFERENCE (2025-12-01):\n";
    std::cout << "  Position: (" << JPL_X << ", " << JPL_Y << ", " << JPL_Z << ") AU\n\n";
    
    // Risultato calcolato diretto (equatoriale, senza rotazione)
    double calc_x = 1.08272;
    double calc_y = 3.08657;
    double calc_z = -0.09158;
    
    std::cout << "CALCULATED (equatorial, no rotation):\n";
    std::cout << "  Position: (" << calc_x << ", " << calc_y << ", " << calc_z << ") AU\n";
    
    double err_au = std::sqrt(
        (calc_x - JPL_X)*(calc_x - JPL_X) +
        (calc_y - JPL_Y)*(calc_y - JPL_Y) +
        (calc_z - JPL_Z)*(calc_z - JPL_Z)
    );
    std::cout << "  Error: " << err_au * AU_TO_KM << " km\n";
    std::cout << "  Status: " << (err_au * AU_TO_KM < 200000 ? "EXCELLENT" : "BAD") << "\n\n";
    
    // Try rotation
    CartesianState rotated = eclipticToEquatorial(calc_x, calc_y, calc_z);
    std::cout << "AFTER ecliptic-to-equatorial rotation:\n";
    std::cout << "  Position: (" << rotated.x << ", " << rotated.y << ", " << rotated.z << ") AU\n";
    
    double err_rot = std::sqrt(
        (rotated.x - JPL_X)*(rotated.x - JPL_X) +
        (rotated.y - JPL_Y)*(rotated.y - JPL_Y) +
        (rotated.z - JPL_Z)*(rotated.z - JPL_Z)
    );
    std::cout << "  Error: " << err_rot * AU_TO_KM << " km\n";
    std::cout << "  Status: " << (err_rot * AU_TO_KM < 200000 ? "EXCELLENT" : "BAD") << "\n\n";
    
    std::cout << "CONCLUSION:\n";
    if (err_au < err_rot) {
        std::cout << "  Elements are in EQUATORIAL frame\n";
        std::cout << "  NO rotation needed\n";
    } else {
        std::cout << "  Elements are in ECLIPTIC frame\n";
        std::cout << "  Rotation to equatorial NEEDED\n";
    }
    std::cout << "\n";
    
    return 0;
}
