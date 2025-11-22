#ifndef IOCCULTCALC_TYPES_H
#define IOCCULTCALC_TYPES_H

#include <string>
#include <cmath>

namespace ioccultcalc {

// Costanti astronomiche
constexpr double AU = 149597870.7;        // km
constexpr double C_LIGHT = 299792.458;    // km/s
constexpr double DEG_TO_RAD = M_PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / M_PI;
constexpr double EARTH_RADIUS = 6378.137; // km

// Struttura per coordinate equatoriali
struct EquatorialCoordinates {
    double ra;   // Right Ascension in radianti
    double dec;  // Declination in radianti
    double distance; // Distanza in km (opzionale)
    
    EquatorialCoordinates() : ra(0), dec(0), distance(0) {}
    EquatorialCoordinates(double r, double d, double dist = 0) 
        : ra(r), dec(d), distance(dist) {}
};

// Struttura per coordinate geografiche
struct GeographicCoordinates {
    double longitude; // in gradi
    double latitude;  // in gradi
    double altitude;  // in metri
    
    GeographicCoordinates() : longitude(0), latitude(0), altitude(0) {}
    GeographicCoordinates(double lon, double lat, double alt = 0)
        : longitude(lon), latitude(lat), altitude(alt) {}
};

// Struttura per vettori 3D
struct Vector3D {
    double x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    
    double magnitude() const {
        return std::sqrt(x*x + y*y + z*z);
    }
    
    Vector3D normalize() const {
        double mag = magnitude();
        if (mag > 0) {
            return Vector3D(x/mag, y/mag, z/mag);
        }
        return *this;
    }
    
    Vector3D operator+(const Vector3D& other) const {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }
    
    Vector3D operator-(const Vector3D& other) const {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }
    
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    Vector3D operator/(double scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }
    
    double operator[](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    
    double& operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    
    double dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    
    Vector3D cross(const Vector3D& other) const {
        return Vector3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }
};

// Struttura per data/ora giuliana
struct JulianDate {
    double jd;
    
    JulianDate() : jd(0) {}
    explicit JulianDate(double j) : jd(j) {}
    
    double toMJD() const { return jd - 2400000.5; }
    static JulianDate fromMJD(double mjd) { return JulianDate(mjd + 2400000.5); }
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_TYPES_H
