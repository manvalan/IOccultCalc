#include "ioccultcalc/coordinates.h"
#include "ioccultcalc/ephemeris.h"
#include <cmath>

namespace ioccultcalc {

Vector3D Coordinates::equatorialToCartesian(const EquatorialCoordinates& eq) {
    double r = eq.distance > 0 ? eq.distance : 1.0;
    return Vector3D(
        r * cos(eq.dec) * cos(eq.ra),
        r * cos(eq.dec) * sin(eq.ra),
        r * sin(eq.dec)
    );
}

EquatorialCoordinates Coordinates::cartesianToEquatorial(const Vector3D& vec) {
    double r = vec.magnitude();
    double dec = asin(vec.z / r);
    double ra = atan2(vec.y, vec.x);
    
    if (ra < 0) ra += 2.0 * M_PI;
    
    return EquatorialCoordinates(ra, dec, r);
}

Vector3D Coordinates::geographicToECEF(const GeographicCoordinates& geo) {
    // WGS84 ellipsoid parameters
    const double a = 6378.137; // km, equatorial radius
    const double f = 1.0 / 298.257223563; // flattening
    const double b = a * (1.0 - f); // polar radius
    const double e2 = 2.0 * f - f * f; // first eccentricity squared
    
    double lat = geo.latitude * DEG_TO_RAD;
    double lon = geo.longitude * DEG_TO_RAD;
    double h = geo.altitude / 1000.0; // converti in km
    
    double N = a / sqrt(1.0 - e2 * sin(lat) * sin(lat));
    
    return Vector3D(
        (N + h) * cos(lat) * cos(lon),
        (N + h) * cos(lat) * sin(lon),
        (N * (1.0 - e2) + h) * sin(lat)
    );
}

GeographicCoordinates Coordinates::ecefToGeographic(const Vector3D& ecef) {
    // WGS84 parameters
    const double a = 6378.137;
    const double f = 1.0 / 298.257223563;
    const double b = a * (1.0 - f);
    const double e2 = 2.0 * f - f * f;
    
    double p = sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    double theta = atan2(ecef.z * a, p * b);
    
    double lat = atan2(ecef.z + e2 * b * pow(sin(theta), 3),
                      p - e2 * a * pow(cos(theta), 3));
    double lon = atan2(ecef.y, ecef.x);
    
    double N = a / sqrt(1.0 - e2 * sin(lat) * sin(lat));
    double h = p / cos(lat) - N;
    
    return GeographicCoordinates(
        lon * RAD_TO_DEG,
        lat * RAD_TO_DEG,
        h * 1000.0 // converti in metri
    );
}

double Coordinates::angularSeparation(const EquatorialCoordinates& pos1,
                                     const EquatorialCoordinates& pos2) {
    // Formula haversine per maggiore precisione a piccole separazioni
    double dra = pos2.ra - pos1.ra;
    double ddec = pos2.dec - pos1.dec;
    
    double a = sin(ddec / 2.0) * sin(ddec / 2.0) +
               cos(pos1.dec) * cos(pos2.dec) *
               sin(dra / 2.0) * sin(dra / 2.0);
    
    return 2.0 * atan2(sqrt(a), sqrt(1.0 - a));
}

EquatorialCoordinates Coordinates::applyPrecession(const EquatorialCoordinates& pos,
                                                   const JulianDate& fromEpoch,
                                                   const JulianDate& toEpoch) {
    // Precessione semplificata (accurata per poche decadi)
    // Per precisione maggiore si dovrebbero usare le matrici IAU 2006
    
    double T0 = (fromEpoch.jd - 2451545.0) / 36525.0;
    double T = (toEpoch.jd - fromEpoch.jd) / 36525.0;
    
    // Parametri di precessione in arcsec
    double zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T * T * T) / 3600.0 * DEG_TO_RAD;
    double z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T * T * T) / 3600.0 * DEG_TO_RAD;
    double theta = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T * T * T) / 3600.0 * DEG_TO_RAD;
    
    // Converti in cartesiano, applica rotazione, riconverti
    Vector3D v = equatorialToCartesian(pos);
    
    // Matrice di precessione (semplificata per angoli piccoli)
    double cosZ = cos(z), sinZ = sin(z);
    double cosTheta = cos(theta), sinTheta = sin(theta);
    double cosZeta = cos(zeta), sinZeta = sin(zeta);
    
    Vector3D vPrec(
        v.x * (cosZ * cosTheta * cosZeta - sinZ * sinZeta) -
        v.y * (cosZ * cosTheta * sinZeta + sinZ * cosZeta) -
        v.z * cosZ * sinTheta,
        
        v.x * (sinZ * cosTheta * cosZeta + cosZ * sinZeta) -
        v.y * (sinZ * cosTheta * sinZeta - cosZ * cosZeta) -
        v.z * sinZ * sinTheta,
        
        v.x * sinTheta * cosZeta -
        v.y * sinTheta * sinZeta +
        v.z * cosTheta
    );
    
    return cartesianToEquatorial(vPrec);
}

EquatorialCoordinates Coordinates::applyAberration(const EquatorialCoordinates& pos,
                                                   const Vector3D& earthVelocity) {
    // Aberrazione annuale (stella)
    // earthVelocity in AU/day
    
    Vector3D dir = equatorialToCartesian(pos).normalize();
    double v = earthVelocity.magnitude() * AU / 86400.0; // km/s
    
    // Approssimazione al primo ordine
    double k = v / C_LIGHT;
    
    Vector3D vNorm = earthVelocity.normalize();
    double cosPhi = dir.dot(vNorm);
    
    // Correzione di aberrazione
    Vector3D dirCorr = dir + vNorm * k * (1.0 - cosPhi);
    dirCorr = dirCorr.normalize();
    
    return cartesianToEquatorial(dirCorr);
}

double Coordinates::positionAngle(const EquatorialCoordinates& from,
                                 const EquatorialCoordinates& to) {
    // Angolo di posizione da "from" verso "to"
    double dra = to.ra - from.ra;
    
    double pa = atan2(sin(dra),
                     cos(from.dec) * tan(to.dec) - sin(from.dec) * cos(dra));
    
    // Normalizza a [0, 2π)
    if (pa < 0) pa += 2.0 * M_PI;
    
    return pa;
}

// ===== NUOVE FUNZIONI PER ORBIT DETERMINATION =====

Vector3D Coordinates::raDecToUnitVector(double ra_rad, double dec_rad) {
    // Converte RA/Dec in vettore unitario
    // Sistema equatoriale J2000: x verso vernal equinox, z verso polo nord
    return Vector3D(
        cos(dec_rad) * cos(ra_rad),
        cos(dec_rad) * sin(ra_rad),
        sin(dec_rad)
    );
}

void Coordinates::vectorToRaDec(const Vector3D& vec, double& ra_rad, double& dec_rad) {
    // Normalizza il vettore
    Vector3D unit = vec.normalize();
    
    // Calcola Dec
    dec_rad = asin(unit.z);
    
    // Calcola RA
    ra_rad = atan2(unit.y, unit.x);
    
    // Normalizza RA a [0, 2π)
    if (ra_rad < 0) {
        ra_rad += 2.0 * M_PI;
    }
}

Vector3D Coordinates::observerPosition(const std::string& mpc_code,
                                      const JulianDate& jd) {
    // Database osservatori MPC (subset dei più usati)
    // Formato: {longitude, latitude, altitude_meters}
    // Coordinate dal Minor Planet Center
    
    struct ObservatoryData {
        std::string code;
        double lon;  // gradi Est
        double lat;  // gradi Nord
        double alt;  // metri
    };
    
    static const ObservatoryData observatories[] = {
        // Codice speciale per geocentro
        {"500", 0.0, 0.0, 0.0},
        
        // Osservatori professionali principali
        {"C51", -110.73443, 31.96103, 2096.0},  // WISE (Mt. Lemmon)
        {"D0", -70.8049, -30.2417, 2282.0},      // Cerro Tololo (LSST site)
        {"T0", -70.7494, -30.2444, 2207.0},      // Cerro Tololo 0.9m
        {"W6", -155.4681, 19.8259, 3055.0},      // Mauna Kea (Pan-STARRS)
        {"M2", -155.5762, 19.8263, 4213.0},      // Mauna Kea (Subaru)
        
        // Osservatori amatoriali attivi
        {"C4", 13.7319, 44.8847, 52.0},          // Farra d'Isonzo (Italia)
        {"204", 14.2958, 40.8644, 458.0},        // Caserta (Italia)
        {"A77", 11.9586, 43.7514, 184.0},        // San Marcello (Italia)
        
        // Osservatori storici
        {"000", 0.0, 51.477, 46.0},              // Greenwich
        {"675", -17.8792, 28.7606, 2326.0},      // La Palma
    };
    
    // Cerca osservatorio nel database
    GeographicCoordinates geo;
    bool found = false;
    
    for (const auto& obs : observatories) {
        if (obs.code == mpc_code) {
            geo.longitude = obs.lon;
            geo.latitude = obs.lat;
            geo.altitude = obs.alt;
            found = true;
            break;
        }
    }
    
    if (!found) {
        // Default: geocentro
        geo.longitude = 0.0;
        geo.latitude = 0.0;
        geo.altitude = 0.0;
    }
    
    return observerPositionFromGeo(geo, jd);
}

Vector3D Coordinates::observerPositionFromGeo(const GeographicCoordinates& geo,
                                             const JulianDate& jd) {
    // Calcola posizione osservatore in coordinate heliocentriche
    // 1. Posizione geocentrica dell'osservatore (rotazione Terra)
    // 2. Posizione eliocentrica della Terra
    // 3. Somma per posizione eliocentrica osservatore
    
    // GMST (Greenwich Mean Sidereal Time)
    double T = (jd.jd - 2451545.0) / 36525.0;
    double gmst = 280.46061837 + 360.98564736629 * (jd.jd - 2451545.0) +
                  0.000387933 * T * T - T * T * T / 38710000.0;
    
    // Normalizza a [0, 360)
    gmst = fmod(gmst, 360.0);
    if (gmst < 0) gmst += 360.0;
    
    // LAST (Local Apparent Sidereal Time)
    double last = gmst + geo.longitude;
    last = fmod(last, 360.0);
    if (last < 0) last += 360.0;
    
    double last_rad = last * DEG_TO_RAD;
    
    // Converti coordinate geografiche in ECEF
    Vector3D ecef = geographicToECEF(geo);
    
    // Ruota da sistema ECEF a sistema equatoriale
    // Rotazione attorno asse Z di LAST
    double cos_last = cos(last_rad);
    double sin_last = sin(last_rad);
    
    Vector3D eq_geocentric(
        ecef.x * cos_last - ecef.y * sin_last,
        ecef.x * sin_last + ecef.y * cos_last,
        ecef.z
    );
    
    // Converti da km ad AU (posizione geocentrica)
    eq_geocentric = eq_geocentric * (1.0 / AU);
    
    // Ottieni posizione eliocentrica della Terra (restituisce equatoriali)
    Vector3D earth_helio_eq = Ephemeris::getEarthPosition(jd);
    
    // Posizione eliocentrica osservatore = Terra + offset geocentrico
    // Entrambi in coordinate EQUATORIALI J2000
    return earth_helio_eq + eq_geocentric;
}

Vector3D Coordinates::equatorialToEcliptic(const Vector3D& eq) {
    // Obliquità dell'eclittica (ε) per J2000.0
    // IAU 2000: ε₀ = 23° 26' 21.406" = 23.4392794444°
    constexpr double epsilon = 23.4392794444 * DEG_TO_RAD;
    
    double cos_eps = cos(epsilon);
    double sin_eps = sin(epsilon);
    
    // Rotazione attorno all'asse X di -ε
    // | x' |   |  1      0          0      | | x |
    // | y' | = |  0   cos(ε)    sin(ε)    | | y |
    // | z' |   |  0  -sin(ε)    cos(ε)    | | z |
    
    return Vector3D(
        eq.x,
        eq.y * cos_eps + eq.z * sin_eps,
        -eq.y * sin_eps + eq.z * cos_eps
    );
}

Vector3D Coordinates::eclipticToEquatorial(const Vector3D& ecl) {
    // Obliquità dell'eclittica (ε) per J2000.0
    constexpr double epsilon = 23.4392794444 * DEG_TO_RAD;
    
    double cos_eps = cos(epsilon);
    double sin_eps = sin(epsilon);
    
    // Rotazione attorno all'asse X di +ε (inversa della precedente)
    // | x' |   |  1      0          0      | | x |
    // | y' | = |  0   cos(ε)   -sin(ε)    | | y |
    // | z' |   |  0   sin(ε)    cos(ε)    | | z |
    
    return Vector3D(
        ecl.x,
        ecl.y * cos_eps - ecl.z * sin_eps,
        ecl.y * sin_eps + ecl.z * cos_eps
    );
}

Vector3D Coordinates::raDecToEclipticUnitVector(double ra_rad, double dec_rad) {
    // Prima converti RA/Dec in vettore equatoriale
    Vector3D eq = raDecToUnitVector(ra_rad, dec_rad);
    
    // Poi converti in eclittico
    return equatorialToEcliptic(eq);
}

} // namespace ioccultcalc
