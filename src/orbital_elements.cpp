#include "ioccultcalc/orbital_elements.h"
#include <cmath>

namespace ioccultcalc {

void EquinoctialElements::toKeplerian(double& ecc, double& inc, double& omega,
                                     double& Omega, double& M) const {
    // Converti elementi equinoziali in Kepleriani
    
    // Eccentricità
    ecc = sqrt(h * h + k * k);
    
    // Inclinazione
    inc = 2.0 * atan(sqrt(p * p + q * q));
    
    // Longitudine del nodo ascendente
    Omega = atan2(p, q);
    if (Omega < 0) Omega += 2.0 * M_PI;
    
    // Argomento del periapside
    double omega_plus_Omega = atan2(h, k);
    if (omega_plus_Omega < 0) omega_plus_Omega += 2.0 * M_PI;
    
    omega = omega_plus_Omega - Omega;
    if (omega < 0) omega += 2.0 * M_PI;
    
    // Anomalia media
    M = lambda - omega_plus_Omega;
    while (M < 0) M += 2.0 * M_PI;
    while (M >= 2.0 * M_PI) M -= 2.0 * M_PI;
}

EquinoctialElements EquinoctialElements::fromKeplerian(double a, double ecc, double inc,
                                                      double omega, double Omega, double M,
                                                      const JulianDate& epoch) {
    EquinoctialElements elem;
    elem.a = a;
    elem.epoch = epoch;
    
    // Converti elementi Kepleriani in equinoziali
    double omega_plus_Omega = omega + Omega;
    
    elem.h = ecc * sin(omega_plus_Omega);
    elem.k = ecc * cos(omega_plus_Omega);
    elem.p = tan(inc / 2.0) * sin(Omega);
    elem.q = tan(inc / 2.0) * cos(Omega);
    elem.lambda = M + omega_plus_Omega;
    
    // Normalizza lambda a [0, 2π)
    while (elem.lambda < 0) elem.lambda += 2.0 * M_PI;
    while (elem.lambda >= 2.0 * M_PI) elem.lambda -= 2.0 * M_PI;
    
    return elem;
}

// Nuovo metodo che ritorna OrbitalElements
OrbitalElements EquinoctialElements::toKeplerian() const {
    OrbitalElements orb;
    orb.a = a;
    orb.epoch = epoch;
    orb.H = H;
    orb.G = G;
    orb.diameter = diameter;
    orb.designation = designation;
    orb.name = name;
    
    // Converti elementi
    toKeplerian(orb.e, orb.i, orb.omega, orb.Omega, orb.M);
    
    return orb;
}

EquinoctialElements EquinoctialElements::fromKeplerian(const OrbitalElements& orb) {
    return fromKeplerian(orb.a, orb.e, orb.i, orb.omega, orb.Omega, orb.M, orb.epoch);
}

EquinoctialElements OrbitalElements::toEquinoctial() const {
    return EquinoctialElements::fromKeplerian(*this);
}

OrbitalElements OrbitalElements::fromEquinoctial(const EquinoctialElements& eq) {
    return eq.toKeplerian();
}

// ===== FUNZIONI PER ORBIT DETERMINATION =====

/**
 * @brief Risolve l'equazione di Keplero per anomalia eccentrica
 * 
 * Usa metodo di Newton-Raphson per risolvere: E - e*sin(E) = M
 */
static double solveKeplerEquation(double M, double e, double tolerance = 1e-12) {
    // Guess iniziale
    double E = M + e * sin(M);
    
    // Newton-Raphson
    for (int iter = 0; iter < 50; iter++) {
        double f = E - e * sin(E) - M;
        double fp = 1.0 - e * cos(E);
        
        if (std::abs(fp) < 1e-15) break;
        
        double dE = -f / fp;
        E += dE;
        
        if (std::abs(dE) < tolerance) break;
    }
    
    return E;
}

OrbitalElements cartesianToOrbitalElements(
    const Vector3D& position,
    const Vector3D& velocity,
    double mu,
    const JulianDate& epoch)
{
    OrbitalElements elem;
    elem.epoch = epoch;
    
    // Vettori fondamentali
    double r = position.magnitude();
    double v = velocity.magnitude();
    
    // Momento angolare specifico: h = r × v
    Vector3D h_vec = position.cross(velocity);
    double h = h_vec.magnitude();
    
    // Vettore del nodo: n = k × h (dove k = [0, 0, 1])
    Vector3D k(0, 0, 1);
    Vector3D n_vec = k.cross(h_vec);
    double n = n_vec.magnitude();
    
    // Vettore eccentricità: e = (v × h)/μ - r/|r|
    Vector3D e_vec = velocity.cross(h_vec) * (1.0 / mu) - position * (1.0 / r);
    elem.e = e_vec.magnitude();
    
    // Energia specifica
    double epsilon = v * v / 2.0 - mu / r;
    
    // Semi-asse maggiore
    if (std::abs(epsilon) > 1e-12) {
        elem.a = -mu / (2.0 * epsilon);
    } else {
        // Orbita parabolica (e ≈ 1)
        elem.a = h * h / mu;
    }
    
    // Inclinazione
    elem.i = acos(h_vec.z / h);
    
    // Longitudine del nodo ascendente
    if (n > 1e-10) {
        elem.Omega = acos(n_vec.x / n);
        if (n_vec.y < 0) {
            elem.Omega = 2.0 * M_PI - elem.Omega;
        }
    } else {
        // Orbita equatoriale
        elem.Omega = 0.0;
    }
    
    // Argomento del periapside
    if (n > 1e-10 && elem.e > 1e-10) {
        elem.omega = acos(n_vec.dot(e_vec) / (n * elem.e));
        if (e_vec.z < 0) {
            elem.omega = 2.0 * M_PI - elem.omega;
        }
    } else if (elem.e > 1e-10) {
        // Orbita equatoriale ma non circolare
        elem.omega = atan2(e_vec.y, e_vec.x);
        if (elem.omega < 0) elem.omega += 2.0 * M_PI;
    } else {
        // Orbita circolare
        elem.omega = 0.0;
    }
    
    // Anomalia vera
    double nu;
    if (elem.e > 1e-10) {
        nu = acos(e_vec.dot(position) / (elem.e * r));
        if (position.dot(velocity) < 0) {
            nu = 2.0 * M_PI - nu;
        }
    } else {
        // Orbita circolare: usa argomento di latitudine
        if (n > 1e-10) {
            nu = acos(n_vec.dot(position) / (n * r));
            if (position.z < 0) {
                nu = 2.0 * M_PI - nu;
            }
        } else {
            nu = atan2(position.y, position.x);
            if (nu < 0) nu += 2.0 * M_PI;
        }
    }
    
    // Anomalia eccentrica
    double E;
    if (elem.e < 1.0) {
        // Ellisse
        double cos_E = (elem.e + cos(nu)) / (1.0 + elem.e * cos(nu));
        double sin_E = sqrt(1.0 - elem.e * elem.e) * sin(nu) / (1.0 + elem.e * cos(nu));
        E = atan2(sin_E, cos_E);
    } else if (elem.e > 1.0) {
        // Iperbole
        double sinh_F = sqrt(elem.e * elem.e - 1.0) * sin(nu) / (1.0 + elem.e * cos(nu));
        double F = asinh(sinh_F);
        elem.M = elem.e * sinh(F) - F;
        return elem;
    } else {
        // Parabola
        double D = sqrt(2.0) * tan(nu / 2.0);
        elem.M = D + D * D * D / 3.0;
        return elem;
    }
    
    // Anomalia media
    elem.M = E - elem.e * sin(E);
    
    // Normalizza a [0, 2π)
    while (elem.M < 0) elem.M += 2.0 * M_PI;
    while (elem.M >= 2.0 * M_PI) elem.M -= 2.0 * M_PI;
    
    return elem;
}

void orbitalElementsToCartesian(
    const OrbitalElements& elements,
    Vector3D& position,
    Vector3D& velocity,
    double mu)
{
    double a = elements.a;
    double e = elements.e;
    double i = elements.i;
    double Omega = elements.Omega;
    double omega = elements.omega;
    double M = elements.M;
    
    // Risolvi equazione di Keplero per E
    double E = solveKeplerEquation(M, e);
    
    // Coordinate nel piano orbitale
    double cos_E = cos(E);
    double sin_E = sin(E);
    
    double x_orb = a * (cos_E - e);
    double y_orb = a * sqrt(1.0 - e * e) * sin_E;
    
    // Velocità nel piano orbitale
    double n = sqrt(mu / (a * a * a));  // Moto medio
    double vx_orb = -a * n * sin_E / (1.0 - e * cos_E);
    double vy_orb = a * n * sqrt(1.0 - e * e) * cos_E / (1.0 - e * cos_E);
    
    // Matrici di rotazione
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_i = cos(i);
    double sin_i = sin(i);
    
    // Rotazione: orbitale → eclittico
    // R_z(Omega) * R_x(i) * R_z(omega)
    
    position.x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * x_orb +
                 (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * y_orb;
    
    position.y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * x_orb +
                 (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * y_orb;
    
    position.z = sin_omega * sin_i * x_orb + cos_omega * sin_i * y_orb;
    
    velocity.x = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * vx_orb +
                 (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * vy_orb;
    
    velocity.y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * vx_orb +
                 (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * vy_orb;
    
    velocity.z = sin_omega * sin_i * vx_orb + cos_omega * sin_i * vy_orb;
}

} // namespace ioccultcalc
