#ifndef IOCCULTCALC_COORDINATES_H
#define IOCCULTCALC_COORDINATES_H

#include "types.h"

namespace ioccultcalc {

class Coordinates {
public:
    // Converte coordinate equatoriali in vettore cartesiano
    static Vector3D equatorialToCartesian(const EquatorialCoordinates& eq);
    
    // Converte vettore cartesiano in coordinate equatoriali
    static EquatorialCoordinates cartesianToEquatorial(const Vector3D& vec);
    
    // Converte coordinate geografiche in vettore cartesiano terrestre (ECEF)
    static Vector3D geographicToECEF(const GeographicCoordinates& geo);
    
    // Converte vettore ECEF in coordinate geografiche
    static GeographicCoordinates ecefToGeographic(const Vector3D& ecef);
    
    // Calcola la distanza angolare tra due punti sulla sfera celeste
    static double angularSeparation(const EquatorialCoordinates& pos1,
                                   const EquatorialCoordinates& pos2);
    
    // Applica precessione delle coordinate (J2000 -> data)
    static EquatorialCoordinates applyPrecession(const EquatorialCoordinates& pos,
                                                 const JulianDate& fromEpoch,
                                                 const JulianDate& toEpoch);
    
    // Applica aberrazione annuale
    static EquatorialCoordinates applyAberration(const EquatorialCoordinates& pos,
                                                const Vector3D& earthVelocity);
    
    // Calcola l'angolo di posizione tra due oggetti celesti
    static double positionAngle(const EquatorialCoordinates& from,
                               const EquatorialCoordinates& to);
    
    // ===== NUOVE FUNZIONI PER ORBIT DETERMINATION =====
    
    // Converte RA/Dec in vettore unitario (per osservazioni)
    static Vector3D raDecToUnitVector(double ra_rad, double dec_rad);
    
    // Converte vettore in RA/Dec (radianti)
    static void vectorToRaDec(const Vector3D& vec, double& ra_rad, double& dec_rad);
    
    // Posizione osservatore eliocentrica (include posizione Terra)
    // mpc_code: codice osservatorio MPC (es. "500" = geocentro, "C51" = WISE)
    // Restituisce posizione in coordinate EQUATORIALI J2000
    static Vector3D observerPosition(const std::string& mpc_code, 
                                    const JulianDate& jd);
    
    // Converte coordinate equatoriali → eclittiche (J2000)
    // Rotazione di -ε (obliquità eclittica) attorno all'asse X
    static Vector3D equatorialToEcliptic(const Vector3D& eq);
    
    // Converte coordinate eclittiche → equatoriali (J2000)
    // Rotazione di +ε attorno all'asse X
    static Vector3D eclipticToEquatorial(const Vector3D& ecl);
    
    // Converte RA/Dec in vettore unitario ECLITTICO
    // Input: RA/Dec in radianti (equatoriali)
    // Output: vettore unitario in coordinate eclittiche
    static Vector3D raDecToEclipticUnitVector(double ra_rad, double dec_rad);
    
    // Posizione osservatore da coordinate geografiche
    static Vector3D observerPositionFromGeo(const GeographicCoordinates& geo,
                                           const JulianDate& jd);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_COORDINATES_H
