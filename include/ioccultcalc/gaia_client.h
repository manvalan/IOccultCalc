#ifndef IOCCULTCALC_GAIA_CLIENT_H
#define IOCCULTCALC_GAIA_CLIENT_H

#include "types.h"
#include <string>
#include <vector>
#include <memory>

namespace ioccultcalc {

// Enum per le versioni del catalogo Gaia
enum class GaiaVersion {
    EDR3,  // Early Data Release 3 (2020)
    DR3    // Data Release 3 (2022) - default
};

// Dati di una stella da Gaia DR3/EDR3
struct GaiaStar {
    std::string sourceId;      // Gaia source_id
    EquatorialCoordinates pos; // RA, Dec (J2016.0)
    double parallax;           // milliarcsec
    double pmra;               // proper motion RA (mas/yr)
    double pmdec;              // proper motion Dec (mas/yr)
    double phot_g_mean_mag;    // G magnitude
    double phot_bp_mean_mag;   // BP magnitude
    double phot_rp_mean_mag;   // RP magnitude
    
    GaiaStar() : parallax(0), pmra(0), pmdec(0), 
                 phot_g_mean_mag(99), phot_bp_mean_mag(99), phot_rp_mean_mag(99) {}
    
    // Propaga la posizione della stella a una data epoca
    EquatorialCoordinates propagateToEpoch(const JulianDate& epoch) const;
};

class GaiaClient {
public:
    GaiaClient();
    ~GaiaClient();
    
    // Query stelle in una regione rettangolare del cielo
    // ra, dec in gradi, radius in gradi
    std::vector<GaiaStar> queryRegion(double raCenterDeg, double decCenterDeg, 
                                     double radiusDeg, double maxMagnitude = 18.0);
    
    // Query stelle in un cono
    std::vector<GaiaStar> queryCone(double raCenterDeg, double decCenterDeg, 
                                   double radiusDeg, double maxMagnitude = 18.0);
    
    // Query stella singola per source_id
    GaiaStar queryStar(const std::string& sourceId);
    
    // Cerca stelle lungo un percorso (path of occultation)
    // path: lista di coordinate RA/Dec lungo il percorso
    std::vector<GaiaStar> queryAlongPath(const std::vector<EquatorialCoordinates>& path,
                                        double widthDeg, double maxMagnitude = 18.0);
    
    // Imposta l'URL del servizio TAP di Gaia
    void setTAPURL(const std::string& url);
    
    // Imposta timeout per le query (secondi)
    void setTimeout(int seconds);
    
    // Imposta il numero massimo di righe restituite
    void setMaxRows(int rows);
    
    // Imposta la versione del catalogo Gaia (EDR3 o DR3)
    void setGaiaVersion(GaiaVersion version);
    
    // Ottieni la versione attualmente configurata
    GaiaVersion getGaiaVersion() const;
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
    // Esegue una query ADQL sul TAP service
    std::string executeADQL(const std::string& adqlQuery);
    
    // Parsing della risposta VOTable
    std::vector<GaiaStar> parseVOTable(const std::string& votable);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_GAIA_CLIENT_H
