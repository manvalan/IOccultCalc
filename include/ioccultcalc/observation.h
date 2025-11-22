#ifndef IOCCULTCALC_OBSERVATION_H
#define IOCCULTCALC_OBSERVATION_H

#include "types.h"
#include <string>
#include <vector>

namespace ioccultcalc {

// Formati di osservazione
enum class ObservationFormat {
    MPC_80COL,      // Minor Planet Center 80-column format
    MPC_ADES,       // Astrometric Data Exchange Standard
    ASTDYS_RWO      // AstDyS .rwo format
};

// Osservazione astrometrica singola
struct AstrometricObservation {
    JulianDate epoch;           // Tempo dell'osservazione
    EquatorialCoordinates obs;  // RA, Dec osservate
    
    // Observatory information
    std::string observatoryCode; // Codice MPC (es. "500" = Geocentro)
    GeographicCoordinates observerLocation;
    
    // Measurement details
    double raError;    // Errore RA (arcsec)
    double decError;   // Errore Dec (arcsec)
    double magnitude;  // Magnitudine osservata
    
    // Metadata
    std::string catalogCode;  // Catalogo stelle usato (es. "V" = Gaia)
    std::string note1;        // Note MPC
    std::string note2;
    std::string discoveryFlag;
    
    // Computed fields
    EquatorialCoordinates computed; // RA, Dec calcolate
    double raResidual;   // O-C RA (arcsec)
    double decResidual;  // O-C Dec (arcsec)
    double totalResidual; // Residuo totale (arcsec)
    bool outlier;        // Flagged as outlier
    bool excluded;       // Esclusa dal fit
    
    AstrometricObservation() 
        : raError(1.0), decError(1.0), magnitude(99.9),
          raResidual(0), decResidual(0), totalResidual(0), 
          outlier(false), excluded(false) {}
    
    // Peso dell'osservazione per il fit
    double getWeight() const {
        if (outlier) return 0.0;
        // Peso inversamente proporzionale all'errore
        return 1.0 / (raError * raError + decError * decError);
    }
};

// Set di osservazioni per un asteroide
struct ObservationSet {
    std::string objectDesignation;
    std::vector<AstrometricObservation> observations;
    
    JulianDate firstObservation;
    JulianDate lastObservation;
    double arcLength; // giorni
    int numberOfObservations;
    int numberOfObservatories;
    
    ObservationSet() : arcLength(0), numberOfObservations(0), numberOfObservatories(0) {}
    
    // Calcola statistiche
    void computeStatistics();
    
    // Filtra outliers (sigma clipping)
    void filterOutliers(double sigmaThreshold = 3.0);
    
    // Ottiene osservazioni in un intervallo temporale
    std::vector<AstrometricObservation> getObservationsInRange(
        const JulianDate& start, const JulianDate& end) const;
    
    // Calcola RMS dei residui
    double getRMSResidual() const;
    double getRAResidualRMS() const;
    double getDecResidualRMS() const;
    
    // Salva osservazioni su file (con residui se presenti)
    void saveToFile(const std::string& filename) const;
};

// Catalogo di osservatori MPC
struct Observatory {
    std::string code;           // Codice MPC
    std::string name;           // Nome
    GeographicCoordinates location;
    bool isSpacecraft;          // Se è un osservatorio spaziale
    
    // Coordinate geocentriche cartesiane (dal file MPC)
    // In unità di raggio equatoriale terrestre (a_earth = 6378.137 km)
    double rho_cos_phi;  // ρ·cos(φ') - componente equatoriale
    double rho_sin_phi;  // ρ·sin(φ') - componente polare
    
    static Observatory fromMPCCode(const std::string& code);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_OBSERVATION_H
