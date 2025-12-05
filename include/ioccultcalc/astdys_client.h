#ifndef IOCCULTCALC_ASTDYS_CLIENT_H
#define IOCCULTCALC_ASTDYS_CLIENT_H

#include "orbital_elements.h"
#include <string>
#include <memory>
#include <vector>
#include <array>

namespace ioccultcalc {

// Forward declarations
struct OrbitState;

/**
 * @struct ChebyshevEphemeris
 * @brief Rappresentazione di effemeridi con polinomi di Chebyshev da AstDyS
 * 
 * Formato: coefficienti di Chebyshev per x, y, z geocentriche
 * Intervallo temporale: tipicamente 1-30 giorni
 * Ordine: 10-20 (sufficiente per precisione sub-arcsec)
 */
struct ChebyshevEphemeris {
    double startJD;                           // Inizio intervallo validità [JD]
    double endJD;                             // Fine intervallo validità [JD]
    int order;                                // Ordine polinomio (numero coefficienti - 1)
    
    std::vector<double> coeff_x;              // Coefficienti Chebyshev per X [AU]
    std::vector<double> coeff_y;              // Coefficienti Chebyshev per Y [AU]
    std::vector<double> coeff_z;              // Coefficienti Chebyshev per Z [AU]
    
    bool geocentric;                          // true = geocentrico, false = eliocentrico
    std::string frame;                        // "ECLIPJ2000" o "ICRF"
    
    /**
     * @brief Valuta posizione a data epoca
     * @param jd Epoca [JD]
     * @return Array {x, y, z} in AU
     */
    std::array<double, 3> evaluate(double jd) const;
    
    /**
     * @brief Valuta posizione e velocità a data epoca
     * @param jd Epoca [JD]
     * @return Array {x, y, z, vx, vy, vz} in AU, AU/day
     */
    std::array<double, 6> evaluateWithVelocity(double jd) const;
    
    /**
     * @brief Converti coordinate cartesiane → RA/Dec
     * @param jd Epoca [JD]
     * @return {RA [deg], Dec [deg], distance [AU]}
     */
    std::array<double, 3> getRADecDist(double jd) const;
    
    /**
     * @brief Verifica se epoca è nell'intervallo di validità
     */
    bool isValid(double jd) const {
        return jd >= startJD && jd <= endJD;
    }
};

class AstDysClient {
public:
    AstDysClient();
    ~AstDysClient();
    
    // Scarica elementi orbitali equinoziali per un asteroide specifico
    // designation può essere numero (es. "433") o designazione (es. "2024 AA")
    EquinoctialElements getElements(const std::string& designation);
    
    // Scarica elementi orbitali KEPLERIAN per un asteroide (epoca recente)
    // Usa il catalogo allnum.cat che ha epoche vicine al presente
    OrbitalElements getRecentElements(const std::string& designation);
    
    // Ottiene elementi OSCULANTI da JPL Horizons per un'epoca specifica
    // Nota: AstDyS fornisce elementi MEDI che richiedono propagazione OrbFit
    // Questo metodo usa Horizons per ottenere elementi osculanti istantanei
    OrbitalElements getOsculatingElements(const std::string& designation, 
                                         const JulianDate& epoch);
    
    // Ottiene stato orbitale (posizione/velocità) da JPL Horizons
    // Metodo PREFERITO per popolare database: più veloce e affidabile
    // degli elementi orbitali, evita problemi di parsing
    OrbitState getStateFromHorizons(const std::string& designation,
                                   const JulianDate& epoch);
    
    // ========== NUOVE API CHEBYSHEV ==========
    
    /**
     * @brief Scarica effemeridi con rappresentazione Chebyshev da AstDyS
     * 
     * Vantaggi rispetto a elementi orbitali:
     * - Precisione sub-arcsec senza propagazione numerica
     * - Velocità: valutazione polinomiale vs integrazione ODE
     * - Coerenza: stessi coefficienti usati da AstDyS per previsioni
     * - Validità: intervallo temporale esplicito (no extrapolazione)
     * 
     * @param designation Numero (es. "433") o designazione (es. "2024 AA")
     * @param startJD Inizio intervallo [JD]
     * @param endJD Fine intervallo [JD]
     * @param order Ordine polinomio (10-20, default: 15)
     * @return Effemeridi Chebyshev per l'intervallo richiesto
     * 
     * Esempio:
     *   auto cheb = client.getChebyshevEphemeris("17030", 
     *                                             2460638.5,  // 2025-11-24
     *                                             2460645.5,  // 2025-12-01
     *                                             15);
     *   auto [ra, dec, dist] = cheb.getRADecDist(2460642.39); // 2025-11-28
     */
    ChebyshevEphemeris getChebyshevEphemeris(const std::string& designation,
                                             double startJD,
                                             double endJD,
                                             int order = 15);
    
    /**
     * @brief Scarica effemeridi Chebyshev multi-segmento
     * 
     * Per intervalli lunghi (> 30 giorni), AstDyS divide in segmenti
     * per mantenere l'accuratezza. Questo metodo gestisce automaticamente
     * i segmenti multipli.
     * 
     * @param designation Asteroide
     * @param startJD Inizio intervallo [JD]
     * @param endJD Fine intervallo [JD]
     * @param segment_days Durata segmento (default: 7 giorni)
     * @param order Ordine polinomio per segmento
     * @return Vector di segmenti Chebyshev contigui
     */
    std::vector<ChebyshevEphemeris> getChebyshevEphemerisMultiSegment(
        const std::string& designation,
        double startJD,
        double endJD,
        double segment_days = 7.0,
        int order = 15);
    
    /**
     * @brief Valuta RA/Dec usando Chebyshev (single call)
     * 
     * Metodo di convenienza per ottenere coordinate senza gestire
     * esplicitamente i coefficienti Chebyshev.
     * 
     * @param designation Asteroide
     * @param jd Epoca [JD]
     * @return {RA [deg], Dec [deg], distance [AU]}
     */
    std::array<double, 3> getRADecChebyshev(const std::string& designation,
                                            double jd);
    
    // ========== FINE NUOVE API CHEBYSHEV ==========
    
    // Scarica elementi per una lista di asteroidi
    std::vector<EquinoctialElements> getElementsBatch(const std::vector<std::string>& designations);
    
    // Cerca asteroidi per nome (restituisce lista di possibili match)
    std::vector<std::string> searchByName(const std::string& name);
    
    // Imposta l'URL base di AstDyS (default: https://newton.spacedys.com/astdys2/)
    void setBaseURL(const std::string& url);
    
    // Imposta timeout per le richieste HTTP (secondi)
    void setTimeout(int seconds);
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
    // Parsing del file .eq (equinoctial elements)
    EquinoctialElements parseEquinoctialFile(const std::string& content, 
                                            const std::string& designation);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ASTDYS_CLIENT_H
