#ifndef IOCCULTCALC_ASTDYS_CLIENT_H
#define IOCCULTCALC_ASTDYS_CLIENT_H

#include "orbital_elements.h"
#include <string>
#include <memory>
#include <vector>

namespace ioccultcalc {

// Forward declarations
struct OrbitState;

class AstDysClient {
public:
    AstDysClient();
    ~AstDysClient();
    
    // Scarica elementi orbitali equinoziali per un asteroide specifico
    // designation può essere numero (es. "433") o designazione (es. "2024 AA")
    EquinoctialElements getElements(const std::string& designation);
    
    // Carica elementi equinoziali da file .eq1 locale
    EquinoctialElements getElementsFromFile(const std::string& filepath);
    
    // Carica osservazioni da file .rwo locale
    std::vector<std::string> getObservationsFromFile(const std::string& filepath);
    
    // Ottiene osservazioni (da file locale o download)
    std::vector<std::string> getObservations(const std::string& designation);
    
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
    
    // Scarica elementi per una lista di asteroidi
    std::vector<EquinoctialElements> getElementsBatch(const std::vector<std::string>& designations);
    
    // Cerca asteroidi per nome (restituisce lista di possibili match)
    std::vector<std::string> searchByName(const std::string& name);
    
    // Imposta l'URL base di AstDyS (default: https://newton.spacedys.com/astdys2/)
    void setBaseURL(const std::string& url);
    
    // Imposta directory locale per file .eq1 (se vuoto, usa download HTTP)
    void setLocalEQ1Directory(const std::string& directory);
    
    // Imposta directory locale per file .rwo (se vuoto, usa download HTTP)
    void setLocalRWODirectory(const std::string& directory);
    
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
