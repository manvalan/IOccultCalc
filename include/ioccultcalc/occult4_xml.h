#ifndef IOCCULTCALC_OCCULT4_XML_H
#define IOCCULTCALC_OCCULT4_XML_H

#include "occultation_predictor.h"
#include <string>
#include <vector>

namespace ioccultcalc {

/**
 * @brief Gestione file XML compatibili con Occult4
 * 
 * Questo modulo permette di importare ed esportare predizioni di occultazioni
 * nel formato XML utilizzato da Occult4 di Dave Herald.
 * 
 * Formato supportato:
 * - Occult4 Prediction XML (versione 4.x)
 * - Steve Preston XML format
 * - IOTA XML standard
 */
class Occult4XMLHandler {
public:
    Occult4XMLHandler();
    ~Occult4XMLHandler();
    
    /**
     * @brief Struttura dati per evento Occult4
     */
    struct Occult4Event {
        // Identificazione
        std::string asteroidNumber;      // Numero MPC
        std::string asteroidName;        // Nome asteroide
        std::string asteroidDesignation; // Designazione completa
        
        // Stella
        std::string starCatalog;         // GAIA, UCAC4, TYC, etc.
        std::string starId;              // ID nel catalogo
        double starRA;                   // RA J2000 (gradi)
        double starDec;                  // Dec J2000 (gradi)
        double starMag;                  // Magnitudine
        
        // Tempo
        double jdEvent;                  // JD dell'evento
        std::string dateTimeUTC;         // Data/ora ISO8601
        
        // Geometria
        double closeApproachDist;        // Distanza minima (arcsec)
        double posAngle;                 // Angolo di posizione (gradi)
        double pathWidth;                // Larghezza percorso (km)
        double maxDuration;              // Durata massima (secondi)
        double uncertainty;              // Incertezza 1-sigma (km)
        
        // Shadow path (linea centrale)
        struct PathPoint {
            double latitude;             // Latitudine (gradi)
            double longitude;            // Longitudine (gradi)
            double jd;                   // JD del punto
            std::string dateTime;        // Data/ora UTC
            double altitude;             // Altezza stella (gradi)
            double sunAltitude;          // Altezza Sole (gradi)
        };
        std::vector<PathPoint> centerLine;
        
        // Limiti nord/sud
        std::vector<PathPoint> northLimit;
        std::vector<PathPoint> southLimit;
        
        // Metadati
        double probability;              // Probabilità evento (0-1)
        double starDistance;             // Distanza stella (parsec)
        double asteroidMag;              // Magnitudine asteroide
        double dropMag;                  // Drop di magnitudine
        std::string eventId;             // ID univoco evento
        
        Occult4Event() 
            : starRA(0), starDec(0), starMag(0),
              jdEvent(0), closeApproachDist(0), posAngle(0),
              pathWidth(0), maxDuration(0), uncertainty(0),
              probability(0), starDistance(0), asteroidMag(0), dropMag(0) {}
    };
    
    /**
     * @brief Opzioni per import/export
     */
    struct XMLOptions {
        bool includeUncertainty;         // Includi bande incertezza
        bool includePathPoints;          // Includi punti percorso
        int pathPointsResolution;        // Punti ogni N km
        bool includeStarData;            // Includi dati stella completi
        bool includeAsteroidData;        // Includi dati asteroide
        bool useGaiaIds;                 // Usa ID Gaia (vs UCAC4)
        std::string organizationName;    // Nome organizzazione
        std::string observerName;        // Nome osservatore
        
        XMLOptions()
            : includeUncertainty(true),
              includePathPoints(true),
              pathPointsResolution(100),
              includeStarData(true),
              includeAsteroidData(true),
              useGaiaIds(true),
              organizationName("IOccultCalc"),
              observerName("") {}
    };
    
    // ============================================================================
    // IMPORT da XML
    // ============================================================================
    
    /**
     * @brief Carica eventi da file XML Occult4
     * @param filename Path del file XML
     * @return Vector di eventi caricati
     * @throw std::runtime_error se il parsing fallisce
     */
    std::vector<Occult4Event> loadFromXML(const std::string& filename);
    
    /**
     * @brief Parsa una stringa XML Occult4
     * @param xmlContent Contenuto XML come stringa
     * @return Vector di eventi parsati
     */
    std::vector<Occult4Event> parseXML(const std::string& xmlContent);
    
    /**
     * @brief Converte evento Occult4 in formato IOccultCalc
     * @param occult4Event Evento formato Occult4
     * @return Evento formato IOccultCalc
     */
    OccultationEvent toIOccultCalcEvent(const Occult4Event& occult4Event);
    
    // ============================================================================
    // EXPORT a XML
    // ============================================================================
    
    /**
     * @brief Esporta evento in formato XML Occult4
     * @param event Evento IOccultCalc da esportare
     * @param filename Path del file XML di output
     * @return true se successo, false altrimenti
     */
    bool exportToXML(const OccultationEvent& event, const std::string& filename);
    
    /**
     * @brief Esporta multipli eventi in un unico file XML
     * @param events Vector di eventi da esportare
     * @param filename Path del file XML di output
     * @return true se successo, false altrimenti
     */
    bool exportMultipleToXML(const std::vector<OccultationEvent>& events, 
                            const std::string& filename);
    
    /**
     * @brief Genera stringa XML per singolo evento
     * @param event Evento da convertire
     * @return Stringa XML formattata
     */
    std::string generateXML(const OccultationEvent& event);
    
    /**
     * @brief Genera stringa XML per multipli eventi
     * @param events Eventi da convertire
     * @return Stringa XML formattata
     */
    std::string generateXML(const std::vector<OccultationEvent>& events);
    
    /**
     * @brief Converte evento IOccultCalc in formato Occult4
     * @param event Evento IOccultCalc
     * @return Evento formato Occult4
     */
    Occult4Event toOccult4Event(const OccultationEvent& event);
    
    // ============================================================================
    // CONFIGURAZIONE
    // ============================================================================
    
    /**
     * @brief Imposta opzioni di import/export
     */
    void setOptions(const XMLOptions& options);
    
    /**
     * @brief Ottieni opzioni correnti
     */
    XMLOptions getOptions() const;
    
    /**
     * @brief Valida file XML Occult4
     * @param filename Path del file da validare
     * @return true se valido, false altrimenti
     */
    bool validateXML(const std::string& filename);
    
    /**
     * @brief Ottieni versione formato XML rilevata
     * @param filename Path del file XML
     * @return Stringa versione (es. "4.2.0")
     */
    std::string detectXMLVersion(const std::string& filename);
    
private:
    XMLOptions options_;
    
    // Parsing helpers
    Occult4Event parseEventNode(void* node);
    std::vector<Occult4Event::PathPoint> parsePathPoints(void* node);
    std::string extractTextContent(void* node);
    double extractDoubleContent(void* node, double defaultValue = 0.0);
    void* findChildNode(void* parent, const std::string& name);
    std::vector<void*> findAllChildNodes(void* parent, const std::string& name);
    
    // Generation helpers
    std::string generateEventXML(const Occult4Event& event);
    std::string escapeXML(const std::string& text);
    std::string formatDouble(double value, int precision = 6);
    std::string formatRA(double raDeg);        // HH:MM:SS.sss
    std::string formatDec(double decDeg);      // +DD:MM:SS.ss
    std::string formatDateTime(double jd);     // ISO8601
    
    // Path point generation
    std::vector<Occult4Event::PathPoint> generatePathPoints(
        const std::vector<ShadowPathPoint>& ioPoints);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_OCCULT4_XML_H
