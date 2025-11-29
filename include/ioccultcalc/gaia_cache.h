/**
 * @file gaia_cache.h
 * @brief Sistema di cache locale per stelle GAIA DR3
 * 
 * Permette di scaricare e mantenere una cache locale delle stelle GAIA
 * per evitare query ripetute al server TAP. La cache è organizzata per
 * regioni di cielo (HEALPix) e può essere filtrata per magnitudine.
 */

#ifndef IOCCULTCALC_GAIA_CACHE_H
#define IOCCULTCALC_GAIA_CACHE_H

#include "gaia_client.h"
#include "data_manager.h"
#include <string>
#include <vector>
#include <map>
#include <memory>

namespace ioccultcalc {

/**
 * @brief Statistiche sulla cache GAIA
 */
struct GaiaCacheStats {
    int total_tiles;           // Numero di tiles scaricate
    int total_stars;           // Numero totale di stelle
    double total_size_mb;      // Dimensione totale (MB)
    double sky_coverage;       // Copertura cielo (percentuale)
    double min_magnitude;      // Magnitudine minima
    double max_magnitude;      // Magnitudine massima
    std::string last_update;   // Data ultimo aggiornamento
    
    GaiaCacheStats() : total_tiles(0), total_stars(0), total_size_mb(0),
                       sky_coverage(0), min_magnitude(99), max_magnitude(0) {}
};

/**
 * @brief Singola tile della cache (regione HEALPix)
 */
struct GaiaTile {
    int healpix_id;            // ID HEALPix (NSIDE=32)
    double ra_center;          // RA centro (gradi)
    double dec_center;         // Dec centro (gradi)
    double radius;             // Raggio (gradi)
    int star_count;            // Numero stelle
    double max_magnitude;      // Magnitudine limite
    std::string last_update;   // Data scaricamento
    std::vector<GaiaStar> stars; // Stelle nella tile
    
    GaiaTile() : healpix_id(-1), ra_center(0), dec_center(0), 
                 radius(0), star_count(0), max_magnitude(18.0) {}
};

/**
 * @brief Gestione cache locale GAIA DR3
 * 
 * La cache è organizzata in tiles HEALPix (NSIDE=32, ~13.4 gradi per tile).
 * Ogni tile contiene tutte le stelle fino alla magnitudine specificata.
 * 
 * Directory structure:
 *   ~/.ioccultcalc/gaia/
 *   ├── cache_index.json      # Indice tiles scaricate
 *   ├── tiles/
 *   │   ├── tile_0001.json    # Tile HEALPix 1
 *   │   ├── tile_0002.json    # Tile HEALPix 2
 *   │   └── ...
 *   └── stats.json            # Statistiche cache
 */
class GaiaCache {
public:
    /**
     * @brief Costruttore
     * @param cacheDir Directory cache (default: DataManager path)
     */
    explicit GaiaCache(const std::string& cacheDir = "");
    
    ~GaiaCache();
    
    // ==================== Query Interface ====================
    
    /**
     * @brief Cerca stelle nella cache per regione
     * @param raDeg RA centro (gradi)
     * @param decDeg Dec centro (gradi)
     * @param radiusDeg Raggio ricerca (gradi)
     * @param maxMagnitude Magnitudine limite
     * @param autoDownload Se true, scarica tiles mancanti
     * @return Liste di stelle trovate
     */
    std::vector<GaiaStar> queryRegion(double raDeg, double decDeg, 
                                     double radiusDeg, double maxMagnitude,
                                     bool autoDownload = true);
    
    /**
     * @brief Cerca stelle lungo un percorso
     * @param path Lista di coordinate lungo il percorso
     * @param widthDeg Larghezza banda (gradi)
     * @param maxMagnitude Magnitudine limite
     * @param autoDownload Se true, scarica tiles mancanti
     * @return Liste di stelle trovate
     */
    std::vector<GaiaStar> queryPath(const std::vector<EquatorialCoordinates>& path,
                                   double widthDeg, double maxMagnitude,
                                   bool autoDownload = true);
    
    /**
     * @brief Verifica se una regione è coperta dalla cache
     * @param raDeg RA centro (gradi)
     * @param decDeg Dec centro (gradi)
     * @param radiusDeg Raggio (gradi)
     * @param maxMagnitude Magnitudine richiesta
     * @return true se la regione è completamente coperta
     */
    bool isCovered(double raDeg, double decDeg, double radiusDeg, 
                   double maxMagnitude) const;
    
    // ==================== Cache Management ====================
    
    /**
     * @brief Scarica tiles per una regione
     * @param raDeg RA centro (gradi)
     * @param decDeg Dec centro (gradi)
     * @param radiusDeg Raggio (gradi)
     * @param maxMagnitude Magnitudine limite
     * @param progressCallback Callback per progresso (opzionale)
     * @return Numero di tiles scaricate
     */
    int downloadRegion(double raDeg, double decDeg, double radiusDeg,
                      double maxMagnitude,
                      std::function<void(int current, int total)> progressCallback = nullptr);
    
    /**
     * @brief Scarica tiles lungo un percorso
     * @param path Lista di coordinate
     * @param widthDeg Larghezza banda (gradi)
     * @param maxMagnitude Magnitudine limite
     * @param progressCallback Callback per progresso
     * @return Numero di tiles scaricate
     */
    int downloadPath(const std::vector<EquatorialCoordinates>& path,
                    double widthDeg, double maxMagnitude,
                    std::function<void(int current, int total)> progressCallback = nullptr);
    
    /**
     * @brief Carica indice cache da disco
     * @return true se successo
     */
    bool loadIndex();
    
    /**
     * @brief Salva indice cache su disco
     * @return true se successo
     */
    bool saveIndex();
    
    /**
     * @brief Ottieni statistiche cache
     */
    GaiaCacheStats getStats() const;
    
    /**
     * @brief Pulisci tiles più vecchie di N giorni
     * @param maxAgeDays Età massima (giorni)
     * @return Numero di tiles rimosse
     */
    int cleanOldTiles(int maxAgeDays);
    
    /**
     * @brief Rimuovi tutte le tiles
     */
    void clearCache();
    
    /**
     * @brief Imposta client GAIA per download
     */
    void setGaiaClient(std::shared_ptr<GaiaClient> client);
    
    /**
     * @brief Imposta versione catalogo Gaia (EDR3 o DR3)
     */
    void setGaiaVersion(GaiaVersion version);
    
    /**
     * @brief Abilita/disabilita auto-download
     */
    void setAutoDownload(bool enable);
    
    /**
     * @brief Imposta verbosità (per debug)
     */
    void setVerbose(bool verbose);
    
    // ==================== Utilities (pubbliche per testing) ====================
    
    /**
     * @brief Converti coordinate → HEALPix index
     */
    int coordsToHealpix(double raDeg, double decDeg, int nside = 32) const;
    
    /**
     * @brief Converti HEALPix index → coordinate
     */
    void healpixToCoords(int healpix, int nside, double& raDeg, double& decDeg) const;
    
    /**
     * @brief Query HEALPix tiles in disco
     */
    std::vector<int> queryDisc(double raDeg, double decDeg, double radiusDeg, int nside = 32) const;
    
private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
    
    // Gestione tiles
    bool loadTile(int healpixId);
    bool saveTile(const GaiaTile& tile);
    std::string getTilePath(int healpixId) const;
    bool tileExists(int healpixId) const;
    
    // Download da GAIA
    GaiaTile downloadTile(int healpixId, double maxMagnitude);
};

/**
 * @brief Builder per configurare la cache
 */
class GaiaCacheBuilder {
public:
    GaiaCacheBuilder();
    
    GaiaCacheBuilder& cacheDir(const std::string& dir);
    GaiaCacheBuilder& autoDownload(bool enable);
    GaiaCacheBuilder& verbose(bool enable);
    GaiaCacheBuilder& gaiaClient(std::shared_ptr<GaiaClient> client);
    
    std::unique_ptr<GaiaCache> build();
    
private:
    std::string cacheDir_;
    bool autoDownload_;
    bool verbose_;
    std::shared_ptr<GaiaClient> client_;
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_GAIA_CACHE_H
