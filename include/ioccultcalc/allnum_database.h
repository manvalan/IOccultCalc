#ifndef IOCCULTCALC_ALLNUM_DATABASE_H
#define IOCCULTCALC_ALLNUM_DATABASE_H

#include "orbital_elements.h"
#include "data_manager.h"
#include <sqlite3.h>
#include <string>
#include <memory>
#include <optional>

namespace ioccultcalc {

/**
 * @class AllnumDatabaseReader
 * @brief Reader per database SQLite allnum.db
 * 
 * Legge elementi orbitali kepleriani dal database SQLite locale
 * creato da build_allnum_database.
 * 
 * Usage:
 *   AllnumDatabaseReader db;
 *   auto elem = db.getElement(17030);
 *   if (elem) {
 *       // Usa elementi
 *   }
 */
class AllnumDatabaseReader {
public:
    AllnumDatabaseReader();
    explicit AllnumDatabaseReader(const std::string& db_path);
    ~AllnumDatabaseReader();
    
    /**
     * @brief Carica elementi orbitali per un asteroide
     * @param number Numero asteroide
     * @return Elementi orbitali se trovati, std::nullopt altrimenti
     */
    std::optional<OrbitalElements> getElement(int number);
    
    /**
     * @brief Carica elementi orbitali per designazione
     * @param designation Designazione (es. "17030")
     * @return Elementi orbitali se trovati, std::nullopt altrimenti
     */
    std::optional<OrbitalElements> getElement(const std::string& designation);
    
    /**
     * @brief Verifica se database è disponibile
     */
    bool isAvailable() const { return db_ != nullptr && is_ready_; }
    
    /**
     * @brief Ottieni numero totale di record nel database
     */
    int getRecordCount() const;
    
    /**
     * @brief Ottieni data ultimo aggiornamento
     */
    std::string getLastUpdateDate() const;
    
private:
    sqlite3* db_;
    std::string db_path_;
    bool is_ready_;
    
    void initializeDatabase(const std::string& db_path);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_ALLNUM_DATABASE_H

