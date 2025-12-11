/**
 * @file allnum_database.cpp
 * @brief Implementazione AllnumDatabaseReader per lettura database SQLite
 */

#include "ioccultcalc/allnum_database.h"
#include "ioccultcalc/types.h"
#include <iostream>
#include <cmath>

namespace ioccultcalc {

AllnumDatabaseReader::AllnumDatabaseReader() {
    auto& dm = DataManager::instance();
    std::string default_path = dm.getDatabaseDir() + "/allnum.db";
    initializeDatabase(default_path);
}

AllnumDatabaseReader::AllnumDatabaseReader(const std::string& db_path) {
    initializeDatabase(db_path);
}

AllnumDatabaseReader::~AllnumDatabaseReader() {
    if (db_) {
        sqlite3_close(db_);
        db_ = nullptr;
    }
}

void AllnumDatabaseReader::initializeDatabase(const std::string& db_path) {
    db_path_ = db_path;
    is_ready_ = false;
    db_ = nullptr;
    
    int rc = sqlite3_open(db_path.c_str(), &db_);
    if (rc != SQLITE_OK) {
        // Database non esiste o non accessibile - non è un errore fatale
        if (db_) {
            sqlite3_close(db_);
            db_ = nullptr;
        }
        return;
    }
    
    // Verifica che le tabelle esistano
    const char* check_sql = "SELECT COUNT(*) FROM allnum_asteroids";
    sqlite3_stmt* stmt;
    rc = sqlite3_prepare_v2(db_, check_sql, -1, &stmt, nullptr);
    if (rc == SQLITE_OK) {
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            is_ready_ = true;
        }
        sqlite3_finalize(stmt);
    }
    
    if (!is_ready_ && db_) {
        sqlite3_close(db_);
        db_ = nullptr;
    }
}

std::optional<OrbitalElements> AllnumDatabaseReader::getElement(int number) {
    return getElement(std::to_string(number));
}

std::optional<OrbitalElements> AllnumDatabaseReader::getElement(const std::string& designation) {
    if (!is_ready_ || !db_) {
        return std::nullopt;
    }
    
    // Converti designation a numero se possibile
    int number = 0;
    try {
        number = std::stoi(designation);
    } catch (...) {
        return std::nullopt;
    }
    
    const char* sql = R"(
        SELECT number, designation, epoch_mjd, a, e, i, node_longitude, perihelion_argument, mean_anomaly, H, G
        FROM allnum_asteroids
        WHERE number = ?
        LIMIT 1
    )";
    
    sqlite3_stmt* stmt;
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        return std::nullopt;
    }
    
    sqlite3_bind_int(stmt, 1, number);
    
    std::optional<OrbitalElements> result;
    
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        OrbitalElements elem;
        elem.designation = designation;
        
        // Epoch (MJD) -> JD
        double epoch_mjd = sqlite3_column_double(stmt, 2);
        elem.epoch.jd = epoch_mjd + 2400000.5;
        
        // Elementi orbitali (già in radianti nel database)
        elem.a = sqlite3_column_double(stmt, 3);
        elem.e = sqlite3_column_double(stmt, 4);
        elem.i = sqlite3_column_double(stmt, 5);  // già in radianti
        elem.Omega = sqlite3_column_double(stmt, 6);  // già in radianti
        elem.omega = sqlite3_column_double(stmt, 7);  // già in radianti
        elem.M = sqlite3_column_double(stmt, 8);  // già in radianti
        
        // H e G
        elem.H = sqlite3_column_double(stmt, 9);
        elem.G = sqlite3_column_double(stmt, 10);
        
        result = elem;
    }
    
    sqlite3_finalize(stmt);
    return result;
}

int AllnumDatabaseReader::getRecordCount() const {
    if (!is_ready_ || !db_) {
        return 0;
    }
    
    const char* sql = "SELECT COUNT(*) FROM allnum_asteroids";
    sqlite3_stmt* stmt;
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        return 0;
    }
    
    int count = 0;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(stmt, 0);
    }
    
    sqlite3_finalize(stmt);
    return count;
}

std::string AllnumDatabaseReader::getLastUpdateDate() const {
    if (!is_ready_ || !db_) {
        return "";
    }
    
    const char* sql = R"(
        SELECT download_date 
        FROM allnum_metadata 
        ORDER BY download_date DESC 
        LIMIT 1
    )";
    
    sqlite3_stmt* stmt;
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        return "";
    }
    
    std::string date;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        const char* date_str = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        if (date_str) {
            date = date_str;
        }
    }
    
    sqlite3_finalize(stmt);
    return date;
}

} // namespace ioccultcalc

