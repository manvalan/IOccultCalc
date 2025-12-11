/**
 * @file build_allnum_database.cpp
 * @brief Tool per creare/aggiornare database SQLite da allnum.cat
 * 
 * Features:
 * - Scarica allnum.cat da AstDyS
 * - Parsa formato OEF2.0 con elementi kepleriani
 * - Inserisce in SQLite con tracciamento data
 * - Verifica se database è vecchio (> 30 giorni) e aggiorna
 * 
 * Usage:
 *   ./build_allnum_database [--force] [--db-path PATH] [--max-age DAYS]
 */

#include "ioccultcalc/orbital_elements.h"
#include "ioccultcalc/data_manager.h"
#include "ioccultcalc/types.h"
#include <sqlite3.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <set>
#include <algorithm>

using namespace ioccultcalc;

class AllnumDatabase {
private:
    sqlite3* db_;
    std::string db_path_;
    
    static int callback(void* data, int argc, char** argv, char** azColName) {
        return 0;
    }
    
public:
    AllnumDatabase(const std::string& db_path) : db_(nullptr), db_path_(db_path) {
        int rc = sqlite3_open(db_path.c_str(), &db_);
        if (rc) {
            throw std::runtime_error("Cannot open database: " + db_path + 
                                   " - " + sqlite3_errmsg(db_));
        }
        initializeSchema();
    }
    
    ~AllnumDatabase() {
        if (db_) {
            sqlite3_close(db_);
        }
    }
    
    void initializeSchema() {
        const char* schema = R"(
            CREATE TABLE IF NOT EXISTS allnum_asteroids (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                number INTEGER UNIQUE NOT NULL,
                designation TEXT UNIQUE NOT NULL,
                epoch_mjd REAL NOT NULL,
                a REAL NOT NULL,
                e REAL NOT NULL,
                i REAL NOT NULL,
                node_longitude REAL NOT NULL,
                perihelion_argument REAL NOT NULL,
                mean_anomaly REAL NOT NULL,
                H REAL,
                G REAL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            
            CREATE INDEX IF NOT EXISTS idx_number ON allnum_asteroids(number);
            CREATE INDEX IF NOT EXISTS idx_designation ON allnum_asteroids(designation);
            CREATE INDEX IF NOT EXISTS idx_epoch ON allnum_asteroids(epoch_mjd);
            CREATE INDEX IF NOT EXISTS idx_magnitude ON allnum_asteroids(H);
            
            CREATE TABLE IF NOT EXISTS allnum_metadata (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                download_date TIMESTAMP NOT NULL,
                total_records INTEGER NOT NULL,
                file_url TEXT,
                file_format TEXT,
                data_epoch_mjd REAL,
                last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
        )";
        
        char* err_msg = nullptr;
        int rc = sqlite3_exec(db_, schema, nullptr, nullptr, &err_msg);
        if (rc != SQLITE_OK) {
            std::string error = err_msg ? std::string(err_msg) : "Unknown error";
            sqlite3_free(err_msg);
            throw std::runtime_error("Schema creation failed: " + error);
        }
    }
    
    bool needsUpdate(int max_age_days = 30) {
        const char* sql = R"(
            SELECT download_date, total_records 
            FROM allnum_metadata 
            ORDER BY download_date DESC 
            LIMIT 1
        )";
        
        sqlite3_stmt* stmt;
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            return true; // Se errore, assume che serve update
        }
        
        if (sqlite3_step(stmt) != SQLITE_ROW) {
            sqlite3_finalize(stmt);
            return true; // Nessun record, serve update
        }
        
        // Leggi data download
        const char* date_str = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        int total_records = sqlite3_column_int(stmt, 1);
        
        sqlite3_finalize(stmt);
        
        if (total_records == 0) {
            return true; // Database vuoto
        }
        
        // Parse data (formato: YYYY-MM-DD HH:MM:SS)
        std::tm tm = {};
        std::istringstream ss(date_str);
        ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
        
        if (ss.fail()) {
            return true; // Data invalida, serve update
        }
        
        std::time_t download_time = std::mktime(&tm);
        std::time_t now = std::time(nullptr);
        double days_diff = std::difftime(now, download_time) / (24.0 * 3600.0);
        
        return days_diff > max_age_days;
    }
    
    void clearDatabase() {
        const char* sql = "DELETE FROM allnum_asteroids; DELETE FROM allnum_metadata;";
        char* err_msg = nullptr;
        int rc = sqlite3_exec(db_, sql, nullptr, nullptr, &err_msg);
        if (rc != SQLITE_OK) {
            std::string error = err_msg ? std::string(err_msg) : "Unknown error";
            sqlite3_free(err_msg);
            throw std::runtime_error("Failed to clear database: " + error);
        }
    }
    
    void insertAsteroid(const OrbitalElements& elem, int number) {
        const char* sql = R"(
            INSERT OR REPLACE INTO allnum_asteroids 
            (number, designation, epoch_mjd, a, e, i, node_longitude, perihelion_argument, mean_anomaly, H, G, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
        )";
        
        sqlite3_stmt* stmt;
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            throw std::runtime_error("Failed to prepare statement: " + 
                                   std::string(sqlite3_errmsg(db_)));
        }
        
        sqlite3_bind_int(stmt, 1, number);
        sqlite3_bind_text(stmt, 2, elem.designation.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_double(stmt, 3, elem.epoch.toMJD());
        sqlite3_bind_double(stmt, 4, elem.a);
        sqlite3_bind_double(stmt, 5, elem.e);
        sqlite3_bind_double(stmt, 6, elem.i);  // già in radianti
        sqlite3_bind_double(stmt, 7, elem.Omega);  // già in radianti
        sqlite3_bind_double(stmt, 8, elem.omega);  // già in radianti
        sqlite3_bind_double(stmt, 9, elem.M);  // già in radianti
        sqlite3_bind_double(stmt, 10, elem.H);
        sqlite3_bind_double(stmt, 11, elem.G);
        
        rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            std::string error = sqlite3_errmsg(db_);
            sqlite3_finalize(stmt);
            throw std::runtime_error("Failed to insert asteroid: " + error);
        }
        
        sqlite3_finalize(stmt);
    }
    
    void insertMetadata(int total_records, const std::string& url, 
                       double data_epoch_mjd) {
        const char* sql = R"(
            INSERT INTO allnum_metadata 
            (download_date, total_records, file_url, file_format, data_epoch_mjd)
            VALUES (datetime('now'), ?, ?, 'OEF2.0', ?)
        )";
        
        sqlite3_stmt* stmt;
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            throw std::runtime_error("Failed to prepare metadata statement");
        }
        
        sqlite3_bind_int(stmt, 1, total_records);
        sqlite3_bind_text(stmt, 2, url.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_double(stmt, 3, data_epoch_mjd);
        
        rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE) {
            std::string error = sqlite3_errmsg(db_);
            sqlite3_finalize(stmt);
            throw std::runtime_error("Failed to insert metadata: " + error);
        }
        
        sqlite3_finalize(stmt);
    }
    
    int getRecordCount() {
        const char* sql = "SELECT COUNT(*) FROM allnum_asteroids";
        sqlite3_stmt* stmt;
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            return 0;
        }
        
        if (sqlite3_step(stmt) == SQLITE_ROW) {
            int count = sqlite3_column_int(stmt, 0);
            sqlite3_finalize(stmt);
            return count;
        }
        
        sqlite3_finalize(stmt);
        return 0;
    }
    
    std::set<int> getExistingNumbers() {
        std::set<int> numbers;
        const char* sql = "SELECT number FROM allnum_asteroids";
        sqlite3_stmt* stmt;
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            return numbers;
        }
        
        while (sqlite3_step(stmt) == SQLITE_ROW) {
            int number = sqlite3_column_int(stmt, 0);
            numbers.insert(number);
        }
        
        sqlite3_finalize(stmt);
        return numbers;
    }
};

// Parse single line from allnum.cat (OEF2.0 format)
OrbitalElements parseAllnumLine(const std::string& line, int& number) {
    if (line.length() < 190) {
        throw std::runtime_error("Line too short");
    }
    
    // Extract number from quotes
    size_t firstQuote = line.find('\'');
    size_t secondQuote = line.find('\'', firstQuote + 1);
    if (secondQuote == std::string::npos) {
        throw std::runtime_error("Invalid line format");
    }
    
    std::string number_str = line.substr(firstQuote + 1, secondQuote - firstQuote - 1);
    number = std::stoi(number_str);
    
    OrbitalElements elem;
    elem.designation = number_str;
    
    // Parse using fixed-width positions (OEF2.0 format)
    // Epoch (MJD) - positions 15-27
    std::string mjd_str = line.substr(15, 13);
    double mjd = std::stod(mjd_str);
    elem.epoch.jd = mjd + 2400000.5;
    
    // a (AU) - positions 30-52 (scientific format)
    std::string a_str = line.substr(30, 23);
    elem.a = std::stod(a_str);
    
    // e - positions 55-77 (scientific format)
    std::string e_str = line.substr(55, 23);
    elem.e = std::stod(e_str);
    
    // i (inclination) - positions 80-102 (scientific format, in DEGREES)
    std::string i_str = line.substr(80, 23);
    double i_deg = std::stod(i_str);
    elem.i = i_deg * M_PI / 180.0; // Convert to radians
    
    // Omega - positions 105-127 (scientific format, in DEGREES)
    std::string Omega_str = line.substr(105, 23);
    double Omega_deg = std::stod(Omega_str);
    elem.Omega = Omega_deg * M_PI / 180.0; // Convert to radians
    
    // omega - positions 130-152 (scientific format, in DEGREES)
    std::string omega_str = line.substr(130, 23);
    double omega_deg = std::stod(omega_str);
    elem.omega = omega_deg * M_PI / 180.0; // Convert to radians
    
    // M (mean anomaly) - positions 155-177 (scientific format, in DEGREES)
    std::string M_str = line.substr(155, 23);
    double M_deg = std::stod(M_str);
    elem.M = M_deg * M_PI / 180.0; // Convert to radians
    
    // H - positions 178-183
    std::string H_str = line.substr(178, 6);
    if (!H_str.empty() && H_str.find_first_not_of(" \t") != std::string::npos) {
        elem.H = std::stod(H_str);
    } else {
        elem.H = 15.0; // default
    }
    
    // G - positions 185-189
    std::string G_str = line.substr(185, 5);
    if (!G_str.empty() && G_str.find_first_not_of(" \t") != std::string::npos) {
        elem.G = std::stod(G_str);
    } else {
        elem.G = 0.15; // default
    }
    
    return elem;
}

double parseAllnumFile(const std::string& filepath, AllnumDatabase& db) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    std::string line;
    bool foundHeader = false;
    int count = 0;
    double data_epoch_mjd = 0.0;
    
    std::cout << "Parsing allnum.cat...\n";
    
    while (std::getline(file, line)) {
        if (line.find("END_OF_HEADER") != std::string::npos) {
            foundHeader = true;
            continue;
        }
        if (!foundHeader) continue;
        
        // Skip comment lines
        if (line.empty() || line[0] == '!') continue;
        
        // Parse line: 'number' epoch a e i Omega omega M H G flag
        if (line.find("'") == 0) {
            try {
                int number;
                OrbitalElements elem = parseAllnumLine(line, number);
                
                // Store epoch for metadata (all records have same epoch)
                if (count == 0) {
                    data_epoch_mjd = elem.epoch.toMJD();
                }
                
                db.insertAsteroid(elem, number);
                count++;
                
                if (count % 10000 == 0) {
                    std::cout << "  Processed " << count << " asteroids...\n";
                }
            } catch (const std::exception& e) {
                // Skip invalid lines
                continue;
            }
        }
    }
    
    std::cout << "✓ Parsed " << count << " asteroids\n";
    return data_epoch_mjd;
}

// Parse single line from MPCORB.DAT format
// Returns true if parsed successfully, false otherwise
bool parseMPCORBLine(const std::string& line, OrbitalElements& elem, int& number) {
    if (line.length() < 103) {
        return false;
    }
    
    // Col 1-7: Designation/Number
    std::string desig = line.substr(0, 7);
    // Remove leading/trailing spaces
    size_t start = desig.find_first_not_of(" \t");
    if (start == std::string::npos) return false;
    size_t end = desig.find_last_not_of(" \t");
    desig = desig.substr(start, end - start + 1);
    
    // Check if it's a numbered asteroid (starts with digit)
    if (desig.empty() || !std::isdigit(desig[0])) {
        return false; // Skip unnumbered asteroids
    }
    
    try {
        number = std::stoi(desig);
        elem.designation = desig;
    } catch (...) {
        return false;
    }
    
    // H magnitude - col 9-13 (0-indexed: 8-13)
    std::string H_str = line.substr(8, 5);
    elem.H = H_str.empty() ? 15.0 : std::stod(H_str);
    
    // G slope parameter - col 15-19 (0-indexed: 14-19)
    std::string G_str = line.substr(14, 5);
    elem.G = G_str.empty() ? 0.15 : std::stod(G_str);
    
    // Epoch (packed MPC format) - col 21-25 (0-indexed: 20-25)
    // Format: e.g., "K257B" = 2025-11-02
    // For now, use approximate epoch (recent MPCORB files use ~2460000 JD)
    // TODO: Implement proper packed date decoder
    std::string epoch_packed = line.substr(20, 5);
    elem.epoch.jd = 2460000.0; // Approximate for recent epochs
    
    // Mean anomaly (M) - col 27-35 (0-indexed: 26-35) - IN DEGREES
    std::string M_str = line.substr(26, 9);
    double M_deg = std::stod(M_str);
    elem.M = M_deg * M_PI / 180.0; // Convert to radians
    
    // Argument of perihelion (omega) - col 38-46 (0-indexed: 37-46) - IN DEGREES
    std::string omega_str = line.substr(37, 9);
    double omega_deg = std::stod(omega_str);
    elem.omega = omega_deg * M_PI / 180.0; // Convert to radians
    
    // Longitude of ascending node (Omega) - col 49-57 (0-indexed: 48-57) - IN DEGREES
    std::string Omega_str = line.substr(48, 9);
    double Omega_deg = std::stod(Omega_str);
    elem.Omega = Omega_deg * M_PI / 180.0; // Convert to radians
    
    // Inclination (i) - col 60-68 (0-indexed: 59-68) - IN DEGREES
    std::string i_str = line.substr(59, 9);
    double i_deg = std::stod(i_str);
    elem.i = i_deg * M_PI / 180.0; // Convert to radians
    
    // Eccentricity (e) - col 71-79 (0-indexed: 70-79)
    std::string e_str = line.substr(70, 9);
    elem.e = std::stod(e_str);
    
    // Semimajor axis (a) - col 81-91 (0-indexed: 80-91) - AU
    std::string a_str = line.substr(80, 11);
    elem.a = std::stod(a_str);
    
    return true;
}

void parseMPCORBFile(const std::string& filepath, AllnumDatabase& db, 
                    const std::set<int>& existing_numbers) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open MPCORB file: " + filepath);
    }
    
    std::string line;
    int count = 0;
    int added_count = 0;
    
    std::cout << "Parsing MPCORB.DAT for missing asteroids...\n";
    std::cout << "  Existing asteroids in database: " << existing_numbers.size() << "\n";
    
    while (std::getline(file, line)) {
        // Skip header lines (usually start with special characters or are short)
        if (line.length() < 103 || line[0] == ' ' && line.find("MPCORB") != std::string::npos) {
            continue;
        }
        
        try {
            int number;
            OrbitalElements elem;
            
            if (parseMPCORBLine(line, elem, number)) {
                count++;
                
                // Only add if not already in database
                if (existing_numbers.find(number) == existing_numbers.end()) {
                    db.insertAsteroid(elem, number);
                    added_count++;
                    
                    if (added_count % 1000 == 0) {
                        std::cout << "  Added " << added_count << " asteroids from MPC...\n";
                    }
                }
            }
        } catch (const std::exception&) {
            // Skip invalid lines
            continue;
        }
    }
    
    std::cout << "✓ Parsed " << count << " asteroids from MPCORB\n";
    std::cout << "✓ Added " << added_count << " missing asteroids to database\n";
}

int main(int argc, char* argv[]) {
    bool force = false;
    std::string db_path;
    int max_age_days = 30;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--force") {
            force = true;
        } else if (arg == "--db-path" && i + 1 < argc) {
            db_path = argv[++i];
        } else if (arg == "--max-age" && i + 1 < argc) {
            max_age_days = std::stoi(argv[++i]);
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  --force          Force update even if database is recent\n"
                      << "  --db-path PATH   Specify database path\n"
                      << "  --max-age DAYS   Maximum age in days before update (default: 30)\n"
                      << "  --help           Show this help\n";
            return 0;
        }
    }
    
    // Default database path
    if (db_path.empty()) {
        auto& dm = DataManager::instance();
        db_path = dm.getDatabaseDir() + "/allnum.db";
    }
    
    try {
        std::cout << "📊 Allnum Database Builder\n";
        std::cout << "Database path: " << db_path << "\n\n";
        
        AllnumDatabase db(db_path);
        
        // Check if update needed
        if (!force && !db.needsUpdate(max_age_days)) {
            int count = db.getRecordCount();
            std::cout << "✓ Database is up-to-date (" << count << " records)\n";
            std::cout << "  Use --force to update anyway\n";
            return 0;
        }
        
        std::cout << "Downloading allnum.cat from AstDyS...\n";
        std::string url = "https://newton.spacedys.com/~astdys2/catalogs/allnum.cat";
        std::string temp_file = "/tmp/allnum.cat";
        
        // Download using curl
        std::string cmd = "curl -s --max-time 300 -o " + temp_file + " \"" + url + "\"";
        int ret = system(cmd.c_str());
        if (ret != 0) {
            throw std::runtime_error("Failed to download allnum.cat from " + url);
        }
        
        // Check if file was downloaded
        std::ifstream check_file(temp_file);
        if (!check_file.good()) {
            throw std::runtime_error("Downloaded file is empty or invalid");
        }
        check_file.close();
        
        std::cout << "✓ Downloaded allnum.cat\n";
        std::cout << "Building database from allnum.cat...\n";
        
        // Clear old data
        db.clearDatabase();
        
        // Parse and insert from allnum.cat
        double data_epoch_mjd = parseAllnumFile(temp_file, db);
        
        // Get existing numbers to find missing ones
        std::set<int> existing_numbers = db.getExistingNumbers();
        int allnum_count = existing_numbers.size();
        std::cout << "\n✓ Loaded " << allnum_count << " asteroids from allnum.cat\n";
        
        // Download MPCORB.DAT to complete missing asteroids
        std::cout << "\nDownloading MPCORB.DAT from MPC to complete database...\n";
        std::string mpc_url = "https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz";
        std::string mpc_gz_file = "/tmp/MPCORB.DAT.gz";
        std::string mpc_file = "/tmp/MPCORB.DAT";
        
        // Download MPCORB.DAT.gz
        std::string mpc_cmd = "curl -s --max-time 600 -o " + mpc_gz_file + " \"" + mpc_url + "\"";
        ret = system(mpc_cmd.c_str());
        if (ret != 0) {
            std::cout << "⚠ Warning: Failed to download MPCORB.DAT.gz, skipping MPC completion\n";
        } else {
            // Decompress
            std::string gunzip_cmd = "gunzip -f " + mpc_gz_file;
            ret = system(gunzip_cmd.c_str());
            if (ret != 0) {
                std::cout << "⚠ Warning: Failed to decompress MPCORB.DAT.gz\n";
            } else {
                // Parse MPCORB and add missing asteroids
                parseMPCORBFile(mpc_file, db, existing_numbers);
                
                // Cleanup MPC files
                std::remove(mpc_file.c_str());
            }
        }
        
        // Insert metadata
        int total_records = db.getRecordCount();
        db.insertMetadata(total_records, url, data_epoch_mjd);
        
        std::cout << "\n✅ Database built successfully!\n";
        std::cout << "   From allnum.cat: " << allnum_count << " asteroids\n";
        std::cout << "   From MPCORB.DAT: " << (total_records - allnum_count) << " asteroids\n";
        std::cout << "   Total records: " << total_records << "\n";
        std::cout << "   Database: " << db_path << "\n";
        
        // Cleanup
        std::remove(temp_file.c_str());
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

