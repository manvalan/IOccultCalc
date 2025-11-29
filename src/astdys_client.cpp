#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/time_utils.h"
#include <curl/curl.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <regex>

namespace ioccultcalc {

// Callback per ricevere dati da libcurl
static size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* output) {
    size_t totalSize = size * nmemb;
    output->append(static_cast<char*>(contents), totalSize);
    return totalSize;
}

class AstDysClient::Impl {
public:
    std::string baseURL;
    std::string localEQ1Dir;
    std::string localRWODir;
    int timeout;
    CURL* curl;
    
    Impl() : baseURL("https://newton.spacedys.com/~astdys2/"), 
             localEQ1Dir(""), 
             localRWODir(""),
             timeout(30) {
        curl_global_init(CURL_GLOBAL_DEFAULT);
        curl = curl_easy_init();
    }
    
    ~Impl() {
        if (curl) {
            curl_easy_cleanup(curl);
        }
        curl_global_cleanup();
    }
    
    std::string httpGet(const std::string& url) {
        if (!curl) {
            throw std::runtime_error("CURL not initialized");
        }
        
        std::string response;
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
        curl_easy_setopt(curl, CURLOPT_TIMEOUT, timeout);
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        
        CURLcode res = curl_easy_perform(curl);
        
        if (res != CURLE_OK) {
            throw std::runtime_error(std::string("HTTP request failed: ") + 
                                   curl_easy_strerror(res));
        }
        
        long responseCode;
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &responseCode);
        
        if (responseCode != 200) {
            throw std::runtime_error("HTTP error code: " + std::to_string(responseCode));
        }
        
        return response;
    }
};

AstDysClient::AstDysClient() : pImpl(new Impl()) {}

AstDysClient::~AstDysClient() = default;

void AstDysClient::setBaseURL(const std::string& url) {
    pImpl->baseURL = url;
    if (pImpl->baseURL.back() != '/') {
        pImpl->baseURL += '/';
    }
}

void AstDysClient::setTimeout(int seconds) {
    pImpl->timeout = seconds;
}

void AstDysClient::setLocalEQ1Directory(const std::string& directory) {
    pImpl->localEQ1Dir = directory;
}

void AstDysClient::setLocalRWODirectory(const std::string& directory) {
    pImpl->localRWODir = directory;
}

EquinoctialElements AstDysClient::getElementsFromFile(const std::string& filepath) {
    // Leggi file locale
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath);
    }
    
    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    file.close();
    
    // Estrai designation dal nome file (es: "433.eq1" -> "433")
    size_t lastSlash = filepath.find_last_of("/\\");
    size_t lastDot = filepath.find_last_of('.');
    std::string designation = filepath.substr(
        lastSlash == std::string::npos ? 0 : lastSlash + 1,
        lastDot == std::string::npos ? std::string::npos : lastDot - (lastSlash == std::string::npos ? 0 : lastSlash + 1)
    );
    
    return parseEquinoctialFile(content, designation);
}

std::vector<std::string> AstDysClient::getObservationsFromFile(const std::string& filepath) {
    // Leggi file .rwo locale
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open .rwo file: " + filepath);
    }
    
    std::vector<std::string> observations;
    std::string line;
    bool inDataSection = false;
    
    while (std::getline(file, line)) {
        // Salta header fino a "END_OF_HEADER"
        if (line.find("END_OF_HEADER") != std::string::npos) {
            inDataSection = true;
            continue;
        }
        
        if (!inDataSection) continue;
        
        // Salta righe vuote e commenti
        if (line.empty() || line[0] == '!' || line[0] == '#') continue;
        
        // Aggiungi riga osservazione
        observations.push_back(line);
    }
    
    file.close();
    return observations;
}

std::vector<std::string> AstDysClient::getObservations(const std::string& designation) {
    // Verifica numero asteroide
    int asteroidNumber = 0;
    if (std::all_of(designation.begin(), designation.end(), ::isdigit)) {
        asteroidNumber = std::stoi(designation);
    } else {
        throw std::runtime_error("AstDyS supports only numbered asteroids");
    }
    
    // Calcola subdirectory
    int dirNumber = asteroidNumber / 1000;
    
    // SE directory locale configurata, prova file locale
    if (!pImpl->localRWODir.empty()) {
        std::string localPath = pImpl->localRWODir + "/" + designation + ".rwo";
        
        std::ifstream file(localPath);
        if (file.good()) {
            std::cout << "📁 Carico osservazioni da file locale: " << localPath << std::endl;
            file.close();
            return getObservationsFromFile(localPath);
        } else {
            std::cout << "⚠️  File .rwo locale non trovato: " << localPath << std::endl;
            std::cout << "   Provo download da AstDyS..." << std::endl;
        }
    }
    
    // ALTRIMENTI scarica da AstDyS
    std::string url = pImpl->baseURL + "mpcobs/numbered/" + std::to_string(dirNumber) + "/" + designation + ".rwo";
    
    std::cout << "🌐 Download osservazioni da AstDyS: " << url << std::endl;
    std::string content = pImpl->httpGet(url);
    
    // Parse contenuto
    std::vector<std::string> observations;
    std::istringstream stream(content);
    std::string line;
    bool inDataSection = false;
    
    while (std::getline(stream, line)) {
        if (line.find("END_OF_HEADER") != std::string::npos) {
            inDataSection = true;
            continue;
        }
        if (!inDataSection) continue;
        if (line.empty() || line[0] == '!' || line[0] == '#') continue;
        observations.push_back(line);
    }
    
    return observations;
}

EquinoctialElements AstDysClient::getElements(const std::string& designation) {
    // Verifica numero asteroide
    int asteroidNumber = 0;
    if (std::all_of(designation.begin(), designation.end(), ::isdigit)) {
        asteroidNumber = std::stoi(designation);
    } else {
        throw std::runtime_error("AstDyS supports only numbered asteroids");
    }
    
    // Calcola il numero della directory (numero / 1000)
    int dirNumber = asteroidNumber / 1000;
    
    // SE directory locale configurata, prova a caricare da file
    if (!pImpl->localEQ1Dir.empty()) {
        std::string localPath = pImpl->localEQ1Dir + "/" + designation + ".eq1";
        
        std::ifstream file(localPath);
        if (file.good()) {
            std::cout << "📁 Carico elementi da file locale: " << localPath << std::endl;
            file.close();
            return getElementsFromFile(localPath);
        } else {
            std::cout << "⚠️  File locale non trovato: " << localPath << std::endl;
            std::cout << "   Provo download da AstDyS..." << std::endl;
        }
    }
    
    // ALTRIMENTI scarica da AstDyS (HTTP)
    // Struttura: https://newton.spacedys.com/~astdys2/epoch/numbered/<num/1000>/<num>.eq1
    // Esempio: 433 -> epoch/numbered/0/433.eq1
    //          11234 -> epoch/numbered/11/11234.eq1
    std::string url = pImpl->baseURL + "epoch/numbered/" + std::to_string(dirNumber) + "/" + designation + ".eq1";
    
    std::cout << "🌐 Download da AstDyS: " << url << std::endl;
    std::string content = pImpl->httpGet(url);
    return parseEquinoctialFile(content, designation);
}

std::vector<EquinoctialElements> AstDysClient::getElementsBatch(
    const std::vector<std::string>& designations) {
    
    std::vector<EquinoctialElements> results;
    results.reserve(designations.size());
    
    for (const auto& desig : designations) {
        try {
            results.push_back(getElements(desig));
        } catch (const std::exception& e) {
            // Log error e continua
            // In una implementazione reale, si potrebbe usare un sistema di logging
        }
    }
    
    return results;
}

std::vector<std::string> AstDysClient::searchByName(const std::string& name) {
    // Query al servizio di ricerca di AstDyS
    std::string url = pImpl->baseURL + "search?name=" + name;
    
    std::string content = pImpl->httpGet(url);
    
    // Parsing del risultato (dipende dal formato della risposta)
    // Questo è un esempio semplificato
    std::vector<std::string> results;
    std::istringstream iss(content);
    std::string line;
    
    while (std::getline(iss, line)) {
        if (!line.empty()) {
            results.push_back(line);
        }
    }
    
    return results;
}

EquinoctialElements AstDysClient::parseEquinoctialFile(const std::string& content,
                                                       const std::string& designation) {
    EquinoctialElements elem;
    elem.designation = designation;
    
    // Parsing del formato .eq di AstDyS
    // NOTA: I valori possono essere spezzati su più righe!
    // Formato: "EQU a h k p q_part1\nq_part2 lambda ... MJD epoch TDT"
    // Cerchiamo la linea EQU e la successiva, le combiniamo, poi estraiamo i numeri
    
    std::istringstream iss(content);
    std::string line;
    std::string equBlock;
    bool inEquBlock = false;
    
    while (std::getline(iss, line)) {
        if (line.empty() || line[0] == '!') continue;
        
        // Rimuovi spazi iniziali
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        line = line.substr(start);
        
        // Header, salta
        if (line.find("format") == 0 || line.find("rectype") == 0 || 
            line.find("refsys") == 0 || line.find("END_OF_HEADER") == 0) {
            continue;
        }
        
        // Nome dell'asteroide
        if (elem.name.empty() && line.find_first_not_of("0123456789") == std::string::npos) {
            elem.name = line;
            continue;
        }
        
        // Inizia blocco EQU
        if (line.find("EQU") == 0) {
            equBlock = line;
            inEquBlock = true;
            continue;
        }
        
        // Continua blocco EQU (riga successiva con resto dei dati e MJD)
        if (inEquBlock && line.find("MAG") != 0 && line.find("LSP") != 0 && 
            line.find("COV") != 0 && line.find("NOR") != 0 && line.find("RMS") != 0 &&
            line.find("EIG") != 0 && line.find("WEA") != 0) {
            equBlock += " " + line;
            
            // Se contiene MJD, il blocco EQU è completo
            if (line.find("MJD") != std::string::npos) {
                inEquBlock = false;
                
                // Parsing del blocco EQU completo
                // Formato: "EQU a h k p q lambda ... MJD epoch TDT"
                size_t mjdPos = equBlock.find("MJD");
                std::string equData = equBlock.substr(3, mjdPos - 3); // Dopo "EQU", prima di "MJD"
                std::string mjdData = equBlock.substr(mjdPos + 3);
                
                // Leggi i 6 elementi equinoziali
                std::istringstream iss_equ(equData);
                iss_equ >> elem.a >> elem.h >> elem.k >> elem.p >> elem.q >> elem.lambda;
                
                // Lambda è in GRADI nel formato AstDyS, converti in radianti
                elem.lambda = elem.lambda * M_PI / 180.0;
                
                // Normalizza mean longitude a [0, 2π]
                while (elem.lambda < 0) elem.lambda += 2.0 * M_PI;
                while (elem.lambda >= 2.0 * M_PI) elem.lambda -= 2.0 * M_PI;
                
                // Leggi epoca MJD
                std::istringstream iss_mjd(mjdData);
                double mjd;
                iss_mjd >> mjd;
                elem.epoch.jd = mjd + 2400000.5;
            }
            continue;
        }
        
        // Linea MAG con H e G
        if (line.find("MAG") == 0) {
            std::istringstream iss_line(line.substr(3));
            iss_line >> elem.H >> elem.G;
            continue;
        }
    }
    
    return elem;
}

OrbitalElements AstDysClient::getRecentElements(const std::string& designation) {
    // Scarica il catalogo completo con epoche recenti
    std::string url = pImpl->baseURL + "catalogs/allnum.cat";
    std::string catalogContent = pImpl->httpGet(url);
    
    // Cerca la linea dell'asteroide
    // Formato: '704'  61000.000000  a e i node argperi M H G flag
    std::string searchPattern = "'" + designation + "'";
    
    std::istringstream stream(catalogContent);
    std::string line;
    bool foundHeader = false;
    
    while (std::getline(stream, line)) {
        // Salta header
        if (line.find("END_OF_HEADER") != std::string::npos) {
            foundHeader = true;
            continue;
        }
        
        if (!foundHeader) continue;
        
        // Cerca la linea con l'asteroide
        if (line.find(searchPattern) == 0) {
            // Parsing: 'num' epoch a e i node argperi M H G flag
            OrbitalElements elem;
            elem.designation = designation;
            
            // Rimuovi i quote dal numero
            size_t firstQuote = line.find('\'');
            size_t secondQuote = line.find('\'', firstQuote + 1);
            std::string dataStr = line.substr(secondQuote + 1);
            
            std::istringstream iss(dataStr);
            double mjd;
            iss >> mjd >> elem.a >> elem.e >> elem.i >> elem.Omega >> elem.omega >> elem.M;
            iss >> elem.H >> elem.G;
            
            // Converti MJD in JD
            elem.epoch.jd = mjd + 2400000.5;
            
            // Converti angoli da gradi a radianti
            // NOTA: Elementi sono in frame ECLM J2000
            elem.i *= M_PI / 180.0;
            elem.Omega *= M_PI / 180.0;
            elem.omega *= M_PI / 180.0;
            elem.M *= M_PI / 180.0;
            
            return elem;
        }
    }
    
    throw std::runtime_error("Asteroid " + designation + " not found in catalog");
}

OrbitalElements AstDysClient::getOsculatingElements(const std::string& designation,
                                                    const JulianDate& epoch) {
    // AstDyS fornisce solo elementi MEDI che richiedono propagazione OrbFit.
    // Per elementi OSCULANTI (istantanei, direttamente convertibili in posizione/velocità)
    // usiamo JPL Horizons.
    
    std::cout << "AstDysClient: Richiesta elementi osculanti da JPL Horizons per " 
              << designation << " all'epoca JD " << epoch.jd << std::endl;
    
    JPLHorizonsClient horizons;
    
    // Horizons richiede target ID: numero per asteroidi numerati
    std::string targetId = designation;
    
    // Se designation è un numero, usa direttamente
    // Altrimenti, prova a estrarre il numero se è del tipo "(704)"
    if (designation.find('(') != std::string::npos) {
        size_t start = designation.find('(');
        size_t end = designation.find(')');
        if (end != std::string::npos) {
            targetId = designation.substr(start + 1, end - start - 1);
        }
    }
    
    std::cout << "  Target ID per Horizons: " << targetId << std::endl;
    
    try {
        OrbitalElements elem = horizons.getOsculatingElements(targetId, epoch, "@sun");
        elem.designation = designation;
        
        std::cout << "  Ricevuti elementi osculanti:" << std::endl;
        std::cout << "    a = " << elem.a << " AU" << std::endl;
        std::cout << "    e = " << elem.e << std::endl;
        std::cout << "    i = " << elem.i * 180.0 / M_PI << "°" << std::endl;
        
        return elem;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Impossibile ottenere elementi osculanti da Horizons: " + 
                                std::string(e.what()));
    }
}

OrbitState AstDysClient::getStateFromHorizons(const std::string& designation,
                                              const JulianDate& epoch) {
    // Metodo PREFERITO: usa vettori di stato direttamente da Horizons
    // Più veloce, più affidabile, nessun problema di parsing elementi
    
    std::cout << "AstDysClient: Richiesta vettori di stato da JPL Horizons per " 
              << designation << " all'epoca MJD " << epoch.toMJD() << std::endl;
    
    JPLHorizonsClient horizons;
    
    // Estrai numero asteroide
    std::string targetId = designation;
    if (designation.find('(') != std::string::npos) {
        size_t start = designation.find('(');
        size_t end = designation.find(')');
        if (end != std::string::npos) {
            targetId = designation.substr(start + 1, end - start - 1);
        }
    }
    
    std::cout << "  Target ID per Horizons: " << targetId << std::endl;
    
    try {
        auto [position, velocity] = horizons.getStateVectors(targetId, epoch, "@sun");
        
        OrbitState state(epoch, position, velocity);
        
        double r = position.magnitude();
        double ra = atan2(position.y, position.x) * 180.0 / M_PI;
        if (ra < 0) ra += 360.0;
        double dec = asin(position.z / r) * 180.0 / M_PI;
        
        std::cout << "  Ricevuto stato orbitale:" << std::endl;
        std::cout << "    r = " << r << " AU" << std::endl;
        std::cout << "    RA = " << ra << "°, Dec = " << dec << "°" << std::endl;
        
        return state;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Impossibile ottenere vettori da Horizons: " + 
                                std::string(e.what()));
    }
}

} // namespace ioccultcalc
