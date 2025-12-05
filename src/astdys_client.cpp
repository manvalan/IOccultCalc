#include "ioccultcalc/astdys_client.h"
#include "ioccultcalc/jpl_horizons_client.h"
#include "ioccultcalc/orbit_propagator.h"
#include "ioccultcalc/time_utils.h"
#include <curl/curl.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <regex>
#include <cmath>
#include <algorithm>

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
    int timeout;
    CURL* curl;
    
    Impl() : baseURL("https://newton.spacedys.com/~astdys2/"), timeout(30) {
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

EquinoctialElements AstDysClient::getElements(const std::string& designation) {
    // Costruisci URL per scaricare file .eq
    // Struttura: https://newton.spacedys.com/~astdys2/epoch/numbered/<num/1000>/<num>.eq0
    // Esempio: 433 -> epoch/numbered/0/433.eq0
    //          15080 -> epoch/numbered/15/15080.eq0
    
    int asteroidNumber = 0;
    
    // Se è un numero, calcola la directory
    if (std::all_of(designation.begin(), designation.end(), ::isdigit)) {
        asteroidNumber = std::stoi(designation);
    } else {
        throw std::runtime_error("AstDyS supports only numbered asteroids");
    }
    
    // Calcola il numero della directory (numero / 1000)
    int dirNumber = asteroidNumber / 1000;
    
    // Usa .eq1 che ha epoche RECENTI invece di .eq0 (vecchio)
    std::string url = pImpl->baseURL + "epoch/numbered/" + std::to_string(dirNumber) + "/" + designation + ".eq1";
    
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
    // Formato tipico (varia, questa è una approssimazione):
    // Nome
    // Epoca MJD
    // a h k p q lambda
    // H G
    
    std::istringstream iss(content);
    std::string line;
    int lineNum = 0;
    
    while (std::getline(iss, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '!') continue;
        
        lineNum++;
        
        // Rimuovi spazi iniziali
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos) continue;
        line = line.substr(start);
        
        // Controlla il tipo di linea
        if (line.find("format") == 0 || line.find("rectype") == 0 || 
            line.find("refsys") == 0 || line.find("END_OF_HEADER") == 0) {
            continue; // Header, salta
        }
        
        // Nome dell'asteroide (solo numero o designazione)
        if (elem.name.empty() && line.find_first_not_of("0123456789") == std::string::npos) {
            elem.name = line;
            continue;
        }
        
        // Linea EQU con elementi
        if (line.find("EQU") == 0) {
            std::istringstream iss_line(line.substr(3)); // Salta "EQU"
            iss_line >> elem.a >> elem.h >> elem.k >> elem.p >> elem.q >> elem.lambda;
            continue;
        }
        
        // Linea MJD con epoca
        if (line.find("MJD") == 0) {
            std::istringstream iss_line(line.substr(3)); // Salta "MJD"
            double mjd;
            iss_line >> mjd;
            elem.epoch.jd = mjd + 2400000.5; // Converti MJD in JD
            continue;
        }
        
        // Linea MAG con H e G
        if (line.find("MAG") == 0) {
            std::istringstream iss_line(line.substr(3)); // Salta "MAG"
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

// ============================================================================
// IMPLEMENTAZIONE NUOVE API CHEBYSHEV
// ============================================================================

std::array<double, 3> ChebyshevEphemeris::evaluate(double jd) const {
    if (!isValid(jd)) {
        throw std::runtime_error("JD " + std::to_string(jd) + 
                               " outside Chebyshev validity range [" +
                               std::to_string(startJD) + ", " + 
                               std::to_string(endJD) + "]");
    }
    
    // Normalizza tempo a [-1, +1]
    double t = 2.0 * (jd - startJD) / (endJD - startJD) - 1.0;
    
    // Valuta polinomi di Chebyshev: T_n(t) con ricorrenza
    // T_0 = 1, T_1 = t, T_{n+1} = 2*t*T_n - T_{n-1}
    
    auto evaluatePoly = [](const std::vector<double>& coeff, double t) -> double {
        if (coeff.empty()) return 0.0;
        if (coeff.size() == 1) return coeff[0];
        
        double T_prev = 1.0;         // T_0
        double T_curr = t;           // T_1
        double result = coeff[0] * T_prev + coeff[1] * T_curr;
        
        for (size_t n = 2; n < coeff.size(); ++n) {
            double T_next = 2.0 * t * T_curr - T_prev;
            result += coeff[n] * T_next;
            T_prev = T_curr;
            T_curr = T_next;
        }
        
        return result;
    };
    
    return {
        evaluatePoly(coeff_x, t),
        evaluatePoly(coeff_y, t),
        evaluatePoly(coeff_z, t)
    };
}

std::array<double, 6> ChebyshevEphemeris::evaluateWithVelocity(double jd) const {
    if (!isValid(jd)) {
        throw std::runtime_error("JD outside Chebyshev validity range");
    }
    
    // Normalizza tempo
    double t = 2.0 * (jd - startJD) / (endJD - startJD) - 1.0;
    double dt_dJD = 2.0 / (endJD - startJD);  // Chain rule: dt/dJD
    
    // Valuta posizione E derivata
    auto evaluatePolyWithDeriv = [dt_dJD](const std::vector<double>& coeff, double t) 
        -> std::pair<double, double> {
        if (coeff.empty()) return {0.0, 0.0};
        if (coeff.size() == 1) return {coeff[0], 0.0};
        
        double T_prev = 1.0, T_curr = t;
        double dT_prev = 0.0, dT_curr = 1.0;
        
        double pos = coeff[0] * T_prev + coeff[1] * T_curr;
        double vel = coeff[0] * dT_prev + coeff[1] * dT_curr;
        
        for (size_t n = 2; n < coeff.size(); ++n) {
            double T_next = 2.0 * t * T_curr - T_prev;
            double dT_next = 2.0 * T_curr + 2.0 * t * dT_curr - dT_prev;
            
            pos += coeff[n] * T_next;
            vel += coeff[n] * dT_next;
            
            T_prev = T_curr; T_curr = T_next;
            dT_prev = dT_curr; dT_curr = dT_next;
        }
        
        vel *= dt_dJD;  // Converti dCoord/dt → dCoord/dJD
        return {pos, vel};
    };
    
    auto [x, vx] = evaluatePolyWithDeriv(coeff_x, t);
    auto [y, vy] = evaluatePolyWithDeriv(coeff_y, t);
    auto [z, vz] = evaluatePolyWithDeriv(coeff_z, t);
    
    return {x, y, z, vx, vy, vz};
}

std::array<double, 3> ChebyshevEphemeris::getRADecDist(double jd) const {
    auto [x, y, z] = evaluate(jd);
    
    // Distanza
    double dist = std::sqrt(x*x + y*y + z*z);
    
    // RA, Dec da coordinate cartesiane
    double ra = std::atan2(y, x) * 180.0 / M_PI;
    if (ra < 0.0) ra += 360.0;
    
    double dec = std::asin(z / dist) * 180.0 / M_PI;
    
    return {ra, dec, dist};
}

ChebyshevEphemeris AstDysClient::getChebyshevEphemeris(const std::string& designation,
                                                       double startJD,
                                                       double endJD,
                                                       int order) {
    // URL per Chebyshev da AstDyS:
    // https://newton.spacedys.com/~astdys2/cgi-bin/astdys?objects=<num>;
    //   epochs=<startJD>:<endJD>:<step>;chebyshev=<order>
    
    // Estrai numero asteroide
    std::string asteroidNum = designation;
    // Rimuovi eventuali parentesi: "(433)" → "433"
    if (asteroidNum.find('(') != std::string::npos) {
        size_t start = asteroidNum.find('(');
        size_t end = asteroidNum.find(')');
        if (end != std::string::npos) {
            asteroidNum = asteroidNum.substr(start + 1, end - start - 1);
        }
    }
    
    // Costruisci URL (NOTA: sintassi ipotetica, da verificare con docs AstDyS)
    std::ostringstream url;
    url << pImpl->baseURL << "cgi-bin/astdys?";
    url << "objects=" << asteroidNum << ";";
    url << "epochs=" << std::fixed << std::setprecision(1) 
        << startJD << ":" << endJD << ":1;";
    url << "chebyshev=" << order << ";";
    url << "format=text";
    
    std::cout << "[ASTDYS CHEBYSHEV] Downloading from: " << url.str() << std::endl;
    
    try {
        std::string response = pImpl->httpGet(url.str());
        
        // Parsing della risposta (formato da definire in base a docs AstDyS)
        // Per ora, implementazione placeholder
        ChebyshevEphemeris cheb;
        cheb.startJD = startJD;
        cheb.endJD = endJD;
        cheb.order = order;
        cheb.geocentric = true;
        cheb.frame = "ICRF";
        
        // TODO: Parse response e popola coeff_x, coeff_y, coeff_z
        // Formato atteso (esempio):
        // CHEBYSHEV COEFFICIENTS
        // Order: 15
        // Interval: [2460638.5, 2460645.5]
        // X: 0.123 -0.456 0.789 ...
        // Y: 0.234 -0.567 0.890 ...
        // Z: 0.345 -0.678 0.901 ...
        
        std::cout << "[ASTDYS CHEBYSHEV] Parsing coefficients..." << std::endl;
        std::cout << "[ASTDYS CHEBYSHEV] Response preview:\n" 
                  << response.substr(0, std::min(size_t(500), response.size())) 
                  << std::endl;
        
        // FALLBACK: Se AstDyS non supporta ancora, usa Horizons per generare dati
        std::cout << "[ASTDYS CHEBYSHEV] ⚠️  API Chebyshev non ancora disponibile su AstDyS" << std::endl;
        std::cout << "[ASTDYS CHEBYSHEV] Uso fallback: genero coefficienti da Horizons" << std::endl;
        
        // Genera coefficienti da sample points via Horizons
        JPLHorizonsClient horizons;
        int nSamples = order + 5;  // Sovracampionamento per stabilità
        std::vector<double> samplesX, samplesY, samplesZ;
        
        for (int i = 0; i < nSamples; ++i) {
            double jd = startJD + (endJD - startJD) * i / (nSamples - 1);
            JulianDate epoch;
            epoch.jd = jd;
            
            auto [pos, vel] = horizons.getStateVectors(asteroidNum, epoch, "@geocenter");
            samplesX.push_back(pos.x);
            samplesY.push_back(pos.y);
            samplesZ.push_back(pos.z);
        }
        
        // Fit Chebyshev (metodo dei minimi quadrati)
        auto fitChebyshev = [nSamples, order](const std::vector<double>& samples) {
            std::vector<double> coeffs(order + 1, 0.0);
            
            // Matrice Chebyshev: T_n(t_i) per t_i distribuiti uniformemente
            for (int n = 0; n <= order; ++n) {
                double sum = 0.0;
                double norm = 0.0;
                
                for (int i = 0; i < nSamples; ++i) {
                    double t = 2.0 * i / (nSamples - 1) - 1.0;  // [-1, +1]
                    
                    // Valuta T_n(t)
                    double T_prev = 1.0, T_curr = t;
                    double T_n = (n == 0) ? T_prev : ((n == 1) ? T_curr : 0.0);
                    
                    for (int k = 2; k <= n; ++k) {
                        double T_next = 2.0 * t * T_curr - T_prev;
                        if (k == n) T_n = T_next;
                        T_prev = T_curr;
                        T_curr = T_next;
                    }
                    
                    sum += samples[i] * T_n;
                    norm += T_n * T_n;
                }
                
                coeffs[n] = (norm > 1e-15) ? (sum / norm) : 0.0;
            }
            
            return coeffs;
        };
        
        cheb.coeff_x = fitChebyshev(samplesX);
        cheb.coeff_y = fitChebyshev(samplesY);
        cheb.coeff_z = fitChebyshev(samplesZ);
        
        std::cout << "[ASTDYS CHEBYSHEV] ✓ Generated " << (order + 1) 
                  << " Chebyshev coefficients from " << nSamples << " samples" << std::endl;
        
        return cheb;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to get Chebyshev ephemeris: " + std::string(e.what()));
    }
}

std::vector<ChebyshevEphemeris> AstDysClient::getChebyshevEphemerisMultiSegment(
    const std::string& designation,
    double startJD,
    double endJD,
    double segment_days,
    int order) {
    
    std::vector<ChebyshevEphemeris> segments;
    
    double currentJD = startJD;
    while (currentJD < endJD) {
        double segmentEnd = std::min(currentJD + segment_days, endJD);
        
        std::cout << "[ASTDYS CHEBYSHEV MULTI] Segment [" 
                  << currentJD << ", " << segmentEnd << "]" << std::endl;
        
        auto segment = getChebyshevEphemeris(designation, currentJD, segmentEnd, order);
        segments.push_back(segment);
        
        currentJD = segmentEnd;
    }
    
    std::cout << "[ASTDYS CHEBYSHEV MULTI] Generated " << segments.size() 
              << " segments" << std::endl;
    
    return segments;
}

std::array<double, 3> AstDysClient::getRADecChebyshev(const std::string& designation,
                                                      double jd) {
    // Window di 1 giorno centrato su jd
    double startJD = jd - 0.5;
    double endJD = jd + 0.5;
    
    auto cheb = getChebyshevEphemeris(designation, startJD, endJD, 10);  // Order 10 sufficiente per 1 giorno
    return cheb.getRADecDist(jd);
}

} // namespace ioccultcalc
