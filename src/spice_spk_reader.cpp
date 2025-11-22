/**
 * @file spice_spk_reader.cpp
 * @brief Implementazione lettore SPK con CSPICE
 */

#include "ioccultcalc/spice_spk_reader.h"
#include <iostream>
#include <cstring>
#include <cstdlib>

// Include CSPICE
extern "C" {
    #include "SpiceUsr.h"
}

namespace ioccultcalc {

class SPICESPKReader::Impl {
public:
    int handle;
    bool loaded;
    std::string filepath;
    
    Impl() : handle(-1), loaded(false) {
        // Inizializza CSPICE error handling
        erract_c("SET", 0, const_cast<char*>("RETURN"));
    }
    
    ~Impl() {
        close();
    }
    
    bool loadFile(const std::string& path) {
        if (loaded) {
            close();
        }
        
        // Carica frames kernel per ECLIPJ2000_DE441 (IAU 2006 obliquity, DE441-compatible)
        const char* home = getenv("HOME");
        if (home) {
            std::string framesPath = std::string(home) + "/.ioccultcalc/ephemerides/eclipj2000_de441.tf";
            FILE* f = fopen(framesPath.c_str(), "r");
            if (f) {
                fclose(f);
                furnsh_c(framesPath.c_str());
                if (failed_c()) {
                    char msg[1841];
                    getmsg_c("SHORT", 1840, msg);
                    std::cerr << "Warning: Could not load frames kernel: " << msg << std::endl;
                    reset_c();
                }
            }
        }
        
        // Carica SPK kernel
        furnsh_c(path.c_str());
        
        if (failed_c()) {
            char msg[1841];
            getmsg_c("SHORT", 1840, msg);
            std::cerr << "SPICE error loading " << path << ": " << msg << std::endl;
            reset_c();
            return false;
        }
        
        filepath = path;
        loaded = true;
        
        std::cerr << "SPICESPKReader: Loaded " << path << std::endl;
        return true;
    }
    
    std::pair<Vector3D, Vector3D> getState(int bodyId, double jd, int centerId) {
        if (!loaded) {
            throw std::runtime_error("SPK file not loaded");
        }
        
        // Converti JD in ET (Ephemeris Time)
        // JD = 2451545.0 corrisponde a J2000 = 0 ET
        double et = (jd - 2451545.0) * 86400.0;  // secondi da J2000
        
        // Buffer per stato (pos + vel, 6 elementi)
        double state[6];
        double lt;  // light time (non usato)
        
        // SPKEZR: calcola stato di bodyId rispetto a centerId
        char targetStr[64];
        char observerStr[64];
        snprintf(targetStr, sizeof(targetStr), "%d", bodyId);
        snprintf(observerStr, sizeof(observerStr), "%d", centerId);
        
        spkezr_c(targetStr, et, "ECLIPJ2000_DE441", "NONE", observerStr, state, &lt);
        
        if (failed_c()) {
            char msg[1841];
            getmsg_c("SHORT", 1840, msg);
            reset_c();
            throw std::runtime_error(std::string("SPICE error for body ") + 
                                   std::to_string(bodyId) + ": " + msg);
        }
        
        // SPICE restituisce in km e km/s, converti in AU e AU/day
        constexpr double KM_TO_AU = 1.0 / 149597870.7;
        constexpr double KMS_TO_AUD = 86400.0 / 149597870.7;
        
        Vector3D pos(state[0] * KM_TO_AU, state[1] * KM_TO_AU, state[2] * KM_TO_AU);
        Vector3D vel(state[3] * KMS_TO_AUD, state[4] * KMS_TO_AUD, state[5] * KMS_TO_AUD);
        
        return {pos, vel};
    }
    
    void close() {
        if (loaded) {
            // CSPICE richiede unload dei kernel
            unload_c(filepath.c_str());
            loaded = false;
            filepath.clear();
        }
    }
};

SPICESPKReader::SPICESPKReader() : pImpl(std::make_unique<Impl>()) {}

SPICESPKReader::~SPICESPKReader() = default;

bool SPICESPKReader::loadFile(const std::string& filepath) {
    return pImpl->loadFile(filepath);
}

bool SPICESPKReader::ensureFileLoaded(const std::string& name) {
    // Cache directory
    const char* home = getenv("HOME");
    if (!home) {
        std::cerr << "HOME not set\n";
        return false;
    }
    
    std::string cacheDir = std::string(home) + "/.ioccultcalc/ephemerides/";
    std::string filepath = cacheDir + name;
    
    // Verifica se esiste
    FILE* f = fopen(filepath.c_str(), "rb");
    if (!f) {
        std::cerr << "SPK file not found: " << filepath << std::endl;
        std::cerr << "Download it from: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/" << std::endl;
        return false;
    }
    fclose(f);
    
    return loadFile(filepath);
}

Vector3D SPICESPKReader::getPosition(int bodyId, double jd, int centerId) {
    auto [pos, vel] = pImpl->getState(bodyId, jd, centerId);
    return pos;
}

std::pair<Vector3D, Vector3D> SPICESPKReader::getState(int bodyId, double jd, int centerId) {
    return pImpl->getState(bodyId, jd, centerId);
}

bool SPICESPKReader::isLoaded() const {
    return pImpl->loaded;
}

std::vector<int> SPICESPKReader::getAvailableBodies() const {
    // TODO: implementare enumerazione corpi
    return {};
}

void SPICESPKReader::close() {
    pImpl->close();
}

} // namespace ioccultcalc
