#ifndef IOCCULTCALC_IOC_SPICE_WRAPPER_H
#define IOCCULTCALC_IOC_SPICE_WRAPPER_H

#include <string>
#include <vector>
#include <memory>
#include <stdexcept>

namespace ioccultcalc {

/**
 * @brief Eccezione per errori SPICE
 */
class SpiceException : public std::runtime_error {
public:
    explicit SpiceException(const std::string& message)
        : std::runtime_error(message) {}
};

/**
 * @brief Struttura per posizione e velocità
 */
struct StateVector {
    double x, y, z;     // Posizione (km)
    double vx, vy, vz;  // Velocità (km/s)
    
    StateVector() : x(0), y(0), z(0), vx(0), vy(0), vz(0) {}
    StateVector(double x_, double y_, double z_, double vx_, double vy_, double vz_)
        : x(x_), y(y_), z(z_), vx(vx_), vy(vy_), vz(vz_) {}
};

/**
 * @brief Wrapper C++ thread-safe per IOC_SPICE
 * 
 * Fornisce un'interfaccia C++ moderna per le funzionalità SPICE
 * con supporto OpenMP per calcoli paralleli ad alte prestazioni.
 */
class IOCSpiceWrapper {
public:
    /**
     * @brief Costruttore
     * @param num_threads Numero di thread OpenMP (0 = auto-detect)
     */
    explicit IOCSpiceWrapper(int num_threads = 0);
    
    /**
     * @brief Distruttore - scarica tutti i kernel
     */
    ~IOCSpiceWrapper();
    
    // Non copiabile
    IOCSpiceWrapper(const IOCSpiceWrapper&) = delete;
    IOCSpiceWrapper& operator=(const IOCSpiceWrapper&) = delete;
    
    /**
     * @brief Carica un kernel SPICE
     * @param kernel_path Path al file kernel (.bsp, .tpc, .tls, etc.)
     * @throws SpiceException se il caricamento fallisce
     */
    void loadKernel(const std::string& kernel_path);
    
    /**
     * @brief Carica multipli kernel
     * @param kernel_paths Vettore di path ai kernel
     */
    void loadKernels(const std::vector<std::string>& kernel_paths);
    
    /**
     * @brief Calcola lo stato di un corpo per un singolo tempo
     * @param target Nome del corpo (es. "CERES", "2000001" per Ceres)
     * @param et Tempo effemeridi (secondi dal J2000)
     * @param observer Osservatore (default: "EARTH")
     * @param frame Frame di riferimento (default: "J2000")
     * @param abcorr Correzione aberrazione (default: "LT+S")
     * @return StateVector con posizione e velocità
     * @throws SpiceException se il calcolo fallisce
     */
    StateVector getState(
        const std::string& target,
        double et,
        const std::string& observer = "EARTH",
        const std::string& frame = "J2000",
        const std::string& abcorr = "LT+S"
    );
    
    /**
     * @brief Calcola lo stato di un corpo per multipli tempi (parallelo)
     * @param target Nome del corpo
     * @param et_times Vettore di tempi effemeridi
     * @param observer Osservatore
     * @param frame Frame di riferimento
     * @param abcorr Correzione aberrazione
     * @return Vettore di StateVector
     * @throws SpiceException se il calcolo fallisce
     */
    std::vector<StateVector> getMultiState(
        const std::string& target,
        const std::vector<double>& et_times,
        const std::string& observer = "EARTH",
        const std::string& frame = "J2000",
        const std::string& abcorr = "LT+S"
    );
    
    /**
     * @brief Calcola gli stati di multipli corpi per un singolo tempo (parallelo)
     * @param targets Vettore di nomi dei corpi
     * @param et Tempo effemeridi
     * @param observer Osservatore
     * @param frame Frame di riferimento
     * @param abcorr Correzione aberrazione
     * @return Vettore di StateVector (uno per corpo)
     * @throws SpiceException se il calcolo fallisce
     */
    std::vector<StateVector> getMultiBodyState(
        const std::vector<std::string>& targets,
        double et,
        const std::string& observer = "EARTH",
        const std::string& frame = "J2000",
        const std::string& abcorr = "LT+S"
    );
    
    /**
     * @brief Converte data UTC in tempo effemeridi
     * @param utc_str Stringa UTC (es. "2025-12-01T00:00:00")
     * @return Tempo effemeridi in secondi dal J2000
     * @throws SpiceException se la conversione fallisce
     */
    double utcToET(const std::string& utc_str);
    
    /**
     * @brief Converte tempo effemeridi in data UTC
     * @param et Tempo effemeridi
     * @return Stringa UTC
     */
    std::string etToUTC(double et);
    
    /**
     * @brief Ottiene il numero di thread configurati
     */
    int getNumThreads() const { return num_threads_; }
    
private:
    int num_threads_;
    bool initialized_;
    
    void checkError(int status, const std::string& operation);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_IOC_SPICE_WRAPPER_H
