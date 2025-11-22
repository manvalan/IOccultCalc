#ifndef IOCCULTCALC_TIME_UTILS_H
#define IOCCULTCALC_TIME_UTILS_H

#include "types.h"
#include <string>
#include <ctime>

namespace ioccultcalc {

class TimeUtils {
public:
    // Converte una data ISO (YYYY-MM-DD o YYYY-MM-DD HH:MM:SS) in Julian Date
    static JulianDate isoToJD(const std::string& isoDate);
    
    // Converte Julian Date in stringa ISO
    static std::string jdToISO(const JulianDate& jd);
    
    // Calcola Julian Date dal calendario gregoriano
    static JulianDate calendarToJD(int year, int month, int day, 
                                    int hour = 0, int minute = 0, double second = 0.0);
    
    // Converte Julian Date in calendario gregoriano
    static void jdToCalendar(const JulianDate& jd, int& year, int& month, int& day,
                            int& hour, int& minute, double& second);
    
    // Ottiene il Julian Date corrente
    static JulianDate now();
    
    // Calcola il Greenwich Mean Sidereal Time
    static double gmst(const JulianDate& jd);
    
    // Calcola il Local Sidereal Time per una data longitudine
    static double lst(const JulianDate& jd, double longitude);
    
    // ========================================================================
    // Time Scale Conversions (TDB, TT, UTC)
    // ========================================================================
    
    /**
     * @brief Converte TT (Terrestrial Time) → TDB (Barycentric Dynamical Time)
     * @param jd_tt Julian Date in scala TT
     * @return Julian Date in scala TDB
     * 
     * TDB è il tempo usato dalle effemeridi JPL (DE441, Horizons).
     * Differenza TDB-TT oscilla ±1.7 ms a causa dell'orbita terrestre.
     */
    static JulianDate ttToTDB(const JulianDate& jd_tt);
    
    /**
     * @brief Converte TDB (Barycentric Dynamical Time) → TT (Terrestrial Time)
     * @param jd_tdb Julian Date in scala TDB
     * @return Julian Date in scala TT
     */
    static JulianDate tdbToTT(const JulianDate& jd_tdb);
    
    /**
     * @brief Converte UTC → TT (Terrestrial Time)
     * @param jd_utc Julian Date in scala UTC
     * @return Julian Date in scala TT
     * 
     * TT = UTC + ΔAT + 32.184 secondi
     * dove ΔAT = TAI - UTC (leap seconds, varia nel tempo)
     */
    static JulianDate utcToTT(const JulianDate& jd_utc);
    
    /**
     * @brief Converte UTC → TDB (conversione diretta per convenienza)
     * @param jd_utc Julian Date in scala UTC
     * @return Julian Date in scala TDB
     */
    static JulianDate utcToTDB(const JulianDate& jd_utc);
    
    /**
     * @brief Stima leap seconds (TAI - UTC) per una data
     * @param jd Julian Date (qualsiasi scala)
     * @return Numero di leap seconds
     * 
     * Nota: Valore esatto richiederebbe tabella IERS.
     * Qui usiamo approssimazione basata su storico leap seconds.
     */
    static int getLeapSeconds(const JulianDate& jd);
};

} // namespace ioccultcalc

#endif // IOCCULTCALC_TIME_UTILS_H
