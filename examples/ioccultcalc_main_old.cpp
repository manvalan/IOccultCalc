/**
 * @file ioccultcalc_main.cpp
 * @brief IOccultCalc - Asteroid occultation prediction tool
 * 
 * Uso:
 *   ioccultcalc <preset_file>
 * 
 * Dove preset_file può essere:
 *   - file.oop  (formato OrbFit)  
 *   - file.json (formato JSON)
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ioccultcalc/astdys_client.h>
#include <ioccultcalc/orbit_propagator.h>
#include <ioccultcalc/orbital_elements.h>
#include <ioccultcalc/coordinates.h>
#include <ioccultcalc/time_utils.h>
#include <ioccultcalc/jpl_ephemeris.h>
#include <ioccultcalc/types.h>
#include <ioccultcalc/config_manager.h>
#include <ioccultcalc/chebyshev_astrometry.h>
#include <ioccultcalc/star_screening.h>

using namespace ioccultcalc;

// Alias per compatibilità con il codice esistente
using OccultationResult = ScreeningResult;

/**
 * @brief Converte RA/Dec (gradi) in radianti
 */
void parseCoordinates(const std::string& raStr, const std::string& decStr,
                     double& ra, double& dec) {
    ra = std::stod(raStr) * M_PI / 180.0;
    dec = std::stod(decStr) * M_PI / 180.0;
}

/**
 * @brief Scarica elementi orbitali da AstDyS per qualsiasi asteroide
 */
EquinoctialElements getOrbitalElements(const std::string& targetId, JulianDate epoch) {
    std::cout << "📡 Downloading orbital elements for asteroid: " << targetId << "\n";
    std::cout << "   Source: AstDyS/Lowell Observatory system\n";
    
    try {
        // Usa AstDyS client per scaricare elementi aggiornati
        ioccultcalc::AstDysClient client;
        
        // AstDyS accetta direttamente il numero come stringa
        auto elements = client.getElements(targetId);
        
        std::cout << "   ✓ Elements downloaded successfully\n";
        std::cout << "   Epoch: " << TimeUtils::jdToISO(elements.epoch) << "\n";
        std::cout << "   Name: " << elements.name << "\n\n";
        
        // Converti a elementi Kepleriani per visualizzazione
        OrbitalElements kepler = elements.toKeplerian();
        
        std::cout << "   Orbital elements:\n";
        std::cout << "     a = " << kepler.a << " AU\n";
        std::cout << "     e = " << kepler.e << "\n";
        std::cout << "     i = " << kepler.i * 180.0 / M_PI << "°\n";
        std::cout << "     Ω = " << kepler.Omega * 180.0 / M_PI << "°\n";
        std::cout << "     ω = " << kepler.omega * 180.0 / M_PI << "°\n";
        std::cout << "     M = " << kepler.M * 180.0 / M_PI << "°\n";
        if (kepler.H > 0) {
            std::cout << "     H = " << kepler.H << " mag\n";
        }
        if (kepler.diameter > 0) {
            std::cout << "     Diameter = " << kepler.diameter << " km\n";
        }
        std::cout << "\n";
        
        return elements;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error downloading elements: " << e.what() << "\n";
        std::cerr << "   Please check:\n";
        std::cerr << "   - Asteroid number/designation is valid\n";
        std::cerr << "   - Network connection is working\n";
        std::cerr << "   - AstDyS service is available\n\n";
        throw;
    }
}

/**
 * @brief Main: ricerca occultazioni in un intervallo temporale
 */
int main(int argc, char* argv[]) {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  IOccultCalc - Occultation Search (Any Asteroid)           ║\n";
    std::cout << "║  Propagazione con AST17 (17 asteroidi massivi)             ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    //
    // Per default la linea di comando può essere:
    //    > ioccultcalc preset.oop  //dove preset.oop è un file .oop con tutto il preset
    //    > ioccultcalc preset.json //dove il preset è in formato JSON
    //
    // NON E' ammessa altra configurazione per il lancio di ioccultcalc
    //
    
    // Verifica argomenti
    if (argc != 2) {
        std::cerr << "❌ Uso errato!\n\n";
        std::cerr << "Uso corretto:\n";
        std::cerr << "  ioccultcalc <preset_file>\n\n";
        std::cerr << "Dove preset_file può essere:\n";
        std::cerr << "  - file.oop  (formato OrbFit)\n";
        std::cerr << "  - file.json (formato JSON)\n\n";
        std::cerr << "Esempi:\n";
        std::cerr << "  ioccultcalc preset_17030_nov2025.oop\n";
        std::cerr << "  ioccultcalc config_survey.json\n\n";
        return 1;
    }
    
    std::string configFile = argv[1];
    
    // Verifica che il file esista
    std::ifstream checkFile(configFile);
    if (!checkFile.good()) {
        std::cerr << "❌ Errore: File '" << configFile << "' non trovato!\n\n";
        return 1;
    }
    checkFile.close();

    // Parametri da preset (OBBLIGATORI)
    std::string asteroidId;
    double starRA = 0.0;
    double starDec = 0.0;
    std::string startDateStr;
    std::string endDateStr;
    double asteroidDiameter = 0.0;
    double stepDays = 0.5;
    bool downloadFromAstDyS = true;
    std::string chebyshevCoeffFile;  // File coefficienti astrometrici (opzionale)
    
    // Setup modello astrometrico Chebyshev (opzionale)
    std::unique_ptr<ChebyshevAstrometry> chebyshevModel;
    
    // Tutto il logic principale in un singolo blocco try
    try {
        // Carica configurazione dal preset
        ConfigManager config;
        
        // Determina formato dal suffisso
        if (configFile.substr(configFile.length() - 4) == ".oop") {
            std::cout << "📄 Loading configuration from OrbFit file: " << configFile << "\n";
            config.loadFromOop(configFile);
        } else {
            std::cout << "📄 Loading configuration from JSON file: " << configFile << "\n";
            config.loadFromJson(configFile);
        }
        
        // Leggi parametri dalla configurazione (OBBLIGATORI)
        auto objSection = config.getSection(ConfigSection::OBJECT);
        if (objSection) {
            auto id = objSection->getParameter("id");
            if (id) {
                asteroidId = id->asString();
            } else {
                throw std::runtime_error("Parametro 'id' mancante nella sezione object");
            }
            
            auto diam = objSection->getParameter("diameter");
            if (diam) asteroidDiameter = diam->asDouble();
        } else {
            throw std::runtime_error("Sezione 'object' mancante nel preset");
        }
        
        auto searchSection = config.getSection(ConfigSection::SEARCH);
        if (searchSection) {
            auto start = searchSection->getParameter("start_jd");
            if (start) {
                JulianDate jd(start->asDouble());
                startDateStr = TimeUtils::jdToISO(jd).substr(0, 10);
            }
            
            auto end = searchSection->getParameter("end_jd");
            if (end) {
                JulianDate jd(end->asDouble());
                endDateStr = TimeUtils::jdToISO(jd).substr(0, 10);
            }
            
            auto step = searchSection->getParameter("step_days");
            if (step) stepDays = step->asDouble();
            
            // Coordinate stella (OBBLIGATORIE)
            auto star_ra = searchSection->getParameter("star_ra");
            auto star_dec = searchSection->getParameter("star_dec");
            if (star_ra && star_dec) {
                starRA = star_ra->asDouble() * M_PI / 180.0;  // gradi -> radianti
                starDec = star_dec->asDouble() * M_PI / 180.0;
        } else {
            throw std::runtime_error("Sezione 'search' mancante nel preset");
        }
        
        auto opsSection = config.getSection(ConfigSection::OPERATIONS);
        if (opsSection) {
            auto download = opsSection->getParameter("download_elements");
            if (download) downloadFromAstDyS = download->asBool();
            
            // Carica file coefficienti Chebyshev se specificato
            auto chebyshev_file = opsSection->getParameter("chebyshev_coefficients");
            if (chebyshev_file) {
                chebyshevCoeffFile = chebyshev_file->asString();
                std::cout << "   Chebyshev coefficients file: " << chebyshevCoeffFile << "\n";
            }
        }
        
        // Verifica che tutti i parametri obbligatori siano stati caricati
        if (asteroidId.empty()) {
            throw std::runtime_error("ID asteroide non specificato nel preset");
        }
        if (startDateStr.empty() || endDateStr.empty()) {
            throw std::runtime_error("Date di inizio/fine non specificate nel preset");
        }
        
        std::cout << "   ✓ Configuration loaded successfully\n";
        
        // Setup modello Chebyshev se specificato
        if (!chebyshevCoeffFile.empty()) {
            std::cout << "\n🔧 Loading Chebyshev astrometric model...\n";
            std::cout << "   File: " << chebyshevCoeffFile << "\n";
            
            chebyshevModel = std::make_unique<ChebyshevAstrometry>(5, 5); // ordine 5x5 di default
            
            if (chebyshevModel->loadCoefficientsFromFile(chebyshevCoeffFile)) {
                if (chebyshevModel->isValid()) {
                    const auto& coeffs = chebyshevModel->getCoefficients();
                    std::cout << "   ✓ Chebyshev model loaded successfully\n";
                    std::cout << "     Order: " << coeffs.max_order_u << "×" << coeffs.max_order_v << "\n";
                    std::cout << "     Reference: RA=" << (coeffs.alpha_base * 180.0 / M_PI) 
                             << "°, Dec=" << (coeffs.delta_base * 180.0 / M_PI) << "°\n";
                    std::cout << "   🎯 High-precision astrometry enabled!\n\n";
                } else {
                    std::cout << "   ⚠️  Invalid coefficients, using standard method\n\n";
                    chebyshevModel.reset();
                }
            } else {
                std::cout << "   ⚠️  Failed to load coefficients, using standard method\n\n";
                chebyshevModel.reset();
            }
        } else {
            std::cout << "   📊 Using standard astrometric method\n\n";
        }
        
        std::cout << "📋 Parameters from preset:\n";
    std::cout << "   Asteroid ID: " << asteroidId << "\n";
    std::cout << "   Star RA:  " << starRA * 180.0 / M_PI << "° (J2000)\n";
    std::cout << "   Star Dec: " << starDec * 180.0 / M_PI << "° (J2000)\n";
    std::cout << "   Start: " << startDateStr << "\n";
    std::cout << "   End:   " << endDateStr << "\n";
    std::cout << "   End:   " << endDateStr << "\n";
    if (asteroidDiameter > 0) {
        std::cout << "   Asteroid diameter (user): " << asteroidDiameter << " km\n";
    }
    std::cout << "\n";
    
        // 1. Carica elementi orbitali da AstDyS
        JulianDate startDate = TimeUtils::isoToJD(startDateStr);
        JulianDate endDate = TimeUtils::isoToJD(endDateStr);
        
        EquinoctialElements elements = getOrbitalElements(asteroidId, startDate);
        
        // Converti a Kepleriani per leggere H e diameter
        OrbitalElements kepler = elements.toKeplerian();
        
        // Se il diametro non è specificato dall'utente, usa quello da AstDyS
        if (asteroidDiameter <= 0 && kepler.diameter > 0) {
            asteroidDiameter = kepler.diameter;
            std::cout << "   Using diameter from AstDyS: " << asteroidDiameter << " km\n\n";
        } else if (asteroidDiameter <= 0 && kepler.H > 0) {
            // Default: stima dal H (magnitudine assoluta)
            // Formula: D (km) ≈ 1329 * 10^(-H/5) / sqrt(albedo)
            // Assumiamo albedo = 0.15 (tipico per asteroidi di tipo S)
            asteroidDiameter = 1329.0 * pow(10.0, -kepler.H / 5.0) / sqrt(0.15);
            std::cout << "   ⚠️  Diameter estimated from H=" << kepler.H 
                     << ": ~" << asteroidDiameter << " km (albedo=0.15)\n\n";
        } else if (asteroidDiameter <= 0) {
            // Fallback: 20 km (tipico NEA)
            asteroidDiameter = 20.0;
            std::cout << "   ⚠️  Diameter not available, using default: 20 km\n\n";
        }
        
        // 2. Setup sistema di screening Fase 1
        std::cout << "🔧 Setting up Phase 1 screening system...\n";
        
        // Configurazione screening
        ScreeningConfig screeningConfig;
        screeningConfig.stepSize = stepDays;
        screeningConfig.candidateThreshold = 60.0;  // 60 arcsec threshold
        screeningConfig.occultationThreshold = 15.0;  // 15 arcsec for occultations
        screeningConfig.probabilityThreshold = 0.01;  // 1% minimum probability
        screeningConfig.useDetailedLogging = true;
        screeningConfig.useProgressCallback = true;
        screeningConfig.useChebyshevIfAvailable = true;
        
        // Opzioni propagatore
        PropagatorOptions opts;
        opts.integrator = IntegratorType::RK4;
        opts.stepSize = 0.05;  // 0.05 giorni = ottimo trade-off
        opts.usePlanetaryPerturbations = true;
        opts.useRelativisticCorrections = false;
        opts.useSolarRadiationPressure = false;
        opts.useNonGravitational = false;
        opts.maxSteps = 1000000;
        
        // Inizializza sistema di screening
        StarScreening screening(screeningConfig);
        screening.setAsteroidElements(elements, opts);
        screening.setStarCoordinates(starRA, starDec);
        screening.setAsteroidProperties(asteroidDiameter, kepler.H);
        
        // Setup modello Chebyshev se disponibile
        if (chebyshevModel && chebyshevModel->isValid()) {
            screening.setChebyshevModel(chebyshevModel.get());
            std::cout << "   ✓ Chebyshev astrometric model integrated\n";
        }
        
        // Callback per progress
        screening.setProgressCallback([](int current, int total, const JulianDate& epoch, int candidates) {
            if (current % 50 == 0 || current == total) {
                std::cout << "\r   Progress: " << std::fixed << std::setprecision(1) 
                         << (100.0 * current / total) << "% [" << current << "/" << total 
                         << "] - Candidates: " << candidates << std::flush;
            }
        });
        
        std::cout << "   ✓ Screening system ready\n\n";
        
        // 3. Esegui screening Fase 1
        std::cout << "🔍 Executing Phase 1 screening...\n";
        std::vector<OccultationResult> candidates = screening.screenTimeInterval(startDate, endDate);
        
        std::cout << "\n\n";
        
        // 4. Mostra statistiche screening
        const auto& stats = screening.getLastStatistics();
        std::cout << "📊 Screening Statistics:\n";
        std::cout << "   Total calculations: " << stats.totalCalculations << "\n";
        std::cout << "   Computation time: " << std::fixed << std::setprecision(2) 
                 << stats.computationTime << " seconds\n";
        std::cout << "   Candidates found: " << stats.candidatesFound << "\n";
        std::cout << "   Predicted occultations: " << stats.occultationsFound << "\n";
        std::cout << "   Average astrometric error: " << std::fixed << std::setprecision(3) 
                 << stats.averageError << " arcsec\n";
        std::cout << "   Chebyshev model used: " << (stats.chebyshevUsed ? "Yes" : "No") << "\n\n";
        
        // 5. Mostra risultati
        std::cout << "═══════════════════════════════════════════════════════════\n";
        std::cout << "RESULTS: " << candidates.size() << " candidate occultations found\n";
        std::cout << "═══════════════════════════════════════════════════════════\n\n";
        
        if (candidates.empty()) {
            std::cout << "No close approaches found in this time span.\n";
            std::cout << "Try:\n";
            std::cout << "  - Longer time span\n";
            std::cout << "  - Different star coordinates\n";
            std::cout << "  - Different asteroid\n\n";
        } else {
            // Ordina per separazione (migliori prima)
            std::sort(candidates.begin(), candidates.end(),
                     [](const OccultationResult& a, const OccultationResult& b) {
                         return a.angularSeparation < b.angularSeparation;
                     });
            
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "┌────────────────────┬──────────┬──────────┬──────────┬──────────┐\n";
            std::cout << "│ Date (UTC)         │ Sep (\")  │ Prob (%) │ Dist(AU) │ Status   │\n";
            std::cout << "├────────────────────┼──────────┼──────────┼──────────┼──────────┤\n";
            
            for (size_t i = 0; i < std::min(candidates.size(), size_t(20)); i++) {
                const auto& c = candidates[i];
                std::string dateStr = TimeUtils::jdToISO(c.epoch);
                std::string status = c.isOccultation ? "✓ OCCULT" : "  Close";
                
                std::cout << "│ " << std::setw(18) << dateStr << " │ "
                         << std::setw(8) << c.angularSeparation << " │ "
                         << std::setw(8) << (c.probability * 100.0) << " │ "
                         << std::setw(8) << c.asteroidDistance << " │ "
                         << std::setw(8) << status << " │\n";
            }
            
            std::cout << "└────────────────────┴──────────┴──────────┴──────────┴──────────┘\n\n";
            
            // Dettagli del migliore candidato
            if (!candidates.empty()) {
                const auto& best = candidates[0];
                std::cout << "🎯 Best candidate:\n";
                std::cout << "   Date: " << TimeUtils::jdToISO(best.epoch) << "\n";
                std::cout << "   Separation: " << best.angularSeparation << " arcsec\n";
                std::cout << "   Probability: " << (best.probability * 100.0) << " %\n";
                std::cout << "   Distance: " << best.asteroidDistance << " AU\n";
                std::cout << "   Shadow width: ~" << best.shadowWidth << " km\n";
                std::cout << "   Asteroid RA/Dec: " << (best.asteroidRA * 180.0 / M_PI) 
                         << "° / " << (best.asteroidDec * 180.0 / M_PI) << "°\n";
                std::cout << "   Star RA/Dec: " << (best.starRA * 180.0 / M_PI) 
                         << "° / " << (best.starDec * 180.0 / M_PI) << "°\n";
                std::cout << "   Method: " << (best.usedChebyshev ? "Chebyshev" : "Standard") << "\n";
                if (best.usedChebyshev) {
                    std::cout << "   Astrometric error: " << best.astrometricError << " arcsec\n";
                }
                std::cout << "\n";
            }
        }
        
        std::cout << "✓ Search complete!\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n\n";
        return 1;
    }
    
    return 0;
}
