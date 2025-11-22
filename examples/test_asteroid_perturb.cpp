/**
 * @file test_asteroid_perturb.cpp
 * @brief Test perturbazioni asteroidi massivi (Ceres, Pallas, Vesta)
 * 
 * Testa che il ForceModel calcoli correttamente le accelerazioni
 * dovute ai grandi asteroidi (Ceres, Pallas, Vesta) usando SB441-N16.bsp
 */

#include "ioccultcalc/force_model.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ioccultcalc;

int main() {
    std::cout << "\n=== Test Perturbazioni Asteroidi Massivi ===\n\n";
    
    try {
        // Epoca di test J2000.0
        double jd = 2451545.0;
        std::cout << "Epoca: JD " << std::fixed << std::setprecision(2) << jd << " (J2000.0)\n\n";
        
        // Posizione test: orbita tipica NEA a ~1.5 AU (simile a Eros)
        Vector3D testPos(1.5, 0.0, 0.1);  // AU, eliocentrica
        Vector3D testVel(0.0, 0.015, 0.0);  // AU/day
        
        std::cout << "Posizione test: (" << std::setprecision(6)
                  << testPos.x << ", " << testPos.y << ", " << testPos.z << ") AU\n";
        std::cout << "Velocità test: (" << testVel.x << ", " 
                  << testVel.y << ", " << testVel.z << ") AU/day\n\n";
    
        // Test 1: ForceModel SENZA asteroidi
        std::cout << "=== Test 1: Standard Config (no asteroidi) ===\n";
        
        ForceModelConfig config1 = ForceModelConfig::standardConfig();
        config1.includeCeres = false;
        config1.includePallas = false;
        config1.includeVesta = false;
        
        ForceModel forceModel1(config1);
        forceModel1.initializeJPL();
        
        std::map<PerturbingBody, Vector3D> contrib1;
        Vector3D accel1 = forceModel1.computeAccelerationWithBreakdown(
            jd, testPos, testVel, contrib1);
        
        std::cout << "Accelerazione totale: " << std::scientific << std::setprecision(6)
                  << accel1.magnitude() << " AU/day²\n";
        std::cout << "Componenti: (" << accel1.x << ", " 
                  << accel1.y << ", " << accel1.z << ")\n\n";
        
        std::cout << "Contributi principali:\n";
        for (const auto& [body, accel] : contrib1) {
            double mag = accel.magnitude();
            if (mag > 1e-15) {
                std::string bodyName;
                switch (body) {
                    case PerturbingBody::SUN: bodyName = "Sun"; break;
                    case PerturbingBody::JUPITER: bodyName = "Jupiter"; break;
                    case PerturbingBody::EARTH: bodyName = "Earth"; break;
                    case PerturbingBody::MARS: bodyName = "Mars"; break;
                    default: bodyName = "Other"; break;
                }
                std::cout << "  " << std::setw(10) << bodyName << ": " 
                          << std::setprecision(3) << mag << " AU/day²\n";
            }
        }
        std::cout << "\n";
        
        // Test 2: ForceModel CON asteroidi
        std::cout << "=== Test 2: High Precision Config (con asteroidi) ===\n";
        
        ForceModelConfig config2 = ForceModelConfig::highPrecisionConfig();
        
        ForceModel forceModel2(config2);
        forceModel2.initializeJPL();
        
        std::map<PerturbingBody, Vector3D> contrib2;
        Vector3D accel2 = forceModel2.computeAccelerationWithBreakdown(
            jd, testPos, testVel, contrib2);
        
        std::cout << "Accelerazione totale: " << std::scientific << std::setprecision(6)
                  << accel2.magnitude() << " AU/day²\n";
        std::cout << "Componenti: (" << accel2.x << ", " 
                  << accel2.y << ", " << accel2.z << ")\n\n";
        
        std::cout << "Contributi principali:\n";
        for (const auto& [body, accel] : contrib2) {
            double mag = accel.magnitude();
            if (mag > 1e-20) {  // Mostra anche asteroidi anche se piccoli
                std::string bodyName;
                switch (body) {
                    case PerturbingBody::SUN: bodyName = "Sun"; break;
                    case PerturbingBody::JUPITER: bodyName = "Jupiter"; break;
                    case PerturbingBody::EARTH: bodyName = "Earth"; break;
                    case PerturbingBody::MARS: bodyName = "Mars"; break;
                    case PerturbingBody::CERES: bodyName = "** CERES **"; break;
                    case PerturbingBody::PALLAS: bodyName = "** PALLAS **"; break;
                    case PerturbingBody::VESTA: bodyName = "** VESTA **"; break;
                    default: bodyName = "Other"; break;
                }
                std::cout << "  " << std::setw(12) << bodyName << ": " 
                          << std::setprecision(3) << mag << " AU/day²\n";
            }
        }
        std::cout << "\n";
    
        // Analizza differenze
        std::cout << "=== Analisi Differenza ===\n";
        
        Vector3D accelDiff = accel2 - accel1;
        double diffMag = accelDiff.magnitude();
        
        std::cout << "Differenza accelerazione: " 
                  << std::scientific << std::setprecision(3)
                  << diffMag << " AU/day²\n";
        
        std::cout << "Componenti: ("
                  << accelDiff.x << ", " 
                  << accelDiff.y << ", " 
                  << accelDiff.z << ")\n\n";
        
        // Stima impatto su 1000 giorni di propagazione
        // Δr ≈ 0.5 * Δa * t² (approssimazione)
        double t = 1000.0;  // giorni
        double posError = 0.5 * diffMag * t * t;  // AU
        double posErrorKm = posError * 1.495978707e8;  // km
        double posErrorArcsec = posError * 206265.0 / testPos.magnitude();  // arcsec
        
        std::cout << "=== Stima Impatto su Propagazione (1000 giorni) ===\n";
        std::cout << "Errore posizione stimato: " << std::fixed << std::setprecision(1)
                  << posErrorKm << " km\n";
        std::cout << "Errore angolare stimato: " << std::setprecision(3)
                  << posErrorArcsec << " arcsec\n\n";
        
        // Verifica se gli asteroidi sono stati caricati
        bool ceresFound = false;
        bool pallasFound = false;
        bool vestaFound = false;
        
        for (const auto& [body, accel] : contrib2) {
            if (body == PerturbingBody::CERES && accel.magnitude() > 0) ceresFound = true;
            if (body == PerturbingBody::PALLAS && accel.magnitude() > 0) pallasFound = true;
            if (body == PerturbingBody::VESTA && accel.magnitude() > 0) vestaFound = true;
        }
        
        std::cout << "=== Diagnostica SPK ===\n";
        std::cout << "Ceres caricato: " << (ceresFound ? "✓ YES" : "✗ NO") << "\n";
        std::cout << "Pallas caricato: " << (pallasFound ? "✓ YES" : "✗ NO") << "\n";
        std::cout << "Vesta caricato: " << (vestaFound ? "✓ YES" : "✗ NO") << "\n\n";
        
        std::cout << "=== Conclusione ===\n";
        if (ceresFound || pallasFound || vestaFound) {
            std::cout << "✓ SPK asteroidi funzionante!\n";
            if (posErrorKm > 100.0) {
                std::cout << "✓ Perturbazioni asteroidi significative (>" 
                          << std::setprecision(0) << posErrorKm << " km stimati in 1000 giorni)\n";
            } else if (posErrorKm > 10.0) {
                std::cout << "✓ Perturbazioni asteroidi rilevabili\n";
            } else {
                std::cout << "○ Perturbazioni asteroidi minori ma misurabili\n";
            }
        } else {
            std::cout << "⚠ SPK asteroidi NON caricato o file sb441-n16.bsp mancante\n";
            std::cout << "  Verificare che sb441-n16.bsp sia nella directory dei dati SPK\n";
        }
        
        std::cout << "\nTest completato.\n\n";
        
    } catch (const std::exception& e) {
        std::cerr << "ERRORE: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
