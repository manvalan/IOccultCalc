/**
 * @file italoccult_starmap_example.cpp
 * @brief Esempio integrazione IOC_StarMap con ItalOccultFormatter
 * 
 * Dimostra come usare i preset approach_chart e final_chart
 * di IOC_StarMap per generare mappe celesti professionali.
 * 
 * @author Michele Bigi
 * @date 2025-12-03
 */

#include "ioc_starmap/generator.h"
#include "ioccultcalc/occultation_predictor.h"
#include "ioccultcalc/orbit_propagator.h"
#include <iostream>
#include <vector>
#include <fstream>

using namespace ioc::starmap;
using namespace ioccultcalc;

/**
 * @brief Calcola traiettoria asteroide per N giorni prima dell'evento
 */
std::vector<AsteroidTrackPoint> calculateAsteroidTrack(
    const OrbitalElements& elements,
    double event_jd,
    int days_before = 10
) {
    std::vector<AsteroidTrackPoint> track;
    
    OrbitPropagator propagator;
    propagator.setElements(elements);
    
    double start_jd = event_jd - days_before;
    
    std::cout << "Calcolo traiettoria asteroide:\n";
    std::cout << "  Da: JD " << start_jd << "\n";
    std::cout << "  A:  JD " << event_jd << "\n";
    std::cout << "  Step: 2 giorni\n\n";
    
    for (int i = 0; i <= days_before; i += 2) {
        double jd = start_jd + i;
        
        propagator.propagateTo(jd);
        EphemerisData eph = propagator.getEphemeris(jd);
        
        AsteroidTrackPoint point;
        point.jd = jd;
        point.pos.ra_deg = eph.ra;
        point.pos.dec_deg = eph.dec;
        
        track.push_back(point);
        
        std::cout << "  T-" << (event_jd - jd) << "d: "
                  << "RA=" << eph.ra << "° "
                  << "Dec=" << eph.dec << "°\n";
    }
    
    return track;
}

/**
 * @brief Genera carta di avvicinamento usando IOC_StarMap
 */
std::string generateApproachChartExample() {
    std::cout << "\n=== APPROACH CHART (6° × 4°) ===\n\n";
    
    // Parametri esempio: (433) Eros occultazione 2026-01-16
    CelestialPosition star;
    star.ra_deg = 185.45;   // 12h 21m 48s
    star.dec_deg = 12.34;   // +12° 20' 24"
    
    double event_jd = 2460330.78945;  // 2026-01-16 06:56:48 UTC
    
    // Elementi orbitali Eros (semplificati per esempio)
    OrbitalElements elements;
    elements.a = 1.458;      // AU
    elements.e = 0.223;
    elements.i = 10.83;      // deg
    elements.Omega = 304.30; // deg
    elements.omega = 178.82; // deg
    elements.M0 = 320.23;    // deg
    elements.epoch = event_jd - 10;
    
    // Calcola traiettoria 10 giorni
    auto track = calculateAsteroidTrack(elements, event_jd, 10);
    
    // Configura IOC_StarMap
    StarMapConfig config;
    config.preset = MapPreset::APPROACH_CHART;
    config.center = star;
    config.event_jd = event_jd;
    config.asteroid_track = track;
    config.mag_limit = 12.0;
    config.show_grid = true;
    config.show_labels = true;
    config.output_format = "svg";
    config.width_px = 600;
    config.height_px = 400;
    
    // Stile personalizzato per ItalOccult
    config.style_options = {
        {"background", "#FFFFFF"},
        {"grid_color", "#CCCCCC"},
        {"star_color", "#000000"},
        {"asteroid_track_color", "#FF0000"},
        {"asteroid_track_width", "2px"},
        {"label_font", "Arial"},
        {"label_size", "10px"}
    };
    
    // Genera SVG
    StarMapGenerator generator;
    std::string svg = generator.getEmbeddedSVG(config);
    
    std::cout << "✓ SVG generato (" << svg.length() << " bytes)\n";
    
    return svg;
}

/**
 * @brief Genera carta finale usando IOC_StarMap
 */
std::string generateFinalChartExample() {
    std::cout << "\n=== FINAL CHART (3° × 2°) ===\n\n";
    
    // Stella occultata
    CelestialPosition star;
    star.ra_deg = 185.45;
    star.dec_deg = 12.34;
    
    // Posizione asteroide al momento dell'evento
    CelestialPosition asteroid;
    asteroid.ra_deg = 185.4501;  // Vicino alla stella
    asteroid.dec_deg = 12.3402;
    
    double event_jd = 2460330.78945;
    double uncertainty_arcsec = 0.15;  // 150 mas
    
    std::cout << "Configurazione:\n";
    std::cout << "  Stella:     RA=" << star.ra_deg << "° Dec=" << star.dec_deg << "°\n";
    std::cout << "  Asteroide:  RA=" << asteroid.ra_deg << "° Dec=" << asteroid.dec_deg << "°\n";
    std::cout << "  Separazione: " << 
        (fabs(star.ra_deg - asteroid.ra_deg) * 3600.0) << " arcsec (RA)\n";
    std::cout << "  Incertezza:  " << uncertainty_arcsec << " arcsec\n\n";
    
    // Configura IOC_StarMap
    StarMapConfig config;
    config.preset = MapPreset::FINAL_CHART;
    config.center = star;
    config.event_jd = event_jd;
    config.asteroid_position = asteroid;
    config.uncertainty_arcsec = uncertainty_arcsec;
    config.mag_limit = 14.0;
    config.show_grid = true;
    config.show_labels = true;
    config.show_uncertainty_ellipse = true;
    config.output_format = "svg";
    config.width_px = 600;
    config.height_px = 400;
    
    // Stile personalizzato
    config.style_options = {
        {"background", "#FFFFFF"},
        {"grid_color", "#DDDDDD"},
        {"star_color", "#000000"},
        {"asteroid_color", "#FF0000"},
        {"uncertainty_color", "#FFC0C0"},
        {"uncertainty_opacity", "0.3"},
        {"compass_enabled", "true"},
        {"scale_bar_enabled", "true"}
    };
    
    // Genera SVG
    StarMapGenerator generator;
    std::string svg = generator.getEmbeddedSVG(config);
    
    std::cout << "✓ SVG generato (" << svg.length() << " bytes)\n";
    
    return svg;
}

/**
 * @brief Genera scheda ItalOccult completa con entrambe le mappe
 */
void generateCompleteItalOccultSheet() {
    std::cout << "\n=== SCHEDA ITALOCCULT COMPLETA ===\n\n";
    
    // Genera approach chart
    std::string approach_svg = generateApproachChartExample();
    
    // Genera final chart
    std::string final_svg = generateFinalChartExample();
    
    // Costruisci HTML
    std::ostringstream html;
    
    html << "<!DOCTYPE html>\n";
    html << "<html>\n<head>\n";
    html << "  <title>ItalOccult - (433) Eros - 2026-01-16</title>\n";
    html << "  <meta charset=\"UTF-8\">\n";
    html << "  <style>\n";
    html << "    @page { size: A4; margin: 10mm; }\n";
    html << "    body { font-family: Arial; margin: 0; padding: 0; }\n";
    html << "    .sheet { width: 210mm; height: 297mm; }\n";
    html << "    .header { height: 20%; border-bottom: 2px solid #000; }\n";
    html << "    .middle { height: 40%; display: flex; }\n";
    html << "    .map-container { width: 50%; padding: 10px; }\n";
    html << "    .bottom { height: 35%; padding: 10px; }\n";
    html << "    .footer { height: 5%; text-align: center; }\n";
    html << "  </style>\n";
    html << "</head>\n<body>\n";
    
    html << "<div class=\"sheet\">\n";
    
    // Header
    html << "  <div class=\"header\">\n";
    html << "    <h1>Scheda Occultazione ItalOccult</h1>\n";
    html << "    <h2>(433) Eros - 2026-01-16 06:56:48 UTC</h2>\n";
    html << "    <p>Stella: HIP 60832 (mag 8.2) | RA 12h 21m 48s | Dec +12° 20' 24\"</p>\n";
    html << "  </div>\n";
    
    // Middle section: mappe
    html << "  <div class=\"middle\">\n";
    html << "    <div class=\"map-container\">\n";
    html << "      <h3>Carta di Avvicinamento (6° × 4°)</h3>\n";
    html << "      " << approach_svg << "\n";
    html << "    </div>\n";
    html << "    <div class=\"map-container\">\n";
    html << "      <h3>Percorso sulla Terra</h3>\n";
    html << "      <p>[IOC_Earth map qui]</p>\n";
    html << "    </div>\n";
    html << "  </div>\n";
    
    // Bottom section: campo dettagliato
    html << "  <div class=\"bottom\">\n";
    html << "    <h3>Campo Dettagliato (3° × 2°)</h3>\n";
    html << "    " << final_svg << "\n";
    html << "  </div>\n";
    
    // Footer
    html << "  <div class=\"footer\">\n";
    html << "    <p>Generato con IOccultCalc + IOC_StarMap | " 
         << "Data: 2025-12-03</p>\n";
    html << "  </div>\n";
    
    html << "</div>\n";
    html << "</body>\n</html>\n";
    
    // Salva file
    std::ofstream outfile("example_starmap_italoccult.html");
    outfile << html.str();
    outfile.close();
    
    std::cout << "\n✓ File salvato: example_starmap_italoccult.html\n";
    std::cout << "  Dimensione: " << html.str().length() << " bytes\n";
    std::cout << "\nPer convertire in PDF:\n";
    std::cout << "  python3 convert_to_pdf.py example_starmap_italoccult.html\n";
}

/**
 * @brief Main
 */
int main() {
    std::cout << "=================================================\n";
    std::cout << "IOC_StarMap Integration Example\n";
    std::cout << "ItalOccult Format - Celestial Charts Generator\n";
    std::cout << "=================================================\n";
    
    try {
        generateCompleteItalOccultSheet();
        
        std::cout << "\n✓ Esempio completato con successo!\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Errore: " << e.what() << "\n";
        return 1;
    }
}
