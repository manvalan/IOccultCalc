# IOC_StarMap - Integrazione per ItalOccult

## 📖 Panoramica

**IOC_StarMap** è una libreria esistente per generare mappe celesti professionali con preset configurati per:
1. **Mappe di Avvicinamento** (Approach Charts) - 6° × 4° FOV
2. **Mappe Finali** (Final Charts) - 3° × 2° FOV con dettagli precisi

Questa documentazione spiega come integrare IOC_StarMap nel formato ItalOccult.

---

## 🎯 Preset Disponibili

### 1. Approach Chart (Carta di Avvicinamento)
**Scopo**: Mostrare il percorso dell'asteroide nei 10 giorni precedenti l'occultazione

**Caratteristiche**:
- Campo visivo: **6° × 4°**
- Centro: Stella occultata
- Griglia equatoriale (RA/Dec)
- Traiettoria asteroide con marche temporali (ogni 2 giorni)
- Stelle fino a magnitudine 12
- Posizione finale asteroide evidenziata

**Preset IOC_StarMap**: `"approach_chart"`

```cpp
// Esempio parametri
{
  "preset": "approach_chart",
  "center_ra": 185.45,      // deg (12h 21m 48s)
  "center_dec": 12.34,      // deg (+12° 20' 24")
  "asteroid_track": [
    {"jd": 2460320.5, "ra": 185.2, "dec": 12.1},
    {"jd": 2460322.5, "ra": 185.3, "dec": 12.2},
    // ... ogni 2 giorni per 10 giorni
  ],
  "event_jd": 2460330.5,
  "width_deg": 6.0,
  "height_deg": 4.0,
  "mag_limit": 12.0,
  "grid": true,
  "labels": true
}
```

---

### 2. Final Chart (Carta Finale / Campo Dettagliato)
**Scopo**: Mostrare il campo stellare preciso al momento dell'occultazione

**Caratteristiche**:
- Campo visivo: **3° × 2°**
- Centro: Stella occultata
- Griglia fine (10' divisioni)
- Stelle fino a magnitudine 14
- Posizione prevista asteroide (con banda incertezza)
- Etichette stelle di riferimento
- Orientamento Nord su

**Preset IOC_StarMap**: `"final_chart"`

```cpp
// Esempio parametri
{
  "preset": "final_chart",
  "center_ra": 185.45,
  "center_dec": 12.34,
  "asteroid_pos": {"ra": 185.4501, "dec": 12.3402},
  "uncertainty_arcsec": 0.15,
  "event_jd": 2460330.5,
  "width_deg": 3.0,
  "height_deg": 2.0,
  "mag_limit": 14.0,
  "grid_spacing_arcmin": 10,
  "reference_stars": true,
  "orientation": "north_up"
}
```

---

## 🔧 API IOC_StarMap

### Funzione Base

```cpp
#include "ioc_starmap/generator.h"

std::string generateStarMap(const StarMapConfig& config);
```

### Strutture Dati

```cpp
namespace ioc::starmap {

struct CelestialPosition {
    double ra_deg;   // Right Ascension [0, 360)
    double dec_deg;  // Declination [-90, +90]
};

struct AsteroidTrackPoint {
    double jd;           // Julian Date
    CelestialPosition pos;
};

enum class MapPreset {
    APPROACH_CHART,  // 6° × 4°, 10-day track
    FINAL_CHART      // 3° × 2°, detailed field
};

struct StarMapConfig {
    MapPreset preset;
    CelestialPosition center;
    double event_jd;
    
    // Approach chart specific
    std::vector<AsteroidTrackPoint> asteroid_track;
    
    // Final chart specific
    CelestialPosition asteroid_position;
    double uncertainty_arcsec;
    
    // Display options
    double mag_limit;
    bool show_grid;
    bool show_labels;
    bool show_uncertainty_ellipse;
    
    // Output
    std::string output_format;  // "svg", "png", "jpg"
    int width_px;
    int height_px;
};

// Main generator class
class StarMapGenerator {
public:
    StarMapGenerator();
    
    // Generate map with preset
    std::string generate(const StarMapConfig& config);
    
    // Quick generation helpers
    std::string generateApproachChart(
        const CelestialPosition& star,
        const std::vector<AsteroidTrackPoint>& track,
        double event_jd,
        double mag_limit = 12.0
    );
    
    std::string generateFinalChart(
        const CelestialPosition& star,
        const CelestialPosition& asteroid,
        double uncertainty_arcsec,
        double event_jd,
        double mag_limit = 14.0
    );
    
    // Export to file
    bool exportToFile(const std::string& svg_content, const std::string& filename);
    
    // Get embedded SVG string (for HTML)
    std::string getEmbeddedSVG(const StarMapConfig& config);
};

} // namespace ioc::starmap
```

---

## 🔗 Integrazione con ItalOccultFormatter

### Header Update (`include/output_formatter.h`)

```cpp
#include "ioc_starmap/generator.h"

class ItalOccultFormatter : public OutputFormatterBase {
private:
    // Genera carta di avvicinamento usando IOC_StarMap
    std::string generateApproachChart(const OWOccultation& occ);
    
    // Genera carta finale usando IOC_StarMap
    std::string generateFinalChart(const OWOccultation& occ);
    
    // Helper per calcolare traiettoria 10 giorni
    std::vector<ioc::starmap::AsteroidTrackPoint> calculateAsteroidTrack(
        const OWOccultation& occ,
        int days_before = 10
    );
    
    ioc::starmap::StarMapGenerator starmap_generator_;
};
```

### Implementation (`src/output_formatter.cpp`)

```cpp
#include "ioc_starmap/generator.h"

std::string ItalOccultFormatter::generateApproachChart(const OWOccultation& occ) {
    using namespace ioc::starmap;
    
    // Calcola traiettoria asteroide (10 giorni prima)
    auto track = calculateAsteroidTrack(occ, 10);
    
    // Configura mappa
    StarMapConfig config;
    config.preset = MapPreset::APPROACH_CHART;
    config.center.ra_deg = occ.star_ra;
    config.center.dec_deg = occ.star_dec;
    config.event_jd = occ.event_jd;
    config.asteroid_track = track;
    config.mag_limit = 12.0;
    config.show_grid = true;
    config.show_labels = true;
    config.output_format = "svg";
    config.width_px = 600;
    config.height_px = 400;
    
    // Genera SVG
    return starmap_generator_.getEmbeddedSVG(config);
}

std::string ItalOccultFormatter::generateFinalChart(const OWOccultation& occ) {
    using namespace ioc::starmap;
    
    StarMapConfig config;
    config.preset = MapPreset::FINAL_CHART;
    config.center.ra_deg = occ.star_ra;
    config.center.dec_deg = occ.star_dec;
    config.event_jd = occ.event_jd;
    config.asteroid_position.ra_deg = occ.asteroid_ra;
    config.asteroid_position.dec_deg = occ.asteroid_dec;
    config.uncertainty_arcsec = occ.uncertainty_mas / 1000.0;  // mas → arcsec
    config.mag_limit = 14.0;
    config.show_grid = true;
    config.show_labels = true;
    config.show_uncertainty_ellipse = true;
    config.output_format = "svg";
    config.width_px = 600;
    config.height_px = 400;
    
    return starmap_generator_.getEmbeddedSVG(config);
}

std::vector<ioc::starmap::AsteroidTrackPoint> 
ItalOccultFormatter::calculateAsteroidTrack(const OWOccultation& occ, int days_before) {
    using namespace ioc::starmap;
    
    std::vector<AsteroidTrackPoint> track;
    
    // Recupera elementi orbitali asteroide
    OrbitalElements elements = occ.asteroid_elements;
    OrbitPropagator propagator;
    propagator.setElements(elements);
    
    // Calcola posizioni ogni 2 giorni
    double start_jd = occ.event_jd - days_before;
    for (int i = 0; i <= days_before; i += 2) {
        double jd = start_jd + i;
        
        // Propaga orbita
        propagator.propagateTo(jd);
        EphemerisData eph = propagator.getEphemeris(jd);
        
        AsteroidTrackPoint point;
        point.jd = jd;
        point.pos.ra_deg = eph.ra;
        point.pos.dec_deg = eph.dec;
        
        track.push_back(point);
    }
    
    return track;
}
```

### HTML Generation Update

```cpp
std::string ItalOccultFormatter::generateSheet(const OWOccultation& occ, size_t index) {
    std::ostringstream sheet;
    
    // ... (header section) ...
    
    // Middle Section: Maps
    sheet << "  <div class=\"middle-section\">\n";
    
    // Left: Approach Chart (IOC_StarMap)
    sheet << "    <div class=\"map-container\" style=\"border-right: 2px solid #333;\">\n";
    sheet << "      <div class=\"map-title\">Carta di Avvicinamento (6° × 4°)</div>\n";
    sheet << "      <div class=\"map-frame\">\n";
    
    std::string approach_svg = generateApproachChart(occ);
    sheet << "        " << approach_svg << "\n";
    
    sheet << "      </div>\n";
    sheet << "    </div>\n";
    
    // Right: Earth Map (IOC_Earth - già implementato)
    sheet << "    <div class=\"map-container\">\n";
    sheet << "      <div class=\"map-title\">Percorso sulla Terra</div>\n";
    sheet << "      <div class=\"map-frame\">\n";
    
    std::string earth_map_url = generateEarthMapURL(occ);
    sheet << "        <img src=\"" << earth_map_url << "\" alt=\"Earth Map\" />\n";
    
    sheet << "      </div>\n";
    sheet << "    </div>\n";
    sheet << "  </div>\n";
    
    // Bottom Section: Detailed Chart
    sheet << "  <div class=\"bottom-section\">\n";
    sheet << "    <div class=\"chart-container\">\n";
    sheet << "      <div class=\"chart-title\">Campo Dettagliato (3° × 2°)</div>\n";
    
    std::string final_svg = generateFinalChart(occ);
    sheet << "      " << final_svg << "\n";
    
    sheet << "    </div>\n";
    sheet << "  </div>\n";
    
    // ... (footer section) ...
    
    return sheet.str();
}
```

---

## 📦 Dipendenze CMakeLists.txt

```cmake
# IOC_StarMap library
find_package(IOC_StarMap REQUIRED)

target_link_libraries(ioccultcalc
    PRIVATE
        IOC_StarMap::Core
        # ... altre librerie ...
)

target_include_directories(ioccultcalc
    PRIVATE
        ${IOC_STARMAP_INCLUDE_DIRS}
)
```

---

## 🎨 Stile SVG per ItalOccult

### Approach Chart Customization

```cpp
// Opzioni aggiuntive per matching stile ItalOccult
config.style_options = {
    {"background", "#FFFFFF"},
    {"grid_color", "#CCCCCC"},
    {"grid_width", "1px"},
    {"star_color", "#000000"},
    {"asteroid_track_color", "#FF0000"},
    {"asteroid_track_width", "2px"},
    {"label_font", "Arial"},
    {"label_size", "10px"},
    {"title_enabled", false}  // Titolo gestito da ItalOccult HTML
};
```

### Final Chart Customization

```cpp
config.style_options = {
    {"background", "#FFFFFF"},
    {"grid_color", "#DDDDDD"},
    {"grid_width", "0.5px"},
    {"star_color", "#000000"},
    {"star_size_factor", 1.2},
    {"asteroid_color", "#FF0000"},
    {"asteroid_size", "8px"},
    {"uncertainty_color", "#FFC0C0"},
    {"uncertainty_opacity", "0.3"},
    {"label_font", "Arial"},
    {"label_size", "9px"},
    {"compass_enabled", true},
    {"scale_bar_enabled", true}
};
```

---

## ✅ Vantaggi IOC_StarMap

1. **Preset Configurati**: Parametri ottimizzati per occultazioni
2. **Precisione Astrometrica**: Riduzione epoche, moti propri, aberrazione
3. **Catalogo Gaia DR3**: Stelle fino a magnitudine 21
4. **Performance**: Rendering ottimizzato con cache
5. **SVG Scalabile**: Qualità vettoriale per stampa
6. **Compatibilità**: Output embeddable in HTML/PDF

---

## 🔄 Workflow Completo

```
ItalOccultFormatter
    ↓
generateSheet()
    ↓
├─→ generateApproachChart()
│       ↓
│   IOC_StarMap::generateApproachChart()
│       ↓ (preset: approach_chart)
│   ├─ Calcola traiettoria 10 giorni
│   ├─ Query Gaia DR3 (6° × 4° FOV)
│   ├─ Applica riduzione astrometrica
│   └─ Genera SVG embedded
│
├─→ generateEarthMapURL()
│       ↓
│   IOC_Earth API (già implementato)
│
└─→ generateFinalChart()
        ↓
    IOC_StarMap::generateFinalChart()
        ↓ (preset: final_chart)
    ├─ Query Gaia DR3 (3° × 2° FOV)
    ├─ Calcola posizione asteroide precisa
    ├─ Genera ellisse incertezza
    └─ Genera SVG embedded
```

---

## 📝 Esempio Output HTML

```html
<div class="middle-section">
  <div class="map-container">
    <div class="map-title">Carta di Avvicinamento (6° × 4°)</div>
    <div class="map-frame">
      <!-- SVG generato da IOC_StarMap approach_chart preset -->
      <svg width="600" height="400" xmlns="http://www.w3.org/2000/svg">
        <rect fill="#FFFFFF" width="600" height="400"/>
        <g id="grid"><!-- Griglia equatoriale --></g>
        <g id="stars"><!-- Stelle Gaia DR3 --></g>
        <g id="asteroid_track">
          <path d="M 100,200 L 150,195 L 200,190 ..." stroke="#FF0000" stroke-width="2"/>
          <circle cx="100" cy="200" r="3" fill="#FF0000"/><!-- T-10d -->
          <circle cx="150" cy="195" r="3" fill="#FF0000"/><!-- T-8d -->
          <!-- ... -->
        </g>
        <g id="labels"><!-- Etichette stelle/date --></g>
      </svg>
    </div>
  </div>
  
  <div class="map-container">
    <div class="map-title">Percorso sulla Terra</div>
    <div class="map-frame">
      <img src="https://ioc-earth.example.com/map?..." />
    </div>
  </div>
</div>

<div class="bottom-section">
  <div class="chart-container">
    <div class="chart-title">Campo Dettagliato (3° × 2°)</div>
    <!-- SVG generato da IOC_StarMap final_chart preset -->
    <svg width="600" height="400" xmlns="http://www.w3.org/2000/svg">
      <rect fill="#FFFFFF" width="600" height="400"/>
      <g id="fine_grid"><!-- Griglia 10' --></g>
      <g id="stars"><!-- Stelle fino a mag 14 --></g>
      <g id="asteroid">
        <circle cx="300" cy="200" r="8" fill="#FF0000"/>
        <ellipse cx="300" cy="200" rx="15" ry="10" fill="#FFC0C0" opacity="0.3"/>
      </g>
      <g id="compass"><!-- Nord up --></g>
      <g id="scale"><!-- Barra scala --></g>
    </svg>
  </div>
</div>
```

---

## 🚀 Test e Validazione

### Test Unit

```cpp
TEST(IOCStarMapIntegration, ApproachChart) {
    OWOccultation occ;
    occ.star_ra = 185.45;
    occ.star_dec = 12.34;
    occ.event_jd = 2460330.5;
    // ... set other fields ...
    
    ItalOccultFormatter formatter;
    std::string svg = formatter.generateApproachChart(occ);
    
    ASSERT_FALSE(svg.empty());
    ASSERT_TRUE(svg.find("<svg") != std::string::npos);
    ASSERT_TRUE(svg.find("id=\"asteroid_track\"") != std::string::npos);
}

TEST(IOCStarMapIntegration, FinalChart) {
    OWOccultation occ;
    // ... setup ...
    
    ItalOccultFormatter formatter;
    std::string svg = formatter.generateFinalChart(occ);
    
    ASSERT_FALSE(svg.empty());
    ASSERT_TRUE(svg.find("<svg") != std::string::npos);
    ASSERT_TRUE(svg.find("id=\"asteroid\"") != std::string::npos);
    ASSERT_TRUE(svg.find("ellipse") != std::string::npos);  // Uncertainty
}
```

### Test Integrazione Completa

```bash
# Genera scheda ItalOccult completa con IOC_StarMap
./build/ioccultcalc \
  --preset preset_test_italoccult.oop \
  --output-format ITALOCCULT \
  --asteroid 433 \
  --start-date 2026-01-01 \
  --end-date 2026-01-31

# Verifica SVG embedded
grep -c '<svg' example_italoccult.html  # Deve essere 2 (approach + final)

# Converti in PDF
python3 convert_to_pdf.py example_italoccult.html example_italoccult.pdf
```

---

## 📚 Riferimenti

- **IOC_StarMap Repository**: [github.com/manvalan/IOC_StarMap](https://github.com/manvalan/IOC_StarMap)
- **Gaia DR3 Documentation**: [gea.esac.esa.int/archive](https://gea.esac.esa.int/archive/)
- **ItalOccult Format Spec**: `ITALOCCULT_FORMAT_DOCUMENTATION.md`
- **IOC_Earth API Spec**: `IOC_EARTH_API_SPEC.md`

---

**Autore**: Michele Bigi (@manvalan)  
**Data**: 3 Dicembre 2025  
**Versione**: 1.0
