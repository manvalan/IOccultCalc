# IOC_StarMap - Riepilogo Integrazione ItalOccult

## 🎯 Obiettivo

Utilizzare la libreria **IOC_StarMap** esistente con i suoi **preset configurati** per generare le mappe celesti nel formato ItalOccult.

---

## 📋 Preset IOC_StarMap Disponibili

### 1. **approach_chart** - Carta di Avvicinamento
- **Campo visivo**: 6° × 4°
- **Scopo**: Traiettoria asteroide 10 giorni prima dell'evento
- **Contenuto**:
  - Griglia equatoriale (RA/Dec)
  - Stelle fino a magnitudine 12 (Gaia DR3)
  - Percorso asteroide con marche temporali (ogni 2 giorni)
  - Posizione finale evidenziata

### 2. **final_chart** - Campo Dettagliato
- **Campo visivo**: 3° × 2°
- **Scopo**: Campo stellare preciso al momento dell'occultazione
- **Contenuto**:
  - Griglia fine (10' divisioni)
  - Stelle fino a magnitudine 14 (Gaia DR3)
  - Posizione prevista asteroide
  - Ellisse incertezza
  - Bussola e scala

---

## 🔧 API Semplificata

```cpp
#include "ioc_starmap/generator.h"

// Genera carta avvicinamento
std::string generateApproachChart(
    const CelestialPosition& star,
    const std::vector<AsteroidTrackPoint>& track,
    double event_jd,
    double mag_limit = 12.0
);

// Genera carta finale
std::string generateFinalChart(
    const CelestialPosition& star,
    const CelestialPosition& asteroid,
    double uncertainty_arcsec,
    double event_jd,
    double mag_limit = 14.0
);
```

---

## 📦 Strutture Dati

```cpp
struct CelestialPosition {
    double ra_deg;   // [0, 360)
    double dec_deg;  // [-90, +90]
};

struct AsteroidTrackPoint {
    double jd;           // Julian Date
    CelestialPosition pos;
};
```

---

## 🔗 Integrazione in ItalOccultFormatter

```cpp
class ItalOccultFormatter : public OutputFormatterBase {
private:
    ioc::starmap::StarMapGenerator starmap_generator_;
    
    std::string generateApproachChart(const OWOccultation& occ);
    std::string generateFinalChart(const OWOccultation& occ);
    
    std::vector<ioc::starmap::AsteroidTrackPoint> calculateAsteroidTrack(
        const OWOccultation& occ, 
        int days_before = 10
    );
};
```

---

## 🎨 Layout ItalOccult con IOC_StarMap

```
┌─────────────────────────────────────────────────────────┐
│  HEADER (20%)                                           │
│  Titolo, asteroide, stella, data, coordinate           │
├──────────────────────────┬──────────────────────────────┤
│  MIDDLE (40%)            │                              │
│  ┌────────────────────┐  │  ┌────────────────────────┐  │
│  │ Approach Chart     │  │  │ Earth Map (IOC_Earth) │  │
│  │ (IOC_StarMap)      │  │  │                        │  │
│  │ 6° × 4° FOV        │  │  │ Ground track           │  │
│  │ 10-day track       │  │  │ Time markers           │  │
│  └────────────────────┘  │  └────────────────────────┘  │
├──────────────────────────┴──────────────────────────────┤
│  BOTTOM (35%)                                           │
│  ┌────────────────────────────────────────────────────┐ │
│  │ Final Chart (IOC_StarMap)                         │ │
│  │ 3° × 2° FOV                                       │ │
│  │ Detailed field + uncertainty ellipse              │ │
│  └────────────────────────────────────────────────────┘ │
├─────────────────────────────────────────────────────────┤
│  FOOTER (5%)                                            │
│  Logo ItalOccult, data generazione                     │
└─────────────────────────────────────────────────────────┘
```

---

## 💻 Esempio Implementazione

```cpp
std::string ItalOccultFormatter::generateApproachChart(const OWOccultation& occ) {
    using namespace ioc::starmap;
    
    // 1. Calcola traiettoria 10 giorni
    auto track = calculateAsteroidTrack(occ, 10);
    
    // 2. Configura parametri
    StarMapConfig config;
    config.preset = MapPreset::APPROACH_CHART;
    config.center = {occ.star_ra, occ.star_dec};
    config.event_jd = occ.event_jd;
    config.asteroid_track = track;
    config.mag_limit = 12.0;
    config.width_px = 600;
    config.height_px = 400;
    
    // 3. Genera SVG embedded
    return starmap_generator_.getEmbeddedSVG(config);
}

std::string ItalOccultFormatter::generateFinalChart(const OWOccultation& occ) {
    using namespace ioc::starmap;
    
    StarMapConfig config;
    config.preset = MapPreset::FINAL_CHART;
    config.center = {occ.star_ra, occ.star_dec};
    config.asteroid_position = {occ.asteroid_ra, occ.asteroid_dec};
    config.uncertainty_arcsec = occ.uncertainty_mas / 1000.0;
    config.event_jd = occ.event_jd;
    config.mag_limit = 14.0;
    config.show_uncertainty_ellipse = true;
    config.width_px = 600;
    config.height_px = 400;
    
    return starmap_generator_.getEmbeddedSVG(config);
}

std::vector<ioc::starmap::AsteroidTrackPoint> 
ItalOccultFormatter::calculateAsteroidTrack(const OWOccultation& occ, int days) {
    std::vector<ioc::starmap::AsteroidTrackPoint> track;
    
    OrbitPropagator propagator;
    propagator.setElements(occ.asteroid_elements);
    
    double start_jd = occ.event_jd - days;
    
    for (int i = 0; i <= days; i += 2) {  // Ogni 2 giorni
        double jd = start_jd + i;
        propagator.propagateTo(jd);
        
        EphemerisData eph = propagator.getEphemeris(jd);
        
        track.push_back({
            jd,
            {eph.ra, eph.dec}
        });
    }
    
    return track;
}
```

---

## 📝 Output HTML con IOC_StarMap

```html
<div class="middle-section">
  <!-- Carta Avvicinamento (IOC_StarMap) -->
  <div class="map-container">
    <div class="map-title">Carta di Avvicinamento (6° × 4°)</div>
    <div class="map-frame">
      <svg width="600" height="400">
        <!-- SVG generato da IOC_StarMap preset approach_chart -->
        <g id="grid">...</g>
        <g id="stars">...</g>
        <g id="asteroid_track">
          <path d="..." stroke="#FF0000"/>
        </g>
      </svg>
    </div>
  </div>
  
  <!-- Mappa Terra (IOC_Earth) -->
  <div class="map-container">
    <div class="map-title">Percorso sulla Terra</div>
    <div class="map-frame">
      <img src="https://ioc-earth.example.com/..." />
    </div>
  </div>
</div>

<div class="bottom-section">
  <!-- Campo Dettagliato (IOC_StarMap) -->
  <div class="chart-container">
    <div class="chart-title">Campo Dettagliato (3° × 2°)</div>
    <svg width="600" height="400">
      <!-- SVG generato da IOC_StarMap preset final_chart -->
      <g id="fine_grid">...</g>
      <g id="stars">...</g>
      <g id="asteroid">
        <circle cx="300" cy="200" r="8" fill="#FF0000"/>
        <ellipse cx="300" cy="200" rx="15" ry="10" fill="#FFC0C0" opacity="0.3"/>
      </g>
    </svg>
  </div>
</div>
```

---

## ✅ Vantaggi

1. **Preset Ottimizzati**: Parametri già configurati per occultazioni
2. **Gaia DR3**: Catalogo stellare preciso fino a mag 21
3. **SVG Embedded**: Qualità vettoriale, nessuna dipendenza esterna
4. **Riduzione Astrometrica**: Moti propri, aberrazione, precessione
5. **Performance**: Cache query, rendering veloce
6. **Compatibilità PDF**: SVG preservato in conversione Playwright

---

## 🔄 Workflow

```
Input: OWOccultation
    ↓
ItalOccultFormatter::generateSheet()
    ↓
├─→ generateApproachChart()
│       ↓
│   IOC_StarMap::approach_chart preset
│       ↓
│   Calcola track 10 giorni → Query Gaia DR3 → SVG
│
├─→ generateEarthMapURL()
│       ↓
│   IOC_Earth API (già implementato)
│
└─→ generateFinalChart()
        ↓
    IOC_StarMap::final_chart preset
        ↓
    Posizione asteroide + uncertainty → Query Gaia DR3 → SVG
```

---

## 📦 CMakeLists.txt

```cmake
find_package(IOC_StarMap REQUIRED)

target_link_libraries(ioccultcalc
    PRIVATE
        IOC_StarMap::Core
)
```

---

## 🚀 Test

```bash
# Genera scheda completa
./build/ioccultcalc --preset preset_italoccult.oop --asteroid 433

# Verifica SVG embedded (deve essere 2: approach + final)
grep -c '<svg' example_italoccult.html

# Converti in PDF
python3 convert_to_pdf.py example_italoccult.html

# Apri PDF
open example_italoccult.pdf
```

---

## 📚 Documentazione Completa

- **IOC_STARMAP_INTEGRATION.md** - Spec completa API e integrazione
- **IOC_EARTH_API_SPEC.md** - Spec IOC_Earth (già implementato)
- **ITALOCCULT_INTEGRATION_GUIDE.md** - Guida generale formato

---

**Pronto per l'implementazione! 🚀**

Usa i preset `approach_chart` e `final_chart` di IOC_StarMap per generare le mappe celesti professionali nelle schede ItalOccult.
