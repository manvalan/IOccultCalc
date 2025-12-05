# 📚 Documentazione IOC_StarMap per ItalOccult

## 📁 File Creati

### 1. **IOC_STARMAP_INTEGRATION.md** (650+ righe)
**Documentazione completa API e integrazione**

- Panoramica libreria IOC_StarMap
- Preset disponibili (approach_chart, final_chart)
- API completa con strutture dati
- Integrazione con ItalOccultFormatter
- Layout HTML con mappe celesti
- Dipendenze CMakeLists.txt
- Personalizzazione stile SVG
- Test e validazione

### 2. **IOC_STARMAP_QUICKSTART.md** (314 righe)
**Guida rapida per iniziare subito**

- Obiettivo e preset disponibili
- API semplificata
- Strutture dati essenziali
- Esempio implementazione
- Layout scheda ItalOccult completo
- Workflow visivo
- Output HTML con SVG embedded
- Test checklist

### 3. **IOC_STARMAP_CODE_UPDATE.sh** (90+ righe)
**Script bash automatico per aggiornamento codice**

- Verifica libreria IOC_StarMap
- Backup file esistenti
- Aggiornamento CMakeLists.txt
- Lista modifiche necessarie al codice
- Test checklist
- Riepilogo documentazione

### 4. **examples/italoccult_starmap_example.cpp** (260+ righe)
**Esempio completo funzionante**

- Calcolo traiettoria asteroide
- Generazione approach chart (6° × 4°)
- Generazione final chart (3° × 2°)
- Output HTML completo
- Configurazione stile personalizzato
- Commenti dettagliati

---

## 🎯 Preset IOC_StarMap

### **approach_chart** - Carta Avvicinamento
```
Campo:     6° × 4°
Mag limit: 12.0
Contenuto: - Griglia equatoriale
           - Stelle Gaia DR3
           - Traiettoria 10 giorni (step 2d)
           - Marche temporali
Output:    SVG embedded
```

### **final_chart** - Campo Dettagliato
```
Campo:     3° × 2°
Mag limit: 14.0
Contenuto: - Griglia fine (10')
           - Stelle Gaia DR3
           - Posizione asteroide
           - Ellisse incertezza
           - Bussola + scala
Output:    SVG embedded
```

---

## 🔧 Integrazione ItalOccultFormatter

### Header (`include/output_formatter.h`)
```cpp
#include "ioc_starmap/generator.h"

class ItalOccultFormatter : public OutputFormatterBase {
private:
    ioc::starmap::StarMapGenerator starmap_generator_;
    
    std::string generateApproachChart(const OWOccultation& occ);
    std::string generateFinalChart(const OWOccultation& occ);
    std::vector<ioc::starmap::AsteroidTrackPoint> 
        calculateAsteroidTrack(const OWOccultation& occ, int days = 10);
};
```

### Implementation (`src/output_formatter.cpp`)
```cpp
std::string ItalOccultFormatter::generateApproachChart(const OWOccultation& occ) {
    // 1. Calcola traiettoria
    auto track = calculateAsteroidTrack(occ, 10);
    
    // 2. Configura IOC_StarMap
    StarMapConfig config;
    config.preset = MapPreset::APPROACH_CHART;
    config.center = {occ.star_ra, occ.star_dec};
    config.asteroid_track = track;
    config.event_jd = occ.event_jd;
    config.mag_limit = 12.0;
    
    // 3. Genera SVG
    return starmap_generator_.getEmbeddedSVG(config);
}
```

---

## 📦 CMakeLists.txt

```cmake
# IOC_StarMap library
add_subdirectory(external/IOC_StarMap)

target_link_libraries(ioccultcalc
    PRIVATE
        IOC_StarMap::Core
)
```

---

## 🎨 Layout ItalOccult

```
┌─────────────────────────────────────────┐
│ HEADER (20%)                            │ 
│ Titolo, dati asteroide/stella           │
├────────────────────┬────────────────────┤
│ APPROACH CHART     │ EARTH MAP          │
│ (IOC_StarMap)      │ (IOC_Earth)        │
│ 6° × 4°            │ Ground track       │
│ Track 10 giorni    │ Time markers       │
├────────────────────┴────────────────────┤
│ FINAL CHART (IOC_StarMap)               │
│ 3° × 2° - Campo dettagliato             │
│ Posizione + incertezza                  │
├─────────────────────────────────────────┤
│ FOOTER (5%) - Logo, data                │
└─────────────────────────────────────────┘
```

---

## ✅ Vantaggi

1. **Preset Ottimizzati** - Parametri già configurati
2. **Gaia DR3** - Catalogo preciso fino a mag 21
3. **SVG Embedded** - Qualità vettoriale, no dipendenze
4. **Astrometria** - Riduzione completa (moti propri, aberrazione)
5. **Performance** - Cache + rendering veloce
6. **PDF Ready** - SVG preservato in Playwright

---

## 🚀 Quick Start

### 1. Setup
```bash
# Clona IOC_StarMap
git submodule add https://github.com/manvalan/IOC_StarMap external/IOC_StarMap

# Aggiorna CMakeLists.txt
# (vedi IOC_STARMAP_CODE_UPDATE.sh)
```

### 2. Compila esempio
```bash
# Aggiungi al CMakeLists.txt:
add_executable(italoccult_starmap_example 
    examples/italoccult_starmap_example.cpp)
target_link_libraries(italoccult_starmap_example 
    IOC_StarMap::Core ioccultcalc)

# Compila
./build.sh

# Esegui
./build/examples/italoccult_starmap_example
```

### 3. Test output
```bash
# Verifica HTML
cat example_starmap_italoccult.html | grep -c '<svg'
# Output: 2 (approach + final)

# Converti PDF
python3 convert_to_pdf.py example_starmap_italoccult.html

# Visualizza
open example_starmap_italoccult.pdf
```

---

## 📊 Output Atteso

**HTML con 2 SVG embedded:**
- Approach Chart: ~8-10 KB (6° × 4° FOV, ~200 stelle)
- Final Chart: ~5-7 KB (3° × 2° FOV, ~80 stelle)

**PDF finale:**
- Dimensione: ~250-300 KB
- Qualità: Vettoriale (no perdita qualità)
- Stampa: Perfetto A4 (210mm × 297mm)

---

## 🔄 Workflow Completo

```
Input: OWOccultation
    ↓
ItalOccultFormatter::generateSheet()
    ↓
├─→ generateApproachChart()
│       ↓
│   calculateAsteroidTrack(10 days)
│       ↓
│   IOC_StarMap::approach_chart preset
│       ↓
│   Query Gaia DR3 (6° × 4°)
│       ↓
│   SVG embedded
│
├─→ generateEarthMapURL()
│       ↓
│   IOC_Earth API (già implementato)
│
└─→ generateFinalChart()
        ↓
    IOC_StarMap::final_chart preset
        ↓
    Query Gaia DR3 (3° × 2°)
        ↓
    SVG embedded + uncertainty ellipse
```

---

## 📝 Checklist Implementazione

- [ ] Clonare IOC_StarMap in `external/`
- [ ] Aggiornare `CMakeLists.txt`
- [ ] Aggiungere `#include "ioc_starmap/generator.h"`
- [ ] Implementare `generateApproachChart()`
- [ ] Implementare `generateFinalChart()`
- [ ] Implementare `calculateAsteroidTrack()`
- [ ] Aggiornare `generateSheet()` con SVG embedded
- [ ] Compilare: `./build.sh`
- [ ] Test: `./build/ioccultcalc --preset preset_italoccult.oop`
- [ ] Verificare SVG: `grep -c '<svg' output.html`
- [ ] Generare PDF: `python3 convert_to_pdf.py`
- [ ] Validare output: `open output.pdf`

---

## 📚 Riferimenti

- **Spec Completa**: `IOC_STARMAP_INTEGRATION.md`
- **Guida Rapida**: `IOC_STARMAP_QUICKSTART.md`
- **Script Update**: `IOC_STARMAP_CODE_UPDATE.sh`
- **Esempio Codice**: `examples/italoccult_starmap_example.cpp`
- **IOC_Earth API**: `IOC_EARTH_API_SPEC.md`
- **ItalOccult Format**: `ITALOCCULT_INTEGRATION_GUIDE.md`

---

## 🎯 Prossimi Step

1. **Implementazione**: Seguire `IOC_STARMAP_CODE_UPDATE.sh`
2. **Test**: Compilare esempio `italoccult_starmap_example.cpp`
3. **Integrazione**: Aggiornare `ItalOccultFormatter`
4. **Validazione**: Generare PDF test con dati reali
5. **Produzione**: Deploy con preset ItalOccult

---

**Tutto pronto per l'integrazione! 🚀**

La libreria IOC_StarMap fornisce preset ottimizzati per mappe celesti professionali.
Basta usare `approach_chart` e `final_chart` per generare SVG di qualità nelle schede ItalOccult.

---

**Autore**: Michele Bigi (@manvalan)  
**Data**: 3 Dicembre 2025  
**Versione**: 1.0
