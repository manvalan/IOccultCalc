# Phase 2: Geometria Precisa Occultazione

## ✅ IMPLEMENTATO

### Architettura Completa
- **Header**: `include/phase2_occultation_geometry.h`
- **Implementazione**: `src/phase2_occultation_geometry.cpp`
- **Compilazione**: ✅ Successo

### Funzionalità Implementate

#### 1. Download Osservazioni RWO da AstDyS
```cpp
std::vector<astdyn::observations::OpticalObservation> downloadRWOObservations(designation);
```
- **URL AstDyS**: `https://newton.spacedys.com/~astdys2/mpcobs/numbered/XX/XXXXX.rwo`
- **Parser**: `astdyn::observations::RWOReader::readFile()`
- **Protocollo**: HTTPS via libcurl
- **Caching**: File temporaneo in `/tmp/`

#### 2. Propagazione Precisa
- **Propagatore**: `astdyn::propagation::Propagator` con RKF78
- **Perturbazioni**: Tutti i pianeti attivi (Giove, Saturno, etc.)
- **Tolleranza**: 1e-12 (massima precisione)
- **Step**: 0.1 giorni per integrazione
- **Coordinate**: Eclittica J2000 → ICRF equatoriale

#### 3. Calcolo Geometria
- **Input**: Lista `CandidateStar` da Phase 1
- **Propagazione densa**: ±5 minuti attorno al CA, step 1 secondo
- **Closest approach**: Ricerca sub-millisecondo del minimo
- **Parametri calcolati**:
  - Distanza minima (milliarcsec)
  - Istante esatto CA
  - Durata massima occultazione
  - Velocità angolare
  - Distanza asteroide (AU)

#### 4. Strutture Dati

**OccultationEvent**: Evento completo con geometria
- Identificazione stella e asteroide
- Geometria globale (CA, durata, chord)
- Shadow path sulla Terra (TODO)
- Predizioni per osservatori (TODO)
- Ellisse incertezza (TODO)

**Phase2Config**: Configurazione completa
- Orbital refinement settings
- Download RWO parameters
- Fitting parameters (TODO: quando OrbitFitter disponibile)
- Propagation settings
- Output options

**Phase2Results**: Risultati completi
- Lista eventi calcolati
- Statistiche fit orbitale
- Performance metrics

### Workflow Completo

```
Phase1CandidateScreening
    ↓
    [Lista candidati con CA approssimati]
    ↓
Phase2OccultationGeometry::calculateGeometry()
    ↓
    ┌─────────────────────────────────────┐
    │ STEP 1: Orbital Refinement         │
    │ - Download RWO da AstDyS            │
    │ - Parse con RWOReader (AstDyn)      │
    │ - Fit orbitale preciso (TODO)      │
    │ - Propagazione a epoca CA           │
    └──────────┬──────────────────────────┘
               ↓
    ┌─────────────────────────────────────┐
    │ STEP 2: Calcolo Geometria          │
    │ - Propagazione densa ±5 min         │
    │ - Trova CA esatto (sub-ms)          │
    │ - Calcola parametri evento          │
    │ - Geometria osservatori (TODO)      │
    └──────────┬──────────────────────────┘
               ↓
    [OccultationEvent con geometria completa]
```

## 🚧 TODO

### 1. OrbitFitter (Alta Priorità)
```cpp
// Da implementare in AstDyn:
astdyn::fitting::OrbitFitter
    ├── fit(ObservationSet) → FitResult
    ├── setInitialOrbit(KeplerianElements)
    ├── setSettings(OrbitFitterSettings)
    └── Algoritmo: Differential Corrections o Gauss-Newton
```

**Requisiti**:
- Least squares fit con osservazioni RWO
- Matrice covarianza 6x6 per incertezze
- Iterazioni fino a convergenza
- RMS residui < 0.5 arcsec tipicamente

### 2. Shadow Path su Terra
```cpp
std::vector<ShadowPathPoint> computeShadowPath(event, num_points);
```
- Proiezione ombra su ellissoide WGS84
- Entry/exit points
- Velocità e direzione ombra
- Limiti nord/sud

### 3. Observer Predictions
```cpp
ObserverGeometry computeForObserver(latitude, longitude, elevation);
```
- CA tempo locale
- Altitudine/azimut target
- Distanza dalla linea centrale
- Visibilità (Sole, Luna)
- Durata prevista

### 4. Uncertainty Ellipse
```cpp
UncertaintyEllipse computeUncertainty(covariance_matrix);
```
- Da matrice covarianza fit
- Semi-assi ellisse 1-sigma
- Position angle
- Incertezza temporale

### 5. Correzioni Astrometriche
- ❌ Parallasse geocentrica
- ❌ Aberrazione stellare/planetaria
- ❌ Proper motion stella (Gaia → epoca evento)
- ❌ Light-time correction
- ❌ Deflessione gravitazionale (Sole)

## 📊 Performance Target

| Operazione | Tempo Atteso |
|-----------|--------------|
| Download RWO | < 2 sec |
| Parse RWO (1000 obs) | < 50 ms |
| Orbital fit | 5-30 sec |
| Propagazione densa (600 pts) | < 100 ms |
| CA optimization | < 10 ms |
| **TOTALE per candidato** | **< 40 sec** |

## 🔧 Uso

```cpp
#include "phase1_candidate_screening.h"
#include "phase2_occultation_geometry.h"

// Phase 1: Screening veloce
Phase1CandidateScreening phase1;
phase1.loadAsteroidFromEQ1("17030_astdys.eq1");
auto candidates = phase1.screenCandidates(config1);

// Phase 2: Geometria precisa
Phase2OccultationGeometry phase2;
phase2.loadAsteroidFromEQ1("17030_astdys.eq1");

Phase2Config config2;
config2.refine_orbit_from_observations = true;
config2.mpc_code = "17030";  // ← Scarica RWO automaticamente
config2.time_window_minutes = 5.0;
config2.time_step_seconds = 1.0;

auto results = phase2.calculateGeometry(candidates.candidates, config2);

// Risultati
for (const auto& event : results.events) {
    std::cout << "Stella " << event.star_source_id << "\n";
    std::cout << "  CA: " << event.closest_approach_mas << " mas\n";
    std::cout << "  Time: MJD " << event.time_ca_mjd_utc << "\n";
    std::cout << "  Duration: " << event.max_duration_sec << " sec\n";
}
```

## 📝 Note Implementative

### Download RWO
- **Numbered asteroids**: `/numbered/17/17030.rwo` (17 = 17030 / 1000)
- **Unnumbered**: `/unnumbered/DESIGNATION.rwo`
- **Fallback**: Se download fallisce, usa elementi nominali

### Fitting (quando implementato)
- **Algoritmo**: Differential corrections o Gauss-Newton
- **Pesi**: Da RWO (peso RA, peso Dec)
- **Outlier rejection**: 3-sigma clipping automatico
- **Convergenza**: χ² variation < 0.1% o max iterations

### Precision Levels
- **Phase 1**: ~1 arcsec (screening)
- **Phase 2 nominal**: ~0.1 arcsec (elementi nominali)
- **Phase 2 refined**: ~0.01 arcsec (con fit RWO)

## ✅ Status Finale

- [x] Architettura completa
- [x] Download RWO da AstDyS
- [x] Parsing con AstDyn
- [x] Propagazione precisa
- [x] Calcolo CA esatto
- [x] Strutture dati complete
- [ ] Orbital fitting (dipende da AstDyn)
- [ ] Shadow path su Terra
- [ ] Observer predictions
- [ ] Correzioni astrometriche complete

**Prossimo Step**: Implementare OrbitFitter in AstDyn, poi integrare in Phase 2.
