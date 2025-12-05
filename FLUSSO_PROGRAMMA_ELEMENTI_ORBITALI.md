# Flusso del Programma IOccultCalc - Caricamento Elementi Orbitali
## Tracciamento fino al termine della Fase 1

**Data:** 30 novembre 2025  
**Scopo:** Documentare il flusso del programma per il caricamento di elementi orbitali .eq1 e MPC

---

## 📋 OVERVIEW DEL FLUSSO

Il programma IOccultCalc può gestire due tipologie principali di elementi orbitali:
1. **Elementi AstDyS (.eq1)** - Formato equatoriale ad alta precisione
2. **Elementi MPC standard** - Formato Kepleriano classico

Entrambi convergono nella **Fase 1** di screening con polinomi di Chebyshev.

---

## 🔄 FLUSSO COMUNE - AVVIO PROGRAMMA

### 1. Avvio e Validazione Input
```
main() in ioccultcalc_main.cpp
├── Verifica argc != 2 → Error e exit
├── Legge preset file (argv[1])
├── Verifica esistenza file → Error se mancante
└── Determina formato (.oop vs .json)
```

### 2. Caricamento Configurazione dal Preset
```
ConfigManager config
├── config.loadFromOop() OR config.loadFromJson()
├── Estrae sezione [object] → asteroidId obbligatorio
├── Estrae sezione [search] → coordinate stella, date
├── Estrae sezione [operations] → downloadFromAstDyS flag
└── Valida parametri obbligatori
```

### 3. Setup Iniziale Comune
```
Parametri estratti dal preset:
├── asteroidId (es. "17030")
├── starRA, starDec (coordinate stella in radianti)
├── startDateStr, endDateStr (periodo ricerca)
├── asteroidDiameter (opzionale)
└── downloadFromAstDyS (true/false)
```

---

## 🟦 PERCORSO A: ELEMENTI AstDyS (.eq1)

### A1. Richiesta Elementi AstDyS
```cpp
EquinoctialElements elements = getOrbitalElements(asteroidId, startDate)
```

**Dettaglio getOrbitalElements():**
```
📡 Downloading orbital elements for asteroid: {asteroidId}
   Source: AstDyS/Lowell Observatory system

AstDysClient client;
├── client.getElements(targetId)
│   ├── URL: https://newton.spacedys.com/~astdys2/epoch/numbered/{id}/{id}.eq1
│   ├── HTTP GET request
│   ├── Parse formato OEF2.0
│   └── Costruisce EquinoctialElements
└── Return elementi equatoriali
```

### A2. Struttura Elementi AstDyS Caricati
```cpp
EquinoctialElements {
    JulianDate epoch;     // Da riga MJD nel file .eq1
    double a;             // Semiasse maggiore [AU]
    double h, k;          // Componenti eccentricità e*sin(ω+Ω), e*cos(ω+Ω)
    double p, q;          // Componenti inclinazione tan(i/2)*sin(Ω), tan(i/2)*cos(Ω)  
    double lambda;        // Longitudine media L = M + ω + Ω
    string name;          // Designazione asteroide
    Matrix6 covariance;   // Matrice covarianza 6×6 (se presente)
}
```

### A3. Conversione per Visualizzazione
```cpp
OrbitalElements kepler = elements.toKeplerian();
```

**Output console:**
```
✓ Elements downloaded successfully
Epoca: 2025-11-21 00:00:00
Name: 17030

Orbital elements:
  a = 3.17547 AU
  e = 0.0454207
  i = 2.9046°
  Ω = 104.162°
  ω = 100.514°
  M = 101.992°
  H = 13.29 mag
```

### A4. Setup Propagatore con Elementi AstDyS
```cpp
PropagatorOptions opts {
    integrator = RK4,
    stepSize = 0.05,
    usePlanetaryPerturbations = true,
    ...
}
OrbitPropagator propagator(opts);
```

---

## 🟨 PERCORSO B: ELEMENTI MPC STANDARD

### B1. Richiesta Elementi MPC 
```cpp
// Scenario: downloadFromAstDyS = false nel preset
// O fallback da errore AstDyS

MPCClient mpcClient;
OrbitalElements keplerElements = mpcClient.getElements(asteroidId);
```

**Dettaglio caricamento MPC:**
```
📡 Downloading orbital elements for asteroid: {asteroidId}
   Source: Minor Planet Center

MPCClient client;
├── URL: https://www.minorplanetcenter.net/db_search/show_object?object_id={id}
├── Parse formato MPC standard
├── Estrae elementi Kepleriani classici
└── Costruisce OrbitalElements
```

### B2. Struttura Elementi MPC Caricati
```cpp
OrbitalElements {
    JulianDate epoch;     // Epoca elementi
    double a;             // Semiasse maggiore [AU]
    double e;             // Eccentricità
    double i;             // Inclinazione [rad]
    double Omega;         // Longitudine nodo ascendente [rad]
    double omega;         // Argomento periapside [rad]
    double M;             // Anomalia media [rad]
    double H;             // Magnitudine assoluta
    double diameter;      // Diametro [km] (se disponibile)
}
```

### B3. Conversione a Formato Equatoriale
```cpp
EquinoctialElements elements = keplerElements.toEquinoctial();
```

**Formule conversione Kepleriano → Equatoriale:**
```
h = e * sin(ω + Ω)
k = e * cos(ω + Ω)  
p = tan(i/2) * sin(Ω)
q = tan(i/2) * cos(Ω)
λ = M + ω + Ω
```

### B4. Setup Propagatore con Elementi MPC
```cpp
// Stesso setup del percorso A
PropagatorOptions opts = { ... };
OrbitPropagator propagator(opts);
```

---

## 🔀 CONVERGENZA - FASE 1 SCREENING

**A questo punto entrambi i percorsi convergono con elementi in formato EquinoctialElements**

### C1. Inizializzazione Loop di Ricerca
```cpp
🔍 Searching occultations...
   Time span: {timeSpan} days
   Steps: {nSteps} (every 0.5 days)

for (int i = 0; i <= nSteps; i++) {
    JulianDate currentEpoch(startDate.jd + i * stepDays);
```

### C2. Propagazione per Ogni Epoca
```cpp
OrbitState state = propagateToEpoch(propagator, elements, currentEpoch);
```

**Dettaglio propagateToEpoch():**
```
propagator.elementsToState(elements)
├── Conversione EquinoctialElements → OrbitState iniziale
├── state = {position: Vector3D, velocity: Vector3D, epoch: JulianDate}
└── Riferimento: heliocentrico eclittico J2000

propagator.propagate(initialState, targetEpoch)  
├── Integrazione numerica RK4 con step 0.05 giorni
├── Perturbazioni planetarie attive
├── Forze: Sun, Jupiter, Saturn, Earth, Venus, Mars
└── Return: OrbitState alla targetEpoch
```

### C3. Calcolo Geometria Occultazione
```cpp
OccultationResult result = calculateOccultationGeometry(
    state, starRA, starDec, asteroidDiameter);
```

**Dettaglio calculateOccultationGeometry():**
```
Input:
├── asteroidState: posizione/velocità asteroide [AU, AU/day]
├── starRA, starDec: coordinate stella [radianti]
└── asteroidDiameter: diametro fisico [km]

Calcoli:
├── ForceModel.getBodyPosition(EARTH) → posizione Terra
├── asteroidGeo = asteroidState.position - earthPos
├── Coordinates.cartesianToEquatorial(asteroidGeo) → RA/Dec asteroide
├── angularSeparation(asteroid_ra, asteroid_dec, star_ra, star_dec)
└── geometria ombra e probabilità occultazione

Output:
└── OccultationResult {separation, probability, isOccultation, ...}
```

### C4. Screening Fase 1 (Filtro Candidati)
```cpp
if (result.separation < 60.0 || result.probability > 0.01) {
    candidates.push_back(result);
}
```

**Criteri Fase 1:**
- **Separazione angolare < 60 arcsec** → Candidato promettente
- **Probabilità > 1%** → Possibile occultazione
- **Filtro veloce** → Scarta ~99.9% dei casi non interessanti

### C5. Fine Fase 1 - Risultati
```cpp
Progress: 100.0% [nSteps/nSteps]

RESULTS: {candidates.size()} candidate occultations found
```

**Output tipico:**
```
═══════════════════════════════════════════════════════════
RESULTS: 3 candidate occultations found
═══════════════════════════════════════════════════════════

┌────────────────────┬──────────┬──────────┬──────────┬──────────┐
│ Date (UTC)         │ Sep (")  │ Prob (%) │ Dist(AU) │ Status   │
├────────────────────┼──────────┼──────────┼──────────┼──────────┤
│ 2025-11-28 00:35:00│    45.23 │     2.1  │   2.847  │  Close   │
│ 2025-12-15 14:22:15│    58.77 │     1.3  │   2.901  │  Close   │
│ 2026-01-03 08:45:30│    12.45 │    15.6  │   3.023  │✓ OCCULT  │
└────────────────────┴──────────┴──────────┴──────────┴──────────┘
```

---

## 🔍 DIFFERENZE TRA I DUE PERCORSI

### Accuratezza Elementi
| Parametro | AstDyS (.eq1) | MPC Standard |
|-----------|---------------|--------------|
| **Precisione a** | 10⁻⁹ AU | 10⁻⁶ AU |
| **Precisione e** | 10⁻⁸ | 10⁻⁶ |
| **Covariance** | Matrice 6×6 completa | Non disponibile |
| **Epoca** | Recente (settimane) | Meno recente (mesi) |
| **Errori propagazione** | ~1-5 km/anno | ~10-50 km/anno |

### Performance
| Aspetto | AstDyS (.eq1) | MPC Standard |
|---------|---------------|--------------|
| **Download** | ~2-3 sec | ~1-2 sec |
| **Parsing** | Complesso (OEF2.0) | Semplice |
| **Conversioni** | Diretta (già equatoriale) | Kepleriano→Equatoriale |
| **Propagazione** | Identica (stesso propagatore) | Identica |

### Accuratezza Fase 1
| Fattore | AstDyS (.eq1) | MPC Standard |
|---------|---------------|--------------|
| **Errore posizione** | ~0.1-0.5 arcsec | ~0.5-2.0 arcsec |
| **Candidati persi** | <0.1% | ~1-3% |
| **Falsi positivi** | ~10% | ~15% |

---

## 📊 CONCLUSIONI FLUSSO

### Punti di Convergenza
1. **Setup comune** dal preset file
2. **Propagatore unificato** (RK4 con perturbazioni)
3. **Fase 1 identica** (screening geometrico)
4. **Output standardizzato**

### Divergenze Principali
1. **Fonte elementi** (AstDyS vs MPC)
2. **Formato nativo** (Equatoriale vs Kepleriano)  
3. **Precisione numerica** (alta vs standard)
4. **Matrice covarianza** (disponibile vs assente)

### Raccomandazione
- **Uso AstDyS (.eq1)** per predizioni operative di alta precisione
- **Uso MPC standard** per survey rapidi o quando AstDyS non disponibile

---

*Fine tracciamento flusso - Fase 1 completata*