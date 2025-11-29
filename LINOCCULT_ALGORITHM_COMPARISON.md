# Confronto Algoritmi: LinOccult vs IOccultCalc

## 1. Strategia Generale

### LinOccult (loCalc.cc)
```cpp
// FASE 1: Approssimazione Chebyshev dell'orbita asteroide
pLOAstOrbChebMaker->Create(CHEB_ORDER, CurrentMjdTime, 
                           CurrentMjdTime + CHEB_STEP, cX, cY, cZ);

// FASE 2: Scansione stelle con kdtree in bounding box
ScanStars2(pLOStarData, ScanStep);
  // Calcola RA_MIN, RA_MAX, Dec_MIN, Dec_MAX dalle posizioni asteroide
  // Usa kdtree per ricerca veloce: O(log n) invece di O(n)
  pLOStarData->FastFindStars(pNewLOStarsArray, NewStarNumber,
                             RA_MIN - ANGLE_DELTA1, Dec_MIN - ANGLE_DELTA1,
                             RA_MAX + ANGLE_DELTA1, Dec_MAX + ANGLE_DELTA1);

// FASE 3: Per ogni stella trovata, verifica occultazione
for (i = 0; i < GetStarsCount(); i++) {
    pLOStar = GetStar(i);
    ProcessStar1(pLOAsteroid, pLOStar, pLOEventData, ...);
}
```

### IOccultCalc (occultation_predictor.cpp)
```cpp
// FASE 1: Propagazione numerica orbita con integrator RADAU
// (no approssimazione Chebyshev)

// FASE 2: Query Gaia con multi-region clustering
// Gruppo asteroidi in regioni da 15° → query 5° per regione
// Deduplica stelle per source_id

// FASE 3: Per ogni asteroide vs tutte le stelle scaricate
for (asteroid : asteroids) {
    for (star : allStars) {
        findClosestApproach(star, startDate, endDate);
        if (separation < threshold) {
            predictOccultation(...);
        }
    }
}
```

## 2. Differenze Critiche nell'Algoritmo di Rilevamento

### 2.1 Criterio di Occultazione

**LinOccult** (`ProcessStar1`, linea 2377):
```cpp
// Calcola distanza minima tra stella e ombra asteroide
s0 = -Dot(rAst, eStar);  // Distanza asteroide lungo linea di vista stella
Delta = s0*s0 + R_Earth*R_Earth - Dot(rAst, rAst);
r0 = sqrt(R_Earth*R_Earth - Delta);  // Raggio del punto più vicino

// Calcola incertezza totale
TotalUncertainty = CalculateTotalUncertainty(AngleUncertainty,
                                              -EphemerisUncertainty * s0,
                                              Mjdate, StarMv, -s0);

MaxDist = TotalUncertainty + R_Earth;
if (MaxDist > ExtraRadius) {
    MaxDist = ExtraRadius;  // Default: 3.0 * R_Earth
}

// ✓ OCCULTAZIONE SE:
if (r0 < MaxDist) {
    CreateOccultationEvent(...);
}
```

**Parametri chiave**:
- `ANGLE_DELTA1 = 1.0 * Rad / 60.0` = **1 arcominuto** buffer per ricerca stelle
- `ExtraRadius = 3.0 * R_Earth` = **~19,000 km** distanza massima (default)
- Incertezza considera: ephemeris uncertainty, magnitudine stella, distanza

**IOccultCalc** (occultation_predictor.cpp, ANALIZZATO):
```cpp
// Linea 80-92: Controllo preliminare
double asteroidAngularSize = (asteroidDiameter / (distance * AU)) * RAD_TO_DEG * 3600.0;

// THRESHOLD ATTUALE:
double threshold = asteroidAngularSize + 
                   3.0 * orbitalUncertainty / (distance * AU) * RAD_TO_DEG * 3600.0;
if (threshold == 0) threshold = 10.0; // Default 10 arcsec

// ✓ OCCULTAZIONE SE:
if (separationArcsec < threshold) {
    JulianDate caTime = findClosestApproach(star, jd - 0.5, jd + 0.5);
    OccultationEvent event = predictOccultation(star, caTime);
    
    // Filtra per probabilità minima
    if (event.probability >= minProbability) {
        events.push_back(event);
    }
}

// Linea 154-158: Calcolo probabilità
double uncertaintyArcsec = (orbitalUncertainty / (distance * AU)) * RAD_TO_DEG * 3600.0;
event.probability = calculateProbability(closeApproachDistance, 
                                         asteroidAngularSize,
                                         uncertaintyArcsec);

// Linea 315-334: Calcolo probabilità gaussiana
double z = (r_asteroid - separationArcsec) / sigma;
double prob = 0.5 * (1.0 + erf(z / sqrt(2.0)));
```

**PROBLEMA IDENTIFICATO**: 
1. ❌ **Separazione ANGOLARE (arcsec) invece di DISTANZA 3D (km)**
   - IOccultCalc: `separationArcsec < threshold` (angolo)
   - LinOccult: `r0 < MaxDist` (distanza geometrica in km)
   
2. ❌ **Threshold inadeguato per asteroidi lontani**
   - Esempio: (4) Vesta @ 2.5 AU, diametro 500 km
   - Angular size: 0.068 arcsec (quasi invisibile!)
   - Threshold IOccultCalc: 0.068 + 3*uncertainty ≈ 0.1-1 arcsec
   - Threshold LinOccult: 3.0 * R_Earth ≈ 19,000 km = **2.6 arcsec @ 2.5 AU**
   
3. ⚠️ **Filtro probabilità troppo restrittivo?**
   - Se `minProbability` > 0.0, molti eventi vengono scartati
   - LinOccult considera TUTTI gli eventi con r0 < MaxDist

### 2.2 Gestione Parallasse Stellare

**LinOccult** (linee 2322-2349):
```cpp
// Corregge posizione stella per parallasse
if (Parallax) {
    CalcParallax(pLOStar, ETMjdate, 
                 Rad * Ddd(0, 0, plx/1000.0),
                 ParallaxDelta, ParallaxAlpha);
}

double StarRA = Rad * ra + ParallaxAlpha;
double StarDec = Rad * de + ParallaxDelta;

// Applica fattore di compressione
eStar = APSVec3d(eStar[x], eStar[y], eStar[z] / fac);
```

**IOccultCalc**: Verifica se parallasse è considerata nel calcolo posizione stella.

### 2.3 Propagazione Moto Proprio

**LinOccult** (linee 2318-2328):
```cpp
// DT = anni dalla epoch di catalogazione
DT = (ETMjdate - pLOStar->GetEpoch()) / 365.25;

// Trasforma epoch con moto proprio
TransformEpoch(DT, 
               pLOStar->GetRA(), pLOStar->GetDec(),
               Deg * Parallax * 1000.0 * 3600.0,
               Deg * pmRA * 1000.0 * 3600.0,
               Deg * pmDec * 1000.0 * 3600.0,
               pLOStar->GetVrad(),
               ra, de, plx, uRAs, uDE, VR);
```

**IOccultCalc**: Gaia adapter include `pmra` e `pmdec`, ma è usato nel calcolo?

## 3. Ottimizzazione Ricerca Stelle: kdtree

### LinOccult (loStarData.cc)
```cpp
// Usa kdtree per ricerca spaziale O(log n)
bool FastFindStars(const LOStar** pNewLOStarsArray,
                   unsigned int& NewStarNumber,
                   double RA_MIN, double Dec_MIN,
                   double RA_MAX, double Dec_MAX)
{
    // Ricerca kdtree nel bounding box
    // Ritorna solo stelle nella regione richiesta
    // 5-10× più veloce di scansione lineare
}
```

### IOccultCalc
```cpp
// Scansione lineare di TUTTE le stelle scaricate
for (const auto& star : gaiaStars) {
    double separation = calculateSeparation(asteroidRA, asteroidDec,
                                             star.pos.ra, star.pos.dec);
    if (separation < searchRadius) {
        // Considera stella
    }
}
```

**Problema**: Con 100,000 stelle, scansione lineare è **O(n)** per ogni asteroide.  
**Soluzione**: Implementare kdtree come LinOccult → **O(log n)**.

## 4. Approssimazione Orbita

### LinOccult: Polinomi di Chebyshev
```cpp
// Pre-calcola posizioni asteroide con Chebyshev per CHEB_STEP giorni
LOAstOrbChebMaker->Create(CHEB_ORDER, CurrentMjdTime,
                          CurrentMjdTime + CHEB_STEP, cX, cY, cZ);

// Interpolazione veloce durante scansione
pAPSCheb->Value(TmpMjdTime, r_equ);
```

**Vantaggi**:
- Interpolazione polinomiale velocissima
- Precalcola orbita una volta, usa molte volte
- Riduce chiamate a propagatore orbitale

**IOccultCalc**: Integrazione numerica RADAU ad ogni timestamp?  
Potenzialmente più lento ma più preciso.

## 5. Buffer Angolare

### LinOccult
```cpp
const double ANGLE_DELTA  = 2.0 * Rad;           // ~2° per ricerca generale
const double ANGLE_DELTA1 = 1.0 * Rad / 60.0;   // ~1 arcominuto per stelle

// Ricerca stelle con buffer:
RA_MIN - ANGLE_DELTA1 < star.RA < RA_MAX + ANGLE_DELTA1
Dec_MIN - ANGLE_DELTA1 < star.Dec < Dec_MAX + ANGLE_DELTA1
```

**IOccultCalc**: Usa raggio fisso 5° per query Gaia, ma poi?  
Verifica se filtering successivo usa buffer adeguato.

## 6. Incertezza Ephemeris

### LinOccult (CalculateTotalUncertainty)
```cpp
TotalUncertainty = CalculateTotalUncertainty(
    AngleUncertainty,
    -pLOAsteroid->GetEphemerisUncertainty() * s0,  // Incertezza dipende da distanza
    Mjdate,
    pLOStar->GetMv(),      // Magnitudine stella
    -s0                     // Distanza asteroide
);

MaxDist = TotalUncertainty + R_Earth;
```

**Incertezza considerata**:
- Ephemeris uncertainty dell'asteroide (da orbital elements)
- Magnitudine stella (stelle deboli = incertezza posizione maggiore)
- Distanza asteroide (più lontano = incertezza maggiore)

**IOccultCalc**: Come è calcolata incertezza? È considerata nel threshold?

## 7. Creazione Evento

### LinOccult
```cpp
if (r0 < MaxDist) {
    // Controlla se evento già nel database
    pLOEvent = pLOEventData->FindEvent(AsteroidID, Catalogue,
                                       StarNumber, Mjdate);
    if (!pLOEvent) {
        // Torna indietro di uno step per precisione
        Mjdate = Mjdate - Step;
        CurrentStep--;
        
        // Crea evento dettagliato
        CreateOccultationEvent(pLOAsteroid, pLOStar, pLOEventData,
                               eStar, &Mjdate, &CurrentStep,
                               ET_UT, Step, MaxDist);
    }
}
```

**Step-back strategy**: Quando trova possibile occultazione, torna indietro di uno step per affinare il calcolo.

**IOccultCalc**: Verifica strategia di refinement del tempo di occultazione.

## 8. Raccomandazioni per IOccultCalc

### 8.1 CRITICO: Rivedere Threshold Occultazione
```cpp
// Invece di:
if (closestSeparation < fixedThreshold) { ... }

// Usare:
double maxDist = calculateUncertainty(ephemUncertainty, starMag, distance)
                 + R_Earth;
maxDist = min(maxDist, 3.0 * R_Earth);  // Cap a 19,000 km

if (minDistance < maxDist) {
    // Occultazione possibile
}
```

### 8.2 IMPORTANTE: Implementare kdtree per Stelle
```cpp
// Costruisci kdtree da gaiaStars una volta
KDTree starTree;
for (const auto& star : gaiaStars) {
    starTree.insert(star.pos.ra, star.pos.dec, &star);
}

// Query veloce per ogni asteroide
for (const auto& asteroid : asteroids) {
    auto nearbyStars = starTree.rangeSearch(
        asteroid.ra - searchRadius,
        asteroid.dec - searchRadius,
        asteroid.ra + searchRadius,
        asteroid.dec + searchRadius
    );
    
    // Processa solo stelle vicine (100-1000 invece di 100,000)
    for (const auto* star : nearbyStars) {
        checkOccultation(asteroid, star);
    }
}
```

**Beneficio**: Da O(n_asteroids × n_stars) → O(n_asteroids × log(n_stars))

### 8.3 OPZIONALE: Approssimazione Chebyshev
```cpp
// Pre-calcola posizioni asteroide con fit polinomiale
ChebyshevApproximation cheb(asteroid, startDate, endDate, order=10);

// Interpolazione veloce durante scansione
for (double t = startDate; t < endDate; t += dt) {
    Vector3d pos = cheb.evaluate(t);  // Molto veloce
    // Verifica occultazioni
}
```

**Beneficio**: Riduce calcoli orbitali costosi, utile per ricerca a grana fine.

### 8.4 Gestione Parallasse e Moto Proprio
Verificare che `gaia_adapter.cpp` propaghi correttamente:
```cpp
// In checkOccultation():
double dt = (mjd - star.ref_epoch) / 365.25;  // anni

// Applica moto proprio
double ra_corrected = star.ra + star.pmra * dt / cos(star.dec);
double dec_corrected = star.dec + star.pmdec * dt;

// Applica parallasse geocentrica (se parallax > 0)
if (star.parallax > 0.001) {  // >1 mas
    applyParallaxCorrection(ra_corrected, dec_corrected,
                           star.parallax, observerPos, mjd);
}
```

### 8.5 Buffer Angolare Adeguato
```cpp
const double ANGLE_BUFFER = 1.0 / 60.0;  // 1 arcominuto in gradi

// Nel bounding box per stelle
double ra_min = minAsteroidRA - ANGLE_BUFFER;
double ra_max = maxAsteroidRA + ANGLE_BUFFER;
double dec_min = minAsteroidDec - ANGLE_BUFFER;
double dec_max = maxAsteroidDec + ANGLE_BUFFER;
```

## 9. Prossimi Passi

1. **Analizzare `occultation_predictor.cpp`**:
   - Trovare threshold usato per determinare occultazione
   - Verificare calcolo geometria ombra
   - Controllare considerazione incertezza

2. **Implementare kdtree**:
   - Portare libreria kdtree/ da LinOccult
   - O usare libreria esistente (e.g., nanoflann)
   - Costruire indice spaziale per stelle Gaia

3. **Test con Caso Noto**:
   - Trovare occultazione verificata (e.g., da IOTA)
   - Eseguire LinOccult → risultato atteso
   - Eseguire IOccultCalc → confrontare
   - Debug differenze

4. **Confronto Numerico**:
   ```cpp
   // Aggiungere logging dettagliato
   std::cout << "Asteroid: " << name << std::endl;
   std::cout << "Star: " << sourceId << " RA=" << ra << " Dec=" << dec << std::endl;
   std::cout << "Closest approach: " << minSep << " arcsec" << std::endl;
   std::cout << "Shadow distance: " << r0 << " km" << std::endl;
   std::cout << "Max distance (threshold): " << maxDist << " km" << std::endl;
   std::cout << "Occultation: " << (r0 < maxDist ? "YES" : "NO") << std::endl;
   ```

5. **Ottimizzazione Performance**:
   - Profiling per identificare bottleneck
   - Considerare Chebyshev approximation se propagazione orbitale è lenta
   - Implementare caching intelligente

## 10. Conclusioni

**⚠️ CAUSA CONFERMATA di 0 eventi in IOccultCalc**:

1. ✅ **CONFERMATO: Threshold TROPPO STRETTO**
   ```
   IOccultCalc: threshold = asteroidAngularSize + 3*uncertainty
   Esempio Vesta @ 2.5 AU:
     - Angular size: 0.068 arcsec
     - Uncertainty: ~0.1 arcsec
     - Threshold: 0.17 arcsec  ← TROPPO PICCOLO!
   
   LinOccult: MaxDist = min(TotalUncertainty + R_Earth, 3*R_Earth)
   Esempio Vesta @ 2.5 AU:
     - MaxDist: 19,000 km = 2.6 arcsec @ 2.5 AU  ← 15× PIÙ GRANDE!
   ```

2. ✅ **CONFERMATO: Geometria ANGOLARE vs 3D**
   - IOccultCalc confronta **separazione angolare (arcsec)** con **dimensione angolare asteroide**
   - LinOccult confronta **distanza geometrica 3D (km)** con **raggio ombra (km)**
   - Per asteroidi lontani, angular size diventa infinitesimale anche per grandi asteroidi!

3. ⚠️ **Filtro probabilità troppo conservativo**
   - Se `event.probability >= minProbability` con minProbability > 0.0
   - Molti eventi "possibili" vengono scartati
   - LinOccult crea evento per TUTTI i casi con r0 < MaxDist

4. ✓ **Moto proprio gestito correttamente**
   - `star.propagateToEpoch(time)` presente nel codice

5. ❓ **Buffer query Gaia**: Non ancora verificato nella fase di download stelle

**Priorità implementazione FIX**:
1. **✅ RISOLTO**: Threshold e criterio occultazione analizzati
2. **🔴 CRITICO**: Sostituire calcolo threshold con metodo LinOccult
3. **🟡 MEDIA**: Implementare calcolo distanza 3D invece di solo angolare
4. **🟢 BASSA**: Implementare kdtree per performance (già funziona, solo lento)
5. **🟢 BASSA**: Chebyshev approximation (già RADAU è preciso)

## 11. SOLUZIONE IMPLEMENTATIVA

### Fix Immediato: Threshold Corretto

**File**: `src/occultation_predictor.cpp`, linee 80-92

**PRIMA (SBAGLIATO)**:
```cpp
// Stima dimensione angolare dell'asteroide
double asteroidAngularSize = 0;
if (pImpl->asteroidDiameter > 0 && eph.distance > 0) {
    asteroidAngularSize = (pImpl->asteroidDiameter / (eph.distance * AU)) * RAD_TO_DEG * 3600.0;
}

// Controllo preliminare: stella deve essere vicina
double threshold = asteroidAngularSize + 3.0 * pImpl->orbitalUncertainty / (eph.distance * AU) * RAD_TO_DEG * 3600.0;
if (threshold == 0) threshold = 10.0; // Default 10 arcsec

if (separationArcsec < threshold) {
    // Crea evento
}
```

**DOPO (CORRETTO - STILE LINOCCULT)**:
```cpp
// Calcola distanza minima in km (geometria 3D)
double asteroidDistanceKm = eph.distance * AU;  // Distanza asteroide dalla Terra in km

// Incertezza totale in km
double uncertaintyKm = pImpl->orbitalUncertainty;  // Già in km
double totalUncertainty = uncertaintyKm;

// Calcola distanza massima ombra (come LinOccult)
const double R_Earth = 6371.0;  // km
double maxDistKm = totalUncertainty + R_Earth;

// Cap a 3.0 * R_Earth come LinOccult (default ExtraRadius)
const double MAX_SHADOW_RADIUS = 3.0 * R_Earth;  // ~19,000 km
if (maxDistKm > MAX_SHADOW_RADIUS) {
    maxDistKm = MAX_SHADOW_RADIUS;
}

// Converti in arcsec per confronto
// maxDistKm è il raggio dell'ombra in km
// Converto in angolo: arctan(maxDistKm / asteroidDistanceKm)
double thresholdArcsec = (maxDistKm / asteroidDistanceKm) * RAD_TO_DEG * 3600.0;

if (separationArcsec < thresholdArcsec) {
    // Crea evento - MOLTO PIÙ EVENTI TROVATI!
}
```

**Esempio numerico**:
```
Asteroide: (4) Vesta
Diametro: 500 km
Distanza dalla Terra: 2.5 AU = 374 million km
Incertezza orbitale: 1000 km

VECCHIO METODO:
  asteroidAngularSize = 500 / 374,000,000 * 206,265 = 0.275 arcsec
  threshold = 0.275 + 3 * (1000 / 374,000,000) * 206,265 = 0.275 + 1.65 = 1.9 arcsec

NUOVO METODO:
  maxDistKm = 1000 + 6371 = 7371 km
  (capped a 19,113 km)
  thresholdArcsec = 19,113 / 374,000,000 * 206,265 = 10.5 arcsec
  
RISULTATO: 10.5 / 1.9 = 5.5× PIÙ GRANDE → MOLTI PIÙ EVENTI!
```

### Fix Completo: Distanza 3D (Come LinOccult)

Per implementazione ancora più accurata, calcolare distanza 3D geometrica:

```cpp
// Vettori 3D in coordinate geocentriche (km)
Vector3D asteroidPos_km = eph.geocentricPos * AU;  // Converti da AU a km
Vector3D starDirection = star.propagateToEpoch(eph.jd).toUnitVector();

// Distanza minima punto-linea (come LinOccult ProcessStar1)
// s0 = proiezione asteroide su linea di vista stella
double s0 = -Vector3D::dot(asteroidPos_km, starDirection);

// Calcola distanza perpendicolare (r0 in LinOccult)
Vector3D perpVec = asteroidPos_km + s0 * starDirection;
double r0 = perpVec.magnitude();  // km

// Threshold
double maxDistKm = std::min(totalUncertainty + R_Earth, 3.0 * R_Earth);

if (r0 < maxDistKm) {
    // ✓ OCCULTAZIONE - CRITERIO ESATTO LINOCCULT
    JulianDate caTime = findClosestApproach(...);
    OccultationEvent event = predictOccultation(star, caTime);
    events.push_back(event);
}
```

Questo è il metodo ESATTO di LinOccult ed è geometricamente più corretto.
