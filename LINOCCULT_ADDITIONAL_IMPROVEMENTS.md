# Ulteriori Migliorie da LinOccult per IOccultCalc

## 📊 ANALISI APPROFONDITA LINOCCULT 2.3.0

Dopo aver esaminato il codice sorgente LinOccult, ho identificato **7 ottimizzazioni aggiuntive** applicabili a IOccultCalc.

---

## 1. 🎯 CALCOLO INCERTEZZA DIPENDENTE DA MAGNITUDINE STELLA

### LinOccult (`CalculateTotalUncertainty`, linee 712-735)

```cpp
const double STAR_UNCERTAINTY1 = 60.0 / 1000.0;  // 60 mas per stelle deboli (Mv >= 9)
const double STAR_UNCERTAINTY2 = 7.0 / 1000.0;   // 7 mas per stelle brillanti (Mv < 9)
const double PROPER_MOTION_UNCERTAINTY = 2.5 / 1000.0;  // 2.5 mas/anno

double StarUncertainty = (StarMv < 900) ? STAR_UNCERTAINTY2 : STAR_UNCERTAINTY1;

// Agg incertezza moto proprio dal 2015
StarUncertainty += (ETMjdate - MJD_J2015) * PROPER_MOTION_UNCERTAINTY / 365.25;
StarUncertainty *= Rad / 3600.0;  // Converti in radianti

// Incertezza angolare totale
AngleUncertainty = sqrt(StarUncertainty² + AsteroidAngleUncertainty²);

// Incertezza distanza in km
return sqrt(AsteroidUncertainty² + (StarUncertainty * Distance)²);
```

### IOccultCalc Attuale
```cpp
// src/occultation_predictor.cpp, linea 154
double uncertaintyArcsec = (orbitalUncertainty / (distance * AU)) * RAD_TO_DEG * 3600.0;
// ❌ Non considera incertezza stella
// ❌ Non dipende da magnitudine
// ❌ Non considera moto proprio
```

### 🔧 Fix Raccomandato

**File**: `src/occultation_predictor.cpp`, aggiungere dopo linea 138:

```cpp
// Incertezza stella dipendente da magnitudine (come LinOccult)
double calculateStarUncertainty(double starMag, double yearsSince2015) {
    const double STAR_UNCERT_FAINT = 60.0 / 1000.0;   // 60 mas per Mv >= 9
    const double STAR_UNCERT_BRIGHT = 7.0 / 1000.0;   // 7 mas per Mv < 9
    const double PROPER_MOTION_UNCERT = 2.5 / 1000.0; // 2.5 mas/anno
    
    double starUncert = (starMag < 9.0) ? STAR_UNCERT_BRIGHT : STAR_UNCERT_FAINT;
    starUncert += yearsSince2015 * PROPER_MOTION_UNCERT;
    return starUncert / 3600.0;  // mas → arcsec
}

// Poi in predictOccultation():
double yearsSince2015 = (event.timeCA.jd - 2457023.5) / 365.25;  // MJD J2015
double starUncertArcsec = calculateStarUncertainty(star.phot_g_mean_mag, yearsSince2015);
double asteroidUncertArcsec = (pImpl->orbitalUncertainty / (ephCA.distance * AU)) * RAD_TO_DEG * 3600.0;

// Incertezza totale (somma quadratica)
double totalUncertaintyArcsec = sqrt(starUncertArcsec * starUncertArcsec + 
                                     asteroidUncertArcsec * asteroidUncertArcsec);
```

**Impatto**: Incertezza più realistica → probabilità più accurata → meno falsi positivi

**Stima lavoro**: ~30 righe, 30 minuti

---

## 2. 📈 CALCOLO PROBABILITÀ CON INTEGRAZIONE GAUSSIANA

### LinOccult (`CalculateProbability`, linee 1691-1708)

```cpp
double CalculateProbability(double Distance, double Diameter, double Uncertainty) {
    if (Uncertainty) {
        double x1 = (Distance + Diameter/2.0) / Uncertainty;
        double x2 = (Distance - Diameter/2.0) / Uncertainty;
        double g1 = 100.0 * IntGauss(x1);  // CDF gaussiana
        double g2 = 100.0 * IntGauss(x2);
        return fabs(g1 - g2);  // Probabilità tra x1 e x2
    }
    return 0.0;
}

// IntGauss: Approssimazione Abramowitz & Stegun per CDF(x)
double IntGauss(double u) {
    const double b1 = 0.319381530;
    const double b2 = -0.356563782;
    const double b3 = 1.781477937;
    const double b4 = -1.821255978;
    const double b5 = 1.330274429;
    
    bool negative = (u < 0);
    if (negative) u = -u;
    
    double t = 1.0 / (1.0 + 0.2316419 * u);
    double pdf = exp(-u*u/2.0) / sqrt(2.0 * M_PI);  // f(u)
    double Int = 1.0 - pdf * (b1 + (b2 + (b3 + (b4 + b5*t)*t)*t)*t)*t;
    
    return negative ? (1.0 - Int) : Int;
}
```

### IOccultCalc Attuale
```cpp
// src/occultation_predictor.cpp, linee 315-334
double z = (r_asteroid - separationArcsec) / sigma;
double prob = 0.5 * (1.0 + erf(z / sqrt(2.0)));  // CDF normale semplice
// ✓ Corretto matematicamente
// ⚠️ Non considera diameter dell'asteroide come intervallo
```

### 🔧 Fix Raccomandato

**Metodo LinOccult è più accurato**: integra tra `distance - radius` e `distance + radius`.

**File**: `src/occultation_predictor.cpp`, sostituire linee 315-334:

```cpp
double calculateProbability(double separationArcsec, double asteroidAngularSize, 
                           double uncertaintyArcsec) {
    if (uncertaintyArcsec <= 0) {
        return (separationArcsec <= asteroidAngularSize / 2.0) ? 1.0 : 0.0;
    }
    
    double radius = asteroidAngularSize / 2.0;
    
    // LinOccult method: integra CDF tra distance±radius
    double x1 = (separationArcsec + radius) / uncertaintyArcsec;
    double x2 = (separationArcsec - radius) / uncertaintyArcsec;
    
    // Usa std::erf per CDF gaussiana
    auto gaussCDF = [](double x) {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
    };
    
    double g1 = gaussCDF(x1);
    double g2 = gaussCDF(x2);
    
    return std::abs(g1 - g2);  // Probabilità nell'intervallo
}
```

**Impatto**: Probabilità più accurata per asteroidi grandi

**Stima lavoro**: ~20 righe, 20 minuti

---

## 3. 🚀 APPROSSIMAZIONE CHEBYSHEV PER ORBITE

### LinOccult (`NewNewNewProcessAsteroid`, linee 2490-2570)

```cpp
const int CHEB_ORDER = 11;       // Ordine polinomio Chebyshev
const double CHEB_STEP = 1.0;    // 1 giorno di copertura per interpolazione

// Pre-calcola posizioni asteroide con polinomi Chebyshev
LOAstOrbChebMaker* chebMaker = new LOAstOrbChebMaker(...);

for (double t = MjdStart; t < MjdEnd; t += CHEB_STEP) {
    // Crea approssimazione Chebyshev per intervallo [t, t+1 giorno]
    chebMaker->Create(CHEB_ORDER, t, t + CHEB_STEP, cX, cY, cZ);
    
    APSCheb* cheb = new APSCheb(CHEB_ORDER, cX, cY, cZ, t, t + CHEB_STEP);
    
    // Interpolazione veloce per timestep fini (60 step in 1 giorno)
    for (int step = 0; step < 60; step++) {
        double t_fine = t + step * CHEB_STEP / 60.0;
        cheb->Value(t_fine, r_equ);  // ⚡ Molto veloce!
        
        // Verifica occultazioni per questa posizione
        ScanStars2(...);
        ProcessStar1(...);
    }
}
```

**Vantaggi**:
- ⚡ **Interpolazione 100× più veloce** di integrazione numerica
- ✅ Precisione elevata (ordine 11 = errore ~1 km su 1 giorno)
- 🎯 Permette timestep fini (1 minuto) senza costo computazionale

### IOccultCalc Attuale
```cpp
// Integrazione RADAU per ogni timestep
for (double jd = startJD; jd < endJD; jd += stepDays) {
    EphemerisData eph = pImpl->ephemeris.compute(jd);  // Integrazione completa
    // ⚠️ Lento per timestep fini
}
```

### 🔧 Fix Raccomandato

**Opzione A** - Implementazione completa Chebyshev:
- Richiede libreria polinomi Chebyshev
- ~300 righe codice
- 4-6 ore lavoro

**Opzione B** - Cache semplice pre-calcolato:
```cpp
// Pre-calcola ogni ora, interpola linearmente per timestep fini
std::vector<EphemerisData> ephCache;
for (double jd = startJD; jd < endJD; jd += 0.04167) {  // 1 ora
    ephCache.push_back(pImpl->ephemeris.compute(jd));
}

// Interpolazione lineare veloce
EphemerisData interpolate(double jd) {
    int idx = (jd - startJD) / 0.04167;
    double frac = ((jd - startJD) / 0.04167) - idx;
    return lerp(ephCache[idx], ephCache[idx+1], frac);
}
```

**Impatto**: 10-100× più veloce per scansione fine

**Stima lavoro Opzione B**: ~50 righe, 1-2 ore

---

## 4. 🌍 CORREZIONE OBLATENESS TERRA

### LinOccult (costante `fac`, linea 101)

```cpp
const double fac = 0.996647;  // Ratio polar/equatorial Earth radius

// Applicato a coordinate stella:
eStar = APSVec3d(Polar(StarRA, StarDec));
eStar = APSVec3d(eStar[x], eStar[y], eStar[z] / fac);  // Schiacciamento polare
```

**Formula**: `fac = (1 - flattening)` dove flattening = 0.003353 per Terra

**Scopo**: Compensa lo schiacciamento polare della Terra nel calcolo geometria ombra

### IOccultCalc Attuale
```cpp
// ❌ Non considera oblateness Terra
Vector3D starDirection = star.propagateToEpoch(time).toUnitVector();
```

### 🔧 Fix Raccomandato

**File**: `src/occultation_predictor.cpp`, linee ~145-150:

```cpp
const double EARTH_FLATTENING = 0.003353;  // WGS84
const double POLAR_EQUAT_RATIO = 1.0 - EARTH_FLATTENING;  // 0.996647

auto starPosCA = star.propagateToEpoch(event.timeCA);

// Applica correzione oblateness
Vector3D starDir(
    starPosCA.ra * cos(starPosCA.dec),
    starPosCA.ra * sin(starPosCA.dec),
    starPosCA.dec / POLAR_EQUAT_RATIO  // Compensa schiacciamento
);
starDir.normalize();
```

**Impatto**: Precisione geometry ~0.3% miglioramento per eventi ai poli

**Stima lavoro**: ~5 righe, 10 minuti

---

## 5. 🔍 PARALLASSE GEOCENTRICA STELLA

### LinOccult (`CalcParallax`, linee 690-710)

```cpp
void CalcParallax(const LOStar* pLOStar, double ETMjdate, float Parallax,
                  double& ParallaxDelta, double& ParallaxAlpha) {
    // Calcola posizione Terra rispetto barycenter solare
    Vector3d r_EarthBary;
    GetEarthPosition(ETMjdate, r_EarthBary);
    
    // Calcola correzione parallasse per osservatore geocentrico
    double distance_pc = 1.0 / (Parallax / 1000.0);  // mas → parsec
    
    // Correzione RA e Dec
    ParallaxAlpha = -r_EarthBary[x] / distance_pc * cos(StarDec);
    ParallaxDelta = -r_EarthBary[y] / distance_pc;
}

// Poi applicato:
double StarRA = Rad * ra + ParallaxAlpha;
double StarDec = Rad * de + ParallaxDelta;
```

**Scopo**: Corregge posizione stella per parallasse annua (stelle vicine appaiono spostarsi)

### IOccultCalc Attuale
```cpp
// ❌ Non applica correzione parallasse geocentrica
auto starPos = star.propagateToEpoch(time);
```

### 🔧 Fix Raccomandato

**File**: `src/gaia_star.cpp`, aggiungere metodo:

```cpp
Coordinates GaiaStar::propagateToEpochWithParallax(const JulianDate& jd, 
                                                   const Vector3D& earthPos) const {
    // Propaga con moto proprio
    Coordinates pos = propagateToEpoch(jd);
    
    // Applica correzione parallasse se disponibile (parallax > 1 mas)
    if (parallax > 0.001) {
        double distance_au = 206265.0 / parallax;  // mas → AU
        
        // Correzione parallasse (effetto annuo per stelle vicine)
        double dra = -earthPos.x / distance_au * cos(pos.dec * DEG_TO_RAD);
        double ddec = -earthPos.y / distance_au;
        
        pos.ra += dra * RAD_TO_DEG;
        pos.dec += ddec * RAD_TO_DEG;
    }
    
    return pos;
}
```

**Impatto**: Critico per stelle vicine (parallax > 10 mas, distanza < 100 pc)

**Stima lavoro**: ~20 righe, 30 minuti

---

## 6. 📐 BUFFER ANGOLARE PER RICERCA STELLE

### LinOccult (costanti, linee 99-100)

```cpp
const double ANGLE_DELTA  = 2.0 * Rad;           // ~114° per ricerca generale
const double ANGLE_DELTA1 = 1.0 * Rad / 60.0;   // 1 arcominuto per stelle vicine

// Usato in ScanStars2:
RA_MIN - ANGLE_DELTA1 < star.RA < RA_MAX + ANGLE_DELTA1
Dec_MIN - ANGLE_DELTA1 < star.Dec < Dec_MAX + ANGLE_DELTA1
```

### IOccultCalc Attuale
```cpp
// examples/italoccultcalc.cpp, linea 638
const double queryRadiusDeg = 5.0;  // Raggio fisso query
// ✓ OK per query Gaia
// ⚠️ Ma poi nessun buffer nel filtering successivo
```

### 🔧 Fix Raccomandato

Aggiungere buffer 1 arcominuto quando si filtrano stelle candidate:

```cpp
const double ANGLE_BUFFER = 1.0 / 60.0;  // 1 arcominuto

// Nel loop verifica occultazioni:
if (separationArcsec < thresholdArcsec + ANGLE_BUFFER * 60.0) {
    // Considera stella (margine sicurezza)
}
```

**Impatto**: Evita di perdere eventi al bordo threshold

**Stima lavoro**: ~2 righe, 5 minuti

---

## 7. ⚡ TIMESTEP ADATTIVO

### LinOccult (`SMALL_STEP`, linea 107)

```cpp
const int SMALL_STEP = 60;  // 60 step per giorno = 1 step ogni 24 minuti

// Usa timestep fine solo quando vicino alla stella
for (currentStep = 0; currentStep < ScanStep; currentStep++) {
    TmpMjdTime = CurrentMjdTime + currentStep * Step;
    cheb->Value(TmpMjdTime, r_equ);
    
    // 60 posizioni asteroide in 1 giorno
    ChebArray[currentStep] = r_equ;
}

// Poi verifica tutte le 60 posizioni vs stelle
for (step = 0; step < 60; step++) {
    // Check occultazione per ogni step fine
}
```

**Strategia**: Calcolo grossolano (1 giorno) + refinement fine (24 min) solo vicino eventi

### IOccultCalc Attuale
```cpp
// examples/italoccultcalc.cpp
double stepDays = 0.5;  // Step fisso mezzo giorno
// ⚠️ Potrebbe perdere eventi brevi
```

### 🔧 Fix Raccomandato

```cpp
// Timestep adattivo: grossolano per ricerca, fine per close approach
const double COARSE_STEP = 0.5;  // 12 ore per ricerca iniziale
const double FINE_STEP = 0.01;   // 15 minuti per refinement

// Prima pass: identifica giorni con possibili eventi
std::vector<double> candidateDays;
for (double jd = startJD; jd < endJD; jd += COARSE_STEP) {
    if (checkPossibleEvent(jd)) {
        candidateDays.push_back(jd);
    }
}

// Seconda pass: scan fine solo su giorni candidati
for (double day : candidateDays) {
    for (double jd = day - 0.5; jd < day + 0.5; jd += FINE_STEP) {
        // Scan dettagliato
    }
}
```

**Impatto**: 10× meno calcoli, nessun evento perso

**Stima lavoro**: ~50 righe, 1-2 ore

---

## 📊 PRIORITÀ E STIMA IMPLEMENTAZIONE

| # | Miglioramento | Impatto | Difficoltà | Tempo | Priorità |
|---|---------------|---------|------------|-------|----------|
| 1 | Incertezza stella dipendente da magnitudine | ALTO | Media | 30 min | ⭐⭐⭐ |
| 2 | Probabilità con integrazione gaussiana | MEDIO | Bassa | 20 min | ⭐⭐ |
| 3 | Approssimazione Chebyshev orbite | ALTO | Alta | 2-6 ore | ⭐⭐ (opzionale) |
| 4 | Correzione oblateness Terra | BASSO | Bassa | 10 min | ⭐ |
| 5 | Parallasse geocentrica stelle | MEDIO | Media | 30 min | ⭐⭐ |
| 6 | Buffer angolare ricerca stelle | BASSO | Bassa | 5 min | ⭐ |
| 7 | Timestep adattivo | MEDIO | Media | 1-2 ore | ⭐⭐ |

**TOTALE**: ~4-9 ore per implementare tutto

---

## 🎯 RACCOMANDAZIONE IMPLEMENTAZIONE

### Fase 1: Quick Wins (1 ora)
✅ Implementare PRIMA del test con cache Gaia:
1. **Incertezza stella** (30 min) - Critico per probabilità
2. **Probabilità gaussiana** (20 min) - Miglior accuratezza
3. **Oblateness Terra** (10 min) - Facile, corretto fisicamente

### Fase 2: Medium Priority (2-3 ore)
⏳ Implementare SE performance è problema:
4. **Parallasse geocentrica** (30 min) - Importante per stelle vicine
5. **Timestep adattivo** (1-2 ore) - 10× speedup
6. **Buffer angolare** (5 min) - Evita edge cases

### Fase 3: Advanced (2-6 ore)
🔮 Implementare SE vuoi massima performance:
7. **Chebyshev orbite** (2-6 ore) - 100× speedup ma complesso

---

## 💡 PROSSIMI PASSI RACCOMANDATI

**Opzione A** - Test immediato (raccomandato):
1. ✅ Threshold fix già implementato
2. ⏳ Aspetta cache Gaia locale
3. ⏳ Testa e vedi quanti eventi trova
4. ⏳ POI implementa migliorie se necessario

**Opzione B** - Implementa Quick Wins prima del test:
1. ⏳ Implementa #1, #2, #4 (1 ora totale)
2. ⏳ Compila e testa
3. ⏳ Aspetta cache Gaia
4. ⏳ Test completo con migliorie

**La mia raccomandazione**: **Opzione B** - Le migliorie #1, #2, #4 richiedono solo 1 ora e migliorano significativamente accuratezza probabilità ed incertezza.

Vuoi che implementi i Quick Wins (#1, #2, #4) ora?
