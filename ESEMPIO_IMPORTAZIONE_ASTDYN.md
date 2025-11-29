# Guida Importazione File AstDyS con ITALOccultLibrary/astdyn

## Panoramica

ITALOccultLibrary/astdyn fornisce parser completi per i file di AstDyS:
- **File .eq1**: Elementi orbitali in formato equinoziale (OEF2.0)
- **File .rwo**: Osservazioni ottiche in formato OrbFit

## File .eq1: Elementi Orbitali Equinoziali

### Formato OEF2.0

```
format  = 'OEF2.0'       ! Versione formato
rectype = 'ML'           ! Multi-line record
refsys  = ECLM J2000     ! Sistema di riferimento eclittico J2000
END_OF_HEADER
11234                    ! Designazione oggetto
 EQU   2.6808535916678031E+00   0.032872036471001   0.036254405825130
       0.103391596538937   -0.042907901689093 235.8395861037268
 MJD     61000.000000000 TDT
 MAG  12.874  0.150
```

### Elementi Equinoziali

Gli **elementi equinoziali** evitano singolarità per orbite circolari e planari:

| Simbolo | Nome | Formula | Descrizione |
|---------|------|---------|-------------|
| **a** | Semiasse maggiore | a | In AU |
| **h** | Componente H | e·sin(ϖ) | ϖ = ω + Ω (longitudine perielio) |
| **k** | Componente K | e·cos(ϖ) | |
| **p** | Componente P | tan(i/2)·sin(Ω) | Ω = longitudine nodo ascendente |
| **q** | Componente Q | tan(i/2)·cos(Ω) | |
| **λ** | Longitudine media | M + ω + Ω | In gradi |

### Conversione a Elementi Kepleriani

```cpp
// Eccentricità
e = sqrt(h² + k²)

// Inclinazione
tan(i/2) = sqrt(p² + q²)
i = 2·atan(sqrt(p² + q²))

// Longitudine nodo ascendente
Ω = atan2(p, q)
if (Ω < 0) Ω += 2π

// Longitudine perielio
ϖ = atan2(h, k)
if (ϖ < 0) ϖ += 2π

// Argomento perielio
ω = ϖ - Ω
if (ω < 0) ω += 2π

// Anomalia media
M = λ - ϖ
// Normalizza M in [0, 2π)
```

### Vantaggi Elementi Equinoziali

1. **No singolarità** per e → 0 (orbite circolari)
2. **No singolarità** per i → 0 (orbite planari)
3. **Migliore convergenza** nei least squares
4. **Usati da OrbFit** e altri software professionali

## Parser EQ1Parser (astdyn)

### Uso Base

```cpp
#include <astdyn/io/EQ1Parser.hpp>

// Parse file .eq1
auto elem = astdyn::io::EQ1Parser::parse("11234.eq1");

std::cout << "Object: " << elem.object_name << "\n";
std::cout << "Epoch: MJD " << elem.epoch_mjd_tdb << " TDB\n";
std::cout << "a = " << elem.a << " AU\n";
std::cout << "h = " << elem.h << "\n";
std::cout << "k = " << elem.k << "\n";
std::cout << "p = " << elem.p << "\n";
std::cout << "q = " << elem.q << "\n";
std::cout << "λ = " << elem.lambda << "°\n";

// Converti a Kepleriani
double a, e, i, Omega, omega, M;
astdyn::io::EQ1Parser::equinoctial_to_keplerian(
    elem, a, e, i, Omega, omega, M
);

// Nota: angoli sono in RADIANTI
std::cout << "e = " << e << "\n";
std::cout << "i = " << i * 180.0/M_PI << "°\n";
std::cout << "Ω = " << Omega * 180.0/M_PI << "°\n";
std::cout << "ω = " << omega * 180.0/M_PI << "°\n";
std::cout << "M = " << M * 180.0/M_PI << "°\n";
```

### Struttura Dati

```cpp
struct EquinoctialElements {
    std::string object_name;    // Designazione (es. "11234")
    double epoch_mjd_tdb;       // Epoca MJD TDB
    double a;                   // Semiasse maggiore (AU)
    double h;                   // e·sin(ϖ)
    double k;                   // e·cos(ϖ)
    double p;                   // tan(i/2)·sin(Ω)
    double q;                   // tan(i/2)·cos(Ω)
    double lambda;              // Longitudine media (gradi)
    
    // Opzionali
    double magnitude;           // Magnitudine assoluta H
    double mag_slope;           // Slope parameter G
};
```

## File .rwo: Osservazioni Ottiche

### Formato RWO

Il formato .rwo è **complesso** e usa posizioni fisse in stile Fortran:

```
version =   2
errmod  = 'fcct14'
RMSast  =   5.17357E-01
END_OF_HEADER
! Object K T N YYYY MM DD.ddddd ... HH MM SS.sss ... sDD MM SS.ss ... Mag ... Cod
11234  O A   1989 01 02.85929 ... 06 09 29.830 ... +17 25 56.40 ... 16.6 ... 046
```

### Campi Principali

| Campo | Colonne | Formato | Descrizione |
|-------|---------|---------|-------------|
| Object | 2-10 | A9 | Designazione oggetto |
| Type | 12 | A1 | 'O' = optical, 'R' = radar, 'S' = satellite |
| Note | 14 | A1 | 'A' = asteroid, 'C' = comet, 'P' = planet |
| Tech | 16 | A1 | Tecnica osservativa |
| Date | 18-38 | I4,1X,I2,1X,F13.10 | YYYY MM DD.ddddddddd |
| RA | ~50-62 | I2,1X,I2,1X,F6.3 | HH MM SS.sss |
| Dec | ~104-117 | A1,I2,1X,I2,1X,F5.2 | sDD MM SS.ss |
| Mag | ~120-125 | F5.2 | Magnitudine |
| Obs Code | Fine | A3 | Codice osservatorio MPC |

**Nota**: Le posizioni esatte variano! Il formato usa **FORTRAN READ** con posizioni fisse.

## Parser RWOReader (astdyn)

### Uso Base

```cpp
#include <astdyn/observations/RWOReader.hpp>

// Parse file .rwo
astdyn::observations::RWOReader reader;
auto observations = reader.readFile("11234.rwo");

std::cout << "Loaded " << observations.size() << " observations\n";

for (const auto& obs : observations) {
    std::cout << "MJD " << obs.mjd_utc 
              << " RA " << obs.ra * 180.0/M_PI << "° "
              << "Dec " << obs.dec * 180.0/M_PI << "° "
              << "Mag " << obs.magnitude.value_or(0.0) << " "
              << "Obs " << obs.observatory_code << "\n";
}
```

### Struttura OpticalObservation

```cpp
struct OpticalObservation {
    std::string object_designation;
    double mjd_utc;              // Tempo osservazione (MJD UTC)
    double ra;                   // Ascensione retta (radianti)
    double dec;                  // Declinazione (radianti)
    double sigma_ra;             // Incertezza RA (radianti)
    double sigma_dec;            // Incertezza Dec (radianti)
    std::optional<double> magnitude;  // Magnitudine
    std::string observatory_code;     // Codice MPC (es. "500", "046")
    
    // Altre info opzionali...
};
```

### Lettura Raggruppata

```cpp
// Raggruppa osservazioni per oggetto
auto grouped = reader.readFileGrouped("multiple_objects.rwo");

for (const auto& [object, obs_set] : grouped) {
    std::cout << "Object " << object << ": " 
              << obs_set.observations.size() << " observations\n";
    std::cout << "  Time span: " << obs_set.time_span_days() << " days\n";
}
```

## Esempio Completo: Orbit Fitting

### Workflow Tipico

```cpp
#include <astdyn/io/EQ1Parser.hpp>
#include <astdyn/observations/RWOReader.hpp>
#include <astdyn/AstDynEngine.hpp>

int main() {
    // 1. Parse elementi iniziali
    auto equ = astdyn::io::EQ1Parser::parse("203.eq1");
    
    // 2. Converti a Kepleriani
    double a, e, i, Omega, omega, M;
    astdyn::io::EQ1Parser::equinoctial_to_keplerian(
        equ, a, e, i, Omega, omega, M
    );
    
    // Crea struttura orbit
    astdyn::orbit::KeplerianOrbit initial_orbit;
    initial_orbit.epoch_mjd_tdb = equ.epoch_mjd_tdb;
    initial_orbit.semi_major_axis = a;
    initial_orbit.eccentricity = e;
    initial_orbit.inclination = i;
    initial_orbit.longitude_ascending_node = Omega;
    initial_orbit.argument_perihelion = omega;
    initial_orbit.mean_anomaly = M;
    
    // 3. Carica osservazioni
    astdyn::observations::RWOReader reader;
    auto observations = reader.readFile("203.rwo");
    
    std::cout << "Loaded " << observations.size() << " observations\n";
    
    // 4. Setup engine
    astdyn::AstDynEngine engine("config.oop");
    engine.set_observations(observations);
    
    // 5. Propaga a epoca media osservazioni
    double mean_epoch = engine.get_mean_observation_epoch();
    auto state_at_mean = engine.propagate_orbit(initial_orbit, mean_epoch);
    engine.set_initial_state(state_at_mean);
    
    // 6. Differential correction
    astdyn::orbit_determination::DifferentialCorrectorSettings settings;
    settings.max_iterations = 20;
    settings.convergence_tolerance = 1e-6;
    settings.reject_outliers = true;
    settings.outlier_sigma = 3.0;
    
    auto result = engine.run_differential_correction(settings);
    
    // 7. Risultati
    std::cout << "Converged: " << result.converged << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "RMS: " << result.rms_arcsec << " arcsec\n";
    std::cout << "Observations used: " << result.num_observations_used 
              << "/" << observations.size() << "\n";
    
    // Orbit fitted
    auto fitted = result.final_orbit;
    std::cout << "\nFitted orbit:\n";
    std::cout << "  a = " << fitted.semi_major_axis << " AU\n";
    std::cout << "  e = " << fitted.eccentricity << "\n";
    std::cout << "  i = " << fitted.inclination * 180.0/M_PI << "°\n";
    
    return 0;
}
```

## Problemi Comuni

### 1. Parser .rwo Fallisce

**Problema**: Il formato .rwo è ESTREMAMENTE dipendente da posizioni fisse Fortran.

**Soluzioni**:
1. Usa `RWOReader` di astdyn (già implementato)
2. Se il file ha varianti, adatta il parser
3. Considera scarica osservazioni da MPC (formato più standard)

### 2. Elementi Non Converge

**Problema**: Elementi equinoziali da .eq1 non convergono nel fitting.

**Possibili cause**:
- Epoca troppo lontana dalle osservazioni
- Propagatore non include perturbazioni
- Osservazioni hanno outlier

**Soluzioni**:
```cpp
// Propaga sempre a epoca media
double mean_epoch = (first_obs + last_obs) / 2.0;
auto state = propagate(orbit, mean_epoch);

// Abilita reject outliers
settings.reject_outliers = true;
settings.outlier_sigma = 3.0;  // 3σ

// Usa RKF78 + perturbazioni
config.propagator = "RKF78";
config.use_planetary_perturbations = true;
```

### 3. Coordinate Errate

**Problema**: RA/Dec sembrano sbagliate.

**Cause comuni**:
- Confusione radianti/gradi
- Sistema di riferimento (J2000 vs altro)
- UTC vs TDB

**Soluzione**:
```cpp
// astdyn usa RADIANTI internamente
double ra_deg = obs.ra * 180.0 / M_PI;
double dec_deg = obs.dec * 180.0 / M_PI;

// Epoca è sempre MJD TDB per elementi
// Osservazioni sono MJD UTC
```

## File di Esempio

### Test con (203) Pompeja

```bash
cd external/ITALOccultLibrary/astdyn

# Elementi
tools/203_astdys_latest.eq1

# Osservazioni (varie versioni)
tools/203_astdys_latest.rwo        # Tutte osservazioni
tools/203_recent100.rwo            # Ultime 100
tools/203_subset.rwo               # Subset test

# Esegui esempio
./build/examples/simple_pompeja_fit \
    tools/203_recent100.rwo \
    tools/203_astdys_latest.eq1 \
    tools/203.oop
```

### Test con (11234) 1999 JS82

```bash
# Elementi
data/11234.eq1

# Osservazioni
data/11234.rwo

# Parse solo elementi
./build/tests/test_eq1_parser data/11234.eq1
```

## Compilazione Esempi astdyn

```bash
cd external/ITALOccultLibrary/astdyn
mkdir build && cd build

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_CSPICE=OFF \
    -DBUILD_TESTS=ON \
    -DBUILD_EXAMPLES=ON

make -j8

# Test parser
./tests/test_eq1_parser ../data/11234.eq1

# Esempio fitting
./examples/simple_pompeja_fit \
    ../tools/203_recent100.rwo \
    ../tools/203_astdys_latest.eq1 \
    ../tools/203.oop
```

## Integrazione con IOccultCalc

### Opzione 1: Usa Elementi Già Fittati ✅ RACCOMANDATO

AstDyS fornisce elementi **già fittati** con migliaia di osservazioni:

```bash
# 1. Scarica .eq1
wget https://newton.spacedys.com/astdys/data/11234.eq1

# 2. Parse con EQ1Parser
./build/tests/test_eq1_parser 11234.eq1

# 3. Copia elementi nel preset IOccultCalc
# preset.oop:
general.
    .asteroid_number = 11234
    .epoch_mjd = 61000.0
    .semimajor_axis_au = 2.680854
    .eccentricity = 0.04893825
    .inclination_deg = 12.7744
    .ascending_node_deg = 112.5386
    .perihelion_arg_deg = 289.6601
    .mean_anomaly_deg = 193.6408
    
    .propagator = 'RKF78'
    
propagator.
    .use_planetary_perturbations = .TRUE.

# 4. Esegui IOccultCalc
./italoccultcalc preset.oop
```

### Opzione 2: Fitting Completo (se necessario)

Solo se elementi AstDyS non disponibili o vuoi migliorarli:

```cpp
// Richiede linking con ITALOccultLibrary/astdyn
#include <astdyn/AstDynEngine.hpp>

// ... usa esempio simple_pompeja_fit.cpp
```

**Pro**: Elementi custom ottimizzati
**Contro**: Complessità integrazione, richiede osservazioni

## Link Utili

- **AstDyS**: https://newton.spacedys.com/astdys/
- **OrbFit Manual**: http://adams.dm.unipi.it/orbfit/
- **MPC Observatory Codes**: https://minorplanetcenter.net/iau/lists/ObsCodes.html
- **JPL Horizons**: https://ssd.jpl.nasa.gov/horizons/

## Conclusioni

### ✅ Per IOccultCalc (produzione)

1. **Usa elementi da AstDyS .eq1** (già fittati)
2. **Parse con EQ1Parser** di astdyn
3. **Converti a Kepleriani**
4. **Usa RKF78 + perturbazioni**

### ⚠️ Orbit Fitting Custom

Solo se necessario:
- Compila astdyn con tutte dipendenze
- Usa `AstDynEngine` completo
- Richiede osservazioni da .rwo o MPC

### 📚 Documentazione Codice

- `astdyn/include/astdyn/io/EQ1Parser.hpp` - Parser .eq1
- `astdyn/src/observations/RWOReader.cpp` - Parser .rwo
- `astdyn/examples/simple_pompeja_fit.cpp` - Esempio completo
- `astdyn/tests/test_eq1_parser.cpp` - Test parser

---

**Autore**: IOccultCalc Development Team  
**Data**: 2025-11-29  
**Basato su**: ITALOccultLibrary/astdyn esempi e documentazione
