# IOccultCalc - Guida alla Configurazione Completa

**Versione**: 2.1.0-rkf78  
**Data**: 29 Novembre 2025  
**Autore**: IOccultCalc Development Team

---

## Indice

1. [Panoramica](#panoramica)
2. [Configurazione Propagatore Orbitale](#configurazione-propagatore-orbitale)
3. [Configurazione Modello di Forze](#configurazione-modello-di-forze)
4. [Configurazione Predictor Occultazioni](#configurazione-predictor-occultazioni)
5. [Configurazione Database Asteroidi](#configurazione-database-asteroidi)
6. [Configurazione Cataloghi Stellari](#configurazione-cataloghi-stellari)
7. [Configurazione Output](#configurazione-output)
8. [File di Configurazione JSON](#file-di-configurazione-json)
9. [Esempi Pratici](#esempi-pratici)
10. [Best Practices](#best-practices)

---

## Panoramica

IOccultCalc supporta configurazioni flessibili tramite:
- **Codice C++**: Oggetti configurazione diretti
- **File JSON**: Configurazione persistente
- **File OOP**: Formato compatibile Occult4
- **Variabili ambiente**: Override temporanei

### Priorità Configurazione

1. Parametri espliciti nel codice
2. File di configurazione specificato
3. File di configurazione di default (`~/.ioccultcalc/config.json`)
4. Valori di default hardcoded

---

## Configurazione Propagatore Orbitale

### Struttura: `PropagatorOptions`

```cpp
#include "ioccultcalc/orbit_propagator.h"

PropagatorOptions opts;
```

### Parametri Obbligatori

Nessuno - tutti i parametri hanno valori di default ragionevoli.

### Parametri Opzionali

#### **integrator** (IntegratorType)
- **Tipo**: Enum `IntegratorType`
- **Default**: `IntegratorType::RKF78`
- **Valori possibili**:
  - `RKF78`: Runge-Kutta-Fehlberg 7(8) adattivo **(RACCOMANDATO)**
  - `RK4`: Runge-Kutta classico 4° ordine (legacy)
  - `RA15`: Everhart RA15 (15° ordine, sperimentale)
- **Descrizione**: Metodo di integrazione numerica
- **Quando usare**:
  - `RKF78`: Default per tutte le applicazioni (273x più veloce di RK4)
  - `RK4`: Solo per compatibilità con vecchi risultati
  - `RA15`: Simulazioni molto lunghe (>100 anni)

**Esempio**:
```cpp
opts.integrator = IntegratorType::RKF78;  // Default
```

---

#### **stepSize** (double)
- **Tipo**: `double`
- **Default**: `0.1` giorni
- **Range**: `0.001` - `100.0` giorni
- **Unità**: Giorni
- **Descrizione**: Step iniziale per l'integrazione
- **Note**:
  - Per RKF78: step iniziale, poi adattato automaticamente
  - Per RK4: step fisso
  - Step più piccoli → maggiore precisione, più lento
  - Step più grandi → minore precisione, più veloce

**Esempio**:
```cpp
opts.stepSize = 0.1;  // Default, ottimo per RKF78
```

---

#### **tolerance** (double)
- **Tipo**: `double`
- **Default**: `1e-12`
- **Range**: `1e-14` - `1e-8`
- **Unità**: Adimensionale (errore relativo)
- **Descrizione**: Tolleranza per controllo errore (solo RKF78)
- **Quando usare**:
  - `1e-12`: Default, eccellente per la maggior parte dei casi
  - `1e-14`: Precisione massima, simulazioni scientifiche critiche
  - `1e-10`: Performance migliorate, precisione comunque alta
  - `1e-8`: Ricerche veloci, precisione ridotta accettabile

**Esempio**:
```cpp
opts.tolerance = 1e-12;  // Default
```

---

#### **usePlanetaryPerturbations** (bool)
- **Tipo**: `bool`
- **Default**: `true`
- **Descrizione**: Include perturbazioni dei pianeti maggiori
- **Impatto**: 
  - `true`: Precisione scientifica (RACCOMANDATO)
  - `false`: Keplero puro, solo per test o oggetti molto distanti
- **Performance**: Trascurabile con DE441 ottimizzato

**Esempio**:
```cpp
opts.usePlanetaryPerturbations = true;  // RACCOMANDATO
```

---

#### **useAsteroidalPerturbations** (bool)
- **Tipo**: `bool`
- **Default**: `true`
- **Descrizione**: Include perturbazioni dei 17 asteroidi massivi (AST17)
- **Lista AST17**:
  - (1) Ceres, (2) Pallas, (3) Juno, (4) Vesta
  - (7) Iris, (10) Hygiea, (15) Eunomia, (16) Psyche
  - (29) Amphitrite, (31) Euphrosyne, (52) Europa, (65) Cybele
  - (87) Sylvia, (88) Thisbe, (511) Davida, (704) Interamnia
  - (951) Gaspra
- **Impatto**:
  - `true`: Precisione migliorata ~5-10% per main-belt
  - `false`: Performance leggermente migliori
- **Quando disabilitare**: NEA, TNO, comete

**Esempio**:
```cpp
opts.useAsteroidalPerturbations = true;  // RACCOMANDATO per main-belt
```

---

#### **useRelativisticCorrections** (bool)
- **Tipo**: `bool`
- **Default**: `false`
- **Descrizione**: Include correzioni relativistiche (EIH - Einstein-Infeld-Hoffmann)
- **Impatto**: 
  - Miglioramento ~1-2 km per asteroidi a <2 AU
  - Overhead computazionale ~5%
- **Quando abilitare**:
  - NEA con q < 1.0 AU
  - Mercury crosser
  - Simulazioni precisione massima
- **Quando disabilitare**:
  - Main-belt standard (effetto trascurabile)
  - Performance critiche

**Esempio**:
```cpp
opts.useRelativisticCorrections = false;  // Default
// Per NEA:
opts.useRelativisticCorrections = true;
```

---

#### **ephemerisType** (EphemerisType)
- **Tipo**: Enum `EphemerisType`
- **Default**: `EphemerisType::JPL_DE441`
- **Valori possibili**:
  - `JPL_DE441`: NASA JPL DE441 (2020-2035, alta precisione) **(RACCOMANDATO)**
  - `JPL_DE440`: NASA JPL DE440 (extended range)
  - `JPL_DE430`: NASA JPL DE430 (legacy)
  - `VSOP87`: VSOP87 analitico (veloce, meno preciso)
- **Descrizione**: Sorgente effemeridi planetarie
- **Storage**: DE441 richiede ~300 MB disco

**Esempio**:
```cpp
opts.ephemerisType = EphemerisType::JPL_DE441;  // Default
```

---

#### **maxIterations** (int)
- **Tipo**: `int`
- **Default**: `10000`
- **Range**: `100` - `1000000`
- **Descrizione**: Numero massimo iterazioni per propagazione
- **Quando modificare**:
  - Aumentare per propagazioni molto lunghe (>10 anni)
  - Diminuire per timeout in ricerche massive

**Esempio**:
```cpp
opts.maxIterations = 10000;  // Default, sufficiente per la maggior parte dei casi
```

---

#### **outputInterval** (double)
- **Tipo**: `double`
- **Default**: `1.0` giorni
- **Range**: `0.001` - `365.0` giorni
- **Unità**: Giorni
- **Descrizione**: Intervallo output per traiettorie
- **Uso**: Quando si genera lista stati orbitali

**Esempio**:
```cpp
opts.outputInterval = 1.0;  // Output giornaliero
```

---

### Esempio Configurazione Completa

```cpp
#include "ioccultcalc/orbit_propagator.h"

// Configurazione alta precisione
PropagatorOptions highPrecision;
highPrecision.integrator = IntegratorType::RKF78;
highPrecision.stepSize = 0.1;
highPrecision.tolerance = 1e-14;
highPrecision.usePlanetaryPerturbations = true;
highPrecision.useAsteroidalPerturbations = true;
highPrecision.useRelativisticCorrections = true;
highPrecision.ephemerisType = EphemerisType::JPL_DE441;

OrbitPropagator propagator(highPrecision);

// Configurazione veloce per survey
PropagatorOptions fastSurvey;
fastSurvey.integrator = IntegratorType::RKF78;
fastSurvey.stepSize = 0.5;
fastSurvey.tolerance = 1e-10;
fastSurvey.usePlanetaryPerturbations = true;
fastSurvey.useAsteroidalPerturbations = false;
fastSurvey.useRelativisticCorrections = false;

OrbitPropagator fastPropagator(fastSurvey);
```

---

## Configurazione Modello di Forze

### Struttura: `ForceModelOptions`

```cpp
#include "ioccultcalc/force_model.h"

ForceModelOptions forceOpts;
```

### Parametri Disponibili

#### **includeSun** (bool)
- **Default**: `true`
- **Descrizione**: Include attrazione gravitazionale solare
- **Note**: **NON DISABILITARE** (richiesto per sistema eliocentrico)

#### **includeMoon** (bool)
- **Default**: `true`
- **Descrizione**: Include perturbazioni lunari
- **Impatto**: Significativo per NEA Earth-approaching

#### **includeMercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune** (bool)
- **Default**: `true` per tutti
- **Descrizione**: Include perturbazioni dei singoli pianeti
- **Performance vs Accuracy**:
  - **Critical**: Jupiter, Saturn (sempre on)
  - **Important**: Earth, Mars (main-belt)
  - **Minor**: Mercury, Venus, Uranus, Neptune

**Esempio configurazione selettiva**:
```cpp
forceOpts.includeMercury = false;  // Trascurabile per main-belt
forceOpts.includeVenus = false;
forceOpts.includeUranus = false;
forceOpts.includeNeptune = false;
// Jupiter, Saturn sempre on
```

#### **includeAsteroids** (bool)
- **Default**: `true`
- **Descrizione**: Include perturbazioni AST17

#### **includeRelativity** (bool)
- **Default**: `false`
- **Descrizione**: Include correzioni relativistiche EIH

---

### Configurazioni Predefinite

```cpp
// FAST: Solo giganti gassosi
ForceModelOptions fast = ForceModelOptions::fastConfig();
// include: Jupiter, Saturn

// STANDARD: Tutti i pianeti
ForceModelOptions standard = ForceModelOptions::standardConfig();
// include: tutti i pianeti + Luna

// HIGH PRECISION: Tutto
ForceModelOptions highPrec = ForceModelOptions::highPrecisionConfig();
// include: tutti i pianeti + Luna + AST17 + relativity
```

---

## Configurazione Predictor Occultazioni

### Struttura: `OccultationPredictorConfig`

```cpp
#include "ioccultcalc/occultation_predictor.h"

OccultationPredictorConfig config;
```

### Parametri Ricerca Occultazioni

#### **startDate** (JulianDate) *OBBLIGATORIO*
- **Tipo**: `JulianDate`
- **Descrizione**: Data inizio ricerca
- **Formato**: JD o MJD

**Esempio**:
```cpp
config.startDate = JulianDate(2460645.5);  // Nov 29, 2025
```

---

#### **endDate** (JulianDate) *OBBLIGATORIO*
- **Tipo**: `JulianDate`
- **Descrizione**: Data fine ricerca
- **Range consigliato**: 1 mese - 1 anno

**Esempio**:
```cpp
config.endDate = JulianDate(2460645.5 + 365);  // 1 anno
```

---

#### **observerLocation** (GeographicCoordinate) *OBBLIGATORIO*
- **Tipo**: `GeographicCoordinate`
- **Campi**:
  - `latitude`: Latitudine (gradi, -90 a +90)
  - `longitude`: Longitudine (gradi, -180 a +180)
  - `altitude`: Altitudine (metri)

**Esempio**:
```cpp
GeographicCoordinate observer;
observer.latitude = 45.4642;    // Torino
observer.longitude = 9.1900;
observer.altitude = 239.0;      // metri
config.observerLocation = observer;
```

---

#### **minAltitude** (double)
- **Tipo**: `double`
- **Default**: `20.0` gradi
- **Range**: `0.0` - `90.0` gradi
- **Descrizione**: Altitudine minima oggetto sull'orizzonte
- **Raccomandazioni**:
  - `30°`: Condizioni ottimali
  - `20°`: Standard
  - `15°`: Condizioni accettabili
  - `<15°`: Refrazione atmosferica problematica

**Esempio**:
```cpp
config.minAltitude = 20.0;  // Default
```

---

#### **maxMagnitude** (double)
- **Tipo**: `double`
- **Default**: `14.0` mag
- **Range**: `6.0` - `18.0` mag
- **Descrizione**: Magnitudine massima stella per occultazione
- **Considerazioni**:
  - Limite visuale: ~6.5 mag
  - CCD amatoriale: ~14-16 mag
  - Survey professionale: ~18 mag

**Esempio**:
```cpp
config.maxMagnitude = 14.0;  // Default CCD amatoriale
```

---

#### **minMagnitudeDrop** (double)
- **Tipo**: `double`
- **Default**: `0.5` mag
- **Range**: `0.1` - `5.0` mag
- **Descrizione**: Drop minimo magnitudine per evento significativo
- **Formula**: `Δm = -2.5 log10(1 - (D_ast/D_star)²)`
- **Raccomandazioni**:
  - `2.0 mag`: Eventi facili (>75% drop)
  - `1.0 mag`: Eventi buoni (>60% drop)
  - `0.5 mag`: Eventi accettabili (>37% drop)
  - `0.2 mag`: Eventi difficili (>18% drop)

**Esempio**:
```cpp
config.minMagnitudeDrop = 0.5;  // Default
```

---

#### **maxShadowWidth** (double)
- **Tipo**: `double`
- **Default**: `500.0` km
- **Range**: `10.0` - `10000.0` km
- **Descrizione**: Larghezza massima ombra per considerare evento
- **Uso**: Filtro efficienza (ombre grandi = più probabili ma meno precise)

**Esempio**:
```cpp
config.maxShadowWidth = 500.0;  // Default
```

---

#### **minDuration** (double)
- **Tipo**: `double`
- **Default**: `0.1` secondi
- **Range**: `0.01` - `60.0` secondi
- **Unità**: Secondi
- **Descrizione**: Durata minima occultazione
- **Considerazioni**:
  - Framerate camera
  - Tempo di integrazione
  - Eventi troppo brevi = difficili da rilevare

**Esempio**:
```cpp
config.minDuration = 0.1;  // Default (100 ms)
```

---

#### **searchRadius** (double)
- **Tipo**: `double`
- **Default**: `3.0` raggi stella
- **Range**: `1.0` - `10.0` raggi stella
- **Descrizione**: Raggio ricerca attorno a stella (in raggi stellari)
- **Note**: 
  - `1.0`: Solo eventi centrali (centrale)
  - `3.0`: Include eventi grazing (raccomandato)
  - `5.0+`: Molti falsi positivi

**Esempio**:
```cpp
config.searchRadius = 3.0;  // Default
```

---

#### **starCatalog** (string)
- **Tipo**: `string`
- **Default**: `"GAIA_DR3"`
- **Valori possibili**:
  - `"GAIA_DR3"`: Gaia Data Release 3 **(RACCOMANDATO)**
  - `"GAIA_EDR3"`: Gaia Early DR3
  - `"UCAC4"`: USNO CCD Astrograph Catalog 4
  - `"TYCHO2"`: Tycho-2
  - `"NOMAD"`: Naval Observatory Merged Astrometric Dataset
- **Descrizione**: Catalogo stellare per ricerca

**Esempio**:
```cpp
config.starCatalog = "GAIA_DR3";  // Default
```

---

#### **useParallelProcessing** (bool)
- **Tipo**: `bool`
- **Default**: `true`
- **Descrizione**: Abilita parallelizzazione OpenMP
- **Performance**: 4-8x più veloce su CPU multi-core

**Esempio**:
```cpp
config.useParallelProcessing = true;  // Default
```

---

#### **numThreads** (int)
- **Tipo**: `int`
- **Default**: `0` (auto-detect)
- **Range**: `0` - numero core CPU
- **Descrizione**: Numero thread per parallelizzazione
- **Note**: `0` = usa tutti i core disponibili

**Esempio**:
```cpp
config.numThreads = 0;  // Auto (RACCOMANDATO)
// O specifico:
config.numThreads = 4;  // 4 thread
```

---

### Esempio Configurazione Predictor

```cpp
OccultationPredictorConfig config;

// Parametri obbligatori
config.startDate = JulianDate(2460645.5);        // Nov 29, 2025
config.endDate = JulianDate(2460645.5 + 31);     // Dicembre 2025
config.observerLocation.latitude = 45.4642;       // Torino
config.observerLocation.longitude = 9.1900;
config.observerLocation.altitude = 239.0;

// Parametri opzionali
config.minAltitude = 30.0;                        // Condizioni ottime
config.maxMagnitude = 12.0;                       // Stelle brillanti
config.minMagnitudeDrop = 1.0;                    // Eventi significativi
config.maxShadowWidth = 300.0;                    // Ombre compatte
config.minDuration = 0.2;                         // >200ms
config.searchRadius = 3.0;                        // Include grazing
config.starCatalog = "GAIA_DR3";                  // Catalogo migliore
config.useParallelProcessing = true;              // Parallelo
config.numThreads = 0;                            // Auto

OccultationPredictor predictor(config);
```

---

## Configurazione Database Asteroidi

### Struttura: `AsteroidDatabaseConfig`

```cpp
#include "ioccultcalc/asteroid_database.h"

AsteroidDatabaseConfig dbConfig;
```

### Parametri

#### **databasePath** (string)
- **Tipo**: `string`
- **Default**: `"~/.ioccultcalc/asteroids.db"`
- **Descrizione**: Path file database SQLite

#### **dataSource** (DataSource)
- **Tipo**: Enum `DataSource`
- **Default**: `DataSource::MPC`
- **Valori**:
  - `MPC`: Minor Planet Center
  - `ASTDYS`: AstDyS-2 (Università di Pisa)
  - `JPL`: JPL Small-Body Database
  - `CUSTOM`: File utente

#### **updateInterval** (int)
- **Tipo**: `int`
- **Default**: `7` giorni
- **Descrizione**: Intervallo aggiornamento automatico database

#### **autoUpdate** (bool)
- **Tipo**: `bool`
- **Default**: `true`
- **Descrizione**: Abilita aggiornamento automatico

**Esempio**:
```cpp
AsteroidDatabaseConfig dbConfig;
dbConfig.databasePath = "~/custom/asteroids.db";
dbConfig.dataSource = DataSource::ASTDYS;
dbConfig.updateInterval = 14;  // Ogni 2 settimane
dbConfig.autoUpdate = true;
```

---

## Configurazione Cataloghi Stellari

### Struttura: `StarCatalogConfig`

```cpp
StarCatalogConfig starConfig;
```

### Parametri Gaia DR3

#### **gaiaDataPath** (string)
- **Default**: `"~/.ioccultcalc/gaia_cache"`
- **Descrizione**: Path cache locale Gaia

#### **gaiaQueryURL** (string)
- **Default**: `"https://gea.esac.esa.int/tap-server/tap"`
- **Descrizione**: URL servizio TAP Gaia

#### **maxQueryRadius** (double)
- **Default**: `5.0` gradi
- **Range**: `0.1` - `10.0` gradi
- **Descrizione**: Raggio massimo query Gaia

#### **maxStarsPerQuery** (int)
- **Default**: `10000`
- **Range**: `100` - `100000`
- **Descrizione**: Numero massimo stelle per query

#### **cacheLifetime** (int)
- **Default**: `30` giorni
- **Descrizione**: Validità cache locale

**Esempio**:
```cpp
StarCatalogConfig starConfig;
starConfig.gaiaDataPath = "~/my_gaia_cache";
starConfig.maxQueryRadius = 3.0;
starConfig.maxStarsPerQuery = 5000;
starConfig.cacheLifetime = 60;  // 2 mesi
```

---

## Configurazione Output

### Struttura: `OutputConfig`

```cpp
#include "ioccultcalc/output_manager.h"

OutputConfig outConfig;
```

### Formati Output

#### **formats** (vector<OutputFormat>)
- **Tipo**: `vector<OutputFormat>`
- **Default**: `{OutputFormat::TEXT}`
- **Valori**:
  - `TEXT`: Testo semplice (human-readable)
  - `JSON`: JSON strutturato (parsing)
  - `XML`: XML Occult4-compatible
  - `KML`: Google Earth
  - `IOTA`: IOTA prediction format
  - `CSV`: Comma-separated values

**Esempio**:
```cpp
outConfig.formats = {
    OutputFormat::TEXT,
    OutputFormat::JSON,
    OutputFormat::KML
};
```

---

#### **outputDirectory** (string)
- **Tipo**: `string`
- **Default**: `"./results"`
- **Descrizione**: Directory output file

#### **filenamePrefix** (string)
- **Tipo**: `string`
- **Default**: `"occultation_"`
- **Descrizione**: Prefisso nome file

#### **includeTimestamp** (bool)
- **Tipo**: `bool`
- **Default**: `true`
- **Descrizione**: Include timestamp nei filename

#### **verbosity** (int)
- **Tipo**: `int`
- **Default**: `1`
- **Range**: `0` - `3`
- **Livelli**:
  - `0`: Solo errori
  - `1`: Info base
  - `2`: Info dettagliate
  - `3`: Debug completo

**Esempio**:
```cpp
OutputConfig outConfig;
outConfig.formats = {OutputFormat::JSON, OutputFormat::KML};
outConfig.outputDirectory = "./predictions_2025";
outConfig.filenamePrefix = "asteroid_occ_";
outConfig.includeTimestamp = true;
outConfig.verbosity = 2;
```

---

## File di Configurazione JSON

### Formato File

```json
{
  "propagator": {
    "integrator": "RKF78",
    "stepSize": 0.1,
    "tolerance": 1e-12,
    "usePlanetaryPerturbations": true,
    "useAsteroidalPerturbations": true,
    "useRelativisticCorrections": false,
    "ephemerisType": "JPL_DE441"
  },
  
  "predictor": {
    "startDate": "2025-11-29",
    "endDate": "2025-12-31",
    "observer": {
      "latitude": 45.4642,
      "longitude": 9.1900,
      "altitude": 239.0,
      "name": "Torino"
    },
    "constraints": {
      "minAltitude": 20.0,
      "maxMagnitude": 14.0,
      "minMagnitudeDrop": 0.5,
      "maxShadowWidth": 500.0,
      "minDuration": 0.1,
      "searchRadius": 3.0
    },
    "starCatalog": "GAIA_DR3",
    "parallel": {
      "enabled": true,
      "numThreads": 0
    }
  },
  
  "database": {
    "path": "~/.ioccultcalc/asteroids.db",
    "source": "MPC",
    "autoUpdate": true,
    "updateInterval": 7
  },
  
  "output": {
    "formats": ["TEXT", "JSON", "KML"],
    "directory": "./results",
    "filenamePrefix": "occultation_",
    "includeTimestamp": true,
    "verbosity": 1
  }
}
```

### Caricamento File

```cpp
#include "ioccultcalc/config_manager.h"

ConfigManager config;
config.loadFromFile("my_config.json");

// Accesso parametri
PropagatorOptions propOpts = config.getPropagatorOptions();
OccultationPredictorConfig predConfig = config.getPredictorConfig();
```

---

## Esempi Pratici

### Esempio 1: Survey Veloce Alta Copertura

```cpp
// Configurazione ottimizzata per survey veloce
PropagatorOptions opts;
opts.integrator = IntegratorType::RKF78;
opts.stepSize = 0.5;                          // Step più grande
opts.tolerance = 1e-10;                        // Tolleranza ridotta
opts.usePlanetaryPerturbations = true;
opts.useAsteroidalPerturbations = false;       // Disabilita per velocità
opts.useRelativisticCorrections = false;

OccultationPredictorConfig config;
config.minAltitude = 15.0;                     // Include basse altitudini
config.maxMagnitude = 16.0;                    // Stelle più deboli
config.minMagnitudeDrop = 0.3;                 // Eventi minori
config.useParallelProcessing = true;
config.numThreads = 0;                         // Usa tutti i core
```

---

### Esempio 2: Precisione Massima Pubblicazione

```cpp
// Configurazione precisione massima
PropagatorOptions opts;
opts.integrator = IntegratorType::RKF78;
opts.stepSize = 0.05;                          // Step piccolo
opts.tolerance = 1e-14;                        // Massima precisione
opts.usePlanetaryPerturbations = true;
opts.useAsteroidalPerturbations = true;        // Includi tutto
opts.useRelativisticCorrections = true;        // Per NEA
opts.ephemerisType = EphemerisType::JPL_DE441;

OccultationPredictorConfig config;
config.minAltitude = 30.0;                     // Solo condizioni ottime
config.maxMagnitude = 12.0;                    // Stelle brillanti
config.minMagnitudeDrop = 1.0;                 // Eventi significativi
config.maxShadowWidth = 200.0;                 // Ombre precise
config.minDuration = 0.5;                      // Eventi lunghi
config.searchRadius = 2.0;                     // Solo centrali/near-miss
```

---

### Esempio 3: Ricerca Target Specifico

```cpp
// Configurazione per target specifico (es. 433 Eros)
PropagatorOptions opts;
opts.integrator = IntegratorType::RKF78;
opts.stepSize = 0.1;
opts.tolerance = 1e-12;
opts.usePlanetaryPerturbations = true;
opts.useAsteroidalPerturbations = false;       // Eros è NEA
opts.useRelativisticCorrections = true;        // NEA, importante!

OccultationPredictorConfig config;
config.startDate = JulianDate(2460645.5);
config.endDate = JulianDate(2460645.5 + 7);    // 1 settimana
config.minAltitude = 25.0;
config.maxMagnitude = 13.0;
config.minMagnitudeDrop = 0.8;                 // Eventi buoni
config.searchRadius = 3.5;                     // Include grazing
```

---

### Esempio 4: Test e Debug

```cpp
// Configurazione test/debug
PropagatorOptions opts;
opts.integrator = IntegratorType::RK4;         // Deterministico
opts.stepSize = 1.0;                           // Step grande (veloce)
opts.usePlanetaryPerturbations = true;
opts.useAsteroidalPerturbations = false;
opts.useRelativisticCorrections = false;

OutputConfig outConfig;
outConfig.verbosity = 3;                       // Debug completo
outConfig.formats = {OutputFormat::TEXT};      // Solo testo
```

---

## Best Practices

### Performance

1. **Usa RKF78**: 273x più veloce di RK4 con precisione superiore
2. **Parallelo sempre ON**: `useParallelProcessing = true`
3. **Threads auto**: `numThreads = 0` per auto-detect
4. **Disabilita AST17 per NEA/TNO**: Non rilevante, risparmia tempo
5. **Aumenta stepSize per survey**: 0.5-1.0 giorni per ricerche massive

### Precisione

1. **RKF78 + tolerance 1e-12**: Standard eccellente per pubblicazioni
2. **Abilita perturbazioni planetarie**: Sempre (minimo overhead)
3. **AST17 per main-belt**: Migliora precisione 5-10%
4. **Relatività per NEA**: Solo se q < 1.0 AU
5. **DE441 sempre**: Meglio di VSOP87, overhead trascurabile

### Ricerca Occultazioni

1. **minAltitude >= 20°**: Sotto sconsigliato (refrazione)
2. **maxMagnitude ragionevole**: 14-16 per CCD, 10-12 per visuale
3. **minMagnitudeDrop 0.5-1.0**: Eventi osservabili
4. **searchRadius 2-3**: Ottimo compromesso copertura/falsi positivi
5. **Periodo 1-3 mesi**: Trade-off copertura/tempo computazione

### Database

1. **Auto-update ON**: Elementi sempre aggiornati
2. **ASTDYS per main-belt**: Più preciso di MPC
3. **Update interval 7-14 giorni**: Ragionevole
4. **Cache Gaia valida**: Riduce query remote

### Output

1. **JSON per parsing**: Strutturato, facile da processare
2. **KML per visualizzazione**: Google Earth immediato
3. **TEXT per lettura umana**: Debug e verifica rapida
4. **Timestamp ON**: Tracciabilità risultati

---

## Checklist Pre-Ricerca

### ✅ Configurazione Minima Richiesta

- [ ] `startDate` e `endDate` impostati
- [ ] `observerLocation` (lat, lon, alt)
- [ ] `minAltitude` ragionevole (≥20°)
- [ ] `maxMagnitude` compatibile con setup
- [ ] Database asteroidi popolato
- [ ] Effemeridi JPL installate

### ✅ Ottimizzazioni Raccomandate

- [ ] `integrator = RKF78` (default)
- [ ] `tolerance = 1e-12` (default)
- [ ] `usePlanetaryPerturbations = true`
- [ ] `useParallelProcessing = true`
- [ ] `numThreads = 0` (auto)
- [ ] Output directory configurata
- [ ] Formati output selezionati

### ✅ Per Pubblicazioni Scientifiche

- [ ] `tolerance = 1e-14`
- [ ] `useAsteroidalPerturbations = true` (main-belt)
- [ ] `useRelativisticCorrections = true` (NEA)
- [ ] `minAltitude >= 30°`
- [ ] `minMagnitudeDrop >= 1.0`
- [ ] Validazione round-trip < 10m
- [ ] Confronto con JPL Horizons

---

## Risoluzione Problemi

### Problema: Propagazione lenta

**Soluzioni**:
1. Verifica `integrator = RKF78` (non RK4)
2. Aumenta `stepSize` (0.5-1.0)
3. Riduci `tolerance` (1e-10)
4. Disabilita `useAsteroidalPerturbations` se non necessario

### Problema: Risultati imprecisi

**Soluzioni**:
1. Diminuisci `stepSize` (0.05-0.1)
2. Aumenta `tolerance` precisione (1e-14)
3. Abilita tutte le perturbazioni
4. Verifica elementi asteroidi aggiornati

### Problema: Pochi eventi trovati

**Soluzioni**:
1. Aumenta `maxMagnitude`
2. Diminuisci `minMagnitudeDrop`
3. Aumenta `searchRadius`
4. Estendi periodo ricerca
5. Abbassa `minAltitude` (con cautela)

### Problema: Troppi falsi positivi

**Soluzioni**:
1. Aumenta `minMagnitudeDrop` (≥1.0)
2. Diminuisci `searchRadius` (≤2.0)
3. Riduci `maxShadowWidth`
4. Aumenta `minDuration`

---

## Supporto e Riferimenti

### Documentazione

- Manual completo: `docs/manual/`
- API Reference: `docs/api/`
- Esempi codice: `examples/`

### File Configurazione Esempio

- `preset_default.json`: Configurazione standard
- `preset_high_precision.json`: Massima precisione
- `preset_fast_search.json`: Survey veloce
- `preset_production_monthly.oop`: Produzione mensile

### Contatti

- GitHub Issues: https://github.com/manvalan/IOccultCalc/issues
- Email: [support email]
- Forum: [forum URL]

---

**Ultima revisione**: 29 Novembre 2025  
**Versione IOccultCalc**: 2.1.0-rkf78
