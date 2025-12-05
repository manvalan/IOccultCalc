# Guida Integrazione ITALOccultLibrary in IOccultCalc

## Panoramica

Questa guida spiega come integrare **ITALOccultLibrary** in **IOccultCalc** per la propagazione degli asteroidi con AstDyn.

## Componenti

### 1. ITALOccultLibrary
Libreria installata in `/usr/local` che fornisce:
- **eq1_parser**: Parser per file .eq1 di OrbFit/AstDys
- **orbital_conversions**: Conversioni tra frame di riferimento
- **astdyn_wrapper**: Wrapper high-level per AstDyn con RKF78

### 2. Layer di Integrazione
File in `include/ioccultcalc/integration/`:
- **italoccult_integration.h**: Header con classe `ITALOccultIntegration`
- Struct `AsteroidState` nel formato IOccultCalc
- Helper function `quickPropagateFromEQ1()`

File in `src/integration/`:
- **italoccult_integration.cpp**: Implementazione della conversione

## Installazione

### Prerequisiti
```bash
# ITALOccultLibrary deve essere installata
ls /usr/local/lib/libitaloccultlib.a
ls /usr/local/include/italoccultlib/

# AstDyn deve essere installato
ls /usr/local/lib/libastdyn.a
ls /usr/local/include/astdyn/

# Boost 1.89+ richiesto da AstDyn
brew install boost
```

### Aggiornamento CMakeLists.txt

Aggiungi dopo `find_package(Eigen3)`:

```cmake
# ITALOccultLibrary per propagazione con AstDyn
find_package(ITALOccultLibrary REQUIRED)
find_package(AstDyn REQUIRED)

# Aggiungi file di integrazione ai sources
set(SOURCES
    ${SOURCES}
    src/integration/italoccult_integration.cpp
)

# Link alla libreria principale
target_link_libraries(ioccultcalc
    ITALOccultLibrary::italoccultlib
    AstDyn::astdyn
    # ... altri link esistenti
)
```

## Utilizzo

### Esempio Base: Propagazione Singola

```cpp
#include <ioccultcalc/integration/italoccult_integration.h>
#include <iostream>

int main() {
    // Crea integrator con alta precisione
    ioccultcalc::ITALOccultIntegration integrator(
        ioccultcalc::PropagationSettings::highAccuracy()
    );
    
    // Carica asteroide da file .eq1
    if (!integrator.loadAsteroidFromEQ1("data/17030.eq1")) {
        std::cerr << "Errore caricamento" << std::endl;
        return 1;
    }
    
    // Propaga a data specifica (MJD TDB)
    double target_mjd = 61007.0;
    auto state = integrator.propagateToEpoch(target_mjd);
    
    // Stampa risultati
    std::cout << "Asteroide: " << state.name << std::endl;
    std::cout << "Epoca: MJD " << state.epoch_mjd_tdb << std::endl;
    std::cout << "Posizione ICRF (AU): ["
              << state.x_icrf_au << ", "
              << state.y_icrf_au << ", "
              << state.z_icrf_au << "]" << std::endl;
    
    double distance_au = std::sqrt(
        state.x_icrf_au * state.x_icrf_au +
        state.y_icrf_au * state.y_icrf_au +
        state.z_icrf_au * state.z_icrf_au
    );
    std::cout << "Distanza: " << distance_au << " AU" << std::endl;
    
    return 0;
}
```

### Esempio Avanzato: Propagazione Multiple Epoche

```cpp
#include <ioccultcalc/integration/italoccult_integration.h>
#include <vector>
#include <iostream>

int main() {
    ioccultcalc::ITALOccultIntegration integrator(
        ioccultcalc::PropagationSettings::highAccuracy()
    );
    
    integrator.loadAsteroidFromEQ1("data/17030.eq1");
    
    // Genera epoche (ogni giorno per 10 giorni)
    std::vector<double> epochs;
    double start_mjd = 61000.0;
    for (int i = 0; i < 10; ++i) {
        epochs.push_back(start_mjd + i);
    }
    
    // Propaga tutte le epoche
    auto states = integrator.propagateToEpochs(epochs);
    
    // Processa risultati
    for (const auto& state : states) {
        double r = std::sqrt(
            state.x_icrf_au * state.x_icrf_au +
            state.y_icrf_au * state.y_icrf_au +
            state.z_icrf_au * state.z_icrf_au
        );
        
        std::cout << "MJD " << state.epoch_mjd_tdb 
                  << ": R = " << r << " AU" << std::endl;
    }
    
    return 0;
}
```

### Helper Function: Quick Propagate

Per uso rapido senza creare oggetti:

```cpp
#include <ioccultcalc/integration/italoccult_integration.h>

// Propagazione rapida one-liner
auto state = ioccultcalc::quickPropagateFromEQ1(
    "data/17030.eq1",
    61007.0,  // MJD target
    ioccultcalc::PropagationSettings::highAccuracy()
);

std::cout << "Posizione: [" 
          << state.x_icrf_au << ", "
          << state.y_icrf_au << ", "
          << state.z_icrf_au << "] AU" << std::endl;
```

## Struct AsteroidState

Il formato output per IOccultCalc:

```cpp
struct AsteroidState {
    std::string name;                    // Nome/numero asteroide
    double epoch_mjd_tdb;                // Epoca in MJD TDB
    
    // Posizione ICRF (AU)
    double x_icrf_au;
    double y_icrf_au;
    double z_icrf_au;
    
    // Velocità ICRF (AU/day)
    double vx_icrf_au_day;
    double vy_icrf_au_day;
    double vz_icrf_au_day;
    
    // Parametri orbitali approssimati
    double semi_major_axis_au;           // Semiasse maggiore (AU)
    double eccentricity;                 // Eccentricità
    double inclination_deg;              // Inclinazione (gradi)
    
    // Metadati
    std::string propagation_info;        // Info sulla propagazione
};
```

## Preset di Configurazione

### High Accuracy (default)
```cpp
PropagationSettings::highAccuracy()
```
- Tolleranza: 1e-12 AU
- Step iniziale: 0.01 giorni
- Tutte le perturbazioni attive (8 pianeti + relatività)
- **Uso**: Calcoli scientifici, previsioni occultazioni

### Fast Mode
```cpp
PropagationSettings::fast()
```
- Tolleranza: 1e-9 AU
- Step iniziale: 0.1 giorni
- Perturbazioni ridotte (4 pianeti principali)
- **Uso**: Survey rapide, studi preliminari

### Custom
```cpp
PropagationSettings custom;
custom.tolerance = 1e-10;
custom.initial_step = 0.05;
custom.include_planets = true;
custom.include_relativity = false;
custom.perturb_mercury = true;
custom.perturb_venus = true;
custom.perturb_earth = true;
custom.perturb_mars = true;
custom.perturb_jupiter = true;
custom.perturb_saturn = false;
custom.perturb_uranus = false;
custom.perturb_neptune = false;
```

## Performance

Test su MacBook Air M1:
- **Propagazione 7 giorni**: < 1 ms (17030 Sierks)
- **Propagazione 10 epoche**: ~ 2-3 ms
- **Accuratezza**: Confronto con JPL Horizons < 10 km (asteroid 17030)

## Frame di Riferimento

- **Input**: ECLM J2000 (dal file .eq1)
- **Conversione interna**: ECLM → ICRF (obliquità 23.4393°)
- **Output**: ICRF J2000 (equatoriale)

**Nota**: IOccultCalc può richiedere ulteriori conversioni per frame observer-centered.

## Troubleshooting

### Errore: "Cannot find ITALOccultLibrary"
```bash
# Verifica installazione
ls /usr/local/lib/cmake/ITALOccultLibrary/

# Se mancante, reinstalla
cd ITALOccultLibrary/italoccultlibrary/build
sudo make install
```

### Errore: "Cannot find AstDyn"
```bash
# Crea config manualmente se necessario
# Vedi /usr/local/lib/cmake/AstDyn/AstDynConfig.cmake
```

### Warning Boost 1.89
Il warning `CMP0167` è innocuo. Per sopprimerlo:
```cmake
cmake_policy(SET CMP0167 NEW)
```

## File di Test

Esempio completo in `test_ioccultcalc_integration.cpp`:
```bash
cd ITALOccultLibrary
g++ -std=c++17 -O2 -o test_ioccultcalc test_ioccultcalc_integration.cpp \
    integration/italoccult_integration.cpp \
    -I/usr/local/include -I/usr/local/include/eigen3 -I./integration \
    -L/usr/local/lib -litaloccultlib -lastdyn

./test_ioccultcalc
```

## Supporto

- **Repository**: https://github.com/manvalan/ITALOccultLibrary
- **AstDyn**: https://github.com/[astdyn-repo]
- **IOccultCalc**: https://github.com/[ioccultcalc-repo]

---

**Data**: 1 Dicembre 2025  
**Versione**: ITALOccultLibrary 1.0.0  
**Autore**: IOccultCalc Integration Team
