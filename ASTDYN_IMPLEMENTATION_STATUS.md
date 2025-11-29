# Stato Implementazione AstDyn/RKF78 per IOccultCalc

## Data: 29 Novembre 2025

## Sommario Implementazione

**STATO**: Implementazione PARZIALE - Infrastruttura creata, stub funzionanti, compilazione da completare

### Cosa È Stato Fatto ✅

1. **CMakeLists.txt Aggiornato**
   - Aggiunta opzione `USE_ASTDYN` (default: ON)
   - Configurazione ITALOccultLibrary/astdyn come subdirectory
   - Linking con Eigen3 e Boost (dipendenze AstDyn)
   - Source file `astdyn_propagator_impl.cpp` aggiunto al build

2. **Header Interfaccia Creato**
   - `include/ioccultcalc/astdyn_interface.h` (400+ righe)
   - Strutture dati: `AstDySElements`, `RWOObservation`, `OrbitFitResult`
   - Classi: `AstDynPropagator`, `AstDynOrbitFitter`, `AstDySClient`
   - Namespace utilità: `AstDySUtil`

3. **Implementazione Stub Creata**
   - `src/astdyn_propagator_impl.cpp` (300+ righe)
   - Versione compilabile con implementazioni placeholder
   - Supporto fallback quando `USE_ASTDYN=OFF`
   - Warning chiari per metodi non implementati

4. **Preset di Test Completo**
   - `preset_11234_rkf78_validation.oop` (200+ righe)
   - Configurazione ottimale per validazione
   - Asteroide (11234) 1999 JS82 come test case
   - Documentazione inline completa

5. **Test e Validazione**
   - `test_propagation_compare.cpp` - Confronto propagazione
   - `compare_jpl_horizons.py` - Validazione automatica JPL
   - `PROPAGATION_VALIDATION_TEST.md` - Documentazione risultati

### Cosa Manca ❌

1. **Implementazione Reale AstDynPropagator**
   - Attualmente è uno stub che ritorna elementi non modificati
   - Necessario: integrazione vera con `astdyn::propagation::Propagator`
   - Conversione formati IOccultCalc ↔ AstDyn
   - Gestione errori e eccezioni

2. **Implementazione AstDynOrbitFitter**
   - Stub placeholder senza funzionalità
   - Necessario: differential correction completo
   - Calcolo residui O-C
   - Outlier detection

3. **Implementazione AstDySClient**
   - Download da newton.spacedys.com
   - Parser file .eq1 e .rwo
   - Gestione cache locale

4. **Fix Errori Compilazione**
   - Alignment signature metodi con header
   - Include paths per ITALOccultLibrary
   - Linking con librerie Eigen/Boost

5. **Integrazione OrbitPropagator**
   - Selezione propagatore basata su config `.oop`
   - Switch tra RK4/RKF78/AstDyn-RKF78
   - Performance benchmarking

## File Creati/Modificati

### Nuovi File:
```
include/ioccultcalc/astdyn_interface.h          # 400+ righe - Interfaccia completa
src/astdyn_propagator_impl.cpp                  # 300+ righe - Stub implementazione
preset_11234_rkf78_validation.oop               # 200+ righe - Preset di test
examples/test_propagation_compare.cpp           # 500+ righe - Test Kepleriano
compare_jpl_horizons.py                         # 300+ righe - Validazione JPL
PROPAGATION_VALIDATION_TEST.md                  # 300+ righe - Documentazione
ASTDYN_IMPLEMENTATION_STATUS.md                 # Questo file
```

### File Modificati:
```
CMakeLists.txt                                  # +50 righe - Supporto AstDyn
```

## Prossimi Passi per Completamento

### Priorità 1: Far Compilare 🔴

1. **Fix signature metodi** in `astdyn_propagator_impl.cpp`:
   ```cpp
   // Correggere:
   void enableRelativity(bool) → void useRelativisticCorrections(bool)
   AstDySElements propagateToEpoch(...) → AstDySElements propagate(...)
   
   // Aggiungere metodi mancanti dall'header
   ```

2. **Fix struct OrbitFitResult**:
   ```cpp
   // Usare nomi corretti:
   converged → (non esiste nell'header, rimuovere)
   num_observations → n_observations
   num_outliers → n_outliers
   improved_elements → fitted_elements
   ```

3. **Rimuovere metodi non nell'header**:
   - `computeResiduals()` in AstDynOrbitFitter
   - `downloadElements/Observations()` in AstDySClient
   - Allineare completamente con dichiarazioni header

### Priorità 2: Implementazione Vera 🟡

4. **Integrare ITALOccultLibrary/AstDyn reale**:
   ```cpp
   #ifdef USE_ASTDYN
   #include "astdyn/propagation/Propagator.hpp"
   #include "astdyn/core/OrbitalElements.hpp"
   // ... altri include necessari
   ```

5. **Implementare conversione formati**:
   - `AstDySElements` (IOccultCalc) ↔ `astdyn::core::OrbitalElements`
   - Gestione sistemi di riferimento (ECLM J2000)
   - Conversione epoche (MJD ↔ astdyn::time)

6. **Implementare propagazione**:
   ```cpp
   AstDySElements AstDynPropagator::propagate(...) {
       // 1. Converti IOccultCalc → AstDyn
       // 2. Chiama astdyn::propagation::Propagator
       // 3. Converti risultato AstDyn → IOccultCalc
       // 4. Gestisci errori
   }
   ```

### Priorità 3: Testing e Validazione 🟢

7. **Compilare con `USE_ASTDYN=ON`**:
   ```bash
   cd build
   cmake -DUSE_ASTDYN=ON -DCMAKE_BUILD_TYPE=Release ..
   make -j8
   ```

8. **Test propagazione**:
   ```bash
   ./test_propagation_compare
   python3 compare_jpl_horizons.py
   ```

9. **Validazione accuratezza**:
   - Confronto con JPL Horizons
   - RMS residui < 1" a ±60 giorni
   - Performance < 0.1s per propagazione

### Priorità 4: Integrazione Sistema 🔵

10. **OrbitPropagator::propagate()** supporta AstDyn:
    ```cpp
    if (options_.propagator == "AstDyn-RKF78") {
        AstDynPropagator astdyn;
        // ... usa AstDyn
    }
    ```

11. **Parser configurazione .oop**:
    ```ini
    general.
        .propagator = 'AstDyn-RKF78'  # Nuovo valore
    
    astdyn.
        .tolerance = 1.0e-12
        .use_planets = .TRUE.
        .use_asteroids = .TRUE.
        .use_relativity = .TRUE.
    ```

12. **Output e logging**:
    - Statistiche propagazione (steps, tempo)
    - Confronto performance RK4 vs RKF78
    - Warning se accuratezza degrada

## Comandi Utili

### Compilazione
```bash
# Full rebuild con AstDyn
rm -rf build && mkdir build && cd build
cmake -DUSE_ASTDYN=ON -DCMAKE_BUILD_TYPE=Release ..
make -j8

# Solo stub (senza ITALOccultLibrary)
cmake -DUSE_ASTDYN=OFF ..
make -j8
```

### Test
```bash
# Test propagazione Kepleriana (baseline)
./test_propagation_compare

# Validazione JPL Horizons
python3 compare_jpl_horizons.py

# Run completo con preset
./italoccultcalc preset_11234_rkf78_validation.oop
```

### Debug
```bash
# Verifica linking
otool -L build/libioccultcalc.a  # macOS
ldd build/libioccultcalc.so      # Linux

# Check simboli
nm build/libioccultcalc.a | grep AstDyn

# Dipendenze Eigen/Boost
cmake --build build --target help | grep astdyn
```

## Note Tecniche

### Dipendenze AstDyn
- **Eigen3** >= 3.4: Algebra lineare
- **Boost** >= 1.70: filesystem, program_options, date_time
- **C++17**: Required

### Performance Attese
| Metodo | Tempo | Accuratezza |
|--------|-------|-------------|
| RK4 puro | 0.01s | 10-50" a ±60gg |
| RKF78 stub | N/A | 0" (nessuna propagazione) |
| **RKF78 reale** | **~0.05s** | **1-3" a ±60gg** |
| JPL Horizons | N/A | < 0.05" (standard) |

### Limitazioni Correnti
1. ❌ Propagazione non implementata (stub ritorna input invariato)
2. ❌ Orbit fitting non funzionante
3. ❌ Download AstDyS non implementato
4. ⚠️ Non compila con errori di signature
5. ⚠️ ITALOccultLibrary/AstDyn API in evoluzione

### Quando Sarà Completo
- ✅ Propagazione RKF78 con perturbazioni complete
- ✅ Accuratezza < 1" a ±30 giorni, 1-3" a ±60 giorni
- ✅ Fitting orbitale con 6000+ osservazioni
- ✅ Download automatico da AstDyS
- ✅ Integrazione seamless con IOccultCalc
- ✅ Preset production-ready

## Riferimenti

- **ITALOccultLibrary**: `external/ITALOccultLibrary/astdyn/`
- **Documentazione RKF78**: `external/ITALOccultLibrary/astdyn/docs/RKF78_INTEGRATOR.md`
- **Test validazione**: `PROPAGATION_VALIDATION_TEST.md`
- **Issue tracking**: Da creare su GitHub

## Autore

Michele Bigi  
Data: 29 Novembre 2025  
Versione: 0.1-alpha (stub)

---

**IMPORTANTE**: Questa è un'implementazione PARZIALE. L'infrastruttura è pronta ma la funzionalità non è ancora operativa. Usare solo a scopo di sviluppo e testing dell'architettura.
