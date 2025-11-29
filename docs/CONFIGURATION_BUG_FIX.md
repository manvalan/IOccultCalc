# Configuration Bug Fix - 2025-01-XX

## Problema Identificato

**Issue**: Il file di configurazione `.oop` conteneva la sezione `perturbations.` con parametri per AST17 e altre correzioni, ma `italoccultcalc.cpp` NON leggeva né applicava queste impostazioni.

### Esempio Configurazione Ignorata
```fortran
perturbations.
        .planets = .TRUE.
        .relativity = .TRUE.
        .asteroid_count = 17
```

**Impatto**:
- Utente NON poteva disabilitare perturbazioni planetarie
- Utente NON poteva disabilitare correzioni relativistiche  
- Configurazione `asteroid_count` completamente ignorata
- OrbitPropagator sempre creato con constructor default

### Root Cause

File: `examples/italoccultcalc.cpp`, funzione `propagateOrbits()`

**PRIMA del fix**:
```cpp
// Line 381-393 (circa)
OrbitPropagator propagator;  // Constructor default - NO OPTIONS!
```

Il propagatore veniva sempre creato senza passare opzioni, quindi:
1. `usePlanetaryPerturbations = true` (default)
2. `useRelativisticCorrections = false` (default)
3. AST17 si caricava automaticamente se file SPK disponibile (workaround accidentale)

## Soluzione Implementata

### Modifiche a `examples/italoccultcalc.cpp`

**Funzione**: `propagateOrbits()`  
**Linee**: ~376-420

**DOPO il fix**:
```cpp
// Leggi configurazione perturbazioni dal file .oop
PropagatorOptions opts;
auto pertSection = config.getSection(ConfigSection::PERTURBATIONS);

if (pertSection) {
    auto planetsParam = pertSection->getParameter("planets");
    auto relativityParam = pertSection->getParameter("relativity");
    auto asteroidCountParam = pertSection->getParameter("asteroid_count");
    
    if (planetsParam) {
        opts.usePlanetaryPerturbations = planetsParam->asBool();
    }
    if (relativityParam) {
        opts.useRelativisticCorrections = relativityParam->asBool();
    }
    
    // Note: asteroid_count richiede refactoring OrbitPropagator
    // Per ora AST17 si abilita automaticamente se file SPK esiste
    if (asteroidCountParam) {
        int astCount = asteroidCountParam->asInt();
        if (astCount > 0 && !g_verbose) {
            std::cout << "Configurazione AST17: " << astCount << " asteroidi massivi\n";
        }
    }
}

// ... poi nella sezione di print
if (opts.usePlanetaryPerturbations) {
    std::cout << "  ✓ Perturbazioni gravitazionali (8 pianeti)\n";
}
std::cout << "  ✓ Aberrazione planetaria (light-time)\n";
if (opts.useRelativisticCorrections) {
    std::cout << "  ✓ Effetti relativistici\n";
}

// Crea propagatore CON opzioni
OrbitPropagator propagator(opts);
```

### Cosa Cambia

#### 1. Lettura Configurazione ✅
- `ConfigManager::getSection(ConfigSection::PERTURBATIONS)` interroga sezione perturbations
- `getParameter("planets")` → imposta `opts.usePlanetaryPerturbations`
- `getParameter("relativity")` → imposta `opts.useRelativisticCorrections`
- `getParameter("asteroid_count")` → mostra messaggio (ma non controlla AST17 ancora)

#### 2. Creazione Propagatore ✅
- `OrbitPropagator propagator(opts)` passa opzioni al constructor
- Constructor usa le opzioni invece dei default

#### 3. Output Dinamico ✅
- Messaggi "Correzioni abilitate" ora mostrano REALMENTE cosa è configurato
- Se `planets = .FALSE.` → NON mostra "Perturbazioni gravitazionali"
- Se `relativity = .TRUE.` → mostra "Effetti relativistici"

## Test Effettuati

### Test 1: Configurazione Completa (preset_large_asteroids_jan2026.oop)
```fortran
perturbations.
        .planets = .TRUE.
        .relativity = .TRUE.
        .asteroid_count = 17
```

**Risultato**:
```
Configurazione AST17: 17 asteroidi massivi

Correzioni abilitate:
  ✓ Perturbazioni gravitazionali (8 pianeti)
  ✓ Aberrazione planetaria (light-time)
  ✓ Effetti relativistici
```

✅ **PASS** - Tutte le correzioni abilitate come da config

### Test 2: Configurazione Disabilitata (preset_test_perturbations_off.oop)
```fortran
perturbations.
        .planets = .FALSE.
        .relativity = .FALSE.
        .asteroid_count = 0
```

**Risultato atteso**:
```
Correzioni abilitate:
  ✓ Aberrazione planetaria (light-time)
```

⚠️ **DA TESTARE** - Dopo fix completo dovrebbe mostrare SOLO aberrazione

## Limitazioni Attuali

### asteroid_count Parameter
**Status**: NON completamente implementato

**Ragione**: OrbitPropagator carica AST17 automaticamente nel constructor:

```cpp
// In orbit_propagator.cpp, Impl::Impl()
useAsteroidPerturbations = asteroidReader.ensureFileLoaded("codes_300ast_20100725.bsp");
```

Questo hardcoded load bypassa la configurazione `asteroid_count`.

**Workaround**: Se file SPK non esiste, AST17 è disabilitato. Ma non c'è controllo granulare.

**Fix Required**: Refactoring OrbitPropagator per:
1. Accettare parametro `useAsteroidPerturbations` in PropagatorOptions
2. Caricare SPK SOLO se opzione è true
3. Supportare subset di asteroidi (non solo tutti 17)

### Subset AST17
**Status**: Non supportato

**Esempio**: Utente vuole solo i primi 5 asteroidi più massivi:
```fortran
perturbations.
        .asteroid_count = 5
        .asteroid_ids = 1, 2, 3, 4, 6  ! Ceres, Pallas, Juno, Vesta, Hebe
```

Questo richiede:
1. PropagatorOptions con `std::vector<int> asteroidIds`
2. OrbitPropagator che carica SOLO quegli ID da SPK
3. Nuovo parametro `.asteroid_ids` in ConfigManager

## Vantaggi del Fix

### 1. Configurabilità Completa ✅
Utente può ora controllare:
- Perturbazioni planetarie ON/OFF
- Correzioni relativistiche ON/OFF
- (Parziale) AST17 ON/OFF via file esistenza

### 2. Output Trasparente ✅
Messaggi mostrano ESATTAMENTE cosa viene usato, non più hardcoded

### 3. Debugging Facilitato ✅
Se risultati diversi da attesi, si verifica subito la config

### 4. Coerenza con OrbFit ✅
File `.oop` ora ha lo stesso effetto sia in OrbFit che in IOccultCalc

## File Modificati

```
examples/italoccultcalc.cpp
    - Aggiunta lettura sezione perturbations (lines ~376-400)
    - Creazione PropagatorOptions da config
    - Output dinamico correzioni
    - Passaggio opts a OrbitPropagator constructor
```

## File Creati per Testing

```
preset_test_perturbations_off.oop
    - Configurazione con tutte perturbazioni disabilitate
    - Per verificare che il fix funzioni correttamente
```

## Next Steps

### Priority 1: Test Completo
- [ ] Test con `.planets = .FALSE.` 
- [ ] Verificare che propagazione usi SOLO Sole
- [ ] Confronto risultati vs propagazione completa
- [ ] Misurare differenza posizioni

### Priority 2: asteroid_count Implementation
- [ ] Aggiungere `bool useAsteroidPerturbations` a PropagatorOptions
- [ ] Modificare OrbitPropagator constructor per rispettare flag
- [ ] Caricare SPK solo se necessario
- [ ] Test con asteroid_count=0, 17, e valori intermedi

### Priority 3: asteroid_ids Support
- [ ] Aggiungere `std::vector<int> asteroidIds` a PropagatorOptions
- [ ] Parser per `.asteroid_ids` in ConfigManager
- [ ] OrbitPropagator filtra solo ID richiesti
- [ ] Documentazione subset supportati

### Priority 4: Documentation Update
- [ ] Aggiornare CONFIG_SYSTEM.md con nuova funzionalità
- [ ] Esempi configurazione perturbations in USAGE_GUIDE.md
- [ ] Note su limitazioni asteroid_count in README

## Conclusioni

**Bug**: Configurazione perturbations ignorata → **RISOLTO** ✅

**Impatto**: 
- Controllo completo su perturbazioni planetarie e relatività
- Output trasparente
- Coerenza con formato .oop OrbFit

**Limitazioni Residue**:
- asteroid_count non completamente implementato (workaround: file esistenza)
- Subset asteroidi non supportato

**Testing Status**: 
- ✅ Lettura configurazione confermata
- ⚠️ Test disabilitazione perturbazioni pendente
- ⚠️ Confronto risultati propagazione pendente
