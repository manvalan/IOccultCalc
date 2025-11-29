# Analisi Frame di Riferimento: IOccultCalc vs OrbFit

## Problema Identificato

### OrbFit
**Frame di Riferimento**: ECLM J2000 (Eclittica Media J2000)
- Tutti gli elementi equinoziali sono in frame eclittico
- File `.oel` header: `refsys = ECLM J2000`
- Propagazione avviene in frame eclittico

### IOccultCalc  
**Frame di Riferimento**: MIXED (problema!)

```cpp
// In orbit_propagator.cpp, elementsToState():

// Step 1: Elementi equinoziali ASSUMONO frame eclittico
position_ecl = ... // Calcolo in eclittica

// Step 2: CONVERSIONE a equatoriale
position.y = position_ecl.y * cos_eps - position_ecl.z * sin_eps;
position.z = position_ecl.y * sin_eps + position_ecl.z * cos_eps;

// Step 3: Propagazione avviene in frame equatoriale (DE441)
// getPlanetPosition() usa SPKReader che restituisce J2000 EQUATORIALE
```

### Conseguenze

1. **Input**: Elementi da OrbFit sono in ECLM J2000 ✓
2. **Conversione iniziale**: Eclittica → Equatoriale ✓ 
3. **Propagazione**: Avviene in frame EQUATORIALE ✓
4. **Output**: Vettori in frame EQUATORIALE (non eclittico!)

**PROBLEMA**: L'output di IOccultCalc è in frame EQUATORIALE, ma lo confrontiamo con elementi OrbFit che sono in frame ECLITTICO!

## Verifica nel Codice

### elementsToState (orbit_propagator.cpp, lines 138-250)

```cpp
// Commento esplicito:
// "Position and velocity are in ECLM J2000 (ecliptic frame)"
// "Convert to EQUATORIAL J2000 for the propagator"

// Rotazione obliquità (eclittica → equatoriale)
constexpr double obliquity = 23.4392911 * M_PI / 180.0;
position.x = position_ecl.x;
position.y = position_ecl.y * cos_eps - position_ecl.z * sin_eps;  // Rotazione
position.z = position_ecl.y * sin_eps + position_ecl.z * cos_eps;
```

### stateToElements (orbit_propagator.cpp, line 252)

```cpp
// Commento:
// "I vettori in input sono in frame EQUATORIALE (da Horizons)"
// "Ma elementsToState assume elementi in frame ECLITTICO"
// "Quindi convertiamo vettori equatoriali → eclittici"
```

Questo conferma che:
- `elementsToState`: ECLITTICO → EQUATORIALE
- `stateToElements`: EQUATORIALE → ECLITTICO  
- Propagazione: avviene in EQUATORIALE

## Soluzione

Per confrontare correttamente IOccultCalc con OrbFit, dobbiamo:

### Opzione 1: Convertire output IOccultCalc
```cpp
// Dopo propagazione
Vector3D pos_io_eq = state_target.position;  // Equatoriale

// Converti a eclittico
double cos_eps = cos(23.4392911 * M_PI / 180.0);
double sin_eps = sin(23.4392911 * M_PI / 180.0);

Vector3D pos_io_ecl;
pos_io_ecl.x = pos_io_eq.x;
pos_io_ecl.y = pos_io_eq.y * cos_eps + pos_io_eq.z * sin_eps;
pos_io_ecl.z = -pos_io_eq.y * sin_eps + pos_io_eq.z * cos_eps;

// Ora pos_io_ecl è confrontabile con OrbFit!
```

### Opzione 2: Convertire elementi OrbFit a equatoriale
Convertire gli elementi OrbFit target da eclittico a equatoriale prima del confronto.

### Opzione 3: Propagare in frame eclittico (meglio!)
Modificare IOccultCalc per propagare direttamente in frame eclittico come OrbFit:
1. Non fare conversione eclittica→equatoriale in `elementsToState`
2. Convertire posizioni pianeti da equatoriale a eclittico
3. Propagare tutto in eclittico
4. Output finale in eclittico (match OrbFit)

## Test Corretto

### test_pompeja.cpp ATTUALE
```cpp
// Output:
X = -2.637101661238 AU    // EQUATORIALE J2000
Y = 0.842822657440 AU     // EQUATORIALE J2000
Z = 0.385823846490 AU     // EQUATORIALE J2000
```

### test_pompeja.cpp CORRETTO
```cpp
// Dopo propagazione, converti a eclittico:
Vector3D pos_ecl = equatorialToEcliptic(state_target.position);

// Output confrontabile con OrbFit:
X_ecl = ... AU    // ECLM J2000 (come OrbFit)
Y_ecl = ... AU    // ECLM J2000
Z_ecl = ... AU    // ECLM J2000
```

## Verifica OrbFit Workflow

Dal log `203.olg`:
```
Starting orbital elements:
   ...elements in ECLM J2000...
   Epoch: MJD 61000.00

Propagated orbital elements:
   ...elements in ECLM J2000...
   Epoch: MJD 61192.00
```

OrbFit mantiene SEMPRE gli elementi in frame ECLM J2000.

## Raccomandazioni

### Priorità 1: Fix Immediato
Aggiungere conversione equatoriale→eclittico in `test_pompeja.cpp` per confronto corretto.

### Priorità 2: Documentazione
Aggiungere commenti espliciti in:
- `orbit_propagator.h`: "Propagation occurs in J2000 EQUATORIAL frame"
- `elementsToState()`: "Converts ECLIPTIC elements to EQUATORIAL state"
- Output functions: Specificare frame di output

### Priorità 3: Consistency Check
Verificare che tutte le funzioni di I/O documentino il frame di riferimento usato.

## Conclusione

La discrepanza di **350 milioni di km** nel confronto precedente era dovuta al mixing di frame:
- IOccultCalc output in EQUATORIALE J2000
- OrbFit output in ECLM J2000
- Confronto diretto senza conversione = risultati completamente sbagliati

Con la conversione corretta, ci aspettiamo differenze dell'ordine di **km o decine di km** per 192 giorni di propagazione.
