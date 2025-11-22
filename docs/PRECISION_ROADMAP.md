# Roadmap: Raggiungere Precisione OrbFit-Level

## Obiettivo
Orbit determination con precisione **sub-arcsecond** (< 0.1") paragonabile a OrbFit professionale.

---

## Stato Attuale vs Target

### Precisione Attuale (Nostri Algoritmi)
- **Herget (3 obs)**: ✅ ALGORITMO CORRETTO (bug era in selezione osservazioni test)
- **Differential Corrector**: ~0.5-1° dopo convergenza
- **Propagazione**: ~1" per anno con RA15

**NOTA**: Analisi dettagliata bug Herget in `docs/HERGET_BUG_ANALYSIS.md`. 
Tutti gli algoritmi fondamentali (Lambert solver, frame conversion, orbital elements) sono corretti.
Il problema era osservazioni con Dec=40° inadatte per IOD.

### Precisione Target (OrbFit-Level)
- **Preliminary orbit**: <1° tutti gli elementi ✅ GIÀ RAGGIUNTO con osservazioni adeguate
- **Differential Corrector**: <0.1" RMS residui
- **Propagazione**: <0.01" per anno

### Gap da Colmare
**2-3 ordini di grandezza** sui residui astrometrici (ridotto da 4 dopo fix IOD)

---

## Componenti per Alta Precisione

### 1. ✅ JPL Ephemerides (FATTO)
- [x] DE441 integrato via SPICE
- [x] Posizioni planetarie precise
- **Impatto**: ~0.001" sulla posizione Terra

### 2. ⚠️ Observatory Positions (PARZIALE)
**Status**: Usiamo codici MPC ma senza correzioni atmosferiche
**Serve**:
- [ ] Earth Orientation Parameters (EOP) da IERS
- [ ] Nutation & precession IAU 2006/2000A
- [ ] Observatory elevations precise
- [ ] Light-time correction iterativa
**Impatto**: ~0.1" sulla posizione osservatorio

### 3. ❌ Asteroid Perturbations (MANCANTE)
**Status**: Non considerate
**Serve**:
- [ ] Masse 16 major asteroids (Ceres, Vesta, Pallas, ...)
- [ ] Effemeridi asteroidi da JPL
- [ ] Force model con perturbazioni
**Impatto**: ~0.01" su decenni per NEAs

### 4. ⚠️ Relativity (PARZIALE)
**Status**: Non implementata
**Serve**:
- [ ] Post-Newtonian corrections (Schwarzschild metric)
- [ ] Light deflection by Sun
- [ ] Shapiro delay
**Impatto**: ~0.01" su grandi distanze

### 5. ❌ Non-Gravitational Forces (MANCANTE)
**Status**: Non considerate
**Serve**:
- [ ] Yarkovsky effect (thermal recoil)
- [ ] YORP effect (spin change)
- [ ] Solar radiation pressure
**Impatto**: ~1" per anno per piccoli NEAs

### 6. ⚠️ Observational Weights (BASICO)
**Status**: Pesi uniformi
**Serve**:
- [ ] Sigma per osservatorio/epoca
- [ ] Outlier rejection robusto
- [ ] Correlation handling
**Impatto**: 2x miglioramento convergenza

### 7. ✅ Preliminary Orbit (ALGORITMO CORRETTO)
**Status**: Herget con Lambert solver funziona correttamente
**Validazione Completata**:
- Lambert solver accurato (<1% errore su test sintetici) ✅
- Conversione coordinate corretta (5 test passati) ✅
- Frame di riferimento corretto (ECLIPJ2000) ✅
- Calcolo orbital elements corretto ✅
**Root Cause Identificata**: Selezione osservazioni inadatte (Dec=40° per Eros i=11°)
**Serve** (bassa priorità):
- [ ] Filtro osservazioni in `selectObservations()` (|Dec|<30°)
- [ ] Selezione intelligente per evitare configurazioni estreme
- [ ] Weighted observations nel differential corrector
**Impatto**: Algoritmo già preciso, serve solo pre-processing migliore

---

## Priority Stack (Massimo Impatto/Sforzo)

### Priority 1: ROBUST DIFFERENTIAL CORRECTOR ⚠️ NUOVA PRIORITÀ
**Problema**: Convergenza lenta/instabile con preliminary orbits buoni
**Impatto**: 10x improvement residui → RMS <0.1"
**Effort**: 1-2 giorni

**Components**:
- [ ] Levenberg-Marquardt damping (non solo Gauss-Newton)
- [ ] Adaptive step size
- [ ] Outlier rejection (3-sigma clipping)
- [ ] Proper observational weights (σ per osservatorio/epoca)
- [ ] Robust convergence criteria

### Priority 2: EARTH ORIENTATION PARAMETERS
**Problema**: Observatory positions approssimate
**Impatto**: 0.1" improvement
**Effort**: 3-5 giorni
**Components**:
- [ ] Download IERS Bulletin A/B
- [ ] Parse EOP (UT1-UTC, pole coordinates)
- [ ] Nutation IAU 2006
- [ ] Precession IAU 2000A

### Priority 3: ASTEROID PERTURBATIONS
**Problema**: Force model incompleto
**Impatto**: 0.01" long-term
**Effort**: 1 settimana
**Components**:
- [ ] Load 16 asteroid masses
- [ ] Query JPL Horizons for asteroid ephemerides
- [ ] Integrate in force_model.cpp
- [ ] Test con Eros (perturbato da Ceres)

### Priority 4: RELATIVITY
**Impatto**: 0.01" su distanze >1 AU
**Effort**: 2-3 giorni
**Components**:
- [ ] Schwarzschild metric terms
- [ ] Light deflection
- [ ] Shapiro delay

---

## Implementation Plan

### Phase 1: FIX CRITICAL BUGS (Settimana 1)
**Goal**: Herget funzionante con <1° errore

Tasks:
1. [ ] Debug findTransferOrbit() con synthetic orbit
2. [ ] Unit tests per ogni step Herget
3. [ ] Validate con 10 asteroids noti
4. [ ] Document algoritmo funzionante

**Success Metric**: i_computed = 11° ± 0.5° per Eros

### Phase 2: ROBUST FITTING (Settimana 2)
**Goal**: Differential corrector con RMS <1"

Tasks:
1. [ ] Implement Levenberg-Marquardt
2. [ ] Adaptive weights per osservatorio
3. [ ] 3-sigma outlier rejection
4. [ ] Covariance matrix output

**Success Metric**: RMS <1" con 100+ observations

### Phase 3: HIGH PRECISION (Settimana 3-4)
**Goal**: Sub-arcsecond precision

Tasks:
1. [ ] Integrate EOP (IERS data)
2. [ ] Add asteroid perturbations
3. [ ] Relativistic corrections
4. [ ] Full validation suite

**Success Metric**: RMS <0.1" matching JPL

---

## Validation Strategy

### Test Cases
1. **433 Eros** (16113 obs, well-known orbit)
   - Reference: JPL Horizons #25
   - Target: RMS <0.15" (OrbFit achieves 0.12")

2. **1 Ceres** (250k+ obs, benchmark)
   - Reference: AstDyS
   - Target: RMS <0.10"

3. **2010 TK7** (Earth Trojan, difficult)
   - Reference: JPL SBDB
   - Target: Convergence + RMS <0.5"

### Metrics
- **Preliminary Orbit Error**: Δi, Δa, Δe vs reference
- **RMS Residuals**: arcsec after differential correction
- **Propagation Error**: arcsec/year vs JPL
- **Computation Time**: seconds per fit

---

## Code Locations

### Files to Modify
```
src/initial_orbit.cpp           # Herget + differential corrector
include/ioccultcalc/initial_orbit.h
src/force_model.cpp             # Add perturbations
src/coordinates.cpp             # EOP integration
tests/test_eros_complete.cpp    # Validation suite
```

### New Files Needed
```
src/eop_reader.cpp              # Earth Orientation Parameters
src/asteroid_masses.cpp         # Load & query masses
src/relativity.cpp              # Post-Newtonian corrections
data/eop/                       # IERS Bulletin A/B
data/masses/                    # Asteroid mass database
```

---

## Resources

### Papers
1. Milani & Gronchi (2010) "Theory of Orbit Determination" - Bible
2. Folkner et al. (2014) "JPL Planetary Ephemerides" - DE430/431
3. Farnocchia et al. (2015) "Near-Earth Objects orbit determination" 
4. Petit & Luzum (2010) "IERS Conventions" - EOP standard

### Software References
- **OrbFit 5.0**: `/Users/michelebigi/Astro/OrbFit/src/`
- **Find_Orb** (Project Pluto): Public domain, similar goals
- **SPICE Toolkit**: NASA/JPL reference

### Data Sources
- **IERS EOP**: https://www.iers.org/IERS/EN/DataProducts/EarthOrientationData
- **JPL Horizons**: https://ssd.jpl.nasa.gov/horizons/
- **MPC Database**: https://minorplanetcenter.net/data
- **AstDyS**: https://newton.spacedys.com/astdys/

---

## Timeline Estimate

| Phase | Duration | Precision Gain |
|-------|----------|----------------|
| Fix Herget | 1 settimana | 5° → 1° |
| Robust Fitting | 1 settimana | 1° → 0.5" |
| EOP Integration | 1 settimana | 0.5" → 0.2" |
| Perturbations | 1 settimana | 0.2" → 0.1" |
| **TOTAL** | **1 mese** | **5° → 0.1"** |

---

## Next Immediate Action

**START HERE** ⚡:
```bash
# 1. Fix Herget inclination bug
cd examples
./debug_herget_inclination

# 2. Compare with OrbFit Gauss
# (requires setup OrbFit environment)

# 3. Iterate until i_error < 1°
```

**Success Criteria**: Quando Herget dà i=11° ± 0.5° per Eros, possiamo passare a Phase 2.

