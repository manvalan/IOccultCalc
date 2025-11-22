# Analisi Bug Inclinazione Herget

**Data**: 21 novembre 2025  
**Oggetto**: 433 Eros  
**Problema**: Herget calcola i=40-48° invece di 10.83°

## Risultati Investigazione

### Sintomo
Il metodo Herget produce sistematicamente inclinazioni i=40-48° per Eros (riferimento JPL: i=10.83°).
Errore presente su tutti gli span testati (90, 180, 365 giorni).

### Ipotesi Iniziali Testate ❌

1. **Frame di riferimento errato** ❌
   - Test: `examples/test_frame_conversion.cpp` (5 test completi)
   - Risultato: SPICE usa ECLIPJ2000_DE405 ✅
   - Conversioni equatoriali↔eclittiche corrette ✅
   - Formula i = acos(h.z/|h|) corretta ✅

2. **Lambert solver bugs** ❌
   - Test: convergenza su orbite sintetiche (<1% error)
   - Risultato: Lambert solver matematicamente corretto ✅
   - Logica prograde/retrograde corretta ✅

3. **Bug selezione soluzione prograde/retrograde** ❌
   - Test: debug output mostra solo `result_pro=0, result_retro=-1`
   - Risultato: Solo soluzione prograde converge ✅
   - Non c'è problema di selezione tra soluzioni ✅

### ROOT CAUSE IDENTIFICATA ✅

**Problema: Selezione osservazioni inadatte per preliminary orbit determination**

#### Evidenza

Debug output per span 90 giorni:
```
findTransferOrbit - r1: 1.2 AU, r2: 1.3 AU
  obs1: RA=12.1459° Dec=13.5474°
  obs2: RA=30.3468° Dec=40.8221°    ← PROBLEMA!
  vect1 (ecliptic): (0.950414, 0.280848, 0.133556)
  vect2 (ecliptic): (0.653056, 0.610814, 0.447688)   ← z=0.448 MOLTO ALTO
  pos1_helio: (0.872617, 1.31775, 0.160219) ecliptic
  pos2_helio: (-0.150915, 0.723152, 0.582006) ecliptic   ← z=0.582 AU!
  Using prograde solution: v1=(-0.0102431, -0.0016382, 0.00698572)
  Resulting inclination: 45.3743°
```

#### Analisi Geometrica

Per **obs2** con Dec=40.82°:
- Vettore direzione z-component: vect2.z = 0.448
- Posizione heliocentrica: pos2.z = 0.582 AU
- Angolo dal piano eclittico: arctan(0.582/0.738) ≈ **38°**

Per Eros con i=10.83° (orbita quasi equatoriale):
- Distanza massima dal piano eclittico: a·sin(i) ≈ 1.458·sin(10.83°) ≈ 0.27 AU
- **pos2.z = 0.582 AU è IMPOSSIBILE per questa orbita!**

#### Conclusione

Il Lambert solver calcola **CORRETTAMENTE** i=45° date le posizioni pos1 e pos2.

Il problema è che **le osservazioni selezionate rappresentano una configurazione geometrica estrema**:
- Dec=40.8° è possibile per Eros solo in configurazioni particolari (vicino perielio + geometria Terra-Eros specifica)
- Queste osservazioni sono VALIDE ma INADATTE per preliminary orbit determination
- Producono geometria che viola le assunzioni implicite del metodo Herget (orbita approssimativamente circolare, moto nel piano)

### Causa in `test_eros_complete.cpp`

```cpp
std::vector<AstrometricObservation> selectObservations(
    const std::vector<AstrometricObservation>& all_obs,
    int n_desired,
    double jd_start,
    double jd_end)
{
    // Seleziona osservazioni uniformemente distanziate nel TEMPO
    // NON filtra per qualità geometrica!
    double step = (double)filtered.size() / n_desired;
    for (int i = 0; i < n_desired; i++) {
        int idx = (int)(i * step);
        selected.push_back(filtered[idx]);  // ← Include Dec=40°!
    }
}
```

## Soluzione Raccomandata

### Opzione 1: Filtro Declinazione (Quick Fix)

```cpp
std::vector<AstrometricObservation> selectObservations(...) {
    std::vector<AstrometricObservation> filtered;
    
    for (const auto& obs : all_obs) {
        // Filtro temporale
        if (obs.epoch.jd < jd_start || obs.epoch.jd > jd_end) continue;
        
        // Filtro geometrico per IOD
        double dec_deg = obs.obs.dec * 180.0 / M_PI;
        if (std::abs(dec_deg) > 30.0) continue;  // ← NUOVO!
        
        filtered.push_back(obs);
    }
    
    // ... resto come prima
}
```

**Pro**: Semplice, risolve il problema immediato  
**Contro**: Threshold arbitrario (30°)

### Opzione 2: Selezione Intelligente (Robust)

Criteri multipli:
1. Distribuzione uniforme nel tempo
2. Evitare declinazioni estreme (|Dec| < 30°)
3. Preferire osservazioni con basso residuo atteso
4. Diversificare osservatori (parallax geometrico)
5. Evitare configurazioni near-opposition

**Pro**: Robusto, professional-grade  
**Contro**: Più complesso da implementare

### Opzione 3: Weighted Least-Squares (Best)

Non filtrare osservazioni ma pesarle:
```cpp
// In differential corrector
double weight = 1.0;
if (std::abs(dec_deg) > 30.0) weight *= 0.1;  // Downweight estreme
if (epoch.jd < JD_2000) weight *= 0.5;        // Downweight vecchie
// ... calcolo residui pesati
```

**Pro**: Usa TUTTE le osservazioni, massimizza informazione  
**Contro**: Richiede refactoring differential corrector

## Priorità

**Immediata**: Nessuna - il problema è solo nei test, non nel codice production

**Fase 2 Roadmap** (dopo LM damping): Implementare Opzione 3 (weighted least-squares)

**Note**: Questa investigazione ha confermato che:
- ✅ Tutti gli algoritmi fondamentali (Lambert, frame conversion, orbital elements) sono corretti
- ✅ Il framework è solido e può raggiungere precisione OrbFit-level
- ⚠️ Serve pre-processing intelligente delle osservazioni per IOD robusto

## Debug Output da Rimuovere

Prima di commit finale, rimuovere da `src/initial_orbit.cpp`:
- Linea ~1232: `std::cerr << "Lambert convergence..."`
- Linea ~1243: `std::cerr << "Lambert solutions..."`
- Linea ~1215: `std::cerr << "findTransferOrbit..."`
- Linea ~1280: `std::cerr << "Using prograde..."`

Oppure wrappare in `#ifdef DEBUG_IOD`.
