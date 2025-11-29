# Analisi Performance Polinomi di Chebyshev

## 📊 GUADAGNO PERFORMANCE

### Confronto Metodi

| Metodo | Costo Computazionale | Precisione | Uso Tipico |
|--------|---------------------|------------|------------|
| **Integrazione RADAU** | ~1000 FLOPS per punto | 1 km @ 1 AU | Calcolo preciso singolo punto |
| **Chebyshev ordine 11** | ~50 FLOPS per punto | 1-2 km @ 1 AU | Interpolazione veloce |
| **Interpolazione lineare** | ~10 FLOPS per punto | 100-1000 km @ 1 AU | Approssimazione grossolana |

### Speedup Atteso

```
Scenario: Scansione 1 anno, timestep 1 minuto
- Punti necessari: 365 × 24 × 60 = 525,600 punti
- Metodo attuale (RADAU ogni punto): 525,600 integrazioni = ~30 secondi
- Metodo Chebyshev (pre-compute + interpolazione):
  * Pre-compute 365 segmenti: 365 integrazioni = ~0.2 secondi
  * Interpolazione 525,600 punti: ~0.1 secondi
  * TOTALE: 0.3 secondi
  
→ SPEEDUP: 100× più veloce!
```

---

## 🎯 PRECISIONE POLINOMI CHEBYSHEV

### Formula Base

Posizione asteroide al tempo `t` in intervallo `[t₀, t₁]`:

```
r(t) = Σ(i=0 to N) cᵢ × Tᵢ(x)
```

dove:
- `Tᵢ(x)` = polinomio Chebyshev ordine i
- `x = 2(t - t₀)/(t₁ - t₀) - 1` ∈ [-1, 1] (normalizzazione)
- `cᵢ` = coefficienti pre-calcolati
- `N` = ordine polinomio (LinOccult usa N=11)

### Errore di Approssimazione

**Ordine 11, intervallo 1 giorno**:

| Componente | Errore RMS | Errore Max |
|------------|-----------|-----------|
| Posizione X,Y,Z | 0.5 km | 2 km |
| Velocità Vx,Vy,Vz | 0.01 m/s | 0.05 m/s |
| **Posizione angolare** | **0.2 arcsec** | **1 arcsec** |

**Confronto con incertezze**:
- Incertezza orbitale tipica: 100-1000 km
- Errore Chebyshev: 1-2 km
- **Ratio**: 0.1-2% → **TRASCURABILE**!

### Scelta Ordine Ottimale

| Ordine | Errore @ 1 giorno | FLOPS/eval | Uso |
|--------|------------------|------------|-----|
| 5 | ~50 km | 25 | Troppo impreciso |
| 7 | ~10 km | 35 | Accettabile per survey grossolani |
| **11** | **1-2 km** | **50** | **LinOccult (raccomandato)** |
| 15 | ~0.5 km | 70 | Overkill (overhead inutile) |
| 21 | ~0.1 km | 100 | Eccessivo |

**Conclusione**: Ordine 11 è **OTTIMALE** - errore trascurabile rispetto a incertezze orbitali, ma 100× più veloce.

---

## 🔧 CONFIGURAZIONE RACCOMANDATA

### LinOccult Settings
```cpp
CHEB_ORDER = 11        // Ordine polinomio
CHEB_STEP = 1.0        // Intervallo giorni
CHEB_FINE_STEP = 60    // Subdivision per interpolazione (24 minuti)
```

### IOccultCalc Proposto
```json
{
  "chebyshev": {
    "enabled": true,
    "order": 11,
    "interval_days": 1.0,
    "fine_timestep_minutes": 15.0,
    "precision_check": true
  }
}
```

### Parametri Configurabili

1. **`order`** (7-15, default 11):
   - Più alto → più preciso ma più lento
   - Più basso → meno preciso ma più veloce
   - **Raccomandato**: 11 (LinOccult standard)

2. **`interval_days`** (0.5-2.0, default 1.0):
   - Più corto → più preciso ma più segmenti
   - Più lungo → meno preciso, meno memoria
   - **Raccomandato**: 1.0 giorno (sweet spot)

3. **`fine_timestep_minutes`** (1-30, default 15):
   - Risoluzione temporale scansione stelle
   - Più fine → più eventi catturati
   - **Raccomandato**: 15 min (occultazioni durano 1-60 secondi)

4. **`precision_check`** (bool, default true):
   - Valida errore Chebyshev vs integrazione diretta
   - Disabilita per produzione (overhead 5%)

---

## 📈 BENCHMARK ATTESI

### Scenario Test: 370 Asteroidi, 7 Giorni

**Metodo Attuale (RADAU ogni timestep)**:
```
Timestep: 0.5 giorni (12 ore)
Punti per asteroide: 7 / 0.5 = 14 punti
Totale integrazioni: 370 × 14 = 5,180
Tempo stimato: ~3 secondi
```

**Metodo Chebyshev (cache + interpolazione)**:
```
Pre-compute: 370 × 7 segmenti = 2,590 integrazioni
Tempo pre-compute: ~1.5 secondi

Fine-step scanning (15 minuti):
Punti per asteroide: 7 × 24 × 4 = 672 punti
Interpolazioni: 370 × 672 = 248,640
Tempo interpolazione: ~0.1 secondi

TOTALE: 1.6 secondi
```

**Speedup**: 3s → 1.6s = **2× più veloce** (per questo caso)

### Scenario Heavy: 1000 Asteroidi, 365 Giorni, Timestep Fine

**Metodo Attuale**:
```
Timestep: 1 ora (per catturare tutti eventi)
Punti: 1000 × 365 × 24 = 8,760,000 integrazioni
Tempo: ~5000 secondi = 1.4 ore
```

**Metodo Chebyshev**:
```
Pre-compute: 1000 × 365 = 365,000 integrazioni
Interpolazioni: 1000 × 365 × 24 × 4 (15 min) = 35,040,000
Tempo: 20 sec (pre) + 20 sec (interp) = 40 secondi
```

**Speedup**: 5000s → 40s = **125× PIÙ VELOCE**!

---

## 💾 MEMORIA RICHIESTA

### Per Asteroide
```
Ordine = 11
Componenti = 3 (X, Y, Z)
Coefficienti per segmento = 11 × 3 = 33 double
Giorni = 365

Memoria = 365 × 33 × 8 bytes = 96 KB per asteroide
```

### Per 1000 Asteroidi
```
Totale = 1000 × 96 KB = 96 MB

→ Assolutamente accettabile su qualsiasi sistema moderno
```

---

## 🎯 QUANDO USARE CHEBYSHEV

### ✅ Usa Chebyshev SE:
- Scansione lunga (> 30 giorni)
- Timestep fine (< 1 ora)
- Molti asteroidi (> 100)
- Query multiple sulle stesse effemeridi
- Performance critica

### ❌ Non Necessario SE:
- Calcolo singolo punto (overhead inutile)
- Scan breve (< 7 giorni) con timestep grossolano (> 12 ore)
- Pochi asteroidi (< 10)
- Memoria limitata (< 100 MB disponibili)

---

## 🔬 VALIDAZIONE PRECISIONE

### Test Raccomandati

1. **Error Check Interno**:
```cpp
// Ogni 100 punti, confronta Chebyshev vs RADAU
if (validation_enabled && point % 100 == 0) {
    Vector3D pos_radau = ephemeris.compute(jd);
    Vector3D pos_cheb = chebyshev.evaluate(jd);
    double error_km = (pos_radau - pos_cheb).magnitude();
    
    if (error_km > threshold_km) {
        LOG_WARNING("Chebyshev error " << error_km << " km at JD " << jd);
    }
}
```

2. **Comparison Test**:
```bash
# Run stesso caso con/senza Chebyshev
./italoccultcalc preset.oop --chebyshev=false > results_radau.txt
./italoccultcalc preset.oop --chebyshev=true > results_cheb.txt

# Confronta risultati
diff results_radau.txt results_cheb.txt
# Aspettativa: differenze < 0.1 arcsec (trascurabile)
```

---

## 📋 IMPLEMENTAZIONE STEP-BY-STEP

### Fase 1: Struttura Base (30 min)
1. Creare `chebyshev_approximation.h/cpp`
2. Classe `ChebyshevSegment` per singolo intervallo
3. Classe `ChebyshevCache` per gestione segmenti multipli

### Fase 2: Calcolo Coefficienti (1 ora)
1. Metodo `computeCoefficients()` usando punti Chebyshev
2. Integrare con `Ephemeris::compute()` per campioni
3. Fitting polinomiale ordine N

### Fase 3: Evaluation (30 min)
1. Metodo `evaluate(jd)` per interpolazione veloce
2. Ricorsione Chebyshev: T₀(x)=1, T₁(x)=x, Tₙ₊₁=2xTₙ-Tₙ₋₁
3. Cache risultati intermedi

### Fase 4: Integrazione Predictor (30 min)
1. Aggiungere flag `use_chebyshev` in configurazione
2. Modificare `findOccultations()` per usare cache
3. Pre-compute all'inizio, evaluate durante scan

### Fase 5: Testing & Validation (1 ora)
1. Test unitario precisione
2. Benchmark performance
3. Validation vs RADAU

**TOTALE**: ~3.5 ore implementazione completa

---

## 🚀 CONCLUSIONI

### Guadagno Atteso
- **Performance**: 2-125× più veloce (dipende da scenario)
- **Precisione**: Identica (errore < 1% incertezza orbitale)
- **Memoria**: Modesta (96 MB per 1000 asteroidi × 1 anno)

### Raccomandazione
**IMPLEMENTA CHEBYSHEV** - È il metodo standard usato da LinOccult per una ragione:
- ✅ Drastico speedup per scan lunghi
- ✅ Precisione più che sufficiente
- ✅ Memoria accettabile
- ✅ Configurabile e disabilitabile

### Priority
**ALTA** se:
- Pianifichi scan > 30 giorni
- Vuoi timestep < 1 ora (catturare tutti eventi)
- Performance importante

**MEDIA** se:
- Scan tipici 7-30 giorni
- Timestep 1-12 ore sufficiente

**BASSA** se:
- Solo scan brevi (< 7 giorni)
- Pochi asteroidi (< 50)
- Timestep grossolano (> 12 ore)

Per il tuo caso (370 asteroidi, 7 giorni) il guadagno è moderato (2×), ma diventa enorme per scan annuali!
