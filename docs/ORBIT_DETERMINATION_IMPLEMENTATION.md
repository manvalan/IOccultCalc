# Implementazione Orbit Determination - Piano di Lavoro

## Stato Attuale (Completato)

✅ **Task 1: Metodo di Gauss**
- Creato `include/ioccultcalc/initial_orbit.h` con API completa
- Implementato `src/initial_orbit.cpp` con:
  - `gauss()`: metodo di Gauss per 3 osservazioni
  - `findBestObsForGauss()`: selezione automatica osservazioni
  - `findRealPolynomialRoots()`: risolutore polinomiale 8° grado
  - Helper functions: `crossThenDot()`, verifiche validità
- Aggiornato CMakeLists.txt con nuovi file

### Note Implementative
- Struttura basata su Find_Orb (gauss.cpp)
- Serie f/g con termini fino al 5° ordine (opzionali)
- Gestione 0-3 soluzioni reali
- Verifica condizioni: d0 != 0, distanze positive, convergenza

### TODO per Completare Gauss
1. Implementare conversione RA/Dec → vettori unitari
2. Implementare conversione cartesiano → elementi orbitali
3. Implementare posizioni osservatore (geocentro/topocentric)
4. Test con osservazioni reali di 433 Eros

---

## Prossimi Task

### Task 2: Metodo di Herget (HIGH PRIORITY)

**Riferimenti:**
- Find_Orb: `orb_func.cpp`, funzione `herget_method()` (linee 1916-2167)
- Documentazione: https://www.projectpluto.com/herget.htm

**Algoritmo:**
```cpp
1. Input: R1, R2 (distanze topocentriche guess)
2. set_distance(obs[0], R1)
3. set_distance(obs[n-1], R2)
4. find_transfer_orbit(orbit, obs[0], obs[n-1]) → stato iniziale
5. set_locs(orbit, t0, obs, n_obs) → propaga e calcola posizioni
6. compute_residuals(obs) → xresid[], yresid[]

7. Per variare R1 e R2:
   delta_R1 = R1 / 10000 (max 0.1 AU)
   delta_R2 = R2 / 10000
   
8. Ricalcola con R1+delta, R2+delta
9. Derivate parziali:
   dx_dr1[i] = (dx_new - dx_old) / delta_R1
   dy_dr1[i] = (dy_new - dy_old) / delta_R1
   dx_dr2[i], dy_dr2[i] analogamente

10. Sistema lineare 2x2:
    a = Σ (dx_dr1² + dy_dr1²)
    b = Σ (dx_dr1*dx_dr2 + dy_dr1*dy_dr2)
    c = Σ (xresid*dx_dr1 + yresid*dy_dr1)
    e = Σ (dx_dr2² + dy_dr2²)
    f = Σ (xresid*dx_dr2 + yresid*dy_dr2)
    
    determ = a*e - b²
    d_r1 = -(e*c - b*f) / determ
    d_r2 = -(a*f - c*b) / determ

11. Return d_r1, d_r2
```

**Funzioni da Implementare:**
```cpp
// src/initial_orbit.cpp

// Setta distanza per osservazione
void setDistance(AstrometricObservation& obs, double distance) {
    // obj_posn = obs_posn + distance * unit_vector
}

// Trova orbita trasferimento
int findTransferOrbit(
    double* orbit,              // state vector [x,y,z,vx,vy,vz]
    const AstrometricObservation& obs1,
    const AstrometricObservation& obs2,
    bool already_have_approximate_orbit);

// Propaga orbita e calcola posizioni per tutte le osservazioni
int setLocs(
    const double* orbit,
    double t0,
    std::vector<AstrometricObservation>& observations);

// Calcola residui RA/Dec per un'osservazione
void getResidualData(
    const AstrometricObservation& obs,
    double* xresid,  // RA residual (arcsec)
    double* yresid); // Dec residual (arcsec)

// Metodo di Herget completo
OrbitalElements herget(
    const std::vector<AstrometricObservation>& observations,
    double r1_guess,
    double r2_guess,
    int max_iterations);
```

**Pseudo-Väisälä:**
Se `r1_guess < 0`, usa pseudo-Väisälä:
```cpp
r2 = find_r_given_solar_r(&obs2, -r1);
r1 = find_r_given_solar_r(&obs1, -r1);
```

**Constrained Orbits:**
Supportare vincoli tipo:
- `e=1` (parabolica)
- `e=0` (circolare)
- `a=2.5` (semiasse fisso)

---

### Task 3: find_transfer_orbit (CRITICAL per Herget)

**Riferimenti:**
- Find_Orb: `orb_func.cpp`, linee 1072-1154

**Algoritmo:**
```cpp
1. Input: orbit[6] (guess iniziale o zero), obs1, obs2
2. Se no guess:
   - Usa metodo semplificato (Euler, sezione conica diretta)
   
3. Loop di raffinamento (max 3 iterazioni):
   a. Integra orbit da t1 a t2
   b. Calcola differenza: delta = orbit_propagato - obs2.obj_posn
   c. Se |delta| < target_diff: CONVERGENZA
   d. Altrimenti aggiusta velocità:
      orbit[3:6] -= delta / delta_t

4. Raffinamento fine con derivate parziali:
   Per ogni componente velocità vx, vy, vz:
   - Perturba di h = 1e-5
   - Integra e calcola delta
   - delta[i][j] = (delta_new[j] - delta_old[j]) / h
   
5. Risolve sistema 3x3:
   discr = determinante(delta[1], delta[2], delta[3])
   orbit[3] -= determinante(delta[0], delta[2], delta[3]) / discr
   orbit[4] -= determinante(delta[1], delta[0], delta[3]) / discr
   orbit[5] -= determinante(delta[1], delta[2], delta[0]) / discr

6. Return 0 se successo, error code altrimenti
```

**Errori:**
- `XFER_NO_CONVERGENCE = -1`
- `XFER_INTEGRATION_FAILED = -2`

---

### Task 4: adjust_herget_results

**Riferimenti:**
- Find_Orb: `orb_func.cpp`, linee 1894-1913

**Algoritmo:**
```cpp
1. Filtra osservazioni valide (non escluse)
2. Verifica: almeno 2 osservazioni
3. Verifica: orbita ragionevole (isUnreasonableOrbit)
4. Verifica: span temporale < maxHergetSpan(r1, r2)
5. Chiama extended_orbit_fit con FIT_FIXED_DISTANCES
6. Chiama set_locs per aggiornare residui
7. Return 0 se successo
```

**maxHergetSpan:**
```cpp
// Span massimo in giorni basato su distanze
// Oggetti vicini cambiano rapidamente
double maxHergetSpan(double r1, double r2) {
    double max_speed = 2. * PI / sqrt((r1 + r2) / 2.);
    return 90.0 / max_speed;  // ~90° di movimento massimo
}
```

---

### Task 5: extended_orbit_fit (Least Squares)

**Riferimenti:**
- Find_Orb: `orb_func.cpp`, linee 1235-1342

**Fit Types:**
```cpp
enum FitType {
    FIT_CLASSIC_HERGET = 2,    // R1, R2
    FIT_HERGET_FULL = 6,       // R1, R2, dRA1, dDec1, dRA2, dDec2
    FIT_FIXED_DISTANCES = 4,   // dRA1, dDec1, dRA2, dDec2 (R1,R2 fissate)
    FIT_FULL_STEP = 6          // x, y, z, vx, vy, vz
};
```

**Algoritmo:**
```cpp
1. Integra orbit da epoch a obs[0].jd
2. Loop: i = -1 a n_params
   a. Se i >= 0: perturba params[i] di delta_val = 1e-6
   b. Integra orbit + perturbazione da obs[0].jd a epoch
   c. set_locs per calcolare posizioni previste
   d. Per ogni osservazione:
      - Calcola residui RA/Dec
      - Se i == -1: memorizza in resids[]
      - Se i >= 0: memorizza slopes = (resids_new - resids_old) / delta_val
      
3. Risolve least squares:
   lsquare = lsquare_init(n_params)
   per ogni osservazione:
       lsquare_add_observation(lsquare, resid, weight, slopes)
   lsquare_solve(lsquare, params)
   
4. Applica correzioni:
   Se FIT_CLASSIC_HERGET:
       R1 += params[0]
       R2 += params[1]
   Se FIT_FIXED_DISTANCES:
       Aggiusta solo RA/Dec residui
   Se FIT_FULL_STEP:
       orbit += params

5. Return 0 se successo
```

---

### Task 6: initial_orbit (Auto-Solve Completo)

**Riferimenti:**
- Find_Orb: `orb_func.cpp`, linee 4208-4480

**Strategia Multi-Metodo:**
```cpp
1. Preprocessing:
   - Ordina osservazioni escluse alla fine
   - Calcola arclen = last_jd - first_jd

2. Loop su lunghezze d'arco decrescenti:
   arclen_start = min(arclen, 90 giorni)
   
   Mentre arclen >= 3 giorni:
   
   a. Seleziona subarco di lunghezza arclen
      start = 0, end = trovato da jd
      n_subarc_obs = end - start + 1
   
   b. Se n_subarc_obs >= 3: Prova Gauss
      Per ogni soluzione (0-3):
          epoch = convenient_gauss(...)
          Se epoch valido:
              set_locs(orbit, epoch, obs, n_subarc_obs)
              score = evaluate_initial_orbit(...)
              Se score < 1000:
                  integrate_orbit a obs[start].jd
                  attempt_improvements(orbit, obs, n_subarc_obs)
                  score_new = evaluate_initial_orbit(...)
                  Se score_new < best_score:
                      best_score = score_new
                      best_orbit = orbit
   
   c. Se best_score >= 50: Prova Herget/Väisälä
      Per vari pseudo_r (0.004, 0.01, 0.1, 1.0, 3.0):
          herget_method(obs, n_subarc_obs, -pseudo_r, -pseudo_r, orbit)
          adjust_herget_results(...)
          score = evaluate_initial_orbit(...)
          Se score < best_score:
              best_score = score
              best_orbit = orbit
          
          Se best_score < acceptable_score_limit (5.0):
              attempt_improvements(...)
              BREAK
   
   d. Se best_score < acceptable_score_limit: BREAK
   
   e. arclen *= 0.7  // Prova arco più breve

3. Return best_score e best_orbit
```

**evaluate_initial_orbit:**
```cpp
double evaluate_initial_orbit(
    const OrbitalElements& orbit,
    const std::vector<AstrometricObservation>& observations,
    const JulianDate& epoch)
{
    // 1. set_locs per calcolare posizioni
    // 2. Calcola RMS residui
    // 3. Penalità per orbite non fisiche:
    //    - e >= 1.0: score *= 10
    //    - a < 0: score *= 100
    //    - Hitting planet: score = 1e10
    // 4. Return score (arcsec)
}
```

**attempt_improvements:**
```cpp
double attempt_improvements(
    OrbitalElements& orbit,
    std::vector<AstrometricObservation>& observations)
{
    double best_score = evaluate_initial_orbit(...)
    
    // Metodo 0: Herget steps (max 5 iter)
    for iter in 0..5:
        r1 = obs[0].r
        r2 = obs[n-1].r
        herget_method(obs, n_obs, r1, r2, orbit, &d_r1, &d_r2)
        r1 += d_r1, r2 += d_r2
        herget_method(obs, n_obs, r1, r2, orbit, NULL, NULL)
        adjust_herget_results(obs, n_obs, orbit)
        score = evaluate_initial_orbit(...)
        Se score >= curr_score: BREAK
        curr_score = score
    
    // Metodo 1: Full steps (max 5 iter)
    for iter in 0..5:
        full_improvement(obs, n_obs, orbit, obs[0].jd, NULL, 
                        NO_ORBIT_SIGMAS_REQUESTED, obs[0].jd)
        score = evaluate_initial_orbit(...)
        Se score >= curr_score: BREAK
        curr_score = score
    
    Return best_score
}
```

---

### Task 7: Sistema Sigmas e Weighting

**Estensione AstrometricObservation:**
```cpp
struct AstrometricObservation {
    // ... campi esistenti ...
    
    // Incertezze (default per CCD moderna)
    double sigma_ra = 0.5;   // arcsec
    double sigma_dec = 0.5;  // arcsec
    double sigma_time = 1.0; // seconds
    double sigma_mag = 0.5;  // magnitudes
    
    // Weighting
    double weight = 1.0;
    bool excluded = false;
    
    // Incertezze ellittiche (opzionale)
    double sigma_major = 0.0;  // Major axis (arcsec)
    double sigma_minor = 0.0;  // Minor axis (arcsec)
    double sigma_pa = 0.0;     // Position angle (deg)
    
    double getWeight() const {
        if (excluded) return 0.0;
        if (weight > 0) return weight;
        
        // Peso da sigma circolare
        if (sigma_major == 0.0) {
            double sigma2 = sigma_ra*sigma_ra + sigma_dec*sigma_dec;
            return 1.0 / sigma2;
        }
        
        // Peso da sigma ellittica
        double area = PI * sigma_major * sigma_minor;
        return 1.0 / area;
    }
};
```

**Assegnazione Automatica Sigmas:**
```cpp
void assignSigmas(std::vector<AstrometricObservation>& obs) {
    for (auto& o : obs) {
        // Basato su catalogo stelle
        if (o.catalogCode == "V") {  // Gaia
            o.sigma_ra = 0.05;  // 50 mas
            o.sigma_dec = 0.05;
        } else if (o.catalogCode == "U") {  // UCAC
            o.sigma_ra = 0.1;
            o.sigma_dec = 0.1;
        } else {  // Default
            o.sigma_ra = 0.5;
            o.sigma_dec = 0.5;
        }
        
        // Basato su precisione temporale
        // 0.01s → sigma_time = 1s
        // 0.001s → sigma_time = 0.1s
        
        // Basato su magnitudine
        if (o.magnitude < 15.0) {
            o.sigma_ra *= 0.5;  // Oggetti brillanti
            o.sigma_dec *= 0.5;
        } else if (o.magnitude > 20.0) {
            o.sigma_ra *= 2.0;  // Oggetti deboli
            o.sigma_dec *= 2.0;
        }
    }
}
```

---

### Task 8: Observation Filtering

**Filtri da Implementare:**

```cpp
class ObservationFilter {
public:
    // 1. Sigma-based rejection
    static void filterByResidual(
        std::vector<AstrometricObservation>& obs,
        double max_sigma = 3.0)
    {
        for (auto& o : obs) {
            if (o.outlier) continue;
            
            double sigma_resid_ra = o.raResidual / o.sigma_ra;
            double sigma_resid_dec = o.decResidual / o.sigma_dec;
            double total_sigma = std::sqrt(
                sigma_resid_ra*sigma_resid_ra + 
                sigma_resid_dec*sigma_resid_dec
            );
            
            if (total_sigma > max_sigma) {
                o.excluded = true;
            }
        }
    }
    
    // 2. Blunder management (probabilistic)
    static void applyBlunderWeighting(
        std::vector<AstrometricObservation>& obs,
        double blunder_probability = 0.02)  // 2%
    {
        for (auto& o : obs) {
            if (o.outlier || o.excluded) continue;
            
            double sigma_resid = o.totalResidual / 
                std::sqrt(o.sigma_ra*o.sigma_ra + o.sigma_dec*o.sigma_dec);
            
            // Peso ridotto per residui grandi
            // w = (1-p) * exp(-σ²/2) / (1-p + p*blunder_weight)
            double gaussian_weight = std::exp(-0.5 * sigma_resid * sigma_resid);
            double blunder_weight = 0.01;  // Peso per blunder
            
            o.weight = (1.0 - blunder_probability) * gaussian_weight /
                      ((1.0 - blunder_probability) * gaussian_weight + 
                       blunder_probability * blunder_weight);
        }
    }
    
    // 3. Over-observing (t-and-k scheme)
    static void applyOverObservingLimits(
        std::vector<AstrometricObservation>& obs,
        double time_range_days = 1.0,
        int max_obs_per_range = 5)
    {
        // Raggruppa osservazioni per intervalli temporali
        std::map<int, std::vector<AstrometricObservation*>> groups;
        
        for (auto& o : obs) {
            int group_id = (int)(o.epoch.jd / time_range_days);
            groups[group_id].push_back(&o);
        }
        
        // Per ogni gruppo con > max_obs_per_range
        for (auto& [id, group] : groups) {
            if (group.size() <= max_obs_per_range) continue;
            
            // Ordina per qualità (residuo/sigma)
            std::sort(group.begin(), group.end(),
                [](auto* a, auto* b) {
                    double qa = a->totalResidual / 
                               std::sqrt(a->sigma_ra*a->sigma_ra + a->sigma_dec*a->sigma_dec);
                    double qb = b->totalResidual / 
                               std::sqrt(b->sigma_ra*b->sigma_ra + b->sigma_dec*b->sigma_dec);
                    return qa < qb;
                });
            
            // Escludi le peggiori
            for (size_t i = max_obs_per_range; i < group.size(); i++) {
                group[i]->excluded = true;
            }
        }
    }
};
```

---

## Timeline Stimata

| Task | Priorità | Tempo Stimato | Dipendenze |
|------|---------|---------------|------------|
| ✅ Task 1: Gauss | HIGH | 4h | Nessuna |
| Task 2: Herget | HIGH | 6h | Task 3, coordinate utils |
| Task 3: find_transfer_orbit | CRITICAL | 4h | OrbitPropagator |
| Task 4: adjust_herget_results | HIGH | 2h | Task 2, Task 5 |
| Task 5: extended_orbit_fit | HIGH | 6h | Least squares library |
| Task 6: initial_orbit | MEDIUM | 4h | Task 2-5 |
| Task 7: Sigmas | MEDIUM | 2h | Nessuna |
| Task 8: Filtering | LOW | 3h | Task 7 |
| Test & Debug | HIGH | 8h | Tutti |
| **TOTALE** | | **39h** | |

---

## File da Creare/Modificare

### Nuovi File
- ✅ `include/ioccultcalc/initial_orbit.h`
- ✅ `src/initial_orbit.cpp`
- `src/coordinate_utils.cpp` (RA/Dec ↔ cartesiano)
- `src/least_squares.cpp` (solver LSQ generico)
- `examples/test_initial_orbit.cpp`

### File da Modificare
- ✅ `CMakeLists.txt`
- `include/ioccultcalc/observation.h` (sigmas, weights)
- `src/orbit_fitter.cpp` (integrazione con InitialOrbit)
- `include/ioccultcalc/orbit_propagator.h` (esporre integrate_orbit)

---

## Testing Strategy

### Test 1: Gauss con 3 Osservazioni
```cpp
// 433 Eros, 3 osservazioni ben separate (2023)
// Verifica: 0-3 soluzioni, soluzioni fisiche, convergenza
```

### Test 2: Herget con Arco Lungo
```cpp
// 433 Eros, 100+ osservazioni su 6 mesi
// Verifica: convergenza da vari R1/R2 guess
```

### Test 3: Auto-Solve
```cpp
// 433 Eros, archi di diverse lunghezze
// Verifica: scelta automatica metodo, robustezza
```

### Test 4: Filtering
```cpp
// Dataset con outlier artificiali
// Verifica: rigetto corretto, weighting appropriato
```

---

## Domande Aperte

1. **Libreria Least Squares:** Usare Eigen, Armadillo, o implementazione custom?
2. **Polynomial Root Finder:** Jenkins-Traub? Laguerre? GSL?
3. **Coordinate Transforms:** Includere nutazione/precessione o solo rotazione semplice?
4. **Perturbatori:** Quali asteroidi includere nel calcolo (Top 300 di massa)?
5. **Ephemeris Backend:** SPICE, JPL DE, o calcolo semplificato?

---

## Riferimenti Completi

- **Find_Orb:** https://github.com/Bill-Gray/find_orb
- **Herget Math:** https://www.projectpluto.com/herget.htm
- **Väisälä Math:** https://www.projectpluto.com/vaisala.htm
- **Boulet Book:** "Methods of Orbit Determination for the Micro Computer"
- **Fundamentals:** Vallado, "Fundamentals of Astrodynamics and Applications"
