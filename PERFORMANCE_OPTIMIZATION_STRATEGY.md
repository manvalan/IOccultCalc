# IOccultCalc - Performance Optimization Strategy

## Situazione Attuale (27 Nov 2025)

**Test in corso**: 1000 asteroidi gennaio 2026
- Stato: 40% completato (query Gaia online in corso)
- Tempo stimato: 15-25 minuti totali
- **Problema**: Sta usando query online invece di catalogo Mag18 locale

## Confronto con LinOccult 2.3.0

| Metrica | LinOccult | IOccultCalc (attuale) | Target |
|---------|-----------|----------------------|--------|
| **Caricamento catalogo** | 30s (9GB→RAM) | 0s (on-demand) | 30s |
| **Ricerca 1000 ast × 31 giorni** | 2-5 min | 15-25 min | 5-8 min |
| **Memoria RAM** | 12 GB | 500 MB | 10 GB |
| **Precisione orbitale** | 2-body Kepler | N-body RA15 | N-body |
| **Ottimizzazioni LinOccult** | 7/7 | 8/8 ✓ | 8/8 ✓ |

## Bottleneck Identificati

### 1. **Query Catalogo Stelle** (70% del tempo)
**Problema**: 
- Query cone-by-cone online/locale
- ~300-400 query per scansione completa
- Ogni query: accesso disco, decompressione, filtraggio

**LinOccult Strategy**:
```
[Startup] Load full catalog → 9GB RAM (30s)
[Search]  In-memory lookups → HEALPix spatial index (istantaneo)
```

**IOccultCalc Current**:
```
[Search] Query cone → decompress chunk → filter → return (50-500ms/query)
```

**Ottimizzazione Proposta**:
```cpp
// OPZIONE A: Preload catalogo completo (come LinOccult)
void preloadFullCatalog() {
    // Carica tutti i chunk HEALPix in RAM
    // Pro: velocità massima
    // Contro: 9-12 GB RAM
}

// OPZIONE B: Preload regioni interessate
void preloadRegions(std::vector<HEALPixRegion> regions) {
    // Carica solo HEALPix necessari per gli asteroidi selezionati
    // Pro: memoria proporzionale agli asteroidi
    // Contro: richiede pre-analisi traiettorie
}

// OPZIONE C: Query batching + caching intelligente
void batchedConeQuery(std::vector<ConeQuery> queries) {
    // Raggruppa query vicine
    // Cache chunks HEALPix già letti
    // Pro: compromesso memoria/velocità
    // Contro: complessità intermedia
}
```

### 2. **Propagazione Orbitale** (25% del tempo)
**Problema**:
- RA15 integrator molto preciso ma lento
- Ogni asteroide: ~31 step giornalieri × calcolo completo

**Current Optimization**:
- ✅ Chebyshev approximation ordine 11 (già implementato)
- Guadagno: 2-10× su scan lunghi
- Precisione: <2 km (trascurabile)

**Ottimizzazioni Aggiuntive**:
```cpp
// OPZIONE D: Adaptive step-size (già implementato ma non attivo)
// Attivare use_adaptive_timestep = true nel preset

// OPZIONE E: Parallel orbit propagation
#pragma omp parallel for
for (int i = 0; i < asteroids.size(); i++) {
    propagate(asteroids[i]);
}
// Guadagno teorico: N_cores × (8 core → 8×)
```

### 3. **Detection Occultazioni** (5% del tempo)
**Non critico** - già ottimizzato con threshold LinOccult

## Piano di Ottimizzazione (3 Fasi)

### **FASE 1: Quick Wins (1-2 ore implementazione)**
**Target: 15-25 min → 8-12 min**

1. ✅ Fix configurazione catalogo Mag18 locale
   - Aggiungere parametro `gaia.catalog_file = 'gaia_mag18.cat.gz'`
   - Eliminare query online

2. 🔄 Attivare adaptive timestep
   - Già implementato, solo configurazione
   - Guadagno: 10-20%

3. 🔧 Query batching semplice
   ```cpp
   // Raggruppa asteroidi in regioni HEALPix
   // 1 query per regione invece di N query per asteroide
   ```
   - Guadagno stimato: 30-40%

### **FASE 2: Catalog Preloading (4-6 ore)**
**Target: 8-12 min → 4-6 min**

Implementare sistema simile a LinOccult:

```cpp
class GaiaCatalogPreloader {
    std::unordered_map<uint64_t, std::vector<GaiaStar>> healpix_cache;
    
    void preloadRegions(std::vector<SkyRegion> regions) {
        // Identifica HEALPix necessari
        auto healpix_ids = computeHEALPixCoverage(regions, NSIDE=64);
        
        // Carica chunks in parallelo
        #pragma omp parallel for
        for (auto hpx_id : healpix_ids) {
            healpix_cache[hpx_id] = loadHEALPixChunk(hpx_id);
        }
    }
    
    std::vector<GaiaStar> queryCone(double ra, double dec, double radius) {
        // Query in-memory - istantaneo
        auto hpx_ids = cone_to_healpix(ra, dec, radius);
        std::vector<GaiaStar> results;
        for (auto id : hpx_ids) {
            if (healpix_cache.count(id)) {
                filter_and_append(healpix_cache[id], results, ra, dec, radius);
            }
        }
        return results;
    }
};
```

**Strategia**:
1. All'avvio: pre-analisi traiettorie approssimate
2. Identifica HEALPix necessari (10-50 su 12288 totali)
3. Carica solo quelli in RAM (~500 MB - 2 GB)
4. Query istantanee durante ricerca

**Guadagno**: 3-5× su query (70% del tempo → 50%)

### **FASE 3: Parallelizzazione (2-3 ore)**
**Target: 4-6 min → 2-3 min**

```cpp
// Parallelize asteroid loop
#pragma omp parallel for schedule(dynamic)
for (size_t i = 0; i < asteroids.size(); i++) {
    // Thread-safe processing
    auto events = processAsteroid(asteroids[i], star_cache);
    
    #pragma omp critical
    {
        all_events.insert(all_events.end(), events.begin(), events.end());
    }
}
```

**Considerazioni**:
- GaiaCatalog deve essere thread-safe
- Ephemeris cache thread-local
- Critical section solo per risultati finali

**Guadagno**: Teorico 8× su 8 core, reale ~4-5× per overhead

## Stima Tempi Finali

| Scenario | Tempo 1000 ast × 31 giorni | RAM | Implementazione |
|----------|---------------------------|-----|-----------------|
| **Attuale** | 15-25 min | 500 MB | ✓ Completo |
| **Fase 1** | 8-12 min | 500 MB | 2 ore |
| **Fase 2** | 4-6 min | 2-4 GB | 6 ore |
| **Fase 3** | 2-3 min | 2-4 GB | 3 ore |
| **LinOccult** | 2-5 min | 12 GB | (riferimento) |

## Trade-off Memoria vs Velocità

```
Strategy A: Full Preload (LinOccult-like)
├─ PRO: Velocità massima (2-3 min)
├─ PRO: Semplice da implementare
├─ CONTRO: 12 GB RAM richiesti
└─ CONTRO: Startup 30s

Strategy B: Selective Preload (IOccultCalc optimized)
├─ PRO: 2-4 GB RAM (acceptable)
├─ PRO: No startup overhead
├─ PRO: Scalabile (memoria ∝ asteroidi)
└─ CONTRO: Richiede pre-analysis

Strategy C: Smart Caching (attuale + batching)
├─ PRO: 500 MB RAM
├─ PRO: Zero modific he architetturali
├─ CONTRO: 8-12 min (2-4× più lento di target)
└─ USO: Sistemi con poca RAM
```

## Prossimi Passi

1. **Immediato**: Attendere completamento test attuale per baseline precisa
2. **Fix catalogo**: Correggere preset per usare Mag18 locale (5 min)
3. **Fase 1**: Implementare quick wins (2 ore)
4. **Benchmarking**: Misurare guadagni reali
5. **Decisione**: Fase 2 vs Fase 3 in base a requisiti utente

## Note Tecniche

### HEALPix Coverage Calculation
```cpp
// Dato un asteroide e periodo, calcola HEALPix necessari
std::set<uint64_t> computeAsteroidHEALPixCoverage(
    const Asteroid& ast, 
    JulianDate start, 
    JulianDate end,
    int nside = 64
) {
    std::set<uint64_t> healpix_ids;
    
    // Propaga con step grossolano per coverage
    for (double jd = start.jd; jd <= end.jd; jd += 1.0) {
        auto pos = propagate(ast, jd);
        
        // Converti RA/Dec in HEALPix
        auto hpx = radec_to_healpix(pos.ra, pos.dec, nside);
        
        // Aggiungi HEALPix + neighbors (per margine sicurezza)
        healpix_ids.insert(hpx);
        for (auto neighbor : healpix_neighbors(hpx)) {
            healpix_ids.insert(neighbor);
        }
    }
    
    return healpix_ids;
}
```

### Memory Requirements
```
Mag18 Full Catalog: 9 GB compressed → ~12 GB decompressed
Per HEALPix (NSIDE=64): ~12 GB / 12288 = ~1 MB ciascuno
100 asteroidi coverage: ~200-500 HEALPix → 200-500 MB
1000 asteroidi coverage: ~2000-5000 HEALPix → 2-5 GB
```

## Conclusione

**Raccomandazione**: Implementare **Fase 1 + Fase 2** per raggiungere 4-6 minuti con 2-4 GB RAM.

Questo mantiene un ottimo compromesso:
- 2-3× più veloce dello stato attuale
- Comparabile con LinOccult
- Requisiti RAM ragionevoli (2-4 GB invece di 12 GB)
- Precisione orbitale superiore (N-body vs 2-body)

