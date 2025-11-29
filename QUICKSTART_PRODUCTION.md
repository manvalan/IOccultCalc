# IOccultCalc - Quick Start Operativo

## 🚀 Installazione Rapida

```bash
# 1. Clone repository
git clone https://github.com/manvalan/IOccultCalc.git
cd IOccultCalc

# 2. Installa dipendenze e compila
./install_production.sh

# 3. Download catalogo Gaia Mag18 (9 GB - ESSENZIALE)
./download_gaia_cache.sh
```

## ⚡ Uso Rapido

### Ricerca Standard (Monthly)

```bash
# Cerca occultazioni gennaio 2026 con tutti gli asteroidi promettenti
./italoccultcalc.sh preset_production_monthly.oop
```

**Output**:
- `output/occultations_YYYYMMDD.txt` - Formato IOTA
- `output/occultations_YYYYMMDD.json` - JSON per API
- `output/occultations_YYYYMMDD.kml` - Google Earth

### Asteroide Singolo

```bash
# Crea preset per (4) Vesta
cat > vesta.oop << EOF
time.
    .start_date = '2026-01-01'
    .end_date = '2026-01-31'

asteroids.
    .selection_mode = 'list'
    .asteroid_list = '4'

search.
    .max_magnitude = 16.0
EOF

./italoccultcalc.sh vesta.oop
```

### Ricerca Annuale (Lunga)

```bash
# Modifica preset per anno intero
sed -i '' 's/2026-01-31/2026-12-31/' preset_production_monthly.oop

# Esegui (tempo stimato: ~1 ora con parallelizzazione)
nohup ./italoccultcalc.sh preset_production_monthly.oop > search_2026.log 2>&1 &

# Monitor progress
tail -f search_2026.log
```

## 📊 Performance

| Scenario | Asteroidi | Periodo | Tempo | Speedup |
|----------|-----------|---------|-------|---------|
| **Monthly** | 375 | 31 giorni | **7 min** | 3.4× |
| **Quarterly** | 375 | 90 giorni | ~20 min | 3.5× |
| **Annual** | 375 | 365 giorni | ~80 min | 3.3× |

*Hardware: MacBook Air M-series, 8 core, 16GB RAM*

## 🎯 Ottimizzazioni Attive

✅ **Query Batching** (93% riduzione query)
✅ **Adaptive Timestep** (coarse 6h → fine 30min)
✅ **Chebyshev Interpolation** (ordine 11, 1-day intervals)
✅ **OpenMP Parallelization** (query + detection su 10 core)
✅ **Catalogo Locale Mag18** (300M stelle, zero query online)
✅ **LinOccult Compatibility** (threshold, uncertainties, CDF)

## 📁 Struttura File

```
IOccultCalc/
├── build/examples/italoccultcalc      # Eseguibile
├── italoccultcalc.sh                  # Launcher
├── preset_production_monthly.oop      # Preset standard
├── output/                            # Risultati
│   ├── occultations_YYYYMMDD.txt
│   ├── occultations_YYYYMMDD.json
│   └── occultations_YYYYMMDD.kml
└── ~/.ioccultcalc/                    # Config
    ├── data/
    │   └── all_numbered_asteroids.json
    ├── ephemerides/
    │   ├── de440s.bsp
    │   └── codes_300ast_20100725.bsp
    └── presets/
```

## 🔧 Configurazione

### Modifica Date

```bash
# Esempio: febbraio 2026
sed -i '' 's/2026-01-01/2026-02-01/' preset.oop
sed -i '' 's/2026-01-31/2026-02-28/' preset.oop
```

### Filtra per Magnitudine

```bash
# Solo stelle luminose (< mag 12)
sed -i '' 's/max_magnitude = 16.0/max_magnitude = 12.0/' preset.oop
```

### Asteroidi Specifici

```bash
# Lista custom (Vesta, Ceres, Pallas)
cat > custom.oop << EOF
asteroids.
    .selection_mode = 'list'
    .asteroid_list = '4,1,2'
EOF
```

## 🐛 Troubleshooting

### Errore: Catalogo non trovato

```bash
✗ Gaia catalog not found: /Users/.../gaia_mag18.cat.gz

# Soluzione: Download catalogo
./download_gaia_cache.sh
# Oppure:
curl -L https://github.com/manvalan/IOC_GaiaLib/releases/download/v2.0.0/gaia_mag18.cat.gz \
     -o ~/catalogs/gaia_mag18.cat.gz
```

### Errore: SPICE kernel non trovato

```bash
SPK file not found: de440s.bsp

# Soluzione: Download ephemeris
cd ~/.ioccultcalc/ephemerides
curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
```

### Performance Lenta

```bash
# Verifica OpenMP attivo
./italoccultcalc.sh preset.oop 2>&1 | grep "Calcolo parallelo"

# Dovrebbe mostrare:
# "Calcolo parallelo su N core..."

# Se non compare, reinstalla libomp:
brew reinstall libomp
./install_production.sh
```

## 📚 Documentazione Avanzata

- **Performance**: `PERFORMANCE_OPTIMIZATION_STRATEGY.md`
- **Chebyshev**: `CHEBYSHEV_PERFORMANCE_ANALYSIS.md`
- **Gaia EDR3**: `GAIA_EDR3_QUICKREF.md`
- **Changelog**: `COMPLETION_CHECKLIST.md`

## 🆘 Supporto

- **Issues**: https://github.com/manvalan/IOccultCalc/issues
- **Discussions**: https://github.com/manvalan/IOccultCalc/discussions
- **Email**: michele.bigi@example.com (sostituire)

## 📈 Confronto con LinOccult

| Feature | LinOccult | IOccultCalc |
|---------|-----------|-------------|
| **Tempo (375 ast, 31d)** | 2-5 min | **7 min** |
| **Catalogo** | RAM preload | Disco + cache |
| **Ottimizzazioni** | Chebyshev | Chebyshev + OpenMP + Batching |
| **Memoria** | 12 GB RAM | 2 GB RAM |
| **Output** | IOTA | IOTA + JSON + KML |
| **LinOccult Compat** | - | ✅ 100% |

## 🎖️ Crediti

- **Steve Preston** (LinOccult) - Algoritmi base
- **IOTA** - Standard formato output
- **Gaia Mission** (ESA) - Catalogo stellare
- **JPL** - Ephemeris DE440/DE441
- **OrbFit Team** - Elementi orbitali

---

**Versione**: 1.0.0  
**Data**: Novembre 2025  
**Autore**: Michele Bigi  
**License**: MIT
