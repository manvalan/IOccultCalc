# 🌟 IOccultCalc - Production Ready v1.0

**Predizione occultazioni asteroidali ottimizzato per performance**

Compatibile 100% con LinOccult | Performance: 3.4× più veloce del baseline | Dati reali MPC

---

## ✨ Features

- ✅ **Dati Reali MPC**: Elementi orbitali da Minor Planet Center (600k+ asteroidi)
- ✅ **Catalogo Gaia DR3 Mag18**: 303 milioni di stelle (G ≤ 18)
- ✅ **Performance Ottimizzate**: Query batching + OpenMP parallelization
- ✅ **LinOccult Compatible**: Stessi algoritmi (threshold, CDF, uncertainties)
- ✅ **Output Multipli**: IOTA, JSON, KML, CSV
- ✅ **Zero Query Online**: Tutto offline dopo setup iniziale

## 🚀 Quick Start (5 minuti)

```bash
# 1. Clone
git clone https://github.com/manvalan/IOccultCalc.git
cd IOccultCalc

# 2. Install + Build
./install_production.sh

# 3. Download dati MPC (172 MB → 1.5 GB)
./download_mpc_data.sh

# 4. Download catalogo Gaia (9 GB)
./download_gaia_cache.sh

# 5. Run!
./italoccultcalc.sh preset_production_monthly.oop
```

**Tempo totale setup**: ~30 minuti (download dipendenti)

## 📊 Performance

### Benchmark (MacBook Air M-series, 8 core)

| Scenario | Asteroidi | Periodo | Stelle | Tempo | vs Baseline |
|----------|-----------|---------|--------|-------|-------------|
| Monthly | 375 | 31 giorni | 448 | **7.0 min** | **3.4× faster** |
| Quarterly | 375 | 90 giorni | ~800 | ~20 min | 3.5× faster |
| Annual | 375 | 365 giorni | ~1500 | ~80 min | 3.3× faster |
| Single (Vesta) | 1 | 31 giorni | 12 | ~15 sec | - |

### Ottimizzazioni Attive

| Ottimizzazione | Tecnologia | Speedup |
|----------------|------------|---------|
| Query Batching | Spatial clustering (375→28 regioni) | 93% riduzione |
| Parallel Query | IOC_GaiaLib OpenMP | 1.5× |
| Parallel Detection | OpenMP asteroid loop (10 core) | 6× |
| Adaptive Timestep | LinOccult algorithm | 10-20% |
| Chebyshev | Polynomial interpolation (order 11) | 2× (opzionale) |

**TOTALE**: 3.4× speedup complessivo

## 📁 Struttura Dati

```
~/.ioccultcalc/
├── data/
│   └── all_numbered_asteroids.json    # MPC: 600k asteroidi (1.5 GB)
├── ephemerides/
│   ├── de440s.bsp                     # Pianeti (32 MB)
│   └── codes_300ast_20100725.bsp      # Asteroidi massivi (9 MB)
└── presets/
    └── standard.oop

~/catalogs/
└── gaia_mag18.cat.gz                  # Gaia DR3 (9 GB, 303M stelle)

IOccultCalc/
├── build/examples/italoccultcalc      # Eseguibile
├── italoccultcalc.sh                  # Launcher
├── output/                            # Risultati
│   ├── occultations_YYYYMMDD.txt      # IOTA format
│   ├── occultations_YYYYMMDD.json     # API/web
│   └── occultations_YYYYMMDD.kml      # Google Earth
└── preset_production_monthly.oop      # Config standard
```

## 🎯 Uso Operativo

### Ricerca Mensile Standard

```bash
./italoccultcalc.sh preset_production_monthly.oop
```

**Output**:
- ~5-20 eventi per mese (dipende da geometria)
- Filtro geografico Italia (36°N-47°N, 6°E-19°E)
- Mag drop > 0.5, Duration > 0.1s

### Personalizzazione Date

```bash
# Febbraio 2026
cat > feb2026.oop << EOF
time.
    .start_date = '2026-02-01'
    .end_date = '2026-02-28'

asteroids.
    .selection_mode = 'auto'
    .max_magnitude = 16.0

search.
    .max_magnitude = 16.0
    .search_radius_deg = 0.1

# Include preset_production_monthly.oop settings
include preset_production_monthly.oop
EOF

./italoccultcalc.sh feb2026.oop
```

### Asteroide Specifico

```bash
# (4) Vesta gennaio 2026
./italoccultcalc.sh preset_vesta_jan2026.oop

# Output: Tutti gli eventi Vesta con dettagli completi
```

### Ricerca Annuale (Background)

```bash
# Anno completo 2026
nohup ./italoccultcalc.sh preset_annual_2026.oop > annual_2026.log 2>&1 &

# Monitor
tail -f annual_2026.log

# Stop se necessario
pkill -f italoccultcalc
```

## 🔧 Configurazione Avanzata

### Preset Anatomy

```ini
# preset_custom.oop

general.
    .propagator = 'RK4'              # RA15, RK4, RK8
    .step_size_days = 0.05           # Precision vs speed

time.
    .start_date = '2026-01-01'
    .end_date = '2026-01-31'
    .interval_days = 1.0

asteroids.
    .selection_mode = 'auto'         # auto, list, file
    .asteroid_list = '4,1,2'         # Se mode=list
    .max_magnitude = 16.0            # H magnitude limit
    .min_diameter_km = 50.0
    .max_diameter_km = 1000.0

search.
    .max_magnitude = 16.0            # Stella limite
    .search_radius_deg = 0.1         # Shadow path width
    .min_probability = 0.01          # 1% minimum

occultation.
    # LinOccult compatibility
    .shadow_threshold_factor = 3.0
    .star_uncertainty_base_mas = 7.0
    .star_pm_degradation_mas_per_year = 2.5
    .use_probability_cdf = .TRUE.
    .earth_oblateness_factor = 0.996647
    
    # Performance
    .use_chebyshev = .TRUE.
    .chebyshev_order = 11
    .use_adaptive_timestep = .TRUE.
    .coarse_timestep_hours = 6.0
    .fine_timestep_hours = 0.5

gaia.
    .use_local_cache = .TRUE.
    .cache_directory = '/Users/USER/catalogs'
    .catalog_file = 'gaia_mag18.cat.gz'

output.
    .format = 'iota'
    .directory = 'output'
    .generate_kml = .TRUE.
    .generate_json = .TRUE.

filters.
    .use_geographic_filter = .TRUE.
    .latitude_min = 36.0             # Italia Sud
    .latitude_max = 47.0             # Italia Nord
    .longitude_min = 6.0
    .longitude_max = 19.0
    .max_distance_km = 2000.0
    .min_mag_drop = 0.5
    .min_duration_sec = 0.1
```

### Selezione Asteroidi

**Automatica** (default):
```ini
asteroids.
    .selection_mode = 'auto'
    .max_magnitude = 16.0
    .min_diameter_km = 50.0
    .max_diameter_km = 1000.0
```
Seleziona automaticamente asteroidi promettenti (H < 16, D = 50-1000 km).

**Lista Esplicita**:
```ini
asteroids.
    .selection_mode = 'list'
    .asteroid_list = '4,1,2,15,704'
```
Vesta, Ceres, Pallas, Eunomia, Interamnia.

**Da File**:
```ini
asteroids.
    .selection_mode = 'file'
    .asteroid_file = 'my_asteroids.txt'
```
File con un numero per riga.

## 📈 Output Formats

### IOTA (Standard)

```
(004) Vesta occults Gaia DR3 5813181197970338560  2026-01-15 03:24:17 UTC

Star: RA 10h 15m 23.45s  Dec +12° 34' 56.78"  Mag 13.2
C/A: 0.123"  PA: 67.8°  Vel: 21.3 km/s
Path width: 525.4 km  Duration: 24.7 sec  Prob: 87%
Uncertainty: ±156 km

Center Line:
  42°15'00"N  012°30'00"E  Rome, Italy
  43°00'00"N  011°52'30"E  Florence, Italy
  ...
```

### JSON (API/Web)

```json
{
  "event_id": "2026011503_004_5813181197970338560",
  "asteroid": {
    "number": 4,
    "name": "Vesta",
    "diameter_km": 525.4,
    "magnitude_H": 3.20
  },
  "star": {
    "gaia_dr3_id": "5813181197970338560",
    "ra_deg": 153.8477083,
    "dec_deg": 12.5824389,
    "magnitude_g": 13.2,
    "epoch": "J2016.0"
  },
  "event": {
    "time_utc": "2026-01-15T03:24:17.000Z",
    "closest_approach_arcsec": 0.123,
    "position_angle_deg": 67.8,
    "velocity_km_s": 21.3,
    "duration_sec": 24.7,
    "mag_drop": 5.8,
    "probability_percent": 87,
    "uncertainty_km": 156
  },
  "path": {
    "width_km": 525.4,
    "center_line": [
      {"lat": 42.25, "lon": 12.5, "location": "Rome, Italy"},
      {"lat": 43.0, "lon": 11.875, "location": "Florence, Italy"}
    ]
  }
}
```

### KML (Google Earth)

Visualizzazione path su mappa interattiva.

## 🐛 Troubleshooting

### "Gaia catalog not found"

```bash
# Download catalogo (9 GB)
./download_gaia_cache.sh

# Oppure manuale
curl -L https://github.com/manvalan/IOC_GaiaLib/releases/download/v2.0.0/gaia_mag18.cat.gz \
     -o ~/catalogs/gaia_mag18.cat.gz
```

### "MPC data not found"

```bash
./download_mpc_data.sh
```

### "SPICE kernel not found"

```bash
cd ~/.ioccultcalc/ephemerides
curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.bsp
```

### Performance Lenta

```bash
# Verifica OpenMP
./italoccultcalc.sh preset.oop 2>&1 | grep "parallel"
# Dovrebbe mostrare: "Calcolo parallelo su N core"

# Se assente, reinstalla
brew reinstall libomp
cd build && cmake .. && make -j4
```

### Out of Memory

```bash
# Riduci asteroidi o periodo
asteroids.
    .max_magnitude = 14.0  # Solo asteroidi grandi

time.
    .end_date = '2026-01-15'  # Periodo più breve
```

## 📚 Documentazione

- **Quick Start**: `QUICKSTART_PRODUCTION.md`
- **Performance**: `PERFORMANCE_OPTIMIZATION_STRATEGY.md`
- **Chebyshev**: `CHEBYSHEV_PERFORMANCE_ANALYSIS.md`
- **Gaia**: `GAIA_EDR3_QUICKREF.md`
- **Integration**: `ASTDYN_INTEGRATION_ANALYSIS.md`

## 🆘 Supporto

- **Issues**: https://github.com/manvalan/IOccultCalc/issues
- **Discussions**: https://github.com/manvalan/IOccultCalc/discussions
- **Email**: [da definire]

## 🎖️ Crediti

- **Steve Preston** - LinOccult algorithms
- **Minor Planet Center** - Orbital elements
- **Gaia Mission (ESA)** - Star catalog (1.8 billion stars)
- **JPL/NAIF** - SPICE ephemeris
- **IOTA** - Output format standards
- **OrbFit Team** - Orbit fitting

## 📝 License

MIT License - see LICENSE file

## 🌟 Citation

Se usi IOccultCalc per pubblicazioni scientifiche:

```bibtex
@software{ioccultcalc2025,
  author = {Bigi, Michele},
  title = {IOccultCalc: High-Performance Asteroid Occultation Prediction},
  year = {2025},
  url = {https://github.com/manvalan/IOccultCalc},
  version = {1.0.0}
}
```

---

**Version**: 1.0.0  
**Date**: November 2025  
**Status**: Production Ready ✅  
**Performance**: 3.4× faster than baseline | 7 minutes for 375 asteroids × 31 days
