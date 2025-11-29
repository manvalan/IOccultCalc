# Guida Integrazione Gaia EDR3/DR3

## Panoramica

ITALOccultCalc supporta sia **Gaia EDR3** (Early Data Release 3, Dicembre 2020) che **Gaia DR3** (Data Release 3, Giugno 2022) per le predizioni di occultazioni asteroidali.

## Differenze tra EDR3 e DR3

| Caratteristica | Gaia EDR3 | Gaia DR3 |
|---|---|---|
| **Data rilascio** | Dicembre 2020 | Giugno 2022 |
| **Epoca riferimento** | J2016.0 | J2016.0 |
| **N. sorgenti** | ~1.8 miliardi | ~1.8 miliardi |
| **Astrometria** | Migliorata | Ulteriormente migliorata |
| **Dati astrofisici** | Limitati | Estesi (RVS, XP spectra, ecc.) |
| **Tabella ADQL** | `gaiaedr3.gaia_source` | `gaiadr3.gaia_source` |

### Quale versione usare?

- **EDR3**: Sufficiente per la maggior parte delle predizioni di occultazioni. Astrometria di alta qualità, più stabile.
- **DR3**: Include miglioramenti astrometrici marginali e dati astrofisici estesi (non necessari per occultazioni).

**Raccomandazione**: Per occultazioni, **EDR3** è una scelta eccellente e stabile. DR3 può essere preferito per analisi più recenti o se si necessitano dati astrofisici aggiuntivi.

## Configurazione

### File .oop (OrbFit-style)

Aggiungi nella sezione `gaia.`:

```fortran
gaia.
        .version = 'EDR3'               ! Oppure 'DR3' (default se omesso)
        .use_local_cache = .TRUE.
        .cache_directory = '~/.ioccultcalc/gaia_edr3'  ! Directory separata consigliata
        .max_magnitude = 15.0
```

### File JSON

```json
{
  "gaia": {
    "version": "EDR3",
    "use_local_cache": true,
    "cache_directory": "~/.ioccultcalc/gaia_edr3",
    "max_magnitude": 15.0
  }
}
```

## Download Cache Locale

### Cache Separata per Versione

È consigliabile mantenere cache separate per EDR3 e DR3:

```bash
# Cache per DR3 (default)
~/.ioccultcalc/gaia/

# Cache per EDR3 (specificare in config)
~/.ioccultcalc/gaia_edr3/
```

### Download Dati EDR3

Il tool `gaia_cache_downloader` usa automaticamente la versione specificata nel `GaiaClient`.

**NOTA IMPORTANTE**: Attualmente il downloader non supporta ancora la selezione esplicita della versione da riga di comando. Per scaricare dati EDR3, è necessario:

1. Creare un preset con `gaia.version = 'EDR3'`
2. Eseguire `italoccultcalc` con il preset, che scaricherà automaticamente i dati EDR3 se `auto_download = .TRUE.`

**Esempio**:

```bash
# 1. Crea preset con EDR3
cat > preset_download_edr3.oop << EOF
gaia.
        .version = 'EDR3'
        .use_local_cache = .TRUE.
        .cache_directory = '~/.ioccultcalc/gaia_edr3'
        .auto_download = .TRUE.
        .max_magnitude = 15.0
EOF

# 2. Esegui italoccultcalc (scaricherà dati EDR3 al bisogno)
./italoccultcalc preset_gaia_edr3_example.oop
```

### Verifica Cache

```bash
# Verifica statistiche cache EDR3
ls -lh ~/.ioccultcalc/gaia_edr3/tiles/

# Output esempio:
# tile_0001.json
# tile_0002.json
# ...
```

## Query Online

Se non usi cache locale, le query al servizio TAP ESA Gaia useranno automaticamente la tabella corretta:

- **EDR3**: `gaiaedr3.gaia_source`
- **DR3**: `gaiadr3.gaia_source`

```fortran
gaia.
        .version = 'EDR3'               ! Usa gaiaedr3.gaia_source
        .use_local_cache = .FALSE.      ! Query online
```

## Esempi Completi

### Esempio 1: Ricerca con EDR3 (Cache Locale)

```fortran
! File: search_edr3.oop

time.
        .start_jd = 2460676.500000      ! 2025-01-01
        .end_jd = 2460706.500000        ! 2025-01-31

gaia.
        .version = 'EDR3'
        .use_local_cache = .TRUE.
        .cache_directory = '~/.ioccultcalc/gaia_edr3'
        .auto_download = .TRUE.
        .max_magnitude = 15.0

asteroid_db.
        .source = 'MPCORB'
        .filter_by = 'H_magnitude'
        .max_H = 9.0
```

Esecuzione:

```bash
./italoccultcalc search_edr3.oop
```

### Esempio 2: Query Online DR3 (No Cache)

```fortran
! File: search_dr3_online.oop

time.
        .start_jd = 2460676.500000
        .end_jd = 2460706.500000

gaia.
        .version = 'DR3'                ! Usa DR3 (default)
        .use_local_cache = .FALSE.      ! Query online
        .max_magnitude = 16.0

asteroid_db.
        .source = 'MPCORB'
        .filter_by = 'H_magnitude'
        .max_H = 10.0
```

### Esempio 3: Confronto EDR3 vs DR3

Puoi eseguire predizioni con entrambe le versioni e confrontare:

```bash
# EDR3
./italoccultcalc preset_test_edr3.oop > results_edr3.txt

# DR3
./italoccultcalc preset_test_dr3.oop > results_dr3.txt

# Confronta
diff results_edr3.txt results_dr3.txt
```

## Troubleshooting

### Errore: "Gaia version not recognized"

Assicurati di usare esattamente `'EDR3'` o `'DR3'` (case-insensitive):

```fortran
gaia.
        .version = 'EDR3'   ! ✓ Corretto
        .version = 'edr3'   ! ✓ Corretto (case-insensitive)
        .version = 'GAIADR3' ! ✗ ERRATO
```

### Cache Mista EDR3/DR3

Se hai cache mista, pulisci e riscarica:

```bash
# Pulisci cache esistente
rm -rf ~/.ioccultcalc/gaia_edr3/tiles/*

# Riscarica con preset corretto
./italoccultcalc preset_gaia_edr3_example.oop
```

### Differenze nei Risultati

Le differenze tra EDR3 e DR3 nelle predizioni sono tipicamente **< 1 mas** in posizione stellare, corrispondenti a **< 50 km** sulla Terra per la maggior parte degli asteroidi. Queste differenze sono trascurabili rispetto alle incertezze orbitali asteroidali.

## API Programmatica (C++)

### Impostare Versione Gaia

```cpp
#include "ioccultcalc/gaia_client.h"

// Crea client e imposta versione
GaiaClient client;
client.setGaiaVersion(GaiaVersion::EDR3);  // O GaiaVersion::DR3

// Query stella
auto stars = client.queryCone(ra, dec, radius, 15.0);
```

### Con GaiaCache

```cpp
#include "ioccultcalc/gaia_cache.h"

// Crea cache e imposta versione
GaiaCache cache("~/.ioccultcalc/gaia_edr3");
cache.setGaiaVersion(GaiaVersion::EDR3);

// Query regione
auto stars = cache.queryRegion(ra, dec, radius, 15.0, true);
```

## Risorse Esterne

- **Gaia EDR3**: https://www.cosmos.esa.int/web/gaia/earlydr3
- **Gaia DR3**: https://www.cosmos.esa.int/web/gaia/dr3
- **TAP Service**: https://gea.esac.esa.int/tap-server/tap
- **ADQL Documentation**: http://www.ivoa.net/documents/ADQL/

## Note Tecniche

### Epoca di Riferimento

Entrambe EDR3 e DR3 usano **J2016.0** come epoca di riferimento. ITALOccultCalc propaga automaticamente le posizioni stellari all'epoca di osservazione usando proper motion e parallasse.

### Precisione Astrometrica

| Magnitudine G | Precisione EDR3 | Precisione DR3 |
|---|---|---|
| G < 15 | ~0.02-0.1 mas | ~0.02-0.08 mas |
| 15 < G < 17 | ~0.1-0.5 mas | ~0.08-0.4 mas |
| G > 17 | ~0.5-1.0 mas | ~0.4-0.8 mas |

Per occultazioni asteroidali, le incertezze orbitali (~100-1000 km) dominano rispetto alle incertezze astrometriche stellari (~1-10 km).

## Supporto

Per problemi o domande:
- GitHub Issues: https://github.com/manvalan/IOccultCalc/issues
- Email: michele.bigi@example.com
