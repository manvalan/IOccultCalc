# Gaia EDR3/DR3 - Quick Reference

## Configurazione Rapida

### File .oop (OrbFit-style)

```fortran
gaia.
    .version = 'EDR3'               ! Oppure 'DR3' (default)
    .use_local_cache = .TRUE.
    .cache_directory = '~/.ioccultcalc/gaia_edr3'
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

## Comandi Utili

### Test Integrazione

```bash
# Test query EDR3 vs DR3
cd /path/to/IOccultCalc/build/examples
./test_gaia_edr3
```

### Esempio Completo

```bash
# Usa preset EDR3 incluso
cd /path/to/IOccultCalc
./italoccultcalc preset_gaia_edr3_example.oop
```

### Cache Locale

```bash
# Directory cache EDR3
~/.ioccultcalc/gaia_edr3/

# Directory cache DR3 (default)
~/.ioccultcalc/gaia/
```

## Tabelle ADQL

| Versione | Tabella ADQL |
|---|---|
| EDR3 | `gaiaedr3.gaia_source` |
| DR3 | `gaiadr3.gaia_source` |

## API C++

```cpp
#include "ioccultcalc/gaia_client.h"

// Seleziona versione
GaiaClient client;
client.setGaiaVersion(GaiaVersion::EDR3);  // O GaiaVersion::DR3

// Query
auto stars = client.queryCone(ra, dec, radius, mag_limit);
```

## FAQ

**Q: Quale versione usare?**  
A: EDR3 è perfetto per occultazioni. DR3 ha miglioramenti marginali.

**Q: Posso avere entrambe le cache?**  
A: Sì! Usa directory separate (`gaia_edr3` vs `gaia`).

**Q: Come cambio versione?**  
A: Modifica `gaia.version` nel file config e riavvia.

## Documentazione Completa

- **Guida Dettagliata**: `docs/GAIA_EDR3_INTEGRATION.md`
- **Riepilogo Implementazione**: `GAIA_EDR3_IMPLEMENTATION_SUMMARY.md`
- **Esempio Preset**: `preset_gaia_edr3_example.oop`
