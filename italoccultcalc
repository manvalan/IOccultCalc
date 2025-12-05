#!/bin/bash
# IOccultCalc Launcher Script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXEC="${SCRIPT_DIR}/build/examples/italoccultcalc"

if [ ! -f "$EXEC" ]; then
    echo "✗ Errore: italoccultcalc non trovato"
    echo "  Esegui: ./install_production.sh"
    exit 1
fi

# Check catalogo
if [ ! -f "${HOME}/catalogs/gaia_mag18.cat.gz" ]; then
    echo "⚠ WARNING: Catalogo Gaia Mag18 non trovato"
    echo "  Download: ./download_gaia_cache.sh"
fi

# Esegui
"$EXEC" "$@"
