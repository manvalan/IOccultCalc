#!/bin/bash
# Download MPC Orbital Elements (Real Data)
# Source: Minor Planet Center Extended Files

set -e

echo "╔══════════════════════════════════════════════════════════╗"
echo "║     MPC Orbital Elements Download                        ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""

DATA_DIR="${HOME}/.ioccultcalc/data"
MPC_URL="https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz"
TEMP_FILE="${DATA_DIR}/mpcorb_extended.json.gz"
OUTPUT_FILE="${DATA_DIR}/all_numbered_asteroids.json"

# Create directory
mkdir -p "${DATA_DIR}"

# Download
echo "📥 Download MPC orbital elements..."
echo "   Source: ${MPC_URL}"
echo "   Size: ~2 GB compressed"
echo ""

if [ -f "${OUTPUT_FILE}" ]; then
    read -p "⚠️  File esistente trovato. Sovrascrivere? (y/N) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "✓ Mantenuto file esistente"
        exit 0
    fi
fi

echo "⏳ Download in corso (può richiedere diversi minuti)..."
if ! curl -L --progress-bar "${MPC_URL}" -o "${TEMP_FILE}"; then
    echo "✗ Errore download"
    exit 1
fi

echo "✓ Download completato"
echo ""

# Decompress and convert
echo "📦 Decompressione..."
gunzip -c "${TEMP_FILE}" > "${OUTPUT_FILE}"
rm "${TEMP_FILE}"

# Stats
SIZE=$(du -h "${OUTPUT_FILE}" | cut -f1)
LINES=$(wc -l < "${OUTPUT_FILE}")

echo "✓ Decompressione completata"
echo ""
echo "📊 Statistiche:"
echo "   File: ${OUTPUT_FILE}"
echo "   Size: ${SIZE}"
echo "   Records: ${LINES}"
echo ""

# Preview
echo "📄 Preview primi 5 asteroidi:"
if command -v jq &> /dev/null; then
    head -20 "${OUTPUT_FILE}" | jq -c 'select(.Number != null) | {Number, Name, H, Epoch}' 2>/dev/null | head -5
else
    head -10 "${OUTPUT_FILE}"
fi

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║           ✓ DATI MPC PRONTI                             ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "Elementi orbitali per $(grep -c '"Number"' "${OUTPUT_FILE}") asteroidi numerati"
echo ""
