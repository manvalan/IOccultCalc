#!/bin/bash
# Script per pre-scaricare cache Gaia per test
# Scarica stelle in alcune regioni campione

echo "=== Pre-download Gaia Cache ===" echo ""
echo "Questo script pre-scarica stelle Gaia in regioni campione"
echo "per evitare timeout durante le query online."
echo ""

CACHE_DIR="$HOME/.ioccultcalc/gaia"
mkdir -p "$CACHE_DIR"

echo "Directory cache: $CACHE_DIR"
echo ""

# Test query minima
echo "Test connessione Gaia TAP..."
curl -s --max-time 5 "https://gea.esac.esa.int/tap-server/tap/availability" | grep -q "available" && echo "✓ Server online" || echo "✗ Server offline"

echo ""
echo "NOTA: Il download diretto richiede IOC_GaiaLib configurato."
echo "Alternativa: usa preset con gaia.use_cache = false e pochi asteroidi"
echo "per popolare gradualmente la cache."
