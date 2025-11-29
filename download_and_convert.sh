#!/bin/bash
# Script per workflow completo: download da AstDyS + conversione preset
# Uso: ./download_and_convert.sh asteroids_list.txt [output_dir]

set -e

# Colori
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parametri
INPUT_FILE="${1}"
OUTPUT_DIR="${2:-astdys_data}"
PRESET_DIR="${3:-presets}"

# Banner
echo -e "${BLUE}╔═══════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   Workflow AstDyS: Download + Conversione Preset         ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════════════════════════╝${NC}"
echo

# Verifica parametri
if [ -z "$INPUT_FILE" ]; then
    echo -e "${RED}Errore: specificare file input${NC}"
    echo "Uso: $0 asteroids_list.txt [output_dir] [preset_dir]"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo -e "${RED}Errore: file $INPUT_FILE non trovato${NC}"
    exit 1
fi

# Step 1: Download da AstDyS
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}STEP 1: Download file da AstDyS${NC}"
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo

python3 download_astdys_data.py "$INPUT_FILE" --output-dir "$OUTPUT_DIR"

if [ $? -ne 0 ]; then
    echo -e "${RED}✗ Errore durante download${NC}"
    exit 1
fi

echo
echo -e "${GREEN}✓ Download completato${NC}"
echo

# Step 2: Conversione preset
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}STEP 2: Conversione file .eq1 in preset .oop${NC}"
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo

# Crea directory preset
mkdir -p "$PRESET_DIR"

# Contatori
total=0
success=0
failed=0

# Converti ogni file .eq1
for eq1_file in "$OUTPUT_DIR"/*.eq1; do
    if [ ! -f "$eq1_file" ]; then
        continue
    fi
    
    total=$((total + 1))
    asteroid_num=$(basename "$eq1_file" .eq1)
    preset_file="$PRESET_DIR/preset_${asteroid_num}.oop"
    
    echo -e "${BLUE}[${total}] Conversione ($asteroid_num)...${NC}"
    
    if python3 eq1_to_preset.py "$eq1_file" --output "$preset_file" 2>&1 | grep -q "✓ Preset scritto"; then
        success=$((success + 1))
        echo -e "  ${GREEN}✓ $preset_file${NC}"
    else
        failed=$((failed + 1))
        echo -e "  ${RED}✗ Errore conversione${NC}"
    fi
done

echo
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${YELLOW}RIEPILOGO FINALE${NC}"
echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo
echo "File .eq1 processati: $total"
echo -e "  ${GREEN}✓ Successo: $success${NC}"
echo -e "  ${RED}✗ Falliti:  $failed${NC}"
echo
echo -e "📁 Dati AstDyS:  ${BLUE}$OUTPUT_DIR${NC}"
echo -e "📁 Preset .oop:  ${BLUE}$PRESET_DIR${NC}"
echo

# Step 3: Genera lista preset
preset_list="$PRESET_DIR/all_presets.txt"
ls -1 "$PRESET_DIR"/preset_*.oop > "$preset_list" 2>/dev/null || true

if [ -s "$preset_list" ]; then
    echo -e "${GREEN}✓ Lista preset: $preset_list${NC}"
    echo
    echo -e "${YELLOW}💡 Prossimi step:${NC}"
    echo "   # Esegui calcolo per un singolo asteroide:"
    echo "   ./italoccultcalc $PRESET_DIR/preset_433.oop"
    echo
    echo "   # Esegui batch per tutti:"
    echo "   for preset in $PRESET_DIR/preset_*.oop; do"
    echo "     echo \"Processing \$preset...\""
    echo "     ./italoccultcalc \"\$preset\" > \"results_\$(basename \$preset .oop).txt\""
    echo "   done"
else
    echo -e "${RED}⚠ Nessun preset generato${NC}"
fi

echo
echo -e "${GREEN}╔═══════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║             WORKFLOW COMPLETATO CON SUCCESSO             ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════════════════════════╝${NC}"
