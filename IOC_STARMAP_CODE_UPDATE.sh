#!/bin/bash
# Script per aggiornare ItalOccultFormatter con IOC_StarMap

echo "🔧 IOC_StarMap Integration - Code Update Script"
echo "================================================"
echo ""

# Colori per output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 1. Verifica esistenza libreria IOC_StarMap
echo -e "${BLUE}[1/5]${NC} Verifica libreria IOC_StarMap..."
if [ ! -d "external/IOC_StarMap" ]; then
    echo -e "${YELLOW}⚠${NC}  IOC_StarMap non trovata in external/"
    echo "   Clona il repository:"
    echo "   git submodule add https://github.com/manvalan/IOC_StarMap external/IOC_StarMap"
    echo ""
else
    echo -e "${GREEN}✓${NC}  IOC_StarMap trovata"
fi

# 2. Backup file esistenti
echo -e "${BLUE}[2/5]${NC} Backup file esistenti..."
if [ -f "include/output_formatter.h" ]; then
    cp include/output_formatter.h include/output_formatter.h.bak
    echo -e "${GREEN}✓${NC}  Backup: include/output_formatter.h.bak"
fi

if [ -f "src/output_formatter.cpp" ]; then
    cp src/output_formatter.cpp src/output_formatter.cpp.bak
    echo -e "${GREEN}✓${NC}  Backup: src/output_formatter.cpp.bak"
fi

# 3. Aggiorna CMakeLists.txt
echo -e "${BLUE}[3/5]${NC} Aggiorna CMakeLists.txt..."

if grep -q "IOC_StarMap" CMakeLists.txt; then
    echo -e "${GREEN}✓${NC}  IOC_StarMap già presente in CMakeLists.txt"
else
    echo ""
    echo "Aggiungi al CMakeLists.txt:"
    echo ""
    echo "# IOC_StarMap library"
    echo "add_subdirectory(external/IOC_StarMap)"
    echo ""
    echo "target_link_libraries(ioccultcalc"
    echo "    PRIVATE"
    echo "        IOC_StarMap::Core"
    echo "        # ... altre librerie ..."
    echo ")"
    echo ""
fi

# 4. Mostra modifiche necessarie al codice
echo -e "${BLUE}[4/5]${NC} Modifiche codice necessarie:"
echo ""
echo "📝 include/output_formatter.h"
echo "   Aggiungi:"
echo "   #include \"ioc_starmap/generator.h\""
echo ""
echo "   Nella classe ItalOccultFormatter:"
echo "   private:"
echo "       ioc::starmap::StarMapGenerator starmap_generator_;"
echo "       std::string generateApproachChart(const OWOccultation& occ);"
echo "       std::string generateFinalChart(const OWOccultation& occ);"
echo "       std::vector<ioc::starmap::AsteroidTrackPoint> calculateAsteroidTrack("
echo "           const OWOccultation& occ, int days_before = 10);"
echo ""
echo "📝 src/output_formatter.cpp"
echo "   Sostituisci i placeholder URL con chiamate IOC_StarMap:"
echo ""
echo "   OLD:"
echo "   std::string url = \"https://starmap.example.com/chart?...\";"
echo ""
echo "   NEW:"
echo "   std::string svg = generateApproachChart(occ);"
echo "   sheet << \"        \" << svg << \"\\n\";"
echo ""

# 5. Test checklist
echo -e "${BLUE}[5/5]${NC} Test Checklist:"
echo ""
echo "   □ Compilazione: ./build.sh"
echo "   □ Test unit: ./build/tests/test_output_formatter"
echo "   □ Genera esempio: ./build/ioccultcalc --preset preset_italoccult.oop"
echo "   □ Verifica SVG: grep -c '<svg' example_italoccult.html"
echo "   □ Converti PDF: python3 convert_to_pdf.py example_italoccult.html"
echo "   □ Verifica PDF: open example_italoccult.pdf"
echo ""

# Riepilogo file documentazione
echo "📚 Documentazione:"
echo "   • IOC_STARMAP_INTEGRATION.md    - Spec completa API"
echo "   • IOC_STARMAP_QUICKSTART.md     - Guida rapida"
echo "   • IOC_STARMAP_CODE_UPDATE.sh    - Questo script"
echo ""

echo -e "${GREEN}✓${NC} Preparazione completata!"
echo ""
echo "Prossimi step:"
echo "1. Clona IOC_StarMap se non presente"
echo "2. Aggiorna CMakeLists.txt"
echo "3. Modifica output_formatter.h e output_formatter.cpp"
echo "4. Compila e testa"
echo ""
echo "Per dettagli implementazione, vedi: IOC_STARMAP_INTEGRATION.md"
