#!/bin/bash
#
# Script per aggiornare database allnum.cat
# Verifica se il database è vecchio (> 30 giorni) e lo aggiorna
#
# Usage:
#   ./update_allnum_database.sh [max_age_days]
#
# Default: aggiorna se database è più vecchio di 30 giorni

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"

MAX_AGE_DAYS=${1:-30}

echo "🔄 Allnum Database Updater"
echo "   Max age: $MAX_AGE_DAYS days"
echo ""

# Check if build_allnum_database exists
if [ ! -f "$BUILD_DIR/tools/build_allnum_database" ]; then
    echo "❌ Error: build_allnum_database not found"
    echo "   Please build the project first:"
    echo "   cd $PROJECT_ROOT && mkdir -p build && cd build && cmake .. && make build_allnum_database"
    exit 1
fi

# Run the updater
"$BUILD_DIR/tools/build_allnum_database" --max-age "$MAX_AGE_DAYS"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✅ Allnum database update completed"
else
    echo ""
    echo "❌ Allnum database update failed (exit code: $EXIT_CODE)"
    exit $EXIT_CODE
fi

