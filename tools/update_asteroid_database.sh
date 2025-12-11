#!/bin/bash
# Update Asteroid Database Script
# Runs monthly via cron or manually

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/build_asteroid_database.py"

# Default: update if older than 30 days
MAX_AGE_DAYS=${1:-30}

echo "=========================================="
echo "Asteroid Database Update Script"
echo "=========================================="
echo ""

# Run update script
python3 "${PYTHON_SCRIPT}" --max-age "${MAX_AGE_DAYS}"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✓ Update completed successfully"
else
    echo ""
    echo "✗ Update failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE

