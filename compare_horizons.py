#!/usr/bin/env python3
"""
Confronta posizione 17030 da IOccultCalc con JPL Horizons
"""

print("=== POSIZIONE 17030 DA IOccultCalc ===")
print("Data: 2025-11-28 00:00 UTC (JD 2460642.5)")
print("  RA  = 356.3218° = 23h 45m 17.2s")
print("  Dec = -0.0046° = -0° 00' 17\"")
print()

print("=== QUERY JPL HORIZONS ===")
print("Esegui manualmente su https://ssd.jpl.nasa.gov/horizons.cgi")
print("  Target: 17030 Sierks")
print("  Observer: Geocentric (500@0)")
print("  Date: 2025-11-28 00:00")
print("  Quantities: 1,19 (RA/Dec, Range)")
print()

print("=== STELLA TARGET ===")
print("UCAC4 552-011427 (Gaia DR3 3411546266140512128)")
print("  RA  = 73.4161° = 4h 53m 39.9s")
print("  Dec = +20.3317° = +20° 19' 54\"")
print()

print("=== DISTANZA ANGOLARE ===")
print("  Δ ≈ 283° (opposti nel cielo!)")
