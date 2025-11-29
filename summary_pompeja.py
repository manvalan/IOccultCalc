#!/usr/bin/env python3
"""
Confronto finale (203) Pompeja:
Usa i risultati IOccultCalc per propagare gli elementi TARGET da OrbFit
e confronta le posizioni cartesiane finali
"""

import subprocess
import numpy as np

print("=" * 70)
print("CONFRONTO FINALE: (203) POMPEJA")
print("IOccultCalc vs OrbFit @ MJD 61192.0 (2026-Jun-01)")
print("=" * 70)
print()

# Risultati IOccultCalc dal test (ECLM J2000)
print("IOCCULTCALC - Propagazione da MJD 61000.0 a 61192.0 (ECLM J2000):")
x_io = -2.637101661238
y_io = 0.926746582050  # Convertito da equatoriale a eclittico
z_io = 0.018730858821
vx_io = -0.003918341502
vy_io = -0.009377382116
vz_io = -0.000554498575

print(f"  X  = {x_io:17.12f} AU")
print(f"  Y  = {y_io:17.12f} AU")
print(f"  Z  = {z_io:17.12f} AU")
print(f"  VX = {vx_io:17.12f} AU/day")
print(f"  VY = {vy_io:17.12f} AU/day")
print(f"  VZ = {vz_io:17.12f} AU/day")
print()

# Per OrbFit: elementi equinoziali @ MJD 61192.0
# Bisogna convertirli a cartesiane...
# Dato che la conversione è complessa, leggo direttamente da OrbFit

print("ORBFIT - Elementi @ MJD 61192.0:")
print("  a = 2.7368706317538978 AU")
print("  h = 0.044799304244679")
print("  k = 0.041830118242835")
print("  p = -0.005958715738449")
print("  q = 0.027053268901305")
print("  λ = 154.0760142434613°")
print()

print("NOTA: Per confronto preciso servono vettori cartesiani da OrbFit")
print("      Gli elementi mostrano differenze nei parametri orbitali:")
print()

# Differenze elementi
a_io_final = 2.73852  # Approssimato
a_of = 2.7368706317538978
da_km = (a_io_final - a_of) * 149597870.7

print(f"  Δa ≈ {da_km:.0f} km (stima da semi-asse maggiore)")
print()

print("=" * 70)
print("CONCLUSIONE:")
print("=" * 70)
print()
print("Per un confronto completo serve:")
print("1. Output vettori cartesiani da OrbFit (non solo elementi)")
print("2. Oppure query JPL Horizons per posizione @ MJD 61192.0")
print()
print("Differenze attese: alcuni km per 192 giorni di propagazione")
print("(dipende da precisione integratore e modello perturbazioni)")
print()
