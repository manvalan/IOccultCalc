#!/usr/bin/env python3
"""
Test: Confronto diretto coordinate Horizons vs valori letti
Obiettivo: Capire se c'è trasformazione di coordinate mancante
"""

import requests
from datetime import datetime

# Query Horizons API
def query_horizons(obj_id, jd):
    """Query JPL Horizons per stato vettoriale"""
    base_url = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    params = {
        'format': 'text',
        'COMMAND': f"'{obj_id}'",
        'EPHEM_TYPE': 'VECTORS',
        'CENTER': "'@sun'",
        'START_TIME': f"'JD{jd}'",
        'STOP_TIME': f"'JD{jd+0.1}'",
        'STEP_SIZE': "'1 d'",
        'OUT_UNITS': "'AU-D'",
        'REF_PLANE': "'ECLIPTIC'",
        'REF_SYSTEM': "'ICRF'",
        'VEC_CORR': "'NONE'",
        'CSV_FORMAT': "'YES'"
    }
    
    response = requests.get(base_url, params=params)
    text = response.text
    
    # Parse risposta
    lines = text.split('\n')
    in_data = False
    
    for line in lines:
        if '$$SOE' in line:
            in_data = True
            continue
        if '$$EOE' in line:
            break
        if in_data and ',' in line:
            parts = line.split(',')
            if len(parts) >= 7:
                try:
                    x = float(parts[2].strip())
                    y = float(parts[3].strip())
                    z = float(parts[4].strip())
                    vx = float(parts[5].strip())
                    vy = float(parts[6].strip())
                    vz = float(parts[7].strip())
                    return x, y, z, vx, vy, vz
                except:
                    pass
    
    return None

# Test
print("=" * 70)
print("TEST: Coordinate Horizons per 433 Eros")
print("=" * 70)
print()

jd = 2460615.5  # 2024-Jan-01
print(f"Epoca: JD {jd} (2024-Jan-01 00:00 UTC)")
print()

print("Query JPL Horizons...")
result = query_horizons('433', jd)

if result:
    x, y, z, vx, vy, vz = result
    print(f"✓ Dati ricevuti\n")
    print(f"REF_PLANE: ECLIPTIC")
    print(f"REF_SYSTEM: ICRF")
    print(f"CENTER: @sun (heliocentric)\n")
    print(f"Posizione (AU):")
    print(f"  X = {x:15.10f}")
    print(f"  Y = {y:15.10f}")
    print(f"  Z = {z:15.10f}")
    print(f"  |r| = {(x**2 + y**2 + z**2)**0.5:15.10f}")
    print()
    print(f"Velocità (AU/day):")
    print(f"  VX = {vx:15.10f}")
    print(f"  VY = {vy:15.10f}")
    print(f"  VZ = {vz:15.10f}")
    print()
    
    # Coordinate attese da IOccultCalc (dai test precedenti)
    x_ioc = -0.70265411
    y_ioc = -1.36210444
    z_ioc = -0.25776641
    
    print("=" * 70)
    print("CONFRONTO con valori IOccultCalc")
    print("=" * 70)
    print()
    print(f"Valori IOccultCalc (dai test):")
    print(f"  X = {x_ioc:15.10f}")
    print(f"  Y = {y_ioc:15.10f}")
    print(f"  Z = {z_ioc:15.10f}")
    print()
    
    dx = (x - x_ioc) * 149597870.7  # km
    dy = (y - y_ioc) * 149597870.7
    dz = (z - z_ioc) * 149597870.7
    
    print(f"Differenza:")
    print(f"  ΔX = {dx:10.3f} km")
    print(f"  ΔY = {dy:10.3f} km")
    print(f"  ΔZ = {dz:10.3f} km")
    print(f"  |Δr| = {(dx**2 + dy**2 + dz**2)**0.5:10.3f} km")
    print()
    
    if abs(dz) > 100:
        print(f"⚠️  ATTENZIONE: Differenza significativa su Z ({dz:.1f} km)")
        print(f"   Questo potrebbe indicare problema di frame di riferimento")
    elif (dx**2 + dy**2 + dz**2)**0.5 < 10:
        print(f"✓ Valori concordano (<10 km) → stato iniziale corretto")
else:
    print("✗ Errore nel download da Horizons")
