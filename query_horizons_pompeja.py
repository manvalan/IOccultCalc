#!/usr/bin/env python3
"""
Query JPL Horizons per (203) Pompeja
Epoca: MJD 61192.0 (2027-Aug-24 00:00)
"""

import requests
from datetime import datetime

# Parametri query Horizons
params = {
    'format': 'text',
    'COMMAND': '203',  # (203) Pompeja
    'OBJ_DATA': 'NO',
    'MAKE_EPHEM': 'YES',
    'EPHEM_TYPE': 'VECTORS',
    'CENTER': '@0',  # Solar System Barycenter
    'START_TIME': '2027-Aug-24 00:00',
    'STOP_TIME': '2027-Aug-24 00:01',
    'STEP_SIZE': '1 d',
    'VEC_TABLE': '2',  # Position and velocity
    'REF_PLANE': 'ECLIPTIC',
    'REF_SYSTEM': 'J2000',
    'VEC_CORR': 'NONE',
    'OUT_UNITS': 'AU-D',
    'CSV_FORMAT': 'NO',
}

print("Querying JPL Horizons for (203) Pompeja...")
print("Epoch: 2027-Aug-24 00:00 (MJD 61192.0)")
print()

url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
response = requests.get(url, params=params)

if response.status_code == 200:
    result = response.text
    
    # Estrai i vettori di stato
    lines = result.split('\n')
    in_data = False
    
    for i, line in enumerate(lines):
        if '$$SOE' in line:
            in_data = True
            continue
        if '$$EOE' in line:
            in_data = False
            break
        if in_data and line.strip():
            print(line)
    
    # Cerca e stampa i valori numerici
    print("\n" + "="*70)
    print("RISULTATI HORIZONS:")
    print("="*70)
    
    for i, line in enumerate(lines):
        if '$$SOE' in line:
            # I dati sono nelle righe successive
            data_lines = lines[i+1:i+10]
            for dl in data_lines:
                if 'X =' in dl or 'VX=' in dl:
                    print(dl.strip())
                    
else:
    print(f"Error: {response.status_code}")
    print(response.text[:500])
