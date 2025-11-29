#!/usr/bin/env python3
"""
Scarica elementi orbitali di Ceres da JPL Horizons
con epoca 29 Nov 2025 (vicina alla data target)
"""

import requests

# Query Horizons per elementi osculanti
url = "https://ssd.jpl.nasa.gov/api/horizons.api"

params = {
    'format': 'text',
    'COMMAND': '1',  # Ceres
    'OBJ_DATA': 'YES',
    'MAKE_EPHEM': 'YES',
    'EPHEM_TYPE': 'ELEMENTS',
    'CENTER': '500@10',  # Sole
    'START_TIME': '2025-Nov-29',
    'STOP_TIME': '2025-Nov-30',
    'STEP_SIZE': '1d',
    'REF_SYSTEM': 'ICRF',
    'REF_PLANE': 'ECLIPTIC',
    'OUT_UNITS': 'AU-D',
    'CSV_FORMAT': 'NO'
}

print("Query JPL Horizons per elementi Ceres (epoca 29 Nov 2025)...")
print("URL:", url)
print()

response = requests.get(url, params=params)

if response.status_code == 200:
    content = response.text
    
    # Estrai sezione elementi
    in_elements = False
    for line in content.split('\n'):
        if '$$SOE' in line:
            in_elements = True
            continue
        if '$$EOE' in line:
            in_elements = False
            
        if in_elements or 'EC=' in line or 'QR=' in line or 'IN=' in line:
            print(line)
    
    # Salva risposta completa
    with open('/tmp/ceres_elements_nov2025.txt', 'w') as f:
        f.write(content)
    print("\n\nRisposta completa salvata in /tmp/ceres_elements_nov2025.txt")
else:
    print(f"Errore HTTP {response.status_code}")
    print(response.text)
