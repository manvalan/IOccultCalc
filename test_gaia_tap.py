#!/usr/bin/env python3
"""
Test connessione Gaia TAP con query semplice
Per verificare se il problema è nel server o nel client
"""

import requests
import time

TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"

# Query ADQL minima: 10 stelle più vicine al Sole
QUERY = """
SELECT TOP 10 
    source_id, ra, dec, parallax, phot_g_mean_mag
FROM gaiadr3.gaia_source
WHERE parallax > 100
ORDER BY parallax DESC
"""

print("Test query Gaia TAP...")
print(f"Server: {TAP_URL}")
print(f"Query: {QUERY.strip()}")
print("")

start = time.time()

try:
    response = requests.post(
        TAP_URL,
        data={
            'REQUEST': 'doQuery',
            'LANG': 'ADQL',
            'FORMAT': 'json',
            'QUERY': QUERY
        },
        timeout=30
    )
    
    elapsed = time.time() - start
    
    if response.status_code == 200:
        print(f"✓ Query OK in {elapsed:.1f} secondi")
        data = response.json()
        print(f"  Stelle ricevute: {len(data.get('data', []))}")
    else:
        print(f"✗ Errore HTTP {response.status_code}")
        print(f"  Response: {response.text[:200]}")
        
except requests.Timeout:
    print(f"✗ Timeout dopo 30 secondi")
except Exception as e:
    print(f"✗ Errore: {e}")
