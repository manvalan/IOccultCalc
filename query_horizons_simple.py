#!/usr/bin/env python3
"""
JPL Horizons Query Script SEMPLIFICATO per IOccultCalc
======================================================

Versione semplificata che usa solo i parametri essenziali.
"""

import sys
import requests
import json
from datetime import datetime, timedelta

def query_horizons_simple(target_id: str, mjd: float) -> tuple:
    """Query semplificata JPL Horizons"""
    
    try:
        # Converti MJD in JD
        jd = mjd + 2400000.5
        
        # Converti JD in data (approximazione semplice)
        # JD 2451545.0 = 2000-01-01 12:00:00 TT
        days_since_j2000 = jd - 2451545.0
        j2000_date = datetime(2000, 1, 1, 12, 0, 0)
        target_date = j2000_date + timedelta(days=days_since_j2000)
        
        print(f"🔍 JPL Horizons Query Semplificata:")
        print(f"   Target: {target_id}")
        print(f"   MJD: {mjd:.6f}")
        print(f"   JD: {jd:.6f}")
        print(f"   Data: {target_date.strftime('%Y-%m-%d %H:%M')}")
        
        # URL diretto più semplice
        base_url = "https://ssd.jpl.nasa.gov/api/horizons.api"
        
        # Calcola data finale (1 giorno dopo per intervallo minimo)
        end_date = target_date + timedelta(days=1)
        
        # Parametri minimi - solo date senza orari
        params = {
            'format': 'json',
            'COMMAND': target_id,
            'OBJ_DATA': 'NO',
            'MAKE_EPHEM': 'YES',
            'EPHEM_TYPE': 'OBSERVER',
            'CENTER': '500',
            'START_TIME': target_date.strftime('%Y-%m-%d'),
            'STOP_TIME': end_date.strftime('%Y-%m-%d'),
            'STEP_SIZE': '1h',
            'QUANTITIES': '1'
        }
        
        print(f"📡 Query in corso...")
        response = requests.get(base_url, params=params, timeout=30)
        
        if response.status_code != 200:
            return None, None, f"HTTP_{response.status_code}"
        
        data = response.json()
        
        if 'result' not in data:
            return None, None, "NO_RESULT_FIELD"
        
        result_text = data['result']
        print(f"📄 Risposta: {len(result_text)} caratteri")
        
        # Debug: mostra parte della risposta
        lines = result_text.split('\n')
        print(f"📄 Prime 5 linee:")
        for i, line in enumerate(lines[:5]):
            print(f"      {line}")
        
        # Cerca sezione $$SOE ... $$EOE
        in_ephemeris = False
        for line in lines:
            if '$$SOE' in line:
                in_ephemeris = True
                continue
            elif '$$EOE' in line:
                break
            elif in_ephemeris and line.strip():
                print(f"📊 Linea effemeridi: {line[:100]}")
                
                # Formato JPL Horizons tipico:
                # 2023-Mar-03 00:00     19 21 08.56 -28 34 37.1
                # Data                  RA (H M S)  Dec (D M S)
                
                # Cerca pattern con spazi multipli (campo separatore tipico JPL)
                import re
                
                # Pattern per catturare RA e Dec in formato sessagesimale
                # RA: HH MM SS.SS   Dec: ±DD MM SS.S
                pattern = r'(\d{4}-\w{3}-\d{2}\s+\d{2}:\d{2})\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s*([-+]?\d+)\s+(\d+)\s+(\d+\.\d+)'
                match = re.search(pattern, line)
                
                if match:
                    ra_h = int(match.group(2))
                    ra_m = int(match.group(3))  
                    ra_s = float(match.group(4))
                    dec_d = int(match.group(5))
                    dec_m = int(match.group(6))
                    dec_s = float(match.group(7))
                    
                    print(f"📍 RA: {ra_h}h {ra_m}m {ra_s}s")
                    print(f"📍 Dec: {dec_d}° {dec_m}' {dec_s}\"")
                    
                    # Converti in gradi decimali
                    ra_deg = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0  # 1h = 15°
                    dec_deg = abs(dec_d) + dec_m/60.0 + dec_s/3600.0
                    if dec_d < 0:
                        dec_deg = -dec_deg
                    
                    print(f"✅ Coordinate convertite!")
                    print(f"📊 RA: {ra_deg:.6f}°")
                    print(f"📊 Dec: {dec_deg:.6f}°")
                    
                    return ra_deg, dec_deg, "SUCCESS"
                else:
                    print(f"❌ Pattern RA/Dec non riconosciuto")
                    # Prova pattern più semplice senza secondi
                    pattern2 = r'(\d+)\s+(\d+)\s+(\d+\.\d+)\s*([-+]?\d+)\s+(\d+)\s+(\d+\.\d+)'
                    match2 = re.search(pattern2, line)
                    if match2:
                        print(f"🔍 Pattern semplice trovato")
                        ra_h = int(match2.group(1))
                        ra_m = int(match2.group(2))
                        ra_s = float(match2.group(3))
                        dec_d = int(match2.group(4))
                        dec_m = int(match2.group(5))
                        dec_s = float(match2.group(6))
                        
                        ra_deg = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0
                        dec_deg = abs(dec_d) + dec_m/60.0 + dec_s/3600.0
                        if dec_d < 0:
                            dec_deg = -dec_deg
                            
                        return ra_deg, dec_deg, "SUCCESS"
        
        return None, None, "NO_COORDINATES_FOUND"
        
    except Exception as e:
        return None, None, f"ERROR: {e}"

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 query_horizons_simple.py <target> <mjd>")
        print("Esempio: python3 query_horizons_simple.py 433 60006.024306")
        sys.exit(1)
    
    target = sys.argv[1]
    mjd = float(sys.argv[2])
    
    print("🌌 JPL HORIZONS QUERY SEMPLIFICATA")
    print("=" * 50)
    
    ra, dec, status = query_horizons_simple(target, mjd)
    
    print("\n" + "=" * 50)
    print("📋 RISULTATO:")
    
    if ra is not None:
        print(f"✅ SUCCESS: {status}")
        print(f"📊 RA: {ra:.6f}°")
        print(f"📊 Dec: {dec:.6f}°")
        
        result = {"success": True, "ra_deg": ra, "dec_deg": dec, "status": status}
        print(f"\n📄 JSON: {json.dumps(result)}")
    else:
        print(f"❌ FAILED: {status}")
        result = {"success": False, "error": status}
        print(f"\n📄 JSON: {json.dumps(result)}")

if __name__ == "__main__":
    main()