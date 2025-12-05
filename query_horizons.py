#!/usr/bin/env python3
"""
JPL Horizons Query Script per IOccultCalc
=========================================

Script Python per interrogare JPL Horizons Web API e ottenere
coordinate di asteroidi per validazione.

Uso:
    python3 query_horizons.py <target_id> <mjd> [observer_code]

Esempi:
    python3 query_horizons.py 433 60006.024306         # Eros al 28 Nov 2025
    python3 query_horizons.py 17030 61000.0            # 17030 all'epoca elementi
    python3 query_horizons.py 1 60970.0 500           # Ceres con osservatore geocentro
"""

import sys
import requests
import json
from datetime import datetime, timedelta
from typing import Tuple, Optional

class HorizonsQuery:
    """Client per JPL Horizons Web API"""
    
    BASE_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'IOccultCalc/1.0 (https://github.com/manvalan/IOccultCalc)'
        })
    
    def mjd_to_jd(self, mjd: float) -> float:
        """Converti MJD in JD"""
        return mjd + 2400000.5
    
    def jd_to_calendar(self, jd: float) -> str:
        """Converti JD in data calendario per Horizons"""
        # JD 2400000.5 = 1858-11-17 00:00:00 MJD epoch
        delta_days = jd - 2400000.5
        epoch = datetime(1858, 11, 17)
        target_date = epoch + timedelta(days=delta_days)
        # Formato Horizons: 'YYYY-MM-DD HH:MM'
        return target_date.strftime("%Y-%m-%d %H:%M")
    
    def query_coordinates(self, target_id: str, mjd: float, 
                         observer: str = "500") -> Tuple[Optional[float], Optional[float], str]:
        """
        Query coordinate RA/Dec da JPL Horizons
        
        Args:
            target_id: ID asteroide (es: "433", "17030")
            mjd: Modified Julian Date
            observer: Codice osservatore ("500" = geocentro)
            
        Returns:
            (ra_deg, dec_deg, status_message)
        """
        
        try:
            # Converti MJD in data calendario
            jd = self.mjd_to_jd(mjd)
            date_str = self.jd_to_calendar(jd)
            
            # Prova diversi formati target
            target_formats = [
                target_id,                    # "433"
                f"2{target_id.zfill(6)}",    # "2000433" 
                f"{target_id};",              # "433;"
                target_id.zfill(6)            # "000433"
            ]
            
            print(f"🔍 Query JPL Horizons:")
            print(f"   Target: {target_id}")
            print(f"   MJD: {mjd:.6f}")
            print(f"   JD: {jd:.6f}")
            print(f"   Date: {date_str}")
            print(f"   Observer: {observer}")
            
            for i, target_format in enumerate(target_formats, 1):
                print(f"\n📡 Tentativo {i}/4: target='{target_format}'")
                
                try:
                    # Parametri query Horizons
                    params = {
                        'format': 'json',
                        'COMMAND': target_format,
                        'OBJ_DATA': 'NO',
                        'MAKE_EPHEM': 'YES',
                        'EPHEM_TYPE': 'OBSERVER',
                        'CENTER': observer,
                        'START_TIME': date_str,
                        'STOP_TIME': date_str,
                        'STEP_SIZE': '1d',
                        'QUANTITIES': '1',          # Solo RA/Dec apparenti
                        'CSV_FORMAT': 'YES',
                        'CAL_FORMAT': 'CAL'
                    }
                    
                    # Esegui query
                    print(f"   → Query in corso...")
                    response = self.session.get(self.BASE_URL, params=params, timeout=30)
                    
                    if response.status_code != 200:
                        print(f"   ❌ HTTP {response.status_code}")
                        continue
                    
                    # Parse risposta JSON
                    data = response.json()
                    
                    if 'result' not in data:
                        print(f"   ❌ Nessun risultato in risposta")
                        print(f"   📄 Keys disponibili: {list(data.keys())}")
                        continue
                    
                    # Estrai coordinate dalla risposta
                    result_text = data['result']
                    print(f"   📄 Lunghezza risposta: {len(result_text)} caratteri")
                    
                    # Debug: mostra prime 10 linee
                    lines = result_text.split('\n')
                    print(f"   📄 Prime 10 linee:")
                    for i, line in enumerate(lines[:10]):
                        print(f"      {i+1}: {line[:80]}")
                    
                    # Cerca la sezione effemeridi
                    in_ephemeris = False
                    coord_line = None
                    
                    for line in lines:
                        if '$$SOE' in line:  # Start of Ephemeris
                            in_ephemeris = True
                            print(f"   📊 Inizio effemeridi trovato")
                            continue
                        elif '$$EOE' in line:  # End of Ephemeris
                            in_ephemeris = False
                            print(f"   📊 Fine effemeridi")
                            break
                        elif in_ephemeris and line.strip():
                            # Linea di effemeridi
                            coord_line = line.strip()
                            print(f"   📍 Linea effemeridi: {coord_line[:100]}")
                            break
                    
                    if not coord_line:
                        print(f"   ❌ Coordinate non trovate nella risposta")
                        continue
                    
                    # Parse coordinate CSV
                    parts = coord_line.split(',')
                    if len(parts) < 3:
                        print(f"   ❌ Formato CSV invalido: {len(parts)} campi")
                        continue
                    
                    # Estrai RA e Dec (formato: HH MM SS.SS, ±DD MM SS.S)
                    ra_str = parts[1].strip()    # RA
                    dec_str = parts[2].strip()   # Dec
                    
                    print(f"   📍 RA raw: '{ra_str}'")
                    print(f"   📍 Dec raw: '{dec_str}'")
                    
                    # Converti RA (HH MM SS.SS → degrees)
                    ra_parts = ra_str.split()
                    if len(ra_parts) != 3:
                        print(f"   ❌ Formato RA invalido")
                        continue
                    
                    ra_h = float(ra_parts[0])
                    ra_m = float(ra_parts[1])  
                    ra_s = float(ra_parts[2])
                    ra_deg = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0  # 15°/ora
                    
                    # Converti Dec (±DD MM SS.S → degrees)
                    dec_parts = dec_str.split()
                    if len(dec_parts) != 3:
                        print(f"   ❌ Formato Dec invalido")
                        continue
                    
                    dec_d = float(dec_parts[0])
                    dec_m = float(dec_parts[1])
                    dec_s = float(dec_parts[2])
                    
                    dec_deg = abs(dec_d) + dec_m/60.0 + dec_s/3600.0
                    if dec_d < 0:
                        dec_deg = -dec_deg
                    
                    print(f"   ✅ SUCCESS!")
                    print(f"   📊 RA: {ra_deg:.6f}°")
                    print(f"   📊 Dec: {dec_deg:.6f}°")
                    
                    return ra_deg, dec_deg, f"SUCCESS with {target_format}"
                    
                except requests.Timeout:
                    print(f"   ❌ Timeout (30s)")
                    continue
                except requests.RequestException as e:
                    print(f"   ❌ Network error: {e}")
                    continue
                except (ValueError, IndexError, KeyError) as e:
                    print(f"   ❌ Parse error: {e}")
                    continue
                except Exception as e:
                    print(f"   ❌ Unexpected error: {e}")
                    continue
            
            return None, None, "ALL_FORMATS_FAILED"
            
        except Exception as e:
            return None, None, f"GENERAL_ERROR: {e}"

def main():
    """Main function"""
    
    if len(sys.argv) < 3:
        print("Uso: python3 query_horizons.py <target_id> <mjd> [observer_code]")
        print()
        print("Esempi:")
        print("  python3 query_horizons.py 433 60006.024306")
        print("  python3 query_horizons.py 17030 61000.0")
        print("  python3 query_horizons.py 1 60970.0 500")
        sys.exit(1)
    
    target_id = sys.argv[1]
    mjd = float(sys.argv[2])
    observer = sys.argv[3] if len(sys.argv) > 3 else "500"
    
    print("🌌 JPL HORIZONS QUERY per IOccultCalc")
    print("=" * 50)
    
    # Crea client e query
    client = HorizonsQuery()
    ra_deg, dec_deg, status = client.query_coordinates(target_id, mjd, observer)
    
    print("\n" + "=" * 50)
    print("📋 RISULTATO FINALE:")
    
    if ra_deg is not None and dec_deg is not None:
        print(f"✅ SUCCESS: {status}")
        print(f"📊 RA:  {ra_deg:.6f}°")
        print(f"📊 Dec: {dec_deg:.6f}°")
        
        # Output formato per C++
        print(f"\n🔧 Per C++:")
        print(f"   double jpl_ra = {ra_deg:.6f};")
        print(f"   double jpl_dec = {dec_deg:.6f};")
        
        # Output JSON per integrazione
        result = {
            "success": True,
            "target_id": target_id,
            "mjd": mjd,
            "observer": observer,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "status": status
        }
        print(f"\n📄 JSON:")
        print(json.dumps(result, indent=2))
        
        sys.exit(0)
        
    else:
        print(f"❌ FAILED: {status}")
        result = {
            "success": False,
            "target_id": target_id,
            "mjd": mjd,
            "observer": observer,
            "error": status
        }
        print(f"\n📄 JSON:")
        print(json.dumps(result, indent=2))
        
        sys.exit(1)

if __name__ == "__main__":
    main()