#!/usr/bin/env python3
"""
Script per scaricare dati effemeridi da JPL Horizons per 17030
Periodo: 26-29 Novembre 2025, intervallo 6 ore
"""

import requests
from datetime import datetime, timedelta

def query_jpl_horizons(target='17030', start_date='2025-11-26', end_date='2025-11-29', step='6h'):
    """
    Query JPL Horizons API per ottenere effemeridi
    """
    
    # URL API JPL Horizons
    url = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    
    # Parametri query
    params = {
        'format': 'text',
        'COMMAND': f"'{target}'",  # Asteroide 17030
        'OBJ_DATA': 'NO',
        'MAKE_EPHEM': 'YES',
        'EPHEM_TYPE': 'OBSERVER',
        'CENTER': '500@399',  # Geocentrico
        'START_TIME': f"'{start_date}'",
        'STOP_TIME': f"'{end_date}'",
        'STEP_SIZE': f"'{step}'",
        'QUANTITIES': "'1,20'",  # RA/Dec e distanza
        'REF_SYSTEM': 'ICRF',
        'CAL_FORMAT': 'CAL',
        'TIME_DIGITS': 'MINUTES',
        'ANG_FORMAT': 'DEG',
        'APPARENT': 'AIRLESS',
        'RANGE_UNITS': 'AU',
        'SUPPRESS_RANGE_RATE': 'NO',
        'SKIP_DAYLT': 'NO',
        'EXTRA_PREC': 'NO',
        'R_T_S_ONLY': 'NO',
        'CSV_FORMAT': 'YES'
    }
    
    print(f"🌐 Query JPL Horizons per asteroide {target}...")
    print(f"   Periodo: {start_date} to {end_date}, step {step}")
    print(f"   URL: {url}")
    print()
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        
        # Salva risposta raw
        with open('jpl_horizons_17030_raw.txt', 'w') as f:
            f.write(response.text)
        print("✅ Risposta salvata in: jpl_horizons_17030_raw.txt")
        
        # Parse i dati
        lines = response.text.split('\n')
        data_section = False
        ephemeris_data = []
        
        for line in lines:
            # Trova inizio dati
            if '$$SOE' in line:
                data_section = True
                continue
            # Fine dati
            if '$$EOE' in line:
                data_section = False
                break
            
            if data_section and line.strip():
                # Parse CSV line: Date, , , RA, Dec, delta, deldot
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 5:
                    try:
                        date_str = parts[0]  # "2025-Nov-26 00:00"
                        ra_str = parts[3]     # RA in gradi
                        dec_str = parts[4]    # Dec in gradi
                        
                        # Converti data in MJD
                        # Parse: "2025-Nov-26 00:00"
                        from datetime import datetime
                        dt = datetime.strptime(date_str, '%Y-%b-%d %H:%M')
                        # MJD = JD - 2400000.5, JD = giorni da 4713 BC
                        jd = dt.toordinal() + 1721424.5 + dt.hour/24.0 + dt.minute/1440.0
                        mjd = jd - 2400000.5
                        
                        ephemeris_data.append({
                            'date': date_str,
                            'mjd': mjd,
                            'ra': float(ra_str),
                            'dec': float(dec_str)
                        })
                    except (ValueError, IndexError) as e:
                        print(f"⚠️  Errore parsing: {line[:50]}... -> {e}")
                        continue
        
        # Genera codice C++
        print(f"\n✅ Trovate {len(ephemeris_data)} effemeridi")
        print("\n📋 Dati da copiare nel codice C++:\n")
        print("std::vector<JPLData> getJPLData() {")
        print("    std::vector<JPLData> jpl;")
        print("    ")
        
        for data in ephemeris_data:
            print(f"    jpl.push_back({{{data['mjd']:.6f}, {data['ra']:.6f}, {data['dec']:.6f}, 2.85, \"{data['date']}\"}});")
        
        print("    ")
        print("    return jpl;")
        print("}")
        print()
        
        # Salva anche in file
        with open('jpl_horizons_17030_data.txt', 'w') as f:
            f.write("MJD,RA(deg),Dec(deg),Date\n")
            for data in ephemeris_data:
                f.write(f"{data['mjd']:.6f},{data['ra']:.6f},{data['dec']:.6f},{data['date']}\n")
        
        print("✅ Dati salvati in: jpl_horizons_17030_data.txt")
        return ephemeris_data
        
    except requests.exceptions.RequestException as e:
        print(f"❌ Errore nella query JPL: {e}")
        print("\n⚠️  ALTERNATIVA: Usa il sito web manualmente:")
        print("   https://ssd.jpl.nasa.gov/horizons/app.html#/")
        print("   1. Target Body: 17030")
        print("   2. Observer Location: Geocentric [500@399]")
        print("   3. Time Specification: 2025-11-26 to 2025-11-29, step 6h")
        print("   4. Table Settings: RA/Dec + Delta")
        return None

if __name__ == '__main__':
    data = query_jpl_horizons()
    
    if not data:
        print("\n💡 Per query manuale, visita:")
        print("   https://ssd.jpl.nasa.gov/horizons/app.html#/")
