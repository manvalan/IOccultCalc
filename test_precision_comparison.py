#!/usr/bin/env python3
"""
Test confronto preciso IOccultCalc vs JPL Horizons
=================================================

Confronta la posizione dell'asteroide calcolata da IOccultCalc
con i dati di riferimento di JPL Horizons.
"""

import sys
import subprocess
import json
import math
from datetime import datetime

def convert_mjd_to_jd(mjd):
    """Converte MJD in JD"""
    return mjd + 2400000.5

def query_jpl_horizons(asteroid_num, mjd):
    """Query JPL Horizons e ritorna RA, Dec in gradi"""
    
    try:
        result = subprocess.run([
            'python3', 'query_horizons_simple.py', 
            str(asteroid_num), str(mjd)
        ], capture_output=True, text=True)
        
        if result.returncode != 0:
            return None, None, f"JPL Error: {result.stderr}"
        
        # Estrai JSON dalla risposta
        lines = result.stdout.split('\n')
        json_line = None
        for line in lines:
            if line.strip().startswith('{"success":'):
                json_line = line.strip()
                break
        
        if not json_line:
            return None, None, "No JSON found in JPL response"
        
        data = json.loads(json_line)
        
        if not data.get('success', False):
            return None, None, f"JPL failed: {data.get('error', 'Unknown')}"
        
        return data['ra_deg'], data['dec_deg'], "SUCCESS"
        
    except Exception as e:
        return None, None, f"Exception: {e}"

def test_ioccultcalc_position(asteroid_num, jd):
    """Testa la posizione calcolata da IOccultCalc
    
    Nota: questo è un test concettuale. IOccultCalc attuale
    ha problemi con file mancanti, ma mostriamo il framework
    per confronti futuri.
    """
    
    print(f"🔄 Test IOccultCalc per asteroide {asteroid_num} al JD {jd}")
    print("   ⚠️  IOccultCalc ha attualmente problemi con:")
    print("      - File effemeridi JPL mancanti (de441.bsp)")
    print("      - Gaia database in downtime")
    print("      - Dipendenze di cataloghi stellari")
    print()
    print("   📋 Per un test completo servirebbe:")
    print("      1. File de441.bsp scaricato")
    print("      2. Gaia funzionante O catalogo alternativo")
    print("      3. Versione IOccultCalc che calcoli solo posizione asteroide")
    
    # Per ora restituiamo valori fittizi per mostrare il framework
    # In un test reale, qui ci sarebbe la chiamata a IOccultCalc
    
    return None, None, "IOccultCalc_not_ready"

def calculate_angular_separation(ra1, dec1, ra2, dec2):
    """Calcola separazione angolare in arcosecondi"""
    
    # Converte in radianti
    ra1_rad = math.radians(ra1)
    dec1_rad = math.radians(dec1)
    ra2_rad = math.radians(ra2)
    dec2_rad = math.radians(dec2)
    
    # Formula della distanza angolare sulla sfera
    cos_sep = (math.sin(dec1_rad) * math.sin(dec2_rad) +
               math.cos(dec1_rad) * math.cos(dec2_rad) * 
               math.cos(ra1_rad - ra2_rad))
    
    # Previeni errori numerici
    cos_sep = max(-1.0, min(1.0, cos_sep))
    
    separation_rad = math.acos(cos_sep)
    separation_arcsec = math.degrees(separation_rad) * 3600
    
    return separation_arcsec

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 test_precision_comparison.py <asteroid> <mjd>")
        print("Esempio: python3 test_precision_comparison.py 433 61040.5")
        sys.exit(1)
    
    asteroid = int(sys.argv[1])
    mjd = float(sys.argv[2])
    jd = convert_mjd_to_jd(mjd)
    
    print("🌌 CONFRONTO PRECISIONE IOccultCalc vs JPL HORIZONS")
    print("=" * 60)
    print(f"🎯 Asteroide: {asteroid}")
    print(f"📅 MJD: {mjd:.6f}")
    print(f"📅 JD: {jd:.6f}")
    print(f"📅 Data: {datetime.fromtimestamp((jd - 2440587.5) * 86400).strftime('%Y-%m-%d %H:%M:%S')} UTC")
    print()
    
    # Test JPL Horizons (fonte di riferimento)
    print("🔄 Query JPL Horizons (fonte autorevole)...")
    jpl_ra, jpl_dec, jpl_status = query_jpl_horizons(asteroid, mjd)
    
    if jpl_ra is not None:
        print(f"✅ JPL Horizons SUCCESS")
        print(f"   📊 RA:  {jpl_ra:.8f}°")
        print(f"   📊 Dec: {jpl_dec:.8f}°")
        
        # Converte in formato sessagesimale per leggibilità
        ra_h = int(jpl_ra / 15)
        ra_m = int((jpl_ra / 15 - ra_h) * 60)
        ra_s = ((jpl_ra / 15 - ra_h) * 60 - ra_m) * 60
        
        dec_d = int(abs(jpl_dec))
        dec_m = int((abs(jpl_dec) - dec_d) * 60)
        dec_s = ((abs(jpl_dec) - dec_d) * 60 - dec_m) * 60
        dec_sign = "-" if jpl_dec < 0 else "+"
        
        print(f"   📍 RA:  {ra_h:02d}h {ra_m:02d}m {ra_s:05.2f}s")
        print(f"   📍 Dec: {dec_sign}{dec_d:02d}° {dec_m:02d}' {dec_s:04.1f}\"")
        
    else:
        print(f"❌ JPL Horizons FAILED: {jpl_status}")
        return 1
    
    print()
    
    # Test IOccultCalc
    print("🔄 Test IOccultCalc...")
    ioc_ra, ioc_dec, ioc_status = test_ioccultcalc_position(asteroid, jd)
    
    if ioc_ra is not None:
        print(f"✅ IOccultCalc SUCCESS")
        print(f"   📊 RA:  {ioc_ra:.8f}°")
        print(f"   📊 Dec: {ioc_dec:.8f}°")
        
        # Calcola errore
        separation = calculate_angular_separation(jpl_ra, jpl_dec, ioc_ra, ioc_dec)
        ra_error = (ioc_ra - jpl_ra) * 3600 * math.cos(math.radians(jpl_dec))  # arcsec
        dec_error = (ioc_dec - jpl_dec) * 3600  # arcsec
        
        print()
        print("📊 ANALISI PRECISIONE:")
        print(f"   🎯 Separazione totale: {separation:.2f}\"")
        print(f"   📏 Errore RA:  {ra_error:+.2f}\" (cos δ corrected)")
        print(f"   📏 Errore Dec: {dec_error:+.2f}\"")
        
        if separation < 1.0:
            print("   ✅ ECCELLENTE precisione (<1\")")
        elif separation < 10.0:
            print("   ✅ BUONA precisione (<10\")")
        elif separation < 60.0:
            print("   ⚠️  MODERATA precisione (<60\")")
        else:
            print("   ❌ SCARSA precisione (>60\")")
    
    else:
        print(f"❌ IOccultCalc FAILED: {ioc_status}")
    
    print()
    print("=" * 60)
    print("📋 RIFERIMENTO JPL HORIZONS VALIDATO:")
    print(f"   Asteroide {asteroid} al MJD {mjd:.6f}")
    print(f"   RA: {jpl_ra:.8f}° | Dec: {jpl_dec:.8f}°")
    print("   Questo è il valore di riferimento per validazioni future.")

if __name__ == "__main__":
    main()