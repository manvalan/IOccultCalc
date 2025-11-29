#!/usr/bin/env python3
"""
Script per confrontare propagazione IOccultCalc con JPL Horizons.

Uso:
    python3 compare_jpl_horizons.py

Richiede:
    pip install astroquery
"""

import sys
from datetime import datetime, timedelta

try:
    from astroquery.jplhorizons import Horizons
except ImportError:
    print("Errore: astroquery non installato")
    print("Installare con: pip3 install astroquery")
    sys.exit(1)

def format_ra_hms(ra_deg):
    """Converti RA da gradi a ore:min:sec"""
    ra_hours = ra_deg / 15.0
    h = int(ra_hours)
    rem = (ra_hours - h) * 60.0
    m = int(rem)
    s = (rem - m) * 60.0
    return f"{h:02d}:{m:02d}:{s:06.3f}"

def format_dec_dms(dec_deg):
    """Converti Dec da gradi a gra:min:sec"""
    sign = '+' if dec_deg >= 0 else '-'
    abs_dec = abs(dec_deg)
    d = int(abs_dec)
    rem = (abs_dec - d) * 60.0
    m = int(rem)
    s = (rem - m) * 60.0
    return f"{sign}{d:02d}:{m:02d}:{s:05.2f}"

def angular_separation(ra1, dec1, ra2, dec2):
    """Calcola separazione angolare in arcsec"""
    import math
    
    ra1_r = math.radians(ra1)
    dec1_r = math.radians(dec1)
    ra2_r = math.radians(ra2)
    dec2_r = math.radians(dec2)
    
    dra = ra2_r - ra1_r
    ddec = dec2_r - dec1_r
    
    a = math.sin(ddec/2)**2 + math.cos(dec1_r)*math.cos(dec2_r)*math.sin(dra/2)**2
    sep_rad = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    
    return sep_rad * 180.0 / math.pi * 3600.0

def query_jpl_horizons(target_id, dates):
    """
    Query JPL Horizons per posizioni geocentriche.
    
    Args:
        target_id: ID asteroide (es. '11234')
        dates: Lista di date in formato 'YYYY-MM-DD'
    
    Returns:
        Lista di dict con RA, Dec, distanze
    """
    print(f"\n🔍 Query JPL Horizons per asteroide {target_id}...")
    
    # Query Horizons
    obj = Horizons(
        id=target_id,
        location='500@0',  # Geocentro
        epochs={'start': dates[0], 'stop': dates[-1], 'step': '1d'}
    )
    
    # Effemeridi astrometriche
    eph = obj.ephemerides(quantities='1,9,19,20')  # RA, Dec, distanze, mag
    
    results = []
    for date_str in dates:
        # Trova riga corrispondente alla data
        date_obj = datetime.strptime(date_str, '%Y-%m-%d')
        
        for row in eph:
            # Parse datetime_str che è nel formato '2025-Sep-21 00:00'
            row_date_str = str(row['datetime_str']).strip()
            
            # Converti formato JPL (es. '2025-Sep-21') in datetime
            try:
                row_date = datetime.strptime(row_date_str.split()[0], '%Y-%b-%d')
            except:
                # Prova altri formati
                try:
                    row_date = datetime.strptime(row_date_str[:10], '%Y-%m-%d')
                except:
                    continue
            
            if row_date == date_obj:
                results.append({
                    'date': date_str,
                    'ra': float(row['RA']),
                    'dec': float(row['DEC']),
                    'delta': float(row['delta']),
                    'r': float(row['r']),
                    'mag': float(row['V']) if 'V' in row.colnames else 0.0
                })
                break
    
    return results

def main():
    print("╔═══════════════════════════════════════════════════════════════╗")
    print("║  Confronto IOccultCalc vs JPL Horizons                      ║")
    print("║  Asteroide: (11234) 1999 JS82                                ║")
    print("╚═══════════════════════════════════════════════════════════════╝")
    
    # Dati IOccultCalc (dall'output del programma C++)
    ioccult_results = [
        {
            'date': '2025-09-21',
            'label': '-60 giorni',
            'ra': 223.093281,
            'dec': -7.231270,
            'delta': 3.7959,
            'r': 2.8120,
            'mag': 18.33
        },
        {
            'date': '2025-11-20',
            'label': 'Epoca',
            'ra': 246.826162,
            'dec': -13.117994,
            'delta': 3.5587,
            'r': 2.8087,
            'mag': 18.60
        },
        {
            'date': '2026-01-19',
            'label': '+60 giorni',
            'ra': 267.840967,
            'dec': -14.686056,
            'delta': 2.9132,
            'r': 2.7989,
            'mag': 18.43
        }
    ]
    
    # Query JPL Horizons
    dates = [r['date'] for r in ioccult_results]
    try:
        jpl_results = query_jpl_horizons('11234', dates)
    except Exception as e:
        print(f"\n❌ Errore query JPL Horizons: {e}")
        print("\nVerificare connessione internet e riprovare.")
        print("In alternativa, consultare manualmente:")
        print("  https://ssd.jpl.nasa.gov/horizons.cgi")
        return 1
    
    print("\n✓ Query completata!\n")
    
    # Confronto risultati
    print("═══════════════════════════════════════════════════════════════")
    print("  CONFRONTO RISULTATI")
    print("═══════════════════════════════════════════════════════════════\n")
    
    for ioc, jpl in zip(ioccult_results, jpl_results):
        print(f"▸ {ioc['label'].upper()}: {ioc['date']}")
        print(f"  {'':20s} {'IOccultCalc':>20s}  {'JPL Horizons':>20s}  {'Δ (arcsec)':>12s}")
        print(f"  {'-'*75}")
        
        # RA
        ra_diff_arcsec = (jpl['ra'] - ioc['ra']) * 3600.0 * 15.0  # *15 perché RA in ore
        print(f"  {'RA':20s} {format_ra_hms(ioc['ra']):>20s}  {format_ra_hms(jpl['ra']):>20s}  {ra_diff_arcsec:>12.2f}")
        print(f"  {'  (gradi)':20s} {ioc['ra']:>20.6f}  {jpl['ra']:>20.6f}")
        
        # Dec
        dec_diff_arcsec = (jpl['dec'] - ioc['dec']) * 3600.0
        print(f"  {'Dec':20s} {format_dec_dms(ioc['dec']):>20s}  {format_dec_dms(jpl['dec']):>20s}  {dec_diff_arcsec:>12.2f}")
        print(f"  {'  (gradi)':20s} {ioc['dec']:>20.6f}  {jpl['dec']:>20.6f}")
        
        # Separazione angolare totale
        sep_arcsec = angular_separation(ioc['ra'], ioc['dec'], jpl['ra'], jpl['dec'])
        print(f"  {'Separazione totale':20s} {'':<20s}  {'':<20s}  {sep_arcsec:>12.2f}")
        
        # Distanze
        delta_diff = jpl['delta'] - ioc['delta']
        r_diff = jpl['r'] - ioc['r']
        print(f"  {'Δ (AU)':20s} {ioc['delta']:>20.4f}  {jpl['delta']:>20.4f}  {delta_diff:>12.6f}")
        print(f"  {'r (AU)':20s} {ioc['r']:>20.4f}  {jpl['r']:>20.4f}  {r_diff:>12.6f}")
        
        print()
    
    # Riepilogo
    print("═══════════════════════════════════════════════════════════════")
    print("  ANALISI ACCURATEZZA")
    print("═══════════════════════════════════════════════════════════════\n")
    
    separations = []
    for ioc, jpl in zip(ioccult_results, jpl_results):
        sep = angular_separation(ioc['ra'], ioc['dec'], jpl['ra'], jpl['dec'])
        separations.append(sep)
    
    print(f"Separazione media:    {sum(separations)/len(separations):8.2f} arcsec")
    print(f"Separazione massima:  {max(separations):8.2f} arcsec")
    print(f"Separazione minima:   {min(separations):8.2f} arcsec\n")
    
    print("╔═══════════════════════════════════════════════════════════════╗")
    print("║  INTERPRETAZIONE                                             ║")
    print("╠═══════════════════════════════════════════════════════════════╣")
    
    max_sep = max(separations)
    if max_sep < 1.0:
        print("║  ✓ ECCELLENTE: < 1\" - Accuratezza eccezionale               ║")
        print("║    Propagatore equivalente a JPL per questo intervallo       ║")
    elif max_sep < 5.0:
        print("║  ✓ BUONO: < 5\" - Accuratezza adeguata per survey            ║")
        print("║    Sufficiente per predizioni occultazioni ±60 giorni        ║")
    elif max_sep < 20.0:
        print("║  ~ ACCETTABILE: < 20\" - Propagazione Kepleriana semplice    ║")
        print("║    OK per test, ma usare RKF78 con perturbazioni per        ║")
        print("║    predizioni reali di occultazioni                          ║")
    else:
        print("║  ✗ INSUFFICIENTE: > 20\" - Errore eccessivo                  ║")
        print("║    RICHIESTO: Propagatore perturbato (RKF78 + pianeti)      ║")
    
    print("║                                                               ║")
    print("║  NOTA: Questo test usa propagazione Kepleriana semplificata. ║")
    print("║  Per accuratezza < 1\", usare astdyn_propagator con RKF78     ║")
    print("║  e perturbazioni complete (8 pianeti + 16 asteroidi).        ║")
    print("╚═══════════════════════════════════════════════════════════════╝\n")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
