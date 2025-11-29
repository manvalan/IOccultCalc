#!/usr/bin/env python3
"""
Script per convertire file .eq1 (AstDyS) in preset IOccultCalc

Uso:
    python3 eq1_to_preset.py 11234.eq1 > preset_11234.oop
    python3 eq1_to_preset.py 11234.eq1 --output preset_11234.oop
    python3 eq1_to_preset.py 11234.eq1 --period 2026-01-01 2026-02-28
"""

import sys
import math
import argparse
from datetime import datetime

def parse_eq1(filepath):
    """Parse file .eq1 (formato equinoziale OEF2.0)"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    data = {}
    found_object = False
    
    for line in lines:
        line = line.strip()
        
        # Skip commenti e header
        if not line or line.startswith('!'):
            continue
        if 'format' in line or 'rectype' in line or 'refsys' in line:
            continue
        if 'END_OF_HEADER' in line:
            continue
        
        # Oggetto (prima linea non-header)
        if not found_object and 'EQU' not in line and 'MJD' not in line:
            data['object'] = line.strip()
            found_object = True
            continue
        
        # Elementi equinoziali
        if 'EQU' in line:
            parts = line.split()
            idx = parts.index('EQU') + 1
            data['a'] = float(parts[idx])
            data['h'] = float(parts[idx+1])
            data['k'] = float(parts[idx+2])
            data['p'] = float(parts[idx+3])
            data['q'] = float(parts[idx+4])
            data['lambda'] = float(parts[idx+5])
            continue
        
        # Epoca
        if 'MJD' in line:
            parts = line.split()
            idx = parts.index('MJD') + 1
            data['epoch_mjd'] = float(parts[idx])
            continue
        
        # Magnitudine (opzionale)
        if 'MAG' in line:
            parts = line.split()
            idx = parts.index('MAG') + 1
            data['magnitude'] = float(parts[idx])
            continue
    
    return data

def equinoctial_to_keplerian(equ):
    """Converti elementi equinoziali a Kepleriani"""
    a = equ['a']
    h = equ['h']
    k = equ['k']
    p = equ['p']
    q = equ['q']
    lambda_deg = equ['lambda']
    
    # Eccentricità
    e = math.sqrt(h*h + k*k)
    
    # Inclinazione
    tan_half_i = math.sqrt(p*p + q*q)
    i_rad = 2.0 * math.atan(tan_half_i)
    i_deg = math.degrees(i_rad)
    
    # Longitudine nodo ascendente
    Omega_rad = math.atan2(p, q)
    if Omega_rad < 0:
        Omega_rad += 2 * math.pi
    Omega_deg = math.degrees(Omega_rad)
    
    # Longitudine perielio
    omega_plus_Omega_rad = math.atan2(h, k)
    if omega_plus_Omega_rad < 0:
        omega_plus_Omega_rad += 2 * math.pi
    
    # Argomento perielio
    omega_rad = omega_plus_Omega_rad - Omega_rad
    if omega_rad < 0:
        omega_rad += 2 * math.pi
    omega_deg = math.degrees(omega_rad)
    
    # Anomalia media
    lambda_rad = math.radians(lambda_deg)
    M_rad = lambda_rad - omega_plus_Omega_rad
    while M_rad < 0:
        M_rad += 2 * math.pi
    while M_rad >= 2 * math.pi:
        M_rad -= 2 * math.pi
    M_deg = math.degrees(M_rad)
    
    return {
        'a': a,
        'e': e,
        'i': i_deg,
        'Omega': Omega_deg,
        'omega': omega_deg,
        'M': M_deg
    }

def mjd_to_date(mjd):
    """Converti MJD a data calendario"""
    jd = mjd + 2400000.5
    a = int(jd + 32044)
    b = (4*a + 3) // 146097
    c = a - (146097*b) // 4
    d = (4*c + 3) // 1461
    e = c - (1461*d) // 4
    m = (5*e + 2) // 153
    
    day = e - (153*m + 2) // 5 + 1
    month = m + 3 - 12 * (m // 10)
    year = 100*b + d - 4800 + m // 10
    
    return f"{year:04d}-{month:02d}-{day:02d}"

def generate_preset(equ_data, kep_data, args):
    """Genera preset IOccultCalc"""
    
    object_num = equ_data.get('object', 'Unknown')
    epoch_mjd = equ_data.get('epoch_mjd', 60000.0)
    epoch_date = mjd_to_date(epoch_mjd)
    
    # Date periodo
    if args.period:
        start_date, end_date = args.period
    else:
        # Default: ±60 giorni da epoca
        start_mjd = epoch_mjd - 60
        end_mjd = epoch_mjd + 60
        start_date = mjd_to_date(start_mjd)
        end_date = mjd_to_date(end_mjd)
    
    preset = f"""! =============================================================================
! Preset generato da file .eq1 (AstDyS)
! =============================================================================
! Oggetto: {object_num}
! Epoca elementi: MJD {epoch_mjd:.1f} ({epoch_date})
! Fonte: AstDyS equinoctial elements (OEF2.0)
! Generato: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
! =============================================================================

general.
    ! Identificazione asteroide
    .asteroid_number = {object_num}
    .asteroid_name = 'Asteroid {object_num}'
    
    ! Epoca elementi (MJD TDB)
    .epoch_mjd = {epoch_mjd:.6f}
    
    ! Elementi orbitali Kepleriani (convertiti da equinoziali)
    .semimajor_axis_au = {kep_data['a']:.10f}
    .eccentricity = {kep_data['e']:.10f}
    .inclination_deg = {kep_data['i']:.8f}
    .ascending_node_deg = {kep_data['Omega']:.8f}
    .perihelion_arg_deg = {kep_data['omega']:.8f}
    .mean_anomaly_deg = {kep_data['M']:.8f}
    
    ! Propagatore: RKF78 per accuratezza professionale
    .propagator = 'RKF78'
    .step_size_days = 0.1
    .tolerance = 1.0e-12
    
    ! Output
    .output_format = 'iota'
    .verbose = .TRUE.

propagator.
    ! Perturbazioni planetarie ESSENZIALI per accuratezza
    .use_planetary_perturbations = .TRUE.
    
    ! Pianeti inclusi (tutti 8)
    .include_mercury = .TRUE.
    .include_venus = .TRUE.
    .include_earth = .TRUE.
    .include_mars = .TRUE.
    .include_jupiter = .TRUE.
    .include_saturn = .TRUE.
    .include_uranus = .TRUE.
    .include_neptune = .TRUE.
    
    ! Correzioni relativistiche
    .use_relativistic_corrections = .TRUE.
    
    ! Effemeridi JPL
    .ephemeris_file = 'de441_part-1.bsp'

timespan.
    ! Periodo predizioni
    .start_date = '{start_date}'
    .end_date = '{end_date}'
    .time_step_hours = 1.0

geographic_filter.
    ! Italia + Europa
    .enabled = .TRUE.
    .min_latitude = 36.0
    .max_latitude = 71.0
    .min_longitude = -10.0
    .max_longitude = 40.0

magnitude_filter.
    ! Filtro stelle
    .min_star_magnitude = 4.0
    .max_star_magnitude = 16.0
    .min_drop_magnitude = 0.3

occultation_geometry.
    ! Geometria occultazioni
    .min_duration_seconds = 0.1
    .max_shadow_velocity_km_s = 50.0
    .min_elongation_deg = 45.0
    .max_moon_angle_deg = 30.0

star_catalog.
    ! Catalogo stelle Gaia EDR3
    .catalog = 'gaia_edr3'
    .use_proper_motion = .TRUE.
    .max_search_radius_deg = 2.0
    .parallax_correction = .TRUE.

orbit_fitting.
    ! Fitting DISABILITATO
    ! Elementi da .eq1 sono già fittati da AstDyS con migliaia di osservazioni
    .enable_fitting = .FALSE.
    
    ! Se abilitato in futuro:
    .observation_source = 'mpc'
    .min_observations = 50
    .outlier_threshold_sigma = 3.0
    .min_rms_arcsec = 0.3
    .max_rms_arcsec = 2.0

output.
    .output_file = 'occultations_{object_num}.txt'
    .format = 'detailed'
    .include_uncertainty = .TRUE.
    .coordinate_system = 'J2000'

! =============================================================================
! NOTE TECNICHE
! =============================================================================
! 
! Elementi Equinoziali Originali (da .eq1):
!   a      = {equ_data['a']:.10f} AU
!   h      = {equ_data['h']:.10f}  (= e·sin(ϖ))
!   k      = {equ_data['k']:.10f}  (= e·cos(ϖ))
!   p      = {equ_data['p']:.10f}  (= tan(i/2)·sin(Ω))
!   q      = {equ_data['q']:.10f}  (= tan(i/2)·cos(Ω))
!   lambda = {equ_data['lambda']:.8f}°  (longitudine media)
!
! Conversione Equinoziale → Kepleriano:
!   e = sqrt(h² + k²) = {kep_data['e']:.8f}
!   i = 2·atan(sqrt(p² + q²)) = {kep_data['i']:.6f}°
!   Ω = atan2(p, q) = {kep_data['Omega']:.6f}°
!   ϖ = atan2(h, k) = {math.degrees(math.atan2(equ_data['h'], equ_data['k'])):.6f}°
!   ω = ϖ - Ω = {kep_data['omega']:.6f}°
!   M = λ - ϖ = {kep_data['M']:.6f}°
!
! Accuratezza Attesa (RKF78 + perturbazioni):
!   ±30 giorni: 0.1-0.3 arcsec
!   ±60 giorni: 1-3 arcsec
!   ±90 giorni: 3-10 arcsec
!   ±180 giorni: 10-50 arcsec
!
! Per validare:
!   python3 compare_jpl_horizons.py {object_num} {start_date}
!
! =============================================================================
"""
    
    return preset

def main():
    parser = argparse.ArgumentParser(
        description='Converti file .eq1 (AstDyS) in preset IOccultCalc',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Esempi:
  %(prog)s 11234.eq1
  %(prog)s 11234.eq1 --output preset_11234.oop
  %(prog)s 11234.eq1 --period 2026-01-01 2026-12-31
  %(prog)s external/ITALOccultLibrary/astdyn/data/11234.eq1 > preset.oop
        """
    )
    
    parser.add_argument('eq1_file', help='File .eq1 da convertire')
    parser.add_argument('--output', '-o', help='File output (default: stdout)')
    parser.add_argument('--period', '-p', nargs=2, metavar=('START', 'END'),
                       help='Periodo predizioni (YYYY-MM-DD YYYY-MM-DD)')
    
    args = parser.parse_args()
    
    try:
        # Parse .eq1
        equ_data = parse_eq1(args.eq1_file)
        
        # Converti a Kepleriani
        kep_data = equinoctial_to_keplerian(equ_data)
        
        # Genera preset
        preset = generate_preset(equ_data, kep_data, args)
        
        # Output
        if args.output:
            with open(args.output, 'w') as f:
                f.write(preset)
            print(f"✓ Preset scritto in: {args.output}", file=sys.stderr)
        else:
            print(preset)
        
        # Mostra info
        print(f"\n✓ Conversione completata:", file=sys.stderr)
        print(f"  Oggetto: {equ_data['object']}", file=sys.stderr)
        print(f"  Epoca: MJD {equ_data['epoch_mjd']:.1f} ({mjd_to_date(equ_data['epoch_mjd'])})", file=sys.stderr)
        print(f"  a = {kep_data['a']:.6f} AU", file=sys.stderr)
        print(f"  e = {kep_data['e']:.6f}", file=sys.stderr)
        print(f"  i = {kep_data['i']:.4f}°", file=sys.stderr)
        
    except Exception as e:
        print(f"Errore: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
