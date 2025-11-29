#!/usr/bin/env python3
"""
Confronto propagazione (203) Pompeja:
IOccultCalc vs OrbFit

Converte elementi equinoziali OrbFit a vettori cartesiani
e confronta con output IOccultCalc
"""

import numpy as np

def equinoctial_to_cartesian(a, h, k, p, q, lam, mu=0.01720209895**2):
    """
    Converte elementi equinoziali a coordinate cartesiane
    
    Parametri:
    a   - semiasse maggiore (AU)
    h   - e*sin(omega+Omega)
    k   - e*cos(omega+Omega)
    p   - tan(i/2)*sin(Omega)
    q   - tan(i/2)*cos(Omega)
    lam - longitudine media (radianti)
    mu  - parametro gravitazionale Sole (k^2 in unità gaussiane)
    
    Ritorna:
    (x, y, z, vx, vy, vz) in AU e AU/day
    """
    
    # Eccentricità
    e = np.sqrt(h**2 + k**2)
    
    # Calcola anomalia eccentrica da anomalia media
    # lam = M + omega + Omega
    # Approssimazione: usiamo lam come mean longitude
    omega_plus_Omega = np.arctan2(h, k)
    M = lam - omega_plus_Omega
    
    # Risolvi equazione di Keplero: E - e*sin(E) = M
    E = M
    for _ in range(20):
        E = M + e * np.sin(E)
    
    # Anomalia vera
    nu = 2.0 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E/2),
        np.sqrt(1 - e) * np.cos(E/2)
    )
    
    # Raggio orbitale
    r = a * (1 - e * np.cos(E))
    
    # Coordinate nel piano orbitale
    X_orb = r * np.cos(nu)
    Y_orb = r * np.sin(nu)
    
    # Velocità nel piano orbitale
    n = np.sqrt(mu / a**3)  # mean motion
    VX_orb = -n * a / np.sqrt(1 - e**2) * np.sin(E)
    VY_orb = n * a / np.sqrt(1 - e**2) * np.cos(E)
    
    # Matrice di rotazione dal piano orbitale a eclittica
    # usando elementi equinoziali p, q
    f = 1.0 / (1 + p**2 + q**2)
    
    # Matrice di trasformazione
    # Da Broucke & Cefola 1972
    a11 = 1 - 2*p**2 * f
    a12 = 2*p*q * f
    a13 = -2*p * f
    
    a21 = 2*p*q * f
    a22 = 1 - 2*q**2 * f
    a23 = 2*q * f
    
    a31 = 2*p * f
    a32 = -2*q * f
    a33 = 1 - 2*(p**2 + q**2) * f
    
    # Trasforma a coordinate eclittiche J2000
    x = a11 * X_orb + a12 * Y_orb
    y = a21 * X_orb + a22 * Y_orb
    z = a31 * X_orb + a32 * Y_orb
    
    vx = a11 * VX_orb + a12 * VY_orb
    vy = a21 * VX_orb + a22 * VY_orb
    vz = a31 * VX_orb + a32 * VY_orb
    
    return x, y, z, vx, vy, vz


print("=" * 70)
print("CONFRONTO PROPAGAZIONE (203) POMPEJA")
print("=" * 70)
print()

# Elementi da OrbFit @ MJD 61192.0 (file 203.oel)
print("ORBFIT - Elementi equinoziali @ MJD 61192.0:")
a_orbfit = 2.7368706317538978
h_orbfit = 0.044799304244679
k_orbfit = 0.041830118242835
p_orbfit = -0.005958715738449
q_orbfit = 0.027053268901305
lam_orbfit = 154.0760142434613 * np.pi / 180.0

print(f"  a = {a_orbfit} AU")
print(f"  h = {h_orbfit}")
print(f"  k = {k_orbfit}")
print(f"  p = {p_orbfit}")
print(f"  q = {q_orbfit}")
print(f"  λ = {lam_orbfit * 180/np.pi}°")
print()

# Converti a cartesiane
x_of, y_of, z_of, vx_of, vy_of, vz_of = equinoctial_to_cartesian(
    a_orbfit, h_orbfit, k_orbfit, p_orbfit, q_orbfit, lam_orbfit
)

print("ORBFIT - Vettori di stato (convertiti):")
print(f"  X  = {x_of:17.12f} AU")
print(f"  Y  = {y_of:17.12f} AU")
print(f"  Z  = {z_of:17.12f} AU")
print(f"  VX = {vx_of:17.12f} AU/day")
print(f"  VY = {vy_of:17.12f} AU/day")
print(f"  VZ = {vz_of:17.12f} AU/day")
print()

# Risultati IOccultCalc dal test
print("IOCCULTCALC - Vettori di stato @ MJD 61192.0:")
x_io = -2.637101661238
y_io = 0.842822657440
z_io = 0.385823846490
vx_io = -0.003918341502
vy_io = -0.008383013015
vz_io = -0.004238850883

print(f"  X  = {x_io:17.12f} AU")
print(f"  Y  = {y_io:17.12f} AU")
print(f"  Z  = {z_io:17.12f} AU")
print(f"  VX = {vx_io:17.12f} AU/day")
print(f"  VY = {vy_io:17.12f} AU/day")
print(f"  VZ = {vz_io:17.12f} AU/day")
print()

# Calcola differenze
AU_TO_KM = 149597870.7
dx = (x_io - x_of) * AU_TO_KM
dy = (y_io - y_of) * AU_TO_KM
dz = (z_io - z_of) * AU_TO_KM
dr = np.sqrt(dx**2 + dy**2 + dz**2)

dvx = (vx_io - vx_of) * AU_TO_KM  # km/day
dvy = (vy_io - vy_of) * AU_TO_KM
dvz = (vz_io - vz_of) * AU_TO_KM
dv = np.sqrt(dvx**2 + dvy**2 + dvz**2)

print("=" * 70)
print("DIFFERENZE (IOccultCalc - OrbFit):")
print("=" * 70)
print()
print("Posizione:")
print(f"  ΔX = {dx:12.3f} km")
print(f"  ΔY = {dy:12.3f} km")
print(f"  ΔZ = {dz:12.3f} km")
print(f"  |Δr| = {dr:12.3f} km")
print()

print("Velocità:")
print(f"  ΔVX = {dvx:12.6f} km/day")
print(f"  ΔVY = {dvy:12.6f} km/day")
print(f"  ΔVZ = {dvz:12.6f} km/day")
print(f"  |Δv| = {dv:12.6f} km/day")
print()

# Propagazione era di 192 giorni
dt_days = 192.0
error_per_day = dr / dt_days
error_per_year = error_per_day * 365.25

print("Analisi errore:")
print(f"  Periodo propagazione: {dt_days} giorni")
print(f"  Errore totale: {dr:.1f} km")
print(f"  Errore/giorno: {error_per_day:.1f} km/day")
print(f"  Errore/anno: {error_per_year:.1f} km/year")
print()

# Percentuale rispetto alla distanza
r_orbfit = np.sqrt(x_of**2 + y_of**2 + z_of**2) * AU_TO_KM
error_pct = (dr / r_orbfit) * 100

print(f"  Distanza dal Sole: {r_orbfit:.0f} km")
print(f"  Errore relativo: {error_pct:.6f}%")
print()

print("=" * 70)
