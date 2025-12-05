#!/usr/bin/env python3
"""Genera scheda ItalOccult per 17030 - 28 novembre 2025"""
from datetime import datetime

event_data = {
    'asteroid': {'number': 17030, 'name': 'Sierks', 'designation': '1999 FC9', 'diameter_km': 6.2, 'h_mag': 14.7, 'ra_deg': 346.8312, 'dec_deg': -5.5842},
    'star': {'gaia_id': 'Gaia DR3 3411546266140512128', 'ra_deg': 346.8313, 'dec_deg': -5.5842, 'mag': 11.2},
    'event': {'date': '2025-11-28', 'time_utc': '00:35:17', 'jd': 2460642.524271, 'duration_sec': 4.8, 'mag_drop': 3.5, 'separation_mas': 36.0},
    'path': {'center_lat': -32.5, 'center_lon': -15.2, 'width_km': 6.2, 'length_km': 9500}
}

print(f"✓ Dati caricati per ({event_data['asteroid']['number']}) {event_data['asteroid']['name']}")
print(f"  Evento: {event_data['event']['date']} {event_data['event']['time_utc']} UTC")
print(f"  Durata: {event_data['event']['duration_sec']} sec")
print(f"  Stella: {event_data['star']['gaia_id']}")
