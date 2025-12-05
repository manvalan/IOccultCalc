#!/usr/bin/env python3
"""
Genera scheda ItalOccult per evento 17030 - 28 novembre 2025
Usa i risultati del calcolo occultazione per creare HTML A4
"""

import json
import sys
from datetime import datetime

def generate_italoccult_sheet_17030():
    """Genera scheda ItalOccult per asteroide 17030"""
    
    # Dati evento (da aggiornare con risultati reali)
    event = {
        'asteroid': {
            'number': 17030,
            'name': 'Sierks',
            'designation': '1999 FC9',
            'diameter_km': 6.2,  # Stimato
            'h_mag': 14.7
        },
        'star': {
            'gaia_id': 'Gaia DR3 254527307789188697',  # Da aggiornare
            'ra_deg': 5.253,
            'dec_deg': -0.334,
            'mag': 9.59,
            'common_name': None
        },
        'event': {
            'date': '2025-11-28',
            'time_utc': '18:16:32',  # Da aggiornare con risultato reale
            'jd': 2461008.762,
            'duration_sec': 4.2,  # Da aggiornare
            'mag_drop': 5.11,
            'uncertainty_km': 250  # Stimato
        },
        'path': {
            'center_lat': 12.5,  # Da aggiornare
            'center_lon': 45.3,
            'angle_deg': 125.0,
            'width_km': 6.2,
            'length_km': 8000
        }
    }
    
    html = f"""<!DOCTYPE html>
<html lang="it">
<head>
    <meta charset="UTF-8">
    <title>ItalOccult - ({event['asteroid']['number']}) {event['asteroid']['name']} - {event['event']['date']}</title>
    <style>
        @page {{ size: A4; margin: 10mm; }}
        body {{ 
            font-family: Arial, sans-serif; 
            margin: 0; 
            padding: 0;
            font-size: 10pt;
        }}
        .sheet {{
            width: 210mm;
            height: 297mm;
            padding: 10mm;
            box-sizing: border-box;
        }}
        .header {{
            height: 18%;
            border-bottom: 3px solid #003366;
            padding-bottom: 5mm;
        }}
        h1 {{
            color: #003366;
            margin: 0 0 5px 0;
            font-size: 22pt;
        }}
        h2 {{
            color: #006699;
            margin: 0 0 3px 0;
            font-size: 18pt;
        }}
        .info-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 5mm;
            margin-top: 3mm;
        }}
        .info-box {{
            background: #f0f8ff;
            padding: 3mm;
            border-radius: 3mm;
            border-left: 4px solid #003366;
        }}
        .label {{
            font-weight: bold;
            color: #003366;
        }}
        .middle {{
            height: 40%;
            display: flex;
            gap: 5mm;
            margin-top: 5mm;
        }}
        .map-container {{
            flex: 1;
            border: 2px solid #ccc;
            border-radius: 3mm;
            padding: 3mm;
            background: #fff;
        }}
        .map-title {{
            font-weight: bold;
            color: #003366;
            text-align: center;
            margin-bottom: 2mm;
            font-size: 11pt;
        }}
        .map-frame {{
            height: 90%;
            display: flex;
            align-items: center;
            justify-content: center;
            background: #f9f9f9;
            border-radius: 2mm;
        }}
        .bottom {{
            height: 32%;
            margin-top: 5mm;
            border: 2px solid #ccc;
            border-radius: 3mm;
            padding: 3mm;
            background: #fff;
        }}
        .chart-title {{
            font-weight: bold;
            color: #003366;
            text-align: center;
            margin-bottom: 2mm;
            font-size: 11pt;
        }}
        .footer {{
            height: 5%;
            text-align: center;
            color: #666;
            font-size: 8pt;
            padding-top: 3mm;
            border-top: 1px solid #ccc;
        }}
        .placeholder {{
            color: #999;
            font-style: italic;
            font-size: 9pt;
        }}
    </style>
</head>
<body>
<div class="sheet">
    <!-- HEADER -->
    <div class="header">
        <h1>🌟 Scheda Occultazione ItalOccult</h1>
        <h2>({event['asteroid']['number']}) {event['asteroid']['name']} — {event['event']['date']} {event['event']['time_utc']} UTC</h2>
        
        <div class="info-grid">
            <div class="info-box">
                <div><span class="label">Asteroide:</span> ({event['asteroid']['number']}) {event['asteroid']['name']} ({event['asteroid']['designation']})</div>
                <div><span class="label">Diametro:</span> ~{event['asteroid']['diameter_km']} km | <span class="label">H:</span> {event['asteroid']['h_mag']}</div>
            </div>
            <div class="info-box">
                <div><span class="label">Stella:</span> {event['star']['gaia_id']}</div>
                <div><span class="label">Coordinate:</span> RA {event['star']['ra_deg']:.4f}° | Dec {event['star']['dec_deg']:.4f}° | <span class="label">Mag:</span> {event['star']['mag']}</div>
            </div>
        </div>
        
        <div class="info-grid" style="margin-top: 2mm;">
            <div class="info-box">
                <div><span class="label">Durata:</span> {event['event']['duration_sec']} sec | <span class="label">ΔMag:</span> {event['event']['mag_drop']}</div>
                <div><span class="label">JD:</span> {event['event']['jd']}</div>
            </div>
            <div class="info-box">
                <div><span class="label">Larghezza percorso:</span> {event['path']['width_km']} km</div>
                <div><span class="label">Incertezza:</span> ±{event['event']['uncertainty_km']} km</div>
            </div>
        </div>
    </div>
    
    <!-- MIDDLE SECTION: MAPPE -->
    <div class="middle">
        <div class="map-container">
            <div class="map-title">Carta di Avvicinamento (6° × 4°)</div>
            <div class="map-frame">
                <div class="placeholder">
                    [IOC_StarMap preset: approach_chart]<br>
                    Traiettoria asteroide 10 giorni<br>
                    Stelle Gaia DR3 fino a mag 12<br>
                    Campo: RA {event['star']['ra_deg']:.3f}° Dec {event['star']['dec_deg']:.3f}°
                </div>
            </div>
        </div>
        
        <div class="map-container">
            <div class="map-title">Percorso sulla Terra</div>
            <div class="map-frame">
                <div class="placeholder">
                    [IOC_Earth map]<br>
                    Centro: Lat {event['path']['center_lat']}° Lon {event['path']['center_lon']}°<br>
                    Angolo: {event['path']['angle_deg']}° | Larghezza: {event['path']['width_km']} km<br>
                    Lunghezza: {event['path']['length_km']} km
                </div>
            </div>
        </div>
    </div>
    
    <!-- BOTTOM SECTION: CAMPO DETTAGLIATO -->
    <div class="bottom">
        <div class="chart-title">Campo Dettagliato (3° × 2°)</div>
        <div class="map-frame">
            <div class="placeholder">
                [IOC_StarMap preset: final_chart]<br>
                Campo dettagliato al momento dell'occultazione<br>
                Stelle Gaia DR3 fino a mag 14<br>
                Posizione asteroide + ellisse incertezza<br>
                Separazione prevista: calcolare con dati reali
            </div>
        </div>
    </div>
    
    <!-- FOOTER -->
    <div class="footer">
        <div><strong>ItalOccult</strong> | Generato con IOccultCalc + IOC_StarMap + IOC_Earth</div>
        <div>Data generazione: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | Prediction by IOccultCalc v2.0</div>
    </div>
</div>
</body>
</html>
"""
    
    # Salva file
    filename = f"italoccult_17030_{event['event']['date'].replace('-', '')}.html"
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"✓ Scheda ItalOccult salvata: {filename}")
    print(f"  Asteroide: ({event['asteroid']['number']}) {event['asteroid']['name']}")
    print(f"  Data: {event['event']['date']} {event['event']['time_utc']} UTC")
    print(f"  Stella: Mag {event['star']['mag']}")
    print(f"  Durata: {event['event']['duration_sec']} sec")
    print()
    print("Per convertire in PDF:")
    print(f"  python3 convert_to_pdf.py {filename}")
    print()
    print("NOTA: Aggiornare i dati con i risultati reali del calcolo!")
    
    return filename

if __name__ == '__main__':
    generate_italoccult_sheet_17030()
