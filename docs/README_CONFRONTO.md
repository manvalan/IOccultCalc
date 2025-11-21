# Confronto Previsioni Occultazioni Asteroidali

## Documento: confronto_occultazioni.pdf

### Descrizione

Documento tecnico in italiano che presenta un'analisi comparativa dettagliata tra le previsioni di occultazioni asteroidali generate da **IOccultCalc** e quelle di riferimento di **Steve Preston** (asteroidoccultation.com).

### Contenuto (16 pagine)

1. **Introduzione** (pagina 2-3)
   - Obiettivi dell'analisi
   - Metodologia di confronto
   - Effemeridi e parametri utilizzati

2. **Eventi Analizzati** (pagina 3)
   - Selezione di 5 eventi rappresentativi
   - 2 eventi futuri: (433) Eros, (15) Eunomia
   - 3 eventi passati: (16) Psyche, (704) Interamnia, (10) Hygiea

3. **Analisi Dettagliata Evento per Evento** (pagine 4-11)
   
   **Evento 1: (433) Eros - 15 Marzo 2026**
   - Scheda completa formato IOTA (IOccultCalc)
   - Scheda formato Preston compatto
   - Tabella comparativa parametri
   - Agreement Score: 96% (Excellent)
   
   **Evento 2: (15) Eunomia - 8 Maggio 2026**
   - Agreement Score: 98% (Excellent)
   - Diff. temporale: +5.1 secondi
   - RMS path: 3.7 km
   
   **Evento 3: (16) Psyche - 22 Settembre 2025**
   - Agreement Score: 89% (Very Good)
   - Diff. temporale: +9.6 secondi
   - RMS path: 8.4 km
   
   **Evento 4: (704) Interamnia - 14 Luglio 2025**
   - Agreement Score: 97% (Excellent)
   - Diff. temporale: -2.9 secondi
   - RMS path: 4.1 km
   
   **Evento 5: (10) Hygiea - 3 Dicembre 2024**
   - Agreement Score: 82% (Good)
   - Diff. temporale: +16.1 secondi
   - RMS path: 12.8 km

4. **Analisi Statistica Complessiva** (pagine 12-13)
   - Distribuzione differenze temporali
   - Distribuzione agreement scores
   - Deviazioni path RMS

5. **Discussione** (pagine 13-14)
   - Interpretazione risultati
   - Confronto con letteratura e standard IOTA
   - Fattori che influenzano le differenze

6. **Conclusioni** (pagine 14-15)
   - Sintesi risultati
   - Validazione metodologia IOccultCalc
   - Raccomandazioni operative
   - Sviluppi futuri

### Risultati Principali

- **Agreement Score Medio**: 92.4%
- **Precisione Temporale RMS**: 9.2 secondi
- **Deviazione Path Media**: 6.3 km RMS (< 5% del diametro)
- **Coordinate Stellari**: < 0.2 arcsec (compatibile con Gaia DR3)

### Conclusione

**IOccultCalc è validato** come strumento affidabile per previsioni di occultazioni asteroidali, con prestazioni superiori agli standard IOTA e accordo eccellente con le previsioni di Steve Preston.

### Generazione PDF

Il documento è stato generato con LaTeX:

```bash
cd docs
pdflatex confronto_occultazioni.tex
pdflatex confronto_occultazioni.tex  # Seconda passata per TOC
```

### Autore

**Michele Bigi**  
Gruppo Astrofili Massesi  
Email: mikbigi@gmail.com

### Data

21 Novembre 2025

### Licenza

MIT License - Stesso del progetto IOccultCalc
