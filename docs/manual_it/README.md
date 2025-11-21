# Manuale Scientifico IOccultCalc - Versione Italiana

Manuale scientifico completo per la libreria IOccultCalc, documentando i metodi di previsione delle occultazioni asteroidali ad alta precisione.

## Autore

**Michele Bigi**  
Email: mikbigi@gmail.com  
Organizzazione: Gruppo Astrofili Massesi  
Web: http://www.astrofilimassesi.it

## Contenuti

Il manuale è organizzato in 16 capitoli + 3 appendici:

1. **Introduzione** - Occultazioni asteroidali, stato dell'arte, architettura
2. **Sistemi di Coordinate** - ICRS, trasformazioni, precessione
3. **Sistemi Temporali** - UTC, TAI, TT, TDB, conversioni
4. **Effemeridi Planetarie** - JPL DE441, SPICE, Chebyshev
5. **Meccanica Orbitale** - Elementi equinoziali, equazione di Keplero
6. **Integrazione Numerica** - RKF78, DOPRI853
7. **Perturbazioni** - Modello N-body, 14 corpi
8. **Correzioni Relativistiche** - PN formalism, light-time
9. **Precessione e Nutazione** - IAU2000A/2006
10. **Astrometria Stellare** - Gaia DR3, proper motion
11. **Determinazione Orbite** - Least squares, differential correction
12. **Forma Asteroidi** - Modelli, limb fitting
13. **Metodo Besseliano** - Shadow path, geometria
14. **Propagazione Incertezze** - Covarianza, Monte Carlo
15. **Implementazione** - Architettura C++, API
16. **Validazione** - Confronti, test, risultati

### Appendici

- **A**: Costanti fisiche e astronomiche
- **B**: Algoritmi dettagliati
- **C**: Tabelle JPL DE441

## Compilazione

### Requisiti

- LaTeX distribution completa (TeXLive 2020+, MiKTeX, MacTeX)
- Pacchetti: babel-italian, siunitx, tikz, pgfplots, natbib

### Comandi

```bash
cd docs/manual_it
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Oppure con latexmk:
```bash
latexmk -pdf main.tex
```

### Output

Il comando genera `main.pdf` (~100 pagine, ~1 MB).

## Personalizzazione

### Informazioni Autore

Modificare in `main.tex`:
- `\author{...}` - Nome e affiliazione
- `\date{...}` - Data
- Sezione copyright

### Stile

- Geometry: `docs/manual_it/main.tex` (margini, dimensioni)
- Fonts: modificare preambolo LaTeX
- Colori: hyperref colorlinks

## Struttura File

```
docs/manual_it/
├── main.tex                    # Documento principale
├── README.md                   # Questo file
├── chapters/                   # Capitoli
│   ├── 01_introduzione.tex
│   ├── 02_sistemi_coordinate.tex
│   ├── ... (14 capitoli rimanenti)
│   └── appendice_*.tex
└── figures/                    # Figure (se presenti)
```

## Traduzione

Questa è la versione italiana del manuale originale in inglese (`docs/manual`).

Traduzioni disponibili:
- 🇬🇧 Inglese: `docs/manual/`
- 🇮🇹 Italiano: `docs/manual_it/`

## Contributi

Contributi alla documentazione sono benvenuti:

- Correzioni
- Miglioramenti stilistici
- Esempi aggiuntivi
- Traduzioni in altre lingue

Inviare pull request su: https://github.com/manvalan/IOccultCalc

## Licenza

Il manuale è rilasciato sotto licenza MIT, come il software IOccultCalc.

Copyright © 2025 Michele Bigi - Gruppo Astrofili Massesi

## Citazione

Se utilizzi IOccultCalc nella tua ricerca, puoi citare:

```bibtex
@manual{ioccultcalc2025,
  title = {IOccultCalc: High-Precision Asteroid Occultation Prediction},
  author = {Michele Bigi},
  organization = {Gruppo Astrofili Massesi},
  year = {2025},
  url = {https://github.com/manvalan/IOccultCalc}
}
```

## Contatti

Per domande sul manuale o sul software:

- Email: mikbigi@gmail.com
- GitHub Issues: https://github.com/manvalan/IOccultCalc/issues
- Gruppo Astrofili Massesi: http://www.astrofilimassesi.it

---

*Ad astra per aspera*
