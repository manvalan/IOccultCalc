IOccultCalc Scientific Manual - Overleaf Upload Instructions
=============================================================

File: manual.zip (64 KB, 23 files)

CONTENTS
--------
✓ main.tex - Main document (book class, 11pt)
✓ references.bib - 46 scientific references
✓ chapters/ - 16 chapter files (01-16)
✓ appendices/ - 3 appendix files (A, B, C)

UPLOAD TO OVERLEAF
------------------
1. Go to https://www.overleaf.com
2. Click "New Project" → "Upload Project"
3. Select manual.zip
4. Overleaf will automatically extract all files

COMPILATION
-----------
Overleaf uses pdfLaTeX by default - this is correct for this project.

Compile sequence (automatic in Overleaf):
1. pdflatex main.tex
2. bibtex main
3. pdflatex main.tex
4. pdflatex main.tex

First compilation may take 30-60 seconds.

REQUIRED PACKAGES
-----------------
All packages are available in Overleaf's TeX Live distribution:

✓ amsmath, amssymb, amsfonts - Math typesetting
✓ natbib - Bibliography (author-year style)
✓ hyperref - PDF hyperlinks
✓ tikz - Diagrams (15+ figures)
✓ algorithm, algorithmic - Pseudocode
✓ geometry - Page layout
✓ fancyhdr - Headers/footers
✓ booktabs - Professional tables
✓ xcolor - Colored links

No additional installation needed!

EXPECTED OUTPUT
---------------
✓ ~100-150 pages PDF
✓ 16 chapters with complete scientific content
✓ 3 appendices with reference data
✓ 46 bibliography entries
✓ 15+ TikZ diagrams
✓ 40+ tables
✓ Hyperlinked cross-references

STRUCTURE
---------
IOccultCalc Manual/
├── main.tex                    (main document)
├── references.bib              (bibliography)
├── chapters/
│   ├── 01_introduction.tex
│   ├── 02_coordinate_systems.tex
│   ├── 03_time_systems.tex
│   ├── 04_planetary_ephemerides.tex
│   ├── 05_orbital_mechanics.tex
│   ├── 06_numerical_integration.tex
│   ├── 07_perturbations.tex
│   ├── 08_relativistic_corrections.tex
│   ├── 09_precession_nutation.tex
│   ├── 10_stellar_astrometry.tex
│   ├── 11_orbit_determination.tex
│   ├── 12_asteroid_shape.tex
│   ├── 13_besselian_method.tex
│   ├── 14_uncertainty_propagation.tex
│   ├── 15_implementation.tex
│   └── 16_validation.tex
└── appendices/
    ├── appendix_a_constants.tex
    ├── appendix_b_algorithms.tex
    └── appendix_c_vsop_tables.tex

TROUBLESHOOTING
---------------
If compilation fails:

1. Check compilation logs in Overleaf (bottom panel)
2. Ensure pdfLaTeX is selected (Menu → Compiler)
3. Try "Recompile from scratch" option
4. Verify all 23 files were uploaded

Common issues:
- Missing references → Run bibtex (automatic in Overleaf)
- Undefined citations [?] → Compile multiple times
- TikZ errors → Usually resolve after full compilation

CUSTOMIZATION
-------------
Edit in Overleaf interface:
- Chapters: Click chapters/XX_name.tex
- Bibliography: Click references.bib
- Layout: Edit main.tex preamble

Changes save automatically.

DOWNLOAD PDF
------------
After compilation:
- Click "Download PDF" button (top right)
- Or click "Menu" → "Download" → "PDF"

FILE SIZE: 64 KB compressed, ~170 KB uncompressed (source)
PDF SIZE: ~1-2 MB (with diagrams and tables)

SHARING
-------
Overleaf sharing options:
1. "Share" button → Invite collaborators by email
2. Generate read-only link
3. Download source ZIP for offline work

LICENSE
-------
Manual content: CC-BY-SA 4.0
You can modify and redistribute with attribution.

SUPPORT
-------
GitHub: manvalan/IOccultCalc
Documentation: docs/manual/README.md

CITATION
--------
@manual{IOccultCalc2025,
  title = {IOccultCalc: High-Precision Asteroid Occultation Prediction},
  author = {{IOccultCalc Development Team}},
  year = {2025},
  note = {Version 2.0}
}

==============================================================
Manual created: November 21, 2025
IOccultCalc v2.0 - High-Precision Occultation Prediction
Target accuracy: ±0.5-1 km (10× better than Occult4)
==============================================================
