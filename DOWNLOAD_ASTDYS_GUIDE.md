# Guida: Download Automatico Dati AstDyS

## 📚 Overview

Sistema completo per scaricare automaticamente file `.eq1` (elementi orbitali) e `.rwo` (osservazioni) da AstDyS e convertirli in preset IOccultCalc.

## 🛠️ Strumenti Disponibili

### 1. download_astdys_data.py ⭐

Script Python per download automatico da AstDyS.

**Caratteristiche**:
- ✅ Download `.eq1` (elementi orbitali equinoziali)
- ✅ Download `.rwo` (osservazioni ottiche MPC)
- ✅ Gestione errori SSL e HTTP
- ✅ Skip file già presenti
- ✅ Ritardo configurabile tra richieste
- ✅ Genera summary per ogni asteroide
- ✅ Lista file scaricati con successo

**Uso**:
```bash
# Download da lista asteroidi
python3 download_astdys_data.py asteroids.txt

# Output in directory specifica
python3 download_astdys_data.py asteroids.txt -o astdys_data

# Solo file .eq1 (elementi)
python3 download_astdys_data.py asteroids.txt --only-eq1

# Solo file .rwo (osservazioni)
python3 download_astdys_data.py asteroids.txt --only-rwo

# Forza re-download
python3 download_astdys_data.py asteroids.txt --force

# Ritardo tra download (rispetta server)
python3 download_astdys_data.py asteroids.txt --delay 2.0

# Test con pochi asteroidi
python3 download_astdys_data.py asteroids.txt --max-asteroids 5
```

### 2. eq1_to_preset.py

Script Python per conversione `.eq1` → preset `.oop`.

**Uso**:
```bash
# Conversione singola
python3 eq1_to_preset.py 433.eq1 -o preset_433.oop

# Con periodo custom
python3 eq1_to_preset.py 433.eq1 -p 2026-01-01 2026-12-31

# Conversione batch
for f in astdys_data/*.eq1; do
    num=$(basename "$f" .eq1)
    python3 eq1_to_preset.py "$f" -o "preset_${num}.oop"
done
```

### 3. download_and_convert.sh ⭐

Script bash che automatizza tutto il workflow.

**Uso**:
```bash
# Workflow completo
./download_and_convert.sh asteroids_list.txt

# Con directory custom
./download_and_convert.sh asteroids_list.txt astdys_data presets
```

## 📋 Formati Input Supportati

### Lista Semplice

```
# Un numero per riga
433
11234
203
704
```

### Output IOccultCalc

```
Occultation predictions for January 2026:
(433) Eros occults star ...
(11234) 1999 JS82 occults ...
```

Lo script estrae automaticamente i numeri tra parentesi.

### File Misto

```
# Commenti ignorati
1  # Ceres
2  # Pallas
(10) Hygiea  # Formato IOccultCalc
433          # Solo numero
```

## 🌐 URL AstDyS

### Struttura Path

Gli asteroidi sono organizzati in subdirectory basate su `numero // 1000`:

```
Asteroide    Subdir    Path
─────────────────────────────────────────────────────
1            0         .../0/1.eq1
433          0         .../0/433.eq1
11234        11        .../11/11234.eq1
```

### Endpoint

**Elementi orbitali (.eq1)**:
```
https://newton.spacedys.com/~astdys2/epoch/numbered/XX/NUM.eq1
```

**Osservazioni (.rwo)**:
```
https://newton.spacedys.com/~astdys2/mpcobs/numbered/XX/NUM.rwo
```

Dove:
- `XX` = `NUM // 1000`
- `NUM` = numero asteroide

**Esempi**:
```bash
# (433) Eros
https://newton.spacedys.com/~astdys2/epoch/numbered/0/433.eq1
https://newton.spacedys.com/~astdys2/mpcobs/numbered/0/433.rwo

# (11234) 1999 JS82
https://newton.spacedys.com/~astdys2/epoch/numbered/11/11234.eq1
https://newton.spacedys.com/~astdys2/mpcobs/numbered/11/11234.rwo
```

## 🚀 Workflow Completo

### Step 1: Prepara Lista Asteroidi

```bash
# Crea lista asteroidi
cat > my_asteroids.txt << EOF
433   # Eros
11234 # 1999 JS82
203   # Pompeja
10    # Hygiea
EOF
```

### Step 2: Download da AstDyS

```bash
# Download automatico
python3 download_astdys_data.py my_asteroids.txt \
    --output-dir astdys_data \
    --delay 1.0
```

**Output**:
```
astdys_data/
├── 433.eq1
├── 433.rwo
├── 433_summary.txt
├── 11234.eq1
├── 11234.rwo
├── 11234_summary.txt
├── asteroids_with_eq1.txt
...
```

### Step 3: Conversione Preset

#### Opzione A: Manuale (controllo totale)

```bash
# Converti singolo asteroide
python3 eq1_to_preset.py astdys_data/433.eq1 \
    --output presets/preset_433.oop \
    --period 2026-01-01 2026-12-31
```

#### Opzione B: Batch Loop

```bash
# Converti tutti
mkdir -p presets
for eq1 in astdys_data/*.eq1; do
    num=$(basename "$eq1" .eq1)
    python3 eq1_to_preset.py "$eq1" -o "presets/preset_${num}.oop"
done
```

#### Opzione C: Script Automatico ⭐ RACCOMANDATO

```bash
# Un solo comando fa tutto
./download_and_convert.sh my_asteroids.txt astdys_data presets
```

### Step 4: Esegui Calcoli

```bash
# Singolo asteroide
./italoccultcalc presets/preset_433.oop > results_433.txt

# Batch tutti gli asteroidi
for preset in presets/preset_*.oop; do
    num=$(basename "$preset" .oop | sed 's/preset_//')
    echo "Calculating occultations for ($num)..."
    ./italoccultcalc "$preset" > "results_${num}.txt"
done
```

## 📊 Esempio Completo

```bash
# 1. Crea lista asteroidi interessanti
cat > priority_asteroids.txt << EOF
# Asteroidi grandi con buone osservazioni
1     # Ceres - 148,000+ obs
2     # Pallas - 127,000+ obs
4     # Vesta - 183,000+ obs
10    # Hygiea - 98,000+ obs
433   # Eros - 48,000+ obs
704   # Interamnia - 47,000+ obs
EOF

# 2. Download e conversione automatica
./download_and_convert.sh priority_asteroids.txt astdys_priority presets_priority

# 3. Verifica preset generati
ls -lh presets_priority/

# 4. Esegui calcolo per uno
./italoccultcalc presets_priority/preset_433.oop

# 5. Batch per tutti (opzionale)
mkdir -p results_priority
for preset in presets_priority/preset_*.oop; do
    num=$(basename "$preset" | sed 's/preset_//;s/.oop//')
    echo "Processing ($num)..."
    ./italoccultcalc "$preset" > "results_priority/occ_${num}.txt" 2>&1
done

# 6. Conta occultazioni trovate
grep -c "Occultation" results_priority/*.txt
```

## 🔧 Opzioni Avanzate

### Gestione Rate Limiting

Per liste grandi (>100 asteroidi), usa ritardo maggiore:

```bash
python3 download_astdys_data.py large_list.txt \
    --delay 2.0 \
    --verbose
```

### Download Parziale per Test

```bash
# Solo primi 10 asteroidi
python3 download_astdys_data.py asteroids.txt \
    --max-asteroids 10 \
    --output-dir test_download
```

### Re-download Forzato

```bash
# Scarica anche se file esiste
python3 download_astdys_data.py asteroids.txt --force
```

### Solo Elementi (no osservazioni)

```bash
# Più veloce, solo per conversione preset
python3 download_astdys_data.py asteroids.txt --only-eq1
```

## 📁 Struttura Output

```
progetto/
├── asteroids_list.txt          # Input
├── download_astdys_data.py     # Script download
├── eq1_to_preset.py            # Script conversione
├── download_and_convert.sh     # Workflow automatico
│
├── astdys_data/                # Dati scaricati
│   ├── 433.eq1                 # Elementi orbitali
│   ├── 433.rwo                 # Osservazioni
│   ├── 433_summary.txt         # Info file
│   ├── asteroids_with_eq1.txt  # Lista successi
│   └── ...
│
├── presets/                    # Preset generati
│   ├── preset_433.oop
│   ├── preset_11234.oop
│   ├── all_presets.txt         # Lista preset
│   └── ...
│
└── results/                    # Risultati calcoli
    ├── occ_433.txt
    ├── occ_11234.txt
    └── ...
```

## ⚠️ Limitazioni e Note

### Asteroidi Non Numerati

Lo script attualmente supporta solo asteroidi **numerati** (con numero definitivo).

Per asteroidi con designazione provvisoria (es: 1999 JS82 prima della numerazione):
- Non supportato dal path AstDyS
- Usare numero definitivo se disponibile
- Oppure download manuale

### Rate Limiting

AstDyS non ha limiti pubblici, ma:
- Usa `--delay 1.0` di default (1 secondo tra richieste)
- Per liste grandi usa `--delay 2.0` o più
- Rispetta il server!

### File Non Disponibili

Alcuni asteroidi potrebbero non avere:
- **File .eq1 mancante**: Elementi non aggiornati su AstDyS
  - Usa elementi da MPC o JPL Horizons
- **File .rwo vuoto**: Poche osservazioni
  - Normale per asteroidi scoperti recentemente
  - Elementi .eq1 comunque validi

### Certificati SSL

Lo script disabilita verifica SSL per compatibilità.  
Se preferisci verificare certificati:
1. Installa certificati Python: `pip install certifi`
2. Rimuovi `ssl_context` dallo script

## 🎯 Best Practices

### 1. Organizzazione File

```bash
# Struttura consigliata
mkdir -p {data,presets,results}/{monthly,special}

# Download mensile routine
./download_and_convert.sh monthly_asteroids.txt \
    data/monthly presets/monthly

# Eventi speciali
./download_and_convert.sh special_events.txt \
    data/special presets/special
```

### 2. Backup Dati

```bash
# Backup periodico dati AstDyS
tar -czf astdys_backup_$(date +%Y%m%d).tar.gz astdys_data/

# Backup preset configurati
tar -czf presets_backup_$(date +%Y%m%d).tar.gz presets/
```

### 3. Verifica Qualità

```bash
# Conta osservazioni per asteroide
for rwo in astdys_data/*.rwo; do
    num=$(basename "$rwo" .rwo)
    obs=$(grep -cv "^!" "$rwo" | tail -1)
    echo "($num): $obs observations"
done | sort -t: -k2 -n -r

# Verifica epoca elementi
for eq1 in astdys_data/*.eq1; do
    num=$(basename "$eq1" .eq1)
    epoch=$(grep "MJD" "$eq1" | awk '{print $2}')
    echo "($num): epoch MJD $epoch"
done
```

### 4. Aggiornamento Periodico

```bash
# Script per aggiornamento mensile
#!/bin/bash
DATE=$(date +%Y-%m)
./download_and_convert.sh asteroids_routine.txt \
    "data/update_$DATE" \
    "presets/update_$DATE"
```

## 📚 Riferimenti

- **AstDyS**: https://newton.spacedys.com/~astdys2/
- **Formato OEF2.0**: Vedi `ESEMPIO_IMPORTAZIONE_ASTDYN.md`
- **MPC Observatory Codes**: https://minorplanetcenter.net/iau/lists/ObsCodes.html

## 🆘 Troubleshooting

### Errore SSL

```
✗ Errore URL: [SSL: CERTIFICATE_VERIFY_FAILED]
```

**Soluzione**: Lo script già gestisce questo (usa `ssl_context`)

### File 404

```
✗ Non trovato (404)
```

**Cause**:
- Asteroide non numerato
- Numero errato
- File non disponibile su AstDyS

**Soluzione**: Verifica numero su https://minorplanetcenter.net

### Download Lento

**Soluzione**: Usa `--only-eq1` se non servono osservazioni

### File Corrotto

```
✗ File invalido (no EQU)
```

**Soluzione**: 
```bash
# Re-download con --force
python3 download_astdys_data.py asteroids.txt --force
```

---

**Autore**: IOccultCalc Development Team  
**Data**: 2025-11-29  
**Versione**: 1.0  
**Testato con**: Python 3.8+, macOS/Linux
