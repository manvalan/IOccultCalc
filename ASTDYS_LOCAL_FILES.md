# Configurazione File AstDyS Locali

## 📋 Panoramica

IOccultCalc può utilizzare file `.eq1` (elementi orbitali) e `.rwo` (osservazioni) scaricati localmente da AstDyS, invece di scaricarli ogni volta da internet.

## 🎯 Vantaggi

✅ **Performance**: Nessun download durante l'esecuzione  
✅ **Offline**: Lavora senza connessione internet  
✅ **Batch**: Elabora centinaia di asteroidi senza rate limiting  
✅ **Riproducibilità**: Usa sempre gli stessi dati orbitali  

## 📥 Workflow Completo

### Step 1: Download File da AstDyS

Usa lo script Python fornito:

```bash
# Download lista asteroidi
python3 download_astdys_data.py asteroids.txt --output-dir astdys_data

# Opzioni utili:
python3 download_astdys_data.py asteroids.txt \
    -o astdys_data \
    --only-eq1 \              # Solo elementi (più veloce)
    --delay 1.0 \             # Ritardo tra richieste
    --force                   # Ri-download anche se esistono
```

**Output**:
```
astdys_data/
├── 433.eq1          # Elementi orbitali
├── 433.rwo          # Osservazioni
├── 11234.eq1
├── 11234.rwo
└── ...
```

### Step 2: Configura Preset IOccultCalc

Aggiungi sezione `astdys` nel file preset `.oop`:

```fortran
astdys.
    # Directory con file .eq1 (elementi orbitali)
    .local_eq1_directory = 'astdys_data'
    
    # Directory con file .rwo (osservazioni) - opzionale
    .local_rwo_directory = 'astdys_data'
```

**Esempio completo** (`preset_local_astdys.oop`):

```fortran
general.
    .propagator = 'RKF78'
    .step_size_days = 0.1

time.
    .start_date = '2026-01-01'
    .end_date = '2026-01-31'

astdys.
    .local_eq1_directory = 'astdys_data'
    .local_rwo_directory = 'astdys_data'

asteroids.
    .asteroid_list = '433,11234,203,704,10'

gaia.
    .use_local_cache = .TRUE.
    .cache_directory = '/Users/username/catalogs'

output.
    .format = 'iota'
    .directory = 'output'
```

### Step 3: Esegui IOccultCalc

```bash
./build/examples/italoccultcalc preset_local_astdys.oop
```

**Output atteso**:
```
📂 Configurazione directory AstDyS locali:
   File .eq1: astdys_data
   File .rwo: astdys_data

📁 Carico elementi da file locale: astdys_data/433.eq1
  ✓ (433) Eros - elementi caricati (epoca MJD 60671.0)
📁 Carico elementi da file locale: astdys_data/11234.eq1
  ✓ (11234) 1999 JS82 - elementi caricati
...
```

## 🔄 Workflow Automatico

Usa lo script bash per tutto il processo:

```bash
# 1. Download + Conversione preset
./download_and_convert.sh asteroids.txt astdys_data presets

# 2. Modifica preset generati per usare file locali
# (i preset generati hanno già .local_eq1_directory configurato!)

# 3. Esegui batch
for preset in presets/preset_*.oop; do
    ./build/examples/italoccultcalc "$preset"
done
```

## ⚙️ Opzioni Configurazione

### Directory Relative vs Assolute

```fortran
# Relativa (rispetto a dove esegui italoccultcalc)
.local_eq1_directory = 'astdys_data'

# Assoluta (sempre valida)
.local_eq1_directory = '/Users/username/astro/astdys_data'
```

### Solo File .eq1 (no osservazioni)

```fortran
astdys.
    # Solo elementi orbitali (sufficiente per propagazione)
    .local_eq1_directory = 'astdys_data'
    
    # NON specificare .local_rwo_directory
```

### Fallback ad AstDyS Online

Se file locale non trovato, IOccultCalc scarica automaticamente da internet:

```
⚠️  File locale non trovato: astdys_data/433.eq1
   Provo download da AstDyS...
🌐 Download da AstDyS: https://newton.spacedys.com/~astdys2/epoch/numbered/0/433.eq1
  ✓ (433) Eros - elementi scaricati
```

### Disabilita Directory Locali

```fortran
astdys.
    # Vuoto = usa sempre download HTTP
    .local_eq1_directory = ''
    .local_rwo_directory = ''
```

## 📂 Struttura File Consigliata

```
progetto/
├── italoccultcalc                  # Eseguibile
├── asteroids_jan2026.txt           # Lista asteroidi
│
├── astdys_data/                    # Dati scaricati
│   ├── 433.eq1
│   ├── 433.rwo
│   ├── 11234.eq1
│   ├── 11234.rwo
│   └── ...
│
├── presets/                        # Preset configurati
│   ├── preset_433.oop              # Con astdys.local_eq1_directory
│   ├── preset_11234.oop
│   └── all_presets.txt
│
└── output/                         # Risultati
    ├── occ_433.txt
    ├── occ_11234.txt
    └── ...
```

## 🔍 Verifica Configurazione

```bash
# Verifica file presenti
ls -lh astdys_data/*.eq1

# Conta file scaricati
ls astdys_data/*.eq1 | wc -l

# Test con singolo asteroide
cat > test_local.oop << EOF
general.
    .propagator = 'RKF78'
time.
    .start_date = '2026-01-01'
    .end_date = '2026-01-15'
astdys.
    .local_eq1_directory = 'astdys_data'
asteroids.
    .asteroid_list = '433'
gaia.
    .use_local_cache = .TRUE.
    .cache_directory = '$HOME/catalogs'
output.
    .format = 'iota'
EOF

./build/examples/italoccultcalc test_local.oop
```

## 📊 Performance Comparison

| Modalità | Tempo Setup | Tempo/Asteroide | Offline | Note |
|----------|-------------|-----------------|---------|------|
| **Download HTTP** | 0s | ~2-5s | ❌ | Richiede internet |
| **File Locali** | 1-10min (download batch) | ~0.1s | ✅ | Veloce dopo setup |

**Esempio**: 100 asteroidi
- HTTP: ~200-500 secondi (3-8 minuti)
- Locali: ~10 secondi dopo download iniziale

## 🚨 Troubleshooting

### File non trovato

```
⚠️  File locale non trovato: astdys_data/433.eq1
```

**Soluzione**:
```bash
# Verifica path corretto
ls astdys_data/433.eq1

# Re-download se necessario
python3 download_astdys_data.py <(echo "433") -o astdys_data --force
```

### Directory non valida

```
Error opening file: astdys_data/433.eq1
```

**Cause**:
1. Path relativo sbagliato (esegui da directory diversa)
2. Permessi file
3. Directory non esiste

**Soluzione**:
```bash
# Usa path assoluto
pwd  # Verifica dove sei
ls -la astdys_data/  # Verifica permessi

# Nel preset:
.local_eq1_directory = '/full/path/to/astdys_data'
```

### File corrotto

```
Parse error: invalid .eq1 format
```

**Soluzione**:
```bash
# Re-download con --force
python3 download_astdys_data.py <(echo "433") -o astdys_data --force

# Verifica contenuto file
head -20 astdys_data/433.eq1
# Deve contenere "EQU" e "OEF2.0"
```

## 🔧 Gestione File

### Aggiornamento Elementi Orbitali

Gli elementi orbitali di AstDyS vengono aggiornati periodicamente:

```bash
# Ri-download elementi aggiornati (sovrascrive)
python3 download_astdys_data.py asteroids.txt -o astdys_data --force

# Backup elementi vecchi prima di aggiornare
cp -r astdys_data astdys_data_backup_$(date +%Y%m%d)
python3 download_astdys_data.py asteroids.txt -o astdys_data --force
```

### Pulizia File Vecchi

```bash
# Rimuovi file più vecchi di 30 giorni
find astdys_data -name "*.eq1" -mtime +30 -delete
find astdys_data -name "*.rwo" -mtime +30 -delete

# Rimuovi solo file .rwo (osservazioni, più grandi)
rm astdys_data/*.rwo
```

### Spazio Disco

```bash
# Verifica spazio usato
du -sh astdys_data

# File tipici:
# .eq1: ~1-2 KB ciascuno
# .rwo: ~100 KB - 5 MB ciascuno

# 100 asteroidi:
#   Solo .eq1: ~200 KB
#   Con .rwo: ~50-500 MB
```

## 📝 Best Practices

### 1. Organizzazione File

```bash
# Directory separate per dataset diversi
astdys_data/
├── monthly/        # Dataset mensile routine
├── special/        # Eventi speciali
└── archive/        # Backup storico
```

### 2. Version Control

```bash
# Traccia versione elementi
echo "$(date +%Y-%m-%d)" > astdys_data/download_date.txt

# Nel preset:
# Commento con data download
# Dataset scaricato il 2025-11-29
astdys.
    .local_eq1_directory = 'astdys_data/monthly'
```

### 3. Automatizzazione

```bash
# Script cron per aggiornamento mensile
# crontab -e
0 2 1 * * /path/to/update_astdys.sh

# update_astdys.sh:
#!/bin/bash
DATE=$(date +%Y-%m)
python3 download_astdys_data.py asteroids_routine.txt \
    -o "astdys_data/update_$DATE" --force
```

## 🔗 Riferimenti

- Script download: `download_astdys_data.py`
- Script workflow completo: `download_and_convert.sh`
- Conversione preset: `eq1_to_preset.py`
- Guida download: `DOWNLOAD_ASTDYS_GUIDE.md`
- Formato .eq1: `ESEMPIO_IMPORTAZIONE_ASTDYN.md`

## ✅ Checklist Setup

- [ ] Script `download_astdys_data.py` funzionante
- [ ] Directory `astdys_data` creata
- [ ] File `.eq1` scaricati per tutti gli asteroidi
- [ ] Preset `.oop` con sezione `astdys` configurata
- [ ] Test con singolo asteroide riuscito
- [ ] Batch processing testato
- [ ] Backup dati creato

---

**Documentazione**: IOccultCalc - File AstDyS Locali  
**Versione**: 1.0  
**Data**: 2025-11-29
