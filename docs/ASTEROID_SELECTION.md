# Selezione Asteroidi Specifici

## Nuova Funzionalità: Lista Esplicita di Asteroidi

ITALOccultCalc ora supporta la selezione di asteroidi specifici tramite tre modalità:

### 1. Singolo Numero (`asteroid_number`)

Per analizzare un solo asteroide:

```json
{
  "sections": [
    {
      "type": "object",
      "parameters": [
        {
          "name": "asteroid_number",
          "type": "int",
          "value": "10",
          "comment": "Solo Hygiea"
        }
      ]
    }
  ]
}
```

### 2. Lista Inline (`asteroid_list`)

Per analizzare più asteroidi specifici:

```json
{
  "sections": [
    {
      "type": "object",
      "parameters": [
        {
          "name": "asteroid_list",
          "type": "string",
          "value": "10, 4, 1, 433",
          "comment": "Hygiea, Vesta, Ceres, Eros"
        }
      ]
    }
  ]
}
```

### 3. File Esterno (`asteroid_list_file`)

Per liste lunghe, usa un file di testo:

```json
{
  "sections": [
    {
      "type": "object",
      "parameters": [
        {
          "name": "asteroid_list_file",
          "type": "string",
          "value": "my_asteroids.txt",
          "comment": "File con lista asteroidi"
        }
      ]
    }
  ]
}
```

**Formato file `my_asteroids.txt`:**
```
# Asteroidi per campagna osservativa gennaio 2026
10    # Hygiea - occultazione 9 gennaio
4     # Vesta
1     # Ceres
433   # Eros

# Asteroidi NEA interessanti
99942  # Apophis
```

- Una riga per numero
- Righe vuote ignorate
- Commenti con `#` (inline o su riga separata)
- Spazi bianchi ignorati

## Comportamento

Quando si usa una lista esplicita:
- I filtri `min_diameter`, `max_diameter`, `max_magnitude` vengono **bypassati**
- Vengono cercati **solo** gli asteroidi nella lista
- L'output mostra `★ LISTA ESPLICITA: N asteroidi [...]`

## Esempio Preset Completo

```json
{
  "sections": [
    {
      "type": "object",
      "parameters": [
        {
          "name": "asteroid_list",
          "type": "string",
          "value": "10",
          "comment": "Solo Hygiea per test validazione"
        }
      ]
    },
    {
      "type": "propag",
      "parameters": [
        {"name": "step_size", "type": "double", "value": "0.05"},
        {"name": "type", "type": "string", "value": "RK4"}
      ]
    },
    {
      "type": "ephemeris",
      "parameters": [
        {"name": "jpl_version", "type": "string", "value": "DE441"}
      ]
    },
    {
      "type": "search",
      "parameters": [
        {"name": "start_jd", "type": "double", "value": "2461047.5"},
        {"name": "end_jd", "type": "double", "value": "2461051.5"},
        {"name": "mag_limit", "type": "double", "value": "14.0"}
      ]
    },
    {
      "type": "perturbations",
      "parameters": [
        {"name": "planets", "type": "bool", "value": "true"},
        {"name": "asteroid_count", "type": "int", "value": "0"}
      ]
    },
    {
      "type": "output",
      "parameters": [
        {"name": "verbosity", "type": "int", "value": "2"}
      ]
    }
  ]
}
```

## Uso da Riga di Comando

```bash
./build/examples/italoccultcalc preset_hygiea_test.json -v
```

Output atteso:
```
================================================================
SELEZIONE ASTEROIDI CANDIDATI
================================================================

Criteri selezione:
  ★ LISTA ESPLICITA: 1 asteroidi [10]
  (filtri diametro/magnitudine ignorati)

✓ Trovati 1 asteroidi candidati REALI

Top 5 asteroidi per priorità:
  1. (10) Hygiea - Score: 8.94 ★★★
     Dimensione: 255 km, H=5.64
```
