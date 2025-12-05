# Calcolo Astrometrico con Polinomi di Chebyshev
## Metodo per Determinazione Posizioni e Criteri di Selezione Occultazioni

**Data:** 30 novembre 2025  
**Scopo:** Documentazione matematica del processo astrometrico con Chebyshev

---

## 1. Calcola la Posizione Polinomiale ($\alpha_{\text{pol}}, \delta_{\text{pol}}$)

Devi prima usare i tuoi coefficienti di Chebyshev ($c_{ij}$ che nel modello astrometrico sono $a_{ij}$ per $\alpha$ e $b_{ij}$ per $\delta$) per trovare la posizione celeste calcolata ($\alpha_{\text{pol}}, \delta_{\text{pol}}$) dalla tua misura strumentale ($u_s, v_s$).

Siano $T_i(u)$ e $T_j(v)$ i Polinomi di Chebyshev di primo tipo.

$$\alpha_{\text{pol}} = \alpha_{\text{base}} + \sum_{i} \sum_{j} a_{ij} T_i(u_s) T_j(v_s)$$

$$\delta_{\text{pol}} = \delta_{\text{base}} + \sum_{i} \sum_{j} b_{ij} T_i(u_s) T_j(v_s)$$

### Spiegazione dei Termini

- **$\alpha_{\text{base}}, \delta_{\text{base}}$**: Coordinate di riferimento (centro campo)
- **$a_{ij}, b_{ij}$**: Coefficienti del modello astrometrico di Chebyshev
- **$u_s, v_s$**: Coordinate strumentali normalizzate ($[-1, +1]$)
- **$T_i(u), T_j(v)$**: Polinomi di Chebyshev di ordine $i, j$

### Definizione Polinomi di Chebyshev

$$T_0(x) = 1$$
$$T_1(x) = x$$
$$T_2(x) = 2x^2 - 1$$
$$T_3(x) = 4x^3 - 3x$$
$$T_{n+1}(x) = 2x T_n(x) - T_{n-1}(x)$$

---

## 2. Calcola la Distanza Angolare ($\Delta\theta$)

La distanza angolare tra la posizione calcolata ($\alpha_{\text{pol}}, \delta_{\text{pol}}$) e la posizione catalogata della stella ($\alpha_s, \delta_s$) è data dalla **Legge Sferica del Coseno** (o dalla formula per piccole distanze, che è più stabile in astrometria).

### Formula Completa (Legge Sferica del Coseno)

$$\cos(\Delta\theta) = \sin(\delta_s) \sin(\delta_{\text{pol}}) + \cos(\delta_s) \cos(\delta_{\text{pol}}) \cos(\alpha_s - \alpha_{\text{pol}})$$

### Formula Approssimata per Piccole Distanze

Se la distanza è piccola (meno di pochi minuti d'arco), si può usare l'**approssimazione euclidea**:

$$\Delta\theta \approx \sqrt{(\Delta\alpha \cos(\delta))^2 + (\Delta\delta)^2}$$

Dove:
- $\Delta\alpha = \alpha_s - \alpha_{\text{pol}}$
- $\Delta\delta = \delta_s - \delta_{\text{pol}}$
- $\delta$ = declinazione media per il fattore di correzione

### Quando Usare Quale Formula

| Distanza Angolare | Formula Raccomandata | Errore Relativo |
|-------------------|---------------------|-----------------|
| $< 1'$ | Approssimazione euclidea | $< 0.01\%$ |
| $1' - 10'$ | Legge sferica del coseno | Esatta |
| $> 10'$ | Legge sferica del coseno | Necessaria |

---

## 3. Criterio di Selezione (Closest Approach)

Per selezionare una candidata velocemente, devi confrontare la distanza calcolata $\Delta\theta$ con una **soglia di vicinanza critica** $\epsilon$. 

Questa soglia è determinata da:
- **Incertezze astrometriche** del modello di Chebyshev
- **Dimensioni angolari** dell'oggetto occultante
- **Errori di propagazione orbitale**

### Definizione della Soglia Critica

$\epsilon$ è tipicamente dell'ordine di:
- **Pochi secondi d'arco ($''$)** per astrometria di alta precisione
- **Decine di $''$** per survey rapidi o oggetti piccoli

### Criterio per Closest Approach (Occultazione)

**Se $\Delta\theta \le \epsilon$, la stella è una candidata valida.**

### Valori Tipici di $\epsilon$

| Tipo di Osservazione | $\epsilon$ Tipico | Motivazione |
|---------------------|------------------|-------------|
| **Occultazioni NEA** | $5'' - 15''$ | Oggetti piccoli, alta precisione richiesta |
| **Occultazioni MBA** | $10'' - 30''$ | Oggetti medi, bilancio precisione/efficienza |
| **Survey veloce** | $30'' - 60''$ | Screening preliminare, alta copertura |
| **Follow-up** | $2'' - 8''$ | Verifica finale, massima precisione |

---

## 4. Implementazione Pratica

### Algoritmo Completo

```cpp
// Passo 1: Calcolo posizione polinomiale
double alpha_pol = alpha_base;
double delta_pol = delta_base;

for (int i = 0; i <= max_order_u; i++) {
    for (int j = 0; j <= max_order_v; j++) {
        double Ti = chebyshev_T(i, u_s);
        double Tj = chebyshev_T(j, v_s);
        
        alpha_pol += a_coeff[i][j] * Ti * Tj;
        delta_pol += b_coeff[i][j] * Ti * Tj;
    }
}

// Passo 2: Calcolo distanza angolare
double delta_alpha = alpha_star - alpha_pol;
double delta_delta = delta_star - delta_pol;

// Scelta della formula
double angular_distance;
if (abs(delta_alpha) < 0.1 && abs(delta_delta) < 0.1) { // ~6 arcmin
    // Approssimazione euclidea per piccole distanze
    double cos_delta = cos((delta_star + delta_pol) / 2.0);
    angular_distance = sqrt(pow(delta_alpha * cos_delta, 2) + pow(delta_delta, 2));
} else {
    // Formula sferica completa
    double cos_dist = sin(delta_star) * sin(delta_pol) + 
                     cos(delta_star) * cos(delta_pol) * cos(delta_alpha);
    angular_distance = acos(cos_dist);
}

// Passo 3: Criterio di selezione
double epsilon = 15.0 * ARCSEC_TO_RAD; // 15 arcsec in radianti
bool is_candidate = (angular_distance <= epsilon);
```

### Funzione Polinomi di Chebyshev

```cpp
double chebyshev_T(int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    
    double T_prev2 = 1.0;  // T0
    double T_prev1 = x;    // T1
    double T_current;
    
    for (int i = 2; i <= n; i++) {
        T_current = 2.0 * x * T_prev1 - T_prev2;
        T_prev2 = T_prev1;
        T_prev1 = T_current;
    }
    
    return T_current;
}
```

---

## 5. Considerazioni di Accuratezza

### Errori del Modello di Chebyshev

La precisione finale dipende da:

1. **Ordine del polinomio** ($N \times M$)
2. **Qualità del fitting** astrometrico
3. **Distribuzione delle stelle di riferimento**
4. **Stabilità numerica** dei coefficienti

### Errore Tipico vs Ordine Polinomio

| Ordine | Errore RMS Tipico | Note |
|--------|------------------|------|
| $3 \times 3$ | $0.1'' - 0.3''$ | Campi piccoli, poche distorsioni |
| $5 \times 5$ | $0.05'' - 0.15''$ | Campi medi, distorsioni moderate |
| $7 \times 7$ | $0.02'' - 0.08''$ | Campi grandi, correzione ottica |

### Propagazione Errori

L'errore totale nella selezione candidati è:

$$\sigma_{\text{totale}} = \sqrt{\sigma_{\text{astrom}}^2 + \sigma_{\text{orbit}}^2 + \sigma_{\text{catalog}}^2}$$

Dove:
- $\sigma_{\text{astrom}}$: errore del modello astrometrico
- $\sigma_{\text{orbit}}$: errore di propagazione orbitale  
- $\sigma_{\text{catalog}}$: errore posizione stella nel catalogo

---

## 6. Ottimizzazioni Computazionali

### Pre-calcolo Polinomi

Per efficienza in loop su molte stelle:

```cpp
// Pre-calcola tutti i polinomi T_i(u_s), T_j(v_s)
vector<double> T_u(max_order_u + 1);
vector<double> T_v(max_order_v + 1);

T_u[0] = 1.0; T_v[0] = 1.0;
T_u[1] = u_s; T_v[1] = v_s;

for (int i = 2; i <= max_order_u; i++) {
    T_u[i] = 2.0 * u_s * T_u[i-1] - T_u[i-2];
}
for (int j = 2; j <= max_order_v; j++) {
    T_v[j] = 2.0 * v_s * T_v[j-1] - T_v[j-2];
}

// Poi usa T_u[i] e T_v[j] nel doppio loop
```

### Screening Gerarchico

1. **Primo livello**: controllo grossolano con $\epsilon_{\text{large}}$
2. **Secondo livello**: calcolo preciso solo per candidati promettenti
3. **Terzo livello**: verifica finale con geometria completa

---

*Fine documentazione calcolo astrometrico Chebyshev*