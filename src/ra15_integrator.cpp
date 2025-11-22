/**
 * @file ra15_integrator.cpp
 * @brief Implementazione integratore RA15 (Radau 15° ordine di Everhart)
 * 
 * Basato su OrbFit 5.0.8 (Andrea Milani et al., Università di Pisa)
 * File originale: src/propag/ra15_mod.f90
 * 
 * Riferimenti:
 * - Everhart, E. (1985). "An efficient integrator that uses Gauss-Radau spacings".
 *   IAU Colloq. 83: Dynamics of Comets, pp. 185-202.
 * - Milani, A., Gronchi, G.F. (2010). "Theory of Orbit Determination".
 *   Cambridge University Press, Chapter 7.
 * - OrbFit 5.0.8: https://adams.dm.unipi.it/orbfit/
 */

#include "ioccultcalc/ra15_integrator.hpp"
#include "ioccultcalc/orbit_propagator.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace ioccultcalc {

// Costanti globali per RA15 (come in OrbFit MODULE ra15_mod)
static double w1_global = 0.0;  // Usato in rapred

// Costruttore
RA15Integrator::RA15Integrator(const RA15Options& options)
    : options_(options)
{
    computeRadauConstants();
    resetArrays();
}

// Calcola costanti Gauss-Radau - TRADUZIONE ESATTA da OrbFit ra15_mod.f90
void RA15Integrator::computeRadauConstants() {
    // Spaziature punti Gauss-Radau (h) - ESATTAMENTE da OrbFit linee 717-724
    h_[0] = 0.0;
    h_[1] = 0.05626256053692215;
    h_[2] = 0.18024069173689236;
    h_[3] = 0.35262471711316964;
    h_[4] = 0.54715362633055538;
    h_[5] = 0.73421017721541053;
    h_[6] = 0.88532094683909577;
    h_[7] = 0.97752061356128750;
    
    // Calcola w e u - OrbFit linee 728-733
    // ncl = false per equazioni secondo ordine y" = f(y,t)
    bool ncl = false;
    
    for (int n = 2; n <= 8; ++n) {
        double ww = n + n * n;  // ww = n + n²
        if (ncl) ww = n;
        w_[n - 2] = 1.0 / ww;  // w array indicizzato da 0
        
        ww = n;
        u_[n - 2] = 1.0 / ww;
    }
    
    // w1 per rapred (OrbFit linea 734)
    w1_global = 0.5;  // half per secondo ordine
    if (ncl) w1_global = 1.0;
    
    // Array nw per indicizzazione (OrbFit linea 712)
    int nw[8] = {0, 0, 1, 3, 6, 10, 15, 21};
    
    // Inizializza primi elementi c, d, r (OrbFit linee 735-737)
    c_[0] = -h_[1];  // c(1) = -h(2) in Fortran (1-indexed)
    d_[0] = h_[1];   // d(1) = h(2)
    r_[0] = 1.0 / (h_[2] - h_[1]);  // r(1) = 1/(h(3)-h(2))
    
    // Loop principale (OrbFit linee 738-758)
    int la = 0;  // In C++ array 0-indexed
    int lc = 0;
    
    for (int k = 3; k <= 7; ++k) {
        int lb = la;
        la = lc + 1;
        lc = nw[k];
        
        c_[la] = -h_[k - 1] * c_[lb];  // h_ è 0-indexed
        c_[lc] = c_[la - 1] - h_[k - 1];
        d_[la] = h_[1] * d_[lb];
        d_[lc] = -c_[lc];
        r_[la] = 1.0 / (h_[k] - h_[1]);
        r_[lc] = 1.0 / (h_[k] - h_[k - 1]);
        
        if (k == 3) continue;  // goto 73 in Fortran
        
        for (int l = 4; l <= k; ++l) {
            int ld = la + l - 3;
            int le = lb + l - 4;
            c_[ld] = c_[le] - h_[k - 1] * c_[le + 1];
            d_[ld] = d_[le] + h_[l - 2] * d_[le + 1];
            r_[ld] = 1.0 / (h_[k] - h_[l - 2]);
        }
    }
}

// Reset storage arrays
void RA15Integrator::resetArrays() {
    for (int l = 0; l < 7; ++l) {
        for (int k = 0; k < 6; ++k) {
            b_[l][k] = 0.0;
            g_[l][k] = 0.0;
            e_[l][k] = 0.0;
            bd_[l][k] = 0.0;
        }
    }
}

// Calcola G da B
void RA15Integrator::computeG(
    const std::array<std::array<double, 6>, 7>& b,
    std::array<std::array<double, 6>, 7>& g)
{
    // Formula: g = d * b
    // d è la matrice dei coefficienti di trasformazione
    
    for (int k = 0; k < 6; ++k) {
        g[0][k] = b[0][k] * d_[0];
        
        int idx = 1;
        for (int j = 1; j < 7; ++j) {
            double sum = 0.0;
            for (int l = 0; l <= j; ++l) {
                sum += b[l][k] * d_[idx + l];
            }
            g[j][k] = sum;
            idx += j + 1;
        }
    }
}

// Iterazione metodo implicito - TRADUZIONE ESATTA da OrbFit SUBROUTINE rasust
double RA15Integrator::iterate(
    int iter,
    double t,
    double tm,
    const JulianDate& tini,
    std::array<double, 6>& x_and_v,
    const std::array<double, 6>& f1)
{
    double t2 = t * t;  // Per equazioni secondo ordine
    bool ncl = false;   // false = equazioni secondo ordine y"=f(y,t)
    bool npq = false;   // ← FIX: false per aggiornare anche velocità durante iterazioni
    
    // w1 per equazioni secondo ordine (OrbFit linea 734)
    double w1 = 0.5;  // half
    
    // Storage ck per convergenza (OrbFit linee 452-459)
    static std::array<std::array<double, 6>, 8> ck;  // ck(nv,8)
    
    // Inizializza ck a zero alla prima iterazione
    if (iter == 1) {
        for (int k = 0; k < 6; ++k) {
            for (int j = 0; j < 8; ++j) {
                ck[j][k] = 0.0;
            }
        }
    }
    
    double epsi = 0.0;
    
    // Loop substep j=2..8 (OrbFit linea 463)
    for (int j = 2; j <= 8; ++j) {
        double s = h_[j - 1];  // h_ è 0-indexed, j è 1-indexed in Fortran
        double q = s;
        if (ncl) q = 1.0;
        
        // Predizione posizione e velocità (OrbFit linee 474-488)
        std::array<double, 6> y, z;
        
        for (int k = 0; k < 3; ++k) {  // Solo posizioni (3 componenti)
            // Formula collapsed series eq. (2.9) e (2.10)
            double temp = w_[2] * b_[2][k] + s * (w_[3] * b_[3][k] + 
                         s * (w_[4] * b_[4][k] + s * (w_[5] * b_[5][k] + 
                         s * w_[6] * b_[6][k])));
            
            y[k] = x_and_v[k] + q * (t * x_and_v[k + 3] + 
                   t2 * s * (f1[k] * w1 + s * (w_[0] * b_[0][k] + 
                   s * (w_[1] * b_[1][k] + s * temp))));
            
            // Velocità predictors (per equazioni secondo ordine)
            if (!npq) {
                temp = u_[2] * b_[2][k] + s * (u_[3] * b_[3][k] + 
                       s * (u_[4] * b_[4][k] + s * (u_[5] * b_[5][k] + 
                       s * u_[6] * b_[6][k])));
                
                z[k] = x_and_v[k + 3] + s * t * (f1[k] + s * (u_[0] * b_[0][k] + 
                       s * (u_[1] * b_[1][k] + s * temp)));
            } else {
                z[k] = x_and_v[k + 3];  // Per nclass=-2 usa velocità corrente
            }
        }
        
        // Valuta forza al substep (OrbFit linee 495-504)
        Vector3D pos_j(y[0], y[1], y[2]);
        Vector3D vel_j(z[0], z[1], z[2]);
        JulianDate t_j(tini.jd + tm + s * t);
        
        Vector3D accel_j = accel_func_(t_j, pos_j, vel_j);
        stats_.num_evals++;
        
        std::array<double, 3> fj = {accel_j.x, accel_j.y, accel_j.z};
        
        // Controllo convergenza (OrbFit linee 505-510)
        int nvv = 3;  // Solo posizioni per convergenza
        for (int k = 0; k < nvv; ++k) {
            epsi += std::abs(fj[k] - ck[j - 1][k]);
            ck[j - 1][k] = fj[k];
        }
        
        // Loop su componenti (OrbFit linea 513)
        for (int k = 0; k < 3; ++k) {
            // Calcola g-value (OrbFit linee 518-521)
            double temp = g_[j - 2][k];  // g(j-1,k) in Fortran
            double g_k = (fj[k] - f1[k]) / s;
            
            // Ricorrenza per calcolare g (OrbFit linee 522-558)
            // Ogni caso j ha formula specifica
            if (j <= 2) {
                g_[0][k] = g_k;
            } else if (j == 3) {
                g_[1][k] = (g_k - g_[0][k]) * r_[0];
            } else if (j == 4) {
                g_[2][k] = ((g_k - g_[0][k]) * r_[1] - g_[1][k]) * r_[2];
            } else if (j == 5) {
                g_[3][k] = (((g_k - g_[0][k]) * r_[3] - g_[1][k]) * r_[4] - g_[2][k]) * r_[5];
            } else if (j == 6) {
                g_[4][k] = ((((g_k - g_[0][k]) * r_[6] - g_[1][k]) * r_[7] - g_[2][k]) * r_[8] - g_[3][k]) * r_[9];
            } else if (j == 7) {
                g_[5][k] = (((((g_k - g_[0][k]) * r_[10] - g_[1][k]) * r_[11] - g_[2][k]) * r_[12] - g_[3][k]) * r_[13] - g_[4][k]) * r_[14];
            } else if (j == 8) {
                g_[6][k] = ((((((g_k - g_[0][k]) * r_[15] - g_[1][k]) * r_[16] - g_[2][k]) * r_[17] - g_[3][k]) * r_[18] - g_[4][k]) * r_[19] - g_[5][k]) * r_[20];
            }
            
            // Aggiorna b-values (OrbFit linee 562-564)
            temp = g_[j - 2][k] - temp;
            b_[j - 2][k] = b_[j - 2][k] + temp;
            
            // Propaga aggiornamento agli altri b usando matrice c (OrbFit linee 569-600)
            if (j == 3) {
                b_[0][k] = b_[0][k] + c_[0] * temp;
            } else if (j == 4) {
                b_[0][k] = b_[0][k] + c_[1] * temp;
                b_[1][k] = b_[1][k] + c_[2] * temp;
            } else if (j == 5) {
                b_[0][k] = b_[0][k] + c_[3] * temp;
                b_[1][k] = b_[1][k] + c_[4] * temp;
                b_[2][k] = b_[2][k] + c_[5] * temp;
            } else if (j == 6) {
                b_[0][k] = b_[0][k] + c_[6] * temp;
                b_[1][k] = b_[1][k] + c_[7] * temp;
                b_[2][k] = b_[2][k] + c_[8] * temp;
                b_[3][k] = b_[3][k] + c_[9] * temp;
            } else if (j == 7) {
                b_[0][k] = b_[0][k] + c_[10] * temp;
                b_[1][k] = b_[1][k] + c_[11] * temp;
                b_[2][k] = b_[2][k] + c_[12] * temp;
                b_[3][k] = b_[3][k] + c_[13] * temp;
                b_[4][k] = b_[4][k] + c_[14] * temp;
            } else if (j == 8) {
                b_[0][k] = b_[0][k] + c_[15] * temp;
                b_[1][k] = b_[1][k] + c_[16] * temp;
                b_[2][k] = b_[2][k] + c_[17] * temp;
                b_[3][k] = b_[3][k] + c_[18] * temp;
                b_[4][k] = b_[4][k] + c_[19] * temp;
                b_[5][k] = b_[5][k] + c_[20] * temp;
            }
        }
    }
    
    return epsi;
}

// Predizione stato finale - TRADUZIONE ESATTA da OrbFit SUBROUTINE rapred
void RA15Integrator::predict(
    std::array<double, 6>& x_and_v,
    double t,
    const std::array<double, 6>& f1)
{
    double t2 = t * t;
    bool ncl = false;  // false = equazioni secondo ordine
    double w1 = 0.5;   // Per secondo ordine
    
    // OrbFit linee 637-644: eq. (2.11) e (2.12)
    for (int k = 0; k < 3; ++k) {
        // Posizione (eq. 2.11)
        x_and_v[k] = x_and_v[k] + x_and_v[k + 3] * t + t2 * (
            f1[k] * w1 + 
            b_[0][k] * w_[0] + 
            b_[1][k] * w_[1] + 
            b_[2][k] * w_[2] + 
            b_[3][k] * w_[3] + 
            b_[4][k] * w_[4] + 
            b_[5][k] * w_[5] + 
            b_[6][k] * w_[6]
        );
        
        // Velocità (eq. 2.12) - solo per equazioni secondo ordine
        if (!ncl) {
            x_and_v[k + 3] = x_and_v[k + 3] + t * (
                f1[k] + 
                b_[0][k] * u_[0] + 
                b_[1][k] * u_[1] + 
                b_[2][k] * u_[2] + 
                b_[3][k] * u_[3] + 
                b_[4][k] * u_[4] + 
                b_[5][k] * u_[5] + 
                b_[6][k] * u_[6]
            );
        }
    }
}

// Extrapolazione coefficienti B
void RA15Integrator::extrapolate(double q, int ns) {
    // q = t_new / t_old
    // Estrapola B da step precedente usando rapporto q
    
    if (ns <= 1) return;  // Non extrapolare primo step
    
    // Formula da OrbFit SUBROUTINE bintrp
    // e = extrapolated, bd = derivative of b
    
    for (int k = 0; k < 6; ++k) {
        for (int l = 0; l < 7; ++l) {
            // Predizione semplice: b_new ≈ b_old * (q^(l+1))
            double factor = std::pow(q, l + 1);
            e_[l][k] = b_[l][k] * factor;
        }
    }
    
    // Copia extrapolazione in B per prossimo step
    for (int l = 0; l < 7; ++l) {
        for (int k = 0; k < 6; ++k) {
            b_[l][k] = e_[l][k];
        }
    }
}

// Calcola nuovo stepsize adattivo
double RA15Integrator::computeNewStepsize(double t, double dir) {
    // Trova massimo di |b[6][k]| (coefficiente più alto)
    double hv = 0.0;
    for (int k = 0; k < 6; ++k) {
        hv = std::max(hv, std::abs(b_[6][k]));
    }
    
    // Formula da OrbFit: hv = max(|b[6]|) * w[6] / t^7
    hv = hv * w_[6] / std::pow(std::abs(t), 7.0);
    
    // Tolleranza ss = 10^(-llev)
    double ss = std::pow(10.0, -options_.llev);
    
    // Nuovo step: (ss/hv)^(1/9) (ordine 15 → esponente 1/9)
    double pw = 1.0 / 9.0;
    double t_new = std::pow(ss / hv, pw) * dir;
    
    // Limita crescita/decrescita
    if (t_new / t > 1.4) {
        t_new = t * 1.4;
    }
    
    return t_new;
}

// Helper: indice matrice triangolare
inline int RA15Integrator::getIndex(int j, int l) const {
    // Matrice triangolare superiore: indice = j*(j-1)/2 + l
    if (l >= j) return 0;  // Fuori triangolo
    return (j * (j - 1)) / 2 + l;
}

// Integrazione principale
OrbitState RA15Integrator::integrate(
    const OrbitState& initial_state,
    const JulianDate& target_epoch,
    AccelerationFunction accel_func)
{
    // Salva funzione accelerazione
    accel_func_ = accel_func;
    
    // Reset statistiche
    stats_ = Statistics{};
    
    // Direzione tempo
    double tf = target_epoch.jd - initial_state.epoch.jd;
    double dir = (tf >= 0.0) ? 1.0 : -1.0;
    
    // Stepsize iniziale
    double tp = options_.h_init * dir;
    if (std::abs(tp) > std::abs(tf)) {
        tp = tf;
    }
    
    // Stato corrente
    std::array<double, 6> x_and_v;
    x_and_v[0] = initial_state.position.x;
    x_and_v[1] = initial_state.position.y;
    x_and_v[2] = initial_state.position.z;
    x_and_v[3] = initial_state.velocity.x;
    x_and_v[4] = initial_state.velocity.y;
    x_and_v[5] = initial_state.velocity.z;
    
    double tm = 0.0;  // Tempo dall'inizio
    bool first_step = true;
    bool last_step = false;
    
    int ns = 0;  // Contatore step
    int ncount = 0;  // Contatore riduzioni step
    
    // Accelerazione iniziale
    Vector3D pos0(x_and_v[0], x_and_v[1], x_and_v[2]);
    Vector3D vel0(x_and_v[3], x_and_v[4], x_and_v[5]);
    Vector3D accel0 = accel_func_(initial_state.epoch, pos0, vel0);
    stats_.num_evals++;
    
    std::array<double, 6> f1;
    f1[0] = accel0.x;
    f1[1] = accel0.y;
    f1[2] = accel0.z;
    f1[3] = 0.0;  // Non usato per equazioni secondo ordine
    f1[4] = 0.0;
    f1[5] = 0.0;
    
    if (options_.verbose) {
        std::cout << "RA15: integrazione da JD " << initial_state.epoch.jd 
                  << " a JD " << target_epoch.jd << " (Δt = " << tf << " giorni)\n";
        std::cout << "Stepsize iniziale: " << std::abs(tp) << " giorni\n";
        std::cout << "Tolleranza: ss = 1e-" << options_.llev << ", eprk = " 
                  << options_.eprk << "\n";
    }
    
    // Loop principale
    while (true) {
        // Check limite step
        if (ns >= options_.max_steps) {
            throw std::runtime_error(
                "RA15: raggiunto limite massimo step (" + 
                std::to_string(options_.max_steps) + ")"
            );
        }
        
        // Reset arrays se primo step
        if (first_step || !options_.fixed_step) {
            resetArrays();
        }
        
        double t = tp;  // Stepsize corrente
        
        // Iterazioni per convergenza
        int ni = first_step ? options_.lit1 : options_.lit2;
        double ep0 = 0.0;
        double ep_current = 0.0;
        
        bool converged = false;
        for (int m = 0; m < ni; ++m) {
            ep_current = iterate(m + 1, t, tm, initial_state.epoch, x_and_v, f1);
            
            if (m == 0) {
                ep0 = ep_current;
            } else if (ep0 > 0.0 && ep_current / ep0 < options_.eprk) {
                converged = true;
                stats_.avg_iterations += m + 1;
                if (options_.verbose && ns % 100 == 0) {
                    std::cout << "Step " << ns << ": t=" << (tm + t) 
                              << ", iter=" << (m + 1) 
                              << ", conv=" << (ep_current / ep0) << "\n";
                }
                break;
            }
        }
        
        // Check convergenza
        if (!converged && !options_.fixed_step) {
            // Riduci stepsize e riprova
            tp *= 0.8;
            ncount++;
            stats_.num_rejections++;
            
            if (ncount > 20) {
                throw std::runtime_error(
                    "RA15: troppi fallimenti convergenza (step = " + 
                    std::to_string(tp) + ")"
                );
            }
            
            if (options_.verbose) {
                std::cout << "RA15: convergenza fallita, riduco step a " 
                          << std::abs(tp) << "\n";
            }
            continue;
        }
        
        // Step accettato - predici stato finale
        predict(x_and_v, t, f1);
        
        tm += t;
        ns++;
        stats_.num_steps++;
        first_step = false;
        
        // Check se finito
        if (last_step || std::abs(tm - tf) < 1e-10) {
            break;
        }
        
        // Calcola nuovo stepsize
        if (!options_.fixed_step) {
            tp = computeNewStepsize(t, dir);
            stats_.final_stepsize = std::abs(tp);
        } else {
            tp = options_.h_init * dir;
        }
        
        // Check se prossimo è ultimo step
        if (dir * (tm + tp) > dir * tf - 1e-8) {
            tp = tf - tm;
            last_step = true;
        }
        
        // Ricalcola accelerazione per prossimo step
        Vector3D pos_new(x_and_v[0], x_and_v[1], x_and_v[2]);
        Vector3D vel_new(x_and_v[3], x_and_v[4], x_and_v[5]);
        JulianDate epoch_new(initial_state.epoch.jd + tm);
        Vector3D accel_new = accel_func_(epoch_new, pos_new, vel_new);
        stats_.num_evals++;
        
        f1[0] = accel_new.x;
        f1[1] = accel_new.y;
        f1[2] = accel_new.z;
        
        // Extrapola coefficienti B per prossimo step
        if (!options_.fixed_step && ns > 1) {
            double q = tp / t;
            extrapolate(q, ns);
        }
    }
    
    // Statistiche finali
    if (stats_.num_steps > 0) {
        stats_.avg_iterations /= stats_.num_steps;
    }
    
    if (options_.verbose) {
        std::cout << "\nRA15 completato:\n"
                  << "  Step totali: " << stats_.num_steps << "\n"
                  << "  Valutazioni forza: " << stats_.num_evals << "\n"
                  << "  Step rifiutati: " << stats_.num_rejections << "\n"
                  << "  Media iterazioni/step: " << std::fixed << std::setprecision(2) 
                  << stats_.avg_iterations << "\n"
                  << "  Stepsize finale: " << stats_.final_stepsize << " giorni\n";
    }
    
    // Costruisci stato finale
    OrbitState final;
    final.position = Vector3D(x_and_v[0], x_and_v[1], x_and_v[2]);
    final.velocity = Vector3D(x_and_v[3], x_and_v[4], x_and_v[5]);
    final.epoch = target_epoch;
    
    return final;
}

} // namespace ioccultcalc
