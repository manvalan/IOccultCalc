#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

// Implementazione diretta del metodo accurato
constexpr double PI = 3.14159265358979323846;

double interpolateLinear(const std::vector<double>& x,
                        const std::vector<double>& y, 
                        double target_x) {
    if (x.empty() || y.empty() || x.size() != y.size()) {
        throw std::runtime_error("Dati invalidi per interpolazione lineare");
    }
    
    if (target_x <= x.front()) return y.front();
    if (target_x >= x.back()) return y.back();
    
    for (size_t i = 0; i < x.size() - 1; i++) {
        if (target_x >= x[i] && target_x <= x[i + 1]) {
            double t = (target_x - x[i]) / (x[i + 1] - x[i]);
            return y[i] + t * (y[i + 1] - y[i]);
        }
    }
    
    return y.back();
}

std::vector<double> calculateChebyshevCoeffs_Accurate(
    const std::vector<double>& times, const std::vector<double>& values, int degree) {
    
    int n_coeffs = degree + 1;
    std::vector<double> coeffs(n_coeffs, 0.0);
    
    // Approssimazione con nodi di Chebyshev ottimali
    std::vector<double> cheb_nodes(n_coeffs);
    std::vector<double> cheb_values(n_coeffs);
    
    // Genera nodi di Chebyshev nell'intervallo [times[0], times.back()]
    for (int i = 0; i < n_coeffs; i++) {
        double theta = PI * (i + 0.5) / n_coeffs;
        double x_cheb = std::cos(theta);  // Nodo in [-1, 1]
        cheb_nodes[i] = times[0] + (times.back() - times[0]) * (x_cheb + 1.0) / 2.0;
        cheb_values[i] = interpolateLinear(times, values, cheb_nodes[i]);
    }
    
    std::cout << "🎯 NODI DI CHEBYSHEV GENERATI:\n";
    for (int i = 0; i < n_coeffs; i++) {
        std::cout << "  Nodo " << i+1 << ": MJD " << cheb_nodes[i] 
                  << " → Valore=" << cheb_values[i] << "\n";
    }
    
    // Calcola coefficienti con DCT sui nodi ottimali
    for (int k = 0; k < n_coeffs; k++) {
        double sum = 0.0;
        for (int j = 0; j < n_coeffs; j++) {
            double theta_j = PI * (j + 0.5) / n_coeffs;
            sum += cheb_values[j] * std::cos(k * theta_j);
        }
        coeffs[k] = (2.0 / n_coeffs) * sum;
        if (k == 0) coeffs[k] /= 2.0;
        
        std::cout << "  C_" << k << " = " << std::scientific << coeffs[k] << "\n";
    }
    
    return coeffs;
}

double evaluateChebyshevPolynomial(const std::vector<double>& coeffs, double x) {
    if (coeffs.empty()) return 0.0;
    if (coeffs.size() == 1) return coeffs[0];
    
    // Valutazione ricorsiva stabile
    double result = coeffs[0];
    if (coeffs.size() > 1) {
        double T_prev = 1.0;      // T_0(x)
        double T_curr = x;        // T_1(x)
        result += coeffs[1] * T_curr;
        
        std::cout << "    T_0=" << T_prev << ", T_1=" << T_curr << "\n";
        
        for (size_t n = 2; n < coeffs.size(); n++) {
            double T_next = 2.0 * x * T_curr - T_prev;
            result += coeffs[n] * T_next;
            
            std::cout << "    T_" << n << "=" << T_next << "\n";
            
            T_prev = T_curr;
            T_curr = T_next;
        }
    }
    
    return result;
}

int main() {
    std::cout << "🔧 TEST METODO CHEBYSHEV ACCURATO (DCT)\n";
    std::cout << "=======================================\n\n";
    
    // Dati test da RKF78 (dal debug precedente)
    std::vector<double> times = {
        60001.024306, 60001.858139, 60002.691306, 60003.524639,
        60004.358139, 60005.191639, 60006.024306, 60006.858139,
        60007.691639, 60008.524639, 60009.358139, 60010.191639, 60011.024306
    };
    
    std::vector<double> ra_values = {
        281.196435, 281.476025, 281.754074, 282.030562,
        282.305466, 282.578766, 282.850441, 283.120468,
        283.388827, 283.655496, 283.920452, 284.183674, 284.445140
    };
    
    std::vector<double> dec_values = {
        -22.196055, -22.179831, -22.163237, -22.146283,
        -22.128980, -22.111339, -22.093371, -22.075085,
        -22.056494, -22.037608, -22.018439, -21.998997, -21.979295
    };
    
    double target_mjd = 60006.024306;  // Centro finestra
    double target_ra_exact = 282.850441;   // Valore RKF78 esatto
    double target_dec_exact = -22.093371;
    
    std::cout << "📊 DATI INPUT:\n";
    std::cout << "  Finestra: [" << times[0] << ", " << times.back() << "] MJD\n";
    std::cout << "  Target: " << target_mjd << " MJD\n";
    std::cout << "  RA RKF78: " << target_ra_exact << "°\n";
    std::cout << "  Dec RKF78: " << target_dec_exact << "°\n\n";
    
    // Test metodo accurato per RA
    std::cout << "🔴 CALCOLO COEFFICIENTI RA (METODO ACCURATO):\n";
    auto ra_coeffs = calculateChebyshevCoeffs_Accurate(times, ra_values, 8);
    
    std::cout << "\n🔵 CALCOLO COEFFICIENTI Dec (METODO ACCURATO):\n";
    auto dec_coeffs = calculateChebyshevCoeffs_Accurate(times, dec_values, 8);
    
    // Valutazione al target
    std::cout << "\n🎯 VALUTAZIONE AL TARGET:\n";
    double t_norm = 2.0 * (target_mjd - times[0]) / (times.back() - times[0]) - 1.0;
    std::cout << "  Tempo normalizzato: " << t_norm << "\n\n";
    
    std::cout << "  📊 VALUTAZIONE RA:\n";
    double ra_result = evaluateChebyshevPolynomial(ra_coeffs, t_norm);
    
    std::cout << "\n  📊 VALUTAZIONE Dec:\n";
    double dec_result = evaluateChebyshevPolynomial(dec_coeffs, t_norm);
    
    std::cout << "\n✅ RISULTATI FINALI:\n";
    std::cout << "  Chebyshev RA:  " << std::fixed << std::setprecision(6) << ra_result << "°\n";
    std::cout << "  Chebyshev Dec: " << dec_result << "°\n";
    std::cout << "  RKF78 RA:      " << target_ra_exact << "°\n";
    std::cout << "  RKF78 Dec:     " << target_dec_exact << "°\n";
    
    double ra_error = std::abs(ra_result - target_ra_exact);
    double dec_error = std::abs(dec_result - target_dec_exact);
    double total_error_arcsec = std::sqrt(ra_error*ra_error + dec_error*dec_error) * 3600.0;
    
    std::cout << "\n📏 ERRORI:\n";
    std::cout << "  ΔRA:  " << ra_error << "° = " << ra_error*3600.0 << "\"\n";
    std::cout << "  ΔDec: " << dec_error << "° = " << dec_error*3600.0 << "\"\n";
    std::cout << "  Totale: " << total_error_arcsec << " arcsec\n\n";
    
    if (total_error_arcsec < 1.0) {
        std::cout << "✅ OTTIMO! Errore sub-arcsecondo\n";
    } else if (total_error_arcsec < 10.0) {
        std::cout << "✅ BUONO! Errore accettabile\n";
    } else if (total_error_arcsec < 100.0) {
        std::cout << "⚠️ MEDIO: Miglioramento necessario\n";
    } else {
        std::cout << "❌ CATTIVO: Algoritmo ancora problematico\n";
    }
    
    return 0;
}