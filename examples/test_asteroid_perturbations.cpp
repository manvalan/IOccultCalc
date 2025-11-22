// Test per identificare quali asteroidi hanno maggiore effetto su 433 Eros
// Calcola l'accelerazione gravitazionale di ciascun asteroide su Eros @ MJD 55918.807

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

struct AsteroidData {
    int number;
    std::string name;
    double GM_km3s2;  // km³/s²
    bool inIOccultCalc;
};

int main() {
    // Lista completa AST17 di OrbFit con masse da Hilton 1997
    std::vector<AsteroidData> asteroids = {
        {1, "Ceres", 62.6284, true},
        {2, "Pallas", 14.3, true},
        {3, "Juno", 1.82, false},       // MANCANTE
        {4, "Vesta", 17.8, true},
        {6, "Hebe", 0.85, false},        // MANCANTE
        {7, "Iris", 0.79, false},        // MANCANTE
        {10, "Hygiea", 5.8, true},
        {15, "Eunomia", 0.38, false},    // MANCANTE
        {16, "Psyche", 0.36, false},     // MANCANTE
        {29, "Amphitrite", 0.13, false}, // MANCANTE
        {52, "Europa", 0.34, false},     // MANCANTE
        {65, "Cybele", 0.15, false},     // MANCANTE
        {87, "Sylvia", 0.15, false},     // MANCANTE
        {88, "Thisbe", 0.16, false},     // MANCANTE
        {511, "Davida", 0.39, false},    // MANCANTE
        {704, "Interamnia", 0.34, false},// MANCANTE
        {134340, "Pluto", 871.0, false}  // MANCANTE
    };
    
    // Posizione approssimativa di Eros @ MJD 55918.807 (Eclittica J2000)
    // r = (-0.223, 1.065, -0.368) AU
    double eros_x = -0.223;
    double eros_y = 1.065;
    double eros_z = -0.368;
    
    // Posizioni approssimative degli asteroidi @ MJD 55918.807
    // (valori tipici per main-belt asteroids)
    // Per una stima grossolana, assumiamo orbite circolari a vari semi-assi maggiori
    struct Position {
        int number;
        double a;  // semi-asse maggiore tipico (AU)
        double angle;  // angolo approssimativo nella cintura
    };
    
    std::vector<Position> positions = {
        {1, 2.77, 180.0},    // Ceres
        {2, 2.77, 45.0},     // Pallas
        {3, 2.67, 90.0},     // Juno
        {4, 2.36, 270.0},    // Vesta
        {6, 2.43, 30.0},     // Hebe
        {7, 2.39, 120.0},    // Iris
        {10, 3.14, 200.0},   // Hygiea
        {15, 2.64, 60.0},    // Eunomia
        {16, 2.92, 150.0},   // Psyche
        {29, 2.55, 240.0},   // Amphitrite
        {52, 3.10, 300.0},   // Europa
        {65, 3.43, 330.0},   // Cybele
        {87, 3.49, 15.0},    // Sylvia
        {88, 2.77, 75.0},    // Thisbe
        {511, 3.18, 100.0},  // Davida
        {704, 3.06, 210.0},  // Interamnia
        {134340, 39.5, 0.0}  // Pluto (molto lontano)
    };
    
    const double AU_TO_KM = 1.496e8;
    const double G = 6.67430e-20;  // km³/(kg·s²)
    
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "=== ANALISI PERTURBAZIONI ASTEROIDALI SU 433 EROS ===\n\n";
    std::cout << "Epoca: MJD 55918.807 (Dec 23, 2011)\n";
    std::cout << "Eros @ r = (" << eros_x << ", " << eros_y << ", " << eros_z << ") AU\n\n";
    
    struct Perturbation {
        std::string name;
        int number;
        double accel;  // m/s²
        double distance;  // AU
        bool inIOccultCalc;
    };
    
    std::vector<Perturbation> perturbations;
    
    for (size_t i = 0; i < asteroids.size(); i++) {
        const auto& ast = asteroids[i];
        const auto& pos = positions[i];
        
        // Posizione approssimativa dell'asteroide (orbita circolare)
        double angle_rad = pos.angle * M_PI / 180.0;
        double ast_x = pos.a * std::cos(angle_rad);
        double ast_y = pos.a * std::sin(angle_rad);
        double ast_z = 0.0;  // semplificazione: orbita nel piano
        
        // Distanza Eros-Asteroide
        double dx = ast_x - eros_x;
        double dy = ast_y - eros_y;
        double dz = ast_z - eros_z;
        double dist_au = std::sqrt(dx*dx + dy*dy + dz*dz);
        double dist_km = dist_au * AU_TO_KM;
        
        // Accelerazione gravitazionale: a = GM / r²
        double accel_km_s2 = ast.GM_km3s2 * 1e-12 / (dist_km * dist_km);
        double accel_m_s2 = accel_km_s2 * 1000.0;
        
        perturbations.push_back({ast.name, ast.number, accel_m_s2, dist_au, ast.inIOccultCalc});
    }
    
    // Ordina per accelerazione decrescente
    std::sort(perturbations.begin(), perturbations.end(),
              [](const Perturbation& a, const Perturbation& b) {
                  return a.accel > b.accel;
              });
    
    std::cout << std::setw(5) << "Num" << " "
              << std::setw(12) << "Nome" << " "
              << std::setw(12) << "Accel (m/s²)" << " "
              << std::setw(10) << "Dist (AU)" << " "
              << std::setw(12) << "IOccultCalc" << "\n";
    std::cout << std::string(60, '-') << "\n";
    
    double total_included = 0.0;
    double total_missing = 0.0;
    
    for (const auto& p : perturbations) {
        std::cout << std::setw(5) << p.number << " "
                  << std::setw(12) << p.name << " "
                  << std::scientific << std::setprecision(3)
                  << std::setw(12) << p.accel << " "
                  << std::fixed << std::setprecision(2)
                  << std::setw(10) << p.distance << " "
                  << std::setw(12) << (p.inIOccultCalc ? "✓ Incluso" : "✗ MANCANTE")
                  << "\n";
        
        if (p.inIOccultCalc) {
            total_included += p.accel;
        } else {
            total_missing += p.accel;
        }
    }
    
    std::cout << "\n=== SOMMARIO ===\n";
    std::cout << "Accelerazione totale inclusa in IOccultCalc: " 
              << std::scientific << total_included << " m/s²\n";
    std::cout << "Accelerazione totale MANCANTE: " 
              << std::scientific << total_missing << " m/s²\n";
    std::cout << "Percentuale mancante: " 
              << std::fixed << std::setprecision(1)
              << (total_missing / (total_included + total_missing) * 100.0) << "%\n\n";
    
    std::cout << "NOTA: Questo è un calcolo APPROSSIMATIVO con posizioni stimate.\n";
    std::cout << "      Per risultati precisi occorrerebbe leggere le effemeridi reali.\n";
    
    return 0;
}
