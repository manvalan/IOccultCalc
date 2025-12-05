/**
 * @file test_euler_angles_debug.cpp
 * @brief Debug delle rotazioni euleriane - confronto con documento ITALOccultLibrary
 */

#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    std::cout << "═════════════════════════════════════════════════════════════════\n";
    std::cout << "  DEBUG: Euler angles and rotation matrices\n";
    std::cout << "═════════════════════════════════════════════════════════════════\n";
    
    // Document example data
    double Omega_deg = 94.058;
    double omega_deg = 110.28;
    double i_deg = 2.9046;
    
    double Omega = Omega_deg * M_PI / 180.0;
    double omega = omega_deg * M_PI / 180.0;
    double i = i_deg * M_PI / 180.0;
    
    std::cout << "\nAngles (from document):\n";
    std::cout << "  Ω = " << Omega_deg << "° = " << std::setprecision(6) << Omega << " rad\n";
    std::cout << "  ω = " << omega_deg << "° = " << omega << " rad\n";
    std::cout << "  i = " << i_deg << "° = " << i << " rad\n";
    
    // Calculate trigonometric values
    double cos_Omega = cos(Omega);
    double sin_Omega = sin(Omega);
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double cos_i = cos(i);
    double sin_i = sin(i);
    
    std::cout << "\nTrigonometric values:\n";
    std::cout << "  cos(ω) = " << std::setprecision(6) << cos_omega 
              << ", sin(ω) = " << sin_omega << "\n";
    std::cout << "  cos(i) = " << cos_i << ", sin(i) = " << sin_i << "\n";
    std::cout << "  cos(Ω) = " << cos_Omega << ", sin(Ω) = " << sin_Omega << "\n";
    
    // Document matrix values (from FASE 2.6)
    double P11_doc = -0.9181;
    double P12_doc = 0.3940;
    double P21_doc = -0.3949;
    double P22_doc = -0.9181;
    double P31_doc = 0.0479;
    double P32_doc = -0.0164;
    
    std::cout << "\nMatrix values from document:\n";
    std::cout << "  P11 = " << P11_doc << ", P12 = " << P12_doc << "\n";
    std::cout << "  P21 = " << P21_doc << ", P22 = " << P22_doc << "\n";
    std::cout << "  P31 = " << P31_doc << ", P32 = " << P32_doc << "\n";
    
    // Calculate matrix elements using three separate rotations
    std::cout << "\n=== METHOD 1: Three separate rotations (my method) ===\n";
    
    // First: ω around Z
    std::cout << "\nRotation 1: ω around Z-axis\n";
    std::cout << "  [cos(ω)  -sin(ω)  0]\n";
    std::cout << "  [sin(ω)   cos(ω)  0]\n";
    std::cout << "  [  0       0      1]\n";
    
    // Second: i around X
    std::cout << "\nRotation 2: i around X-axis\n";
    std::cout << "  [1    0      0 ]\n";
    std::cout << "  [0  cos(i) -sin(i)]\n";
    std::cout << "  [0  sin(i)  cos(i)]\n";
    
    // Third: Ω around Z
    std::cout << "\nRotation 3: Ω around Z-axis\n";
    std::cout << "  [cos(Ω)  -sin(Ω)  0]\n";
    std::cout << "  [sin(Ω)   cos(Ω)  0]\n";
    std::cout << "  [  0       0      1]\n";
    
    // Compose: R = R_z(Ω) · R_x(i) · R_z(ω)
    std::cout << "\n=== Composed matrix P = R_z(Ω) · R_x(i) · R_z(ω) ===\n";
    
    // First R_z(ω):
    // [cos_omega  -sin_omega  0]
    // [sin_omega   cos_omega  0]
    // [  0           0        1]
    
    // Then R_x(i):
    // [1      0        0    ]
    // [0   cos_i   -sin_i  ]
    // [0   sin_i    cos_i  ]
    
    // Intermediate: R_x(i) · R_z(ω)
    double r_11 = cos_omega;
    double r_12 = -sin_omega;
    double r_21 = sin_omega * cos_i;
    double r_22 = cos_omega * cos_i;
    double r_31 = sin_omega * sin_i;
    double r_32 = -cos_omega * sin_i;
    
    std::cout << "\nIntermediate R_x(i) · R_z(ω):\n";
    std::cout << "  r11 = " << std::setprecision(6) << r_11 << "\n";
    std::cout << "  r12 = " << r_12 << "\n";
    std::cout << "  r21 = " << r_21 << "\n";
    std::cout << "  r22 = " << r_22 << "\n";
    std::cout << "  r31 = " << r_31 << "\n";
    std::cout << "  r32 = " << r_32 << "\n";
    
    // Finally R_z(Ω):
    // [cos_Omega  -sin_Omega  0]
    // [sin_Omega   cos_Omega  0]
    // [  0           0        1]
    
    double P11 = cos_Omega * r_11 - sin_Omega * r_21;
    double P12 = cos_Omega * r_12 - sin_Omega * r_22;
    double P21 = sin_Omega * r_11 + cos_Omega * r_21;
    double P22 = sin_Omega * r_12 + cos_Omega * r_22;
    double P31 = r_31;
    double P32 = r_32;
    
    std::cout << "\nFinal matrix P (my calculation):\n";
    std::cout << "  P11 = " << P11 << ", P12 = " << P12 << "\n";
    std::cout << "  P21 = " << P21 << ", P22 = " << P22 << "\n";
    std::cout << "  P31 = " << P31 << ", P32 = " << P32 << "\n";
    
    std::cout << "\nComparison with document:\n";
    std::cout << "  P11: calc=" << std::setprecision(5) << P11 << " doc=" << P11_doc 
              << " Δ=" << (P11 - P11_doc) << "\n";
    std::cout << "  P12: calc=" << P12 << " doc=" << P12_doc 
              << " Δ=" << (P12 - P12_doc) << "\n";
    std::cout << "  P21: calc=" << P21 << " doc=" << P21_doc 
              << " Δ=" << (P21 - P21_doc) << "\n";
    std::cout << "  P22: calc=" << P22 << " doc=" << P22_doc 
              << " Δ=" << (P22 - P22_doc) << "\n";
    std::cout << "  P31: calc=" << P31 << " doc=" << P31_doc 
              << " Δ=" << (P31 - P31_doc) << "\n";
    std::cout << "  P32: calc=" << P32 << " doc=" << P32_doc 
              << " Δ=" << (P32 - P32_doc) << "\n";
    
    // Test on orbital plane position
    double x_orb_doc = 2.7393;
    double y_orb_doc = 1.3337;
    
    double x_ecl_calc = P11 * x_orb_doc + P12 * y_orb_doc;
    double y_ecl_calc = P21 * x_orb_doc + P22 * y_orb_doc;
    
    std::cout << "\n=== TEST: Apply matrix to orbital position ===\n";
    std::cout << "Orbital position: (" << x_orb_doc << ", " << y_orb_doc << ")\n";
    std::cout << "Calculated ecliptic: (" << std::setprecision(6) << x_ecl_calc 
              << ", " << y_ecl_calc << ")\n";
    std::cout << "Expected from document: (-1.9899, -2.3073)\n";
    
    double err_x = x_ecl_calc - (-1.9899);
    double err_y = y_ecl_calc - (-2.3073);
    std::cout << "Error: Δx=" << err_x << ", Δy=" << err_y << "\n";
    
    return 0;
}
