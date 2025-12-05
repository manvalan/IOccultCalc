#include <chebyshev_rkf78_propagation.h>
#include <eq1_parser.h>
#include <iostream>

int main() {
    EQ1Parser parser;
    auto result = parser.parse("17030_astdys.eq1");
    if (!result.has_value()) {
        std::cerr << "Errore parse\n";
        return 1;
    }
    
    auto elements = result->elements;
    ChebyshevRKF78Propagator prop(elements, 61000.0);
    
    auto positions = prop.propagateForChebyshev(60642.0, 60643.0, 50);
    
    std::cout << "Posizioni propagate: " << positions.size() << std::endl;
    if (!positions.empty()) {
        std::cout << "Prima: " << positions[0].transpose() << std::endl;
        std::cout << "Ultima: " << positions.back().transpose() << std::endl;
    }
    
    return 0;
}
