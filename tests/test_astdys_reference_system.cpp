/**
 * @file test_astdys_reference_system.cpp
 * @brief Verifica sistema di riferimento elementi AstDyS .eq1
 */

#include "ioccultcalc/astdys_client.h"
#include <iostream>
#include <sstream>
#include <string>

using namespace ioccultcalc;

int main() {
    std::cout << "\n╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║  VERIFICA SISTEMA DI RIFERIMENTO ASTDYS .eq1              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n\n";
    
    try {
        // Scarica file direttamente usando curl
        std::string url = "https://newton.spacedys.com/~astdys2/epoch/numbered/17/17030.eq1";
        
        // Usa curl per scaricare (semplificato per il test)
        std::string cmd = "curl -s \"" + url + "\"";
        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            throw std::runtime_error("Failed to execute curl");
        }
        
        std::string content;
        char buffer[128];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            content += buffer;
        }
        pclose(pipe);
        
        std::cout << "Header del file .eq1:\n";
        std::cout << "─────────────────────────────────────────────────────────\n";
        
        std::istringstream iss(content);
        std::string line;
        int header_lines = 0;
        
        while (std::getline(iss, line) && header_lines < 10) {
            if (line.find("END_OF_HEADER") != std::string::npos) {
                std::cout << line << "\n";
                break;
            }
            std::cout << line << "\n";
            header_lines++;
        }
        
        std::cout << "─────────────────────────────────────────────────────────\n\n";
        
        // Analizza header
        iss.clear();
        iss.str(content);
        std::string refsys;
        std::string format;
        std::string rectype;
        
        while (std::getline(iss, line)) {
            if (line.find("refsys") == 0 || line.find("refsys") == 2) {
                // Estrai valore dopo "="
                size_t eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    refsys = line.substr(eq_pos + 1);
                    // Rimuovi spazi e commenti
                    size_t comment_pos = refsys.find('!');
                    if (comment_pos != std::string::npos) {
                        refsys = refsys.substr(0, comment_pos);
                    }
                    // Trim
                    size_t start = refsys.find_first_not_of(" \t");
                    size_t end = refsys.find_last_not_of(" \t");
                    if (start != std::string::npos && end != std::string::npos) {
                        refsys = refsys.substr(start, end - start + 1);
                    }
                }
            }
            if (line.find("format") == 0 || line.find("format") == 2) {
                size_t eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    format = line.substr(eq_pos + 1);
                    size_t comment_pos = format.find('!');
                    if (comment_pos != std::string::npos) {
                        format = format.substr(0, comment_pos);
                    }
                    size_t start = format.find_first_not_of(" \t");
                    size_t end = format.find_last_not_of(" \t");
                    if (start != std::string::npos && end != std::string::npos) {
                        format = format.substr(start, end - start + 1);
                    }
                }
            }
            if (line.find("rectype") == 0 || line.find("rectype") == 2) {
                size_t eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    rectype = line.substr(eq_pos + 1);
                    size_t comment_pos = rectype.find('!');
                    if (comment_pos != std::string::npos) {
                        rectype = rectype.substr(0, comment_pos);
                    }
                    size_t start = rectype.find_first_not_of(" \t");
                    size_t end = rectype.find_last_not_of(" \t");
                    if (start != std::string::npos && end != std::string::npos) {
                        rectype = rectype.substr(start, end - start + 1);
                    }
                }
            }
            if (line.find("END_OF_HEADER") != std::string::npos) {
                break;
            }
        }
        
        std::cout << "Analisi sistema di riferimento:\n";
        std::cout << "  Format: " << format << "\n";
        std::cout << "  Record Type: " << rectype << "\n";
        std::cout << "  Reference System: " << refsys << "\n\n";
        
        // Interpretazione
        std::cout << "Interpretazione:\n";
        if (refsys.find("ECLM") != std::string::npos || refsys.find("ECLM") != std::string::npos) {
            std::cout << "  ✓ ECLM = Eclittico Medio (Mean Ecliptic)\n";
            std::cout << "  ✓ Sistema: Eclittico Medio J2000\n";
            std::cout << "  ✓ Elementi: Equinoctiali eclittici medi\n";
        } else if (refsys.find("ECL") != std::string::npos) {
            std::cout << "  ✓ ECL = Eclittico\n";
            std::cout << "  ✓ Sistema: Eclittico J2000\n";
        } else if (refsys.find("EQ") != std::string::npos || refsys.find("EQU") != std::string::npos) {
            std::cout << "  ⚠ EQ/EQU = Equatoriale\n";
            std::cout << "  ⚠ Sistema: Equatoriale (ICRF)\n";
            std::cout << "  ⚠ ATTENZIONE: Gli elementi potrebbero essere equatoriali!\n";
        } else {
            std::cout << "  ? Sistema non riconosciuto: " << refsys << "\n";
        }
        
        std::cout << "\nFormato elementi:\n";
        if (rectype.find("ML") != std::string::npos) {
            std::cout << "  ✓ ML = Mean Longitude (Longitudine Media)\n";
            std::cout << "  ✓ Elementi equinoctiali con longitudine media\n";
        } else if (rectype.find("1L") != std::string::npos) {
            std::cout << "  ✓ 1L = Single Line (Una linea)\n";
        }
        
        std::cout << "\nConclusione:\n";
        if (refsys.find("ECLM") != std::string::npos) {
            std::cout << "  Gli elementi AstDyS .eq1 sono EQUINOCTIALI ECLITTICI MEDI J2000\n";
            std::cout << "  Questo è corretto per la propagazione con AstDyn\n";
        } else if (refsys.find("EQ") != std::string::npos || refsys.find("EQU") != std::string::npos) {
            std::cout << "  ⚠ ATTENZIONE: Gli elementi potrebbero essere EQUATORIALI\n";
            std::cout << "  ⚠ Questo richiederebbe una conversione prima della propagazione\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}

