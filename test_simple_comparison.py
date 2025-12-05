#!/usr/bin/env python3
"""
Test semplice di confronto con IOccultCalc
==========================================

Confronta i risultati IOccultCalc con JPL Horizons usando solo i binari esistenti.
"""

import sys
import subprocess
import json

def run_ioccultcalc_basic(asteroid_num: int, mjd: float):
    """Esegue IOccultCalc con il binario basic esistente"""
    
    try:
        # Usa l'esempio basic se esiste
        cmd = ['./build/examples/example_basic', str(asteroid_num)]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd='/Users/michelebigi/VisualStudio Code/GitHub/IOccultCalc')
        
        if result.returncode == 0:
            print(f"✅ IOccultCalc basic output:")
            print(result.stdout)
            return True
        else:
            print(f"❌ IOccultCalc error: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"❌ Errore esecuzione IOccultCalc: {e}")
        return False

def main():
    if len(sys.argv) < 3:
        print("Uso: python3 test_simple_comparison.py <asteroid> <mjd>")
        print("Esempio: python3 test_simple_comparison.py 433 60006.024306")
        sys.exit(1)
    
    asteroid = int(sys.argv[1])
    mjd = float(sys.argv[2])
    
    print("🌌 CONFRONTO SEMPLICE IOccultCalc vs JPL")
    print("=" * 50)
    print(f"🎯 Asteroide: {asteroid}")
    print(f"📅 MJD: {mjd:.6f}")
    print()
    
    # Test IOccultCalc
    print("🔄 Test IOccultCalc...")
    ioccult_ok = run_ioccultcalc_basic(asteroid, mjd)
    
    print()
    print("=" * 30)
    
    # Test JPL Horizons
    print("🔄 Test JPL Horizons...")
    try:
        result = subprocess.run(['python3', 'query_horizons_simple.py', str(asteroid), str(mjd)], 
                              capture_output=True, text=True, 
                              cwd='/Users/michelebigi/VisualStudio Code/GitHub/IOccultCalc')
        
        if result.returncode == 0:
            print("✅ JPL Query completata:")
            print(result.stdout)
        else:
            print("❌ Errore JPL query:")
            print(result.stderr)
            
    except Exception as e:
        print(f"❌ Errore query JPL: {e}")
    
    print()
    print("=" * 50)
    print("📊 CONCLUSIONE: I risultati sono mostrati sopra per confronto manuale")

if __name__ == "__main__":
    main()