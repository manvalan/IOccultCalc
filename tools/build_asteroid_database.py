#!/usr/bin/env python3
"""
Build Asteroid Database from AstDyS and MPC
============================================

Combines orbital elements from AstDyS allnum.cat with physical properties
from MPC Extended JSON to create a unified asteroid database.

Features:
- Downloads allnum.cat from AstDyS
- Downloads MPC Extended JSON
- Merges orbital elements (AstDyS) with physical properties (MPC)
- Tracks parsing date for periodic updates
- Supports monthly auto-update

Usage:
    python3 build_asteroid_database.py [--force] [--output OUTPUT_FILE]
    
Options:
    --force          Force update even if database is recent
    --output FILE    Output JSON file path (default: ~/.ioccultcalc/data/all_numbered_asteroids.json)
    --max-age DAYS   Maximum age in days before forcing update (default: 30)
"""

import json
import gzip
import os
import sys
import argparse
import urllib.request
import urllib.error
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Any
import time

# URLs
ASTDYS_ALLNUM_URL = "https://newton.spacedys.com/astdys2/catalogs/allnum.cat"
MPC_EXTENDED_URL = "https://minorplanetcenter.net/Extended_Files/mpcorb_extended.json.gz"

# Default paths
DEFAULT_DATA_DIR = Path.home() / ".ioccultcalc" / "data"
DEFAULT_OUTPUT_FILE = DEFAULT_DATA_DIR / "all_numbered_asteroids.json"


def download_file(url: str, output_path: Path, decompress: bool = False) -> bool:
    """Download a file from URL to local path."""
    print(f"Downloading {url}...")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"✓ Downloaded to {output_path}")
        
        if decompress and output_path.suffix == '.gz':
            print(f"Decompressing {output_path}...")
            with gzip.open(output_path, 'rb') as f_in:
                decompressed_path = output_path.with_suffix('')
                with open(decompressed_path, 'wb') as f_out:
                    f_out.write(f_in.read())
            output_path.unlink()  # Remove .gz file
            print(f"✓ Decompressed to {decompressed_path}")
            return True
        
        return True
    except urllib.error.URLError as e:
        print(f"✗ Error downloading {url}: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def parse_allnum_cat(file_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Parse AstDyS allnum.cat file.
    
    Format: 'num' epoch a e i node argperi M H G flag
    """
    print(f"Parsing allnum.cat from {file_path}...")
    asteroids = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    found_header_end = False
    parsed_count = 0
    
    for line in lines:
        line = line.strip()
        
        # Skip until END_OF_HEADER
        if not found_header_end:
            if "END_OF_HEADER" in line:
                found_header_end = True
            continue
        
        if not line or line.startswith('#'):
            continue
        
        # Parse line: 'num' epoch a e i node argperi M H G flag
        # Example: '704'  61000.000000  2.793  0.123  10.5  80.3  73.6  77.4  3.5  0.15  0
        
        try:
            # Extract number in quotes
            if line[0] != "'":
                continue
            
            quote_end = line.find("'", 1)
            if quote_end == -1:
                continue
            
            asteroid_num = line[1:quote_end].strip()
            if not asteroid_num.isdigit():
                continue
            
            # Parse rest of line
            data_part = line[quote_end + 1:].strip()
            parts = data_part.split()
            
            if len(parts) < 9:
                continue
            
            epoch_mjd = float(parts[0])
            a = float(parts[1])
            e = float(parts[2])
            i_deg = float(parts[3])
            node_deg = float(parts[4])
            argperi_deg = float(parts[5])
            M_deg = float(parts[6])
            H = float(parts[7])
            G = float(parts[8])
            
            asteroids[asteroid_num] = {
                "number": int(asteroid_num),
                "designation": asteroid_num,
                "epoch_mjd": epoch_mjd,
                "epoch_jd": epoch_mjd + 2400000.5,
                "a": a,
                "e": e,
                "i": i_deg,
                "node": node_deg,
                "argperi": argperi_deg,
                "M": M_deg,
                "H": H,
                "G": G,
                "source": "AstDyS"
            }
            
            parsed_count += 1
            if parsed_count % 10000 == 0:
                print(f"  Parsed {parsed_count} asteroids...")
                
        except (ValueError, IndexError) as e:
            # Skip malformed lines
            continue
    
    print(f"✓ Parsed {parsed_count} asteroids from allnum.cat")
    return asteroids


def parse_mpc_extended(file_path: Path) -> Dict[str, Dict[str, Any]]:
    """
    Parse MPC Extended JSON file.
    
    Returns dict keyed by asteroid number with physical properties.
    """
    print(f"Parsing MPC Extended JSON from {file_path}...")
    
    with open(file_path, 'r') as f:
        mpc_data = json.load(f)
    
    mpc_asteroids = {}
    
    # MPC Extended JSON is an array of asteroid objects
    if isinstance(mpc_data, list):
        for ast in mpc_data:
            num = ast.get("number")
            if num:
                num_str = str(num)
                mpc_asteroids[num_str] = {
                    "name": ast.get("name", ""),
                    "designation": ast.get("designation", ""),
                    "diameter": ast.get("diameter"),
                    "albedo": ast.get("albedo"),
                    "spec_T": ast.get("spec_T"),  # Spectral type
                    "orbit_type": ast.get("orbit_type"),  # MBA, NEA, etc.
                    "H": ast.get("H"),
                    "G": ast.get("G"),
                    "source": "MPC"
                }
    
    print(f"✓ Parsed {len(mpc_asteroids)} asteroids from MPC Extended")
    return mpc_asteroids


def merge_asteroid_data(astdys_data: Dict, mpc_data: Dict) -> List[Dict[str, Any]]:
    """
    Merge AstDyS orbital elements with MPC physical properties.
    
    AstDyS provides high-precision orbital elements.
    MPC provides physical properties (diameter, albedo, etc.).
    """
    print("Merging AstDyS and MPC data...")
    
    merged = []
    merged_count = 0
    mpc_only_count = 0
    
    # Start with all AstDyS asteroids (they have orbital elements)
    for num, astdys_ast in astdys_data.items():
        merged_ast = astdys_ast.copy()
        
        # Add MPC data if available
        if num in mpc_data:
            mpc_ast = mpc_data[num]
            merged_ast.update({
                "name": mpc_ast.get("name", ""),
                "diameter": mpc_ast.get("diameter"),
                "albedo": mpc_ast.get("albedo"),
                "spec_T": mpc_ast.get("spec_T"),
                "orbit_type": mpc_ast.get("orbit_type"),
            })
            # Use MPC H and G if available (may be more recent)
            if mpc_ast.get("H") is not None:
                merged_ast["H"] = mpc_ast["H"]
            if mpc_ast.get("G") is not None:
                merged_ast["G"] = mpc_ast["G"]
            merged_count += 1
        else:
            mpc_only_count += 1
        
        merged.append(merged_ast)
    
    print(f"✓ Merged {merged_count} asteroids with MPC data")
    print(f"  {mpc_only_count} asteroids have only AstDyS data")
    
    return merged


def check_update_needed(output_file: Path, max_age_days: int = 30) -> bool:
    """Check if database needs update based on file age."""
    if not output_file.exists():
        print(f"Database file not found: {output_file}")
        return True
    
    # Check metadata in existing file
    try:
        with open(output_file, 'r') as f:
            data = json.load(f)
        
        metadata = data.get("metadata", {})
        last_update_str = metadata.get("last_update")
        
        if last_update_str:
            try:
                last_update = datetime.fromisoformat(last_update_str.replace('Z', '+00:00'))
                age_days = (datetime.now() - last_update.replace(tzinfo=None)).days
                
                if age_days < max_age_days:
                    print(f"Database is {age_days} days old (max: {max_age_days} days)")
                    print(f"Last update: {last_update_str}")
                    return False
                else:
                    print(f"Database is {age_days} days old, update needed")
                    return True
            except ValueError:
                print("Could not parse last_update date, forcing update")
                return True
    except (json.JSONDecodeError, KeyError):
        print("Could not read metadata from existing file, forcing update")
        return True
    
    return True


def build_database(output_file: Path, force: bool = False, max_age_days: int = 30) -> bool:
    """Main function to build asteroid database."""
    
    # Check if update needed
    if not force and not check_update_needed(output_file, max_age_days):
        print("Database is up to date. Use --force to update anyway.")
        return True
    
    # Create data directory
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Temporary files
    temp_dir = Path.home() / ".ioccultcalc" / "cache"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    allnum_path = temp_dir / "allnum.cat"
    mpc_path = temp_dir / "mpcorb_extended.json"
    
    # Download files
    print("\n=== Downloading data files ===")
    if not download_file(ASTDYS_ALLNUM_URL, allnum_path):
        return False
    
    if not download_file(MPC_EXTENDED_URL, mpc_path, decompress=True):
        return False
    
    # Parse files
    print("\n=== Parsing data files ===")
    astdys_data = parse_allnum_cat(allnum_path)
    mpc_data = parse_mpc_extended(mpc_path)
    
    # Merge data
    print("\n=== Merging data ===")
    merged_asteroids = merge_asteroid_data(astdys_data, mpc_data)
    
    # Create output structure
    print("\n=== Creating output file ===")
    output_data = {
        "metadata": {
            "version": "2.0",
            "source": "AstDyS_allnum_cat + MPC_Extended",
            "astdys_url": ASTDYS_ALLNUM_URL,
            "mpc_url": MPC_EXTENDED_URL,
            "last_update": datetime.now().isoformat(),
            "total_asteroids": len(merged_asteroids),
            "with_mpc_data": sum(1 for a in merged_asteroids if a.get("diameter") is not None),
            "parsing_date": datetime.now().isoformat()
        },
        "asteroids": merged_asteroids
    }
    
    # Write output
    print(f"Writing to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"✓ Database created: {output_file}")
    print(f"  Total asteroids: {len(merged_asteroids)}")
    print(f"  With MPC data: {output_data['metadata']['with_mpc_data']}")
    print(f"  Last update: {output_data['metadata']['last_update']}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Build asteroid database from AstDyS and MPC",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force update even if database is recent"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        default=str(DEFAULT_OUTPUT_FILE),
        help=f"Output JSON file path (default: {DEFAULT_OUTPUT_FILE})"
    )
    
    parser.add_argument(
        "--max-age",
        type=int,
        default=30,
        help="Maximum age in days before forcing update (default: 30)"
    )
    
    args = parser.parse_args()
    
    output_path = Path(args.output)
    
    print("=" * 70)
    print("Asteroid Database Builder")
    print("Combining AstDyS allnum.cat + MPC Extended JSON")
    print("=" * 70)
    print()
    
    success = build_database(output_path, force=args.force, max_age_days=args.max_age)
    
    if success:
        print("\n✓ Database build completed successfully!")
        sys.exit(0)
    else:
        print("\n✗ Database build failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()

