# -*- coding: utf-8 -*-
import sys
import glob

def check_mol2_charge(filename):
    total_charge = 0.0
    atom_count = 0
    in_atom_section = False
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('@<TRIPOS>ATOM'):
                    in_atom_section = True
                    continue
                elif line.startswith('@<TRIPOS>BOND') or line.startswith('@<TRIPOS>SUBSTRUCTURE'):
                    in_atom_section = False
                    continue
                    
                if in_atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 9:
                        charge = float(parts[8])
                        total_charge += charge
                        atom_count += 1
                        
        print(f"[{filename}]")
        print(f"  Atoms:  {atom_count}")
        print(f"  Charge: {total_charge:+.6f}\n")
        
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
    except Exception as e:
        print(f"Error processing '{filename}': {e}")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Check specific files provided as arguments
        for arg in sys.argv[1:]:
            # Handle wildcards manually for Windows/some shells
            for file in glob.glob(arg):
                check_mol2_charge(file)
    else:
        print("Usage: python check_charge.py <file1.mol2> [file2.mol2 ...]")
        print("Example: python check_charge.py lig_part*.mol2")
