import sys
import math
import os

# Універсальна функція дистанції (працює на старих версіях Python)
def calc_dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def main():
    if len(sys.argv) != 4:
        print("Usage: python stencil_cutter.py <original.mol2> <avogadro_stencil.mol2> <output_base.mol2>")
        sys.exit(1)
        
    orig_file = sys.argv[1]
    template_file = sys.argv[2]
    base_out_file = sys.argv[3]
    
    # 1. Читаємо координати з обрізаного файлу (трафарету)
    valid_coords = []
    in_atoms = False
    with open(template_file, 'r') as f:
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atoms = True
                continue
            elif line.startswith('@<TRIPOS>'):
                in_atoms = False
                
            if in_atoms and line.strip():
                parts = line.split()
                if len(parts) >= 5:
                    valid_coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                    
    # 2. Опрацьовуємо оригінальний файл (зберігаємо заряди та типи)
    with open(orig_file, 'r') as f:
        orig_lines = f.readlines()
        
    out_atoms = []
    out_bonds = []
    old_to_new_idx = {}
    
    total_charge = 0.0 
    
    current_section = None
    new_atom_idx = 1
    new_bond_idx = 1
    
    for i, line in enumerate(orig_lines):
        if line.startswith('@<TRIPOS>ATOM'):
            current_section = 'ATOM'
            continue
        elif line.startswith('@<TRIPOS>BOND'):
            current_section = 'BOND'
            continue
        elif line.startswith('@<TRIPOS>'):
            current_section = line.strip()
            continue
            
        if current_section == 'ATOM' and line.strip():
            parts = line.split()
            if len(parts) < 9:
                continue
                
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            
            # Зіставляємо координати (допуск 0.005 Ангстрем)
            keep = False
            for vc in valid_coords:
                if calc_dist((x, y, z), vc) < 0.005:
                    keep = True
                    break
            
            if keep:
                old_idx = int(parts[0])
                old_to_new_idx[old_idx] = new_atom_idx
                
                charge = float(parts[8])
                total_charge += charge
                
                # Зберігаємо формати за допомогою універсального .format()
                new_line = "{:>7} {:<7} {:>9} {:>9} {:>9} {:<7} {:>5} {:<7} {:>10}\n".format(
                    new_atom_idx, parts[1], parts[2], parts[3], parts[4], parts[5], parts[6], parts[7], parts[8]
                )
                out_atoms.append(new_line)
                new_atom_idx += 1
                
        elif current_section == 'BOND' and line.strip():
            parts = line.split()
            if len(parts) < 4:
                continue
            a1 = int(parts[1])
            a2 = int(parts[2])
            
            # Переносимо зв'язки з новою нумерацією
            if a1 in old_to_new_idx and a2 in old_to_new_idx:
                new_a1 = old_to_new_idx[a1]
                new_a2 = old_to_new_idx[a2]
                new_line = "{:>6} {:>5} {:>5} {:>4}\n".format(
                    new_bond_idx, new_a1, new_a2, parts[3]
                )
                out_bonds.append(new_line)
                new_bond_idx += 1

    # Генеруємо ім'я файлу із зарядом
    name, ext = os.path.splitext(base_out_file)
    if not ext:
        ext = ".mol2"
    charge_str = "{:+.4f}".format(total_charge)
    final_out_file = "{}_Q{}{}".format(name, charge_str, ext)

    # 3. Записуємо результат
    with open(final_out_file, 'w') as f:
        for i, line in enumerate(orig_lines):
            if line.startswith('@<TRIPOS>ATOM'):
                break
            if i > 0 and orig_lines[i-2].startswith('@<TRIPOS>MOLECULE'):
                parts = line.split()
                if len(parts) >= 2:
                    parts[0] = str(len(out_atoms))
                    parts[1] = str(len(out_bonds))
                    f.write(" " + "  ".join(parts) + "\n")
            else:
                f.write(line)
                
        f.write('@<TRIPOS>ATOM\n')
        for atom_line in out_atoms:
            f.write(atom_line)
            
        f.write('@<TRIPOS>BOND\n')
        for bond_line in out_bonds:
            f.write(bond_line)
            
        # Переносимо секцію SUBSTRUCTURE
        subst_start = -1
        for i, line in enumerate(orig_lines):
            if line.startswith('@<TRIPOS>SUBSTRUCTURE'):
                subst_start = i
                break
        if subst_start != -1:
            for line in orig_lines[subst_start:]:
                f.write(line)
                
    # Вивід результатів у термінал
    print("\n==========================================")
    print("DONE! Saved atoms: {}, bonds: {}".format(len(out_atoms), len(out_bonds)))
    print("TOTAL FRAGMENT CHARGE: {:+.5f}".format(total_charge))
    print("Saved file as: {}".format(final_out_file))
    print("==========================================\n")

if __name__ == '__main__':
    main()
