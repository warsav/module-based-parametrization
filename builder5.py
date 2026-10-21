# -*- coding: utf-8 -*-
import sys
import math
import os
import glob
import argparse
import re

def calc_dist(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def read_mol2(filename):
    atoms = []
    bonds = []
    with open(filename, 'r') as f:
        section = None
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                section = 'ATOM'; continue
            elif line.startswith('@<TRIPOS>BOND'):
                section = 'BOND'; continue
            elif line.startswith('@<TRIPOS>'):
                section = 'OTHER'; continue
                
            if section == 'ATOM' and line.strip():
                p = line.split()
                if len(p) >= 9:
                    atoms.append({
                        'id': int(p[0]), 'name': p[1],
                        'x': float(p[2]), 'y': float(p[3]), 'z': float(p[4]),
                        'type': p[5], 'res_id': p[6], 'res_name': p[7], 'charge': float(p[8])
                    })
            elif section == 'BOND' and line.strip():
                p = line.split()
                if len(p) >= 4:
                    bonds.append({'id': int(p[0]), 'a1': int(p[1]), 'a2': int(p[2]), 'type': p[3]})
    return atoms, bonds

def write_preview_mol2(filename, atoms, bonds):
    """
    Writes ONLY the temporary preview assembly for PyMOL.
    """
    with open(filename, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\nSYMMETRIC_ASSEMBLY_PREVIEW\n")
        f.write(f"{len(atoms)} {len(bonds)} 1 0 0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for idx, a in enumerate(atoms):
            f.write(f"{idx+1:>7} {a['name']:<7} {a['x']:>9.4f} {a['y']:>9.4f} {a['z']:>9.4f} {a['type']:<7} {a['res_id']:>5} {a['res_name']:<7} {a['charge']:>10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for idx, b in enumerate(bonds):
            f.write(f"{idx+1:>7} {b['a1']:>5} {b['a2']:>5} {b['type']:>4}\n")

def save_updated_mol2(original_file, new_file, updated_atoms):
    """
    SAFE PARSER: Reads the original file line-by-line, replaces ONLY the charge 
    value in the ATOM section, and strictly copies everything else.
    """
    charge_map = {a['id']: a['charge'] for a in updated_atoms}
    
    with open(original_file, 'r') as f_in, open(new_file, 'w') as f_out:
        section = None
        for line in f_in:
            if line.startswith('@<TRIPOS>'):
                section = line.strip()
                f_out.write(line)
                continue
            
            if section == '@<TRIPOS>ATOM' and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    atom_id = int(parts[0])
                    if atom_id in charge_map:
                        new_q = charge_map[atom_id]
                        a_id, a_name, a_x, a_y, a_z, a_type, a_resid, a_resname = parts[:8]
                        # Reconstruct the line using strict antechamber mol2 formatting
                        formatted_line = f"{int(a_id):>7} {a_name:<7} {float(a_x):>9.4f} {float(a_y):>9.4f} {float(a_z):>9.4f} {a_type:<7} {int(a_resid):>5} {a_resname:<7} {new_q:>10.6f}\n"
                        f_out.write(formatted_line)
                        continue
            
            f_out.write(line)

def verify_and_extract_dowry(clean_atoms, full_atoms, full_bonds, frag_num):
    c_to_f_map = {}
    f_to_c_map = {}
    
    for ca in clean_atoms:
        for fa in full_atoms:
            if calc_dist((ca['x'], ca['y'], ca['z']), (fa['x'], fa['y'], fa['z'])) < 0.001:
                c_to_f_map[ca['id']] = fa['id']
                f_to_c_map[fa['id']] = ca['id']
                break
                
    if len(c_to_f_map) != len(clean_atoms):
        raise ValueError(f"CRITICAL ERROR: Failed to map all clean atoms to full topology in Part {frag_num}!")

    all_full_ids = set(a['id'] for a in full_atoms)
    mapped_full_ids = set(f_to_c_map.keys())
    missing_by_coords = all_full_ids - mapped_full_ids

    cut_bonds = []
    for b in full_bonds:
        in_c1 = b['a1'] in f_to_c_map
        in_c2 = b['a2'] in f_to_c_map
        if in_c1 != in_c2:
            clean_side_id = b['a1'] if in_c1 else b['a2']
            buffer_side_id = b['a2'] if in_c1 else b['a1']
            cut_bonds.append((clean_side_id, buffer_side_id))

    adj_list = {a['id']: [] for a in full_atoms}
    for b in full_bonds:
        adj_list[b['a1']].append(b['a2'])
        adj_list[b['a2']].append(b['a1'])

    missing_by_bfs = set()
    dowry_charges = {}

    for clean_side_id, buffer_side_id in cut_bonds:
        visited = set()
        queue = [buffer_side_id]
        buffer_charge_sum = 0.0
        
        while queue:
            curr = queue.pop(0)
            if curr not in visited:
                visited.add(curr)
                atom = next(a for a in full_atoms if a['id'] == curr)
                buffer_charge_sum += atom['charge']
                
                for neighbor in adj_list[curr]:
                    if neighbor not in mapped_full_ids and neighbor not in visited:
                        queue.append(neighbor)
                        
        missing_by_bfs.update(visited)
        clean_id = f_to_c_map[clean_side_id]
        dowry_charges[clean_id] = buffer_charge_sum

    if missing_by_coords != missing_by_bfs:
        print(f"\n[!!!] VERIFICATION FAILED FOR PART {frag_num}")
        raise ValueError(f"Overlap identification methods mismatched for Part {frag_num}.")

    print(f"  [+] Part {frag_num}: Overlap double-verification PASSED. Found cuts: {len(cut_bonds)}.")
    return dowry_charges

def get_fragment_number(filename):
    """Витягує номер фрагмента з назви файлу типу lig_partX_..."""
    match = re.search(r'lig_part(\d+)', filename)
    if match:
        return int(match.group(1))
    return 999 # Fallback для сортування, якщо номер не знайдено

def main():
    parser = argparse.ArgumentParser(description="Symmetric Charge Equilibration Builder (Dynamic Fragments)")
    parser.add_argument("-c", "--clean", nargs="+", help="Explicit list of CLEAN .mol2 files in correct assembly order")
    parser.add_argument("-f", "--full", nargs="+", help="Explicit list of FULL (BCC) .mol2 files in correct assembly order")
    args = parser.parse_args()

    print("==================================================================")
    print("      SYMMETRIC CHARGE EQUILIBRATION BUILDER (ALGORITHM V2.2)     ")
    print("           [DYNAMIC FRAGMENTS & GROMACS OPTIMIZED]                ")
    print("==================================================================\n")

    clean_files_list = []
    full_files_list = []

    # Визначаємо файли: чи беремо з аргументів, чи шукаємо автоматично
    if args.clean and args.full:
        if len(args.clean) != len(args.full):
            print("CRITICAL ERROR: Number of clean files must match number of full files when using flags.")
            sys.exit(1)
        clean_files_list = args.clean
        full_files_list = args.full
        print(f"[*] Using explicitly provided files. Count: {len(clean_files_list)}")
    else:
        # Автоматичне визначення
        found_clean = glob.glob("lig_part*_clean_Q*.mol2")
        if not found_clean:
            print("CRITICAL ERROR: No 'lig_part*_clean_Q*.mol2' files found in directory.")
            sys.exit(1)
            
        # Сортуємо за номером фрагмента у назві
        clean_files_list = sorted(found_clean, key=get_fragment_number)
        
        for clean_file in clean_files_list:
            frag_num = get_fragment_number(clean_file)
            expected_full = f"lig_part{frag_num}_bcc_gaff2.mol2"
            if not os.path.exists(expected_full):
                print(f"CRITICAL ERROR: Found {clean_file} but missing matching full file {expected_full}")
                sys.exit(1)
            full_files_list.append(expected_full)
            
        print(f"[*] Auto-detected {len(clean_files_list)} fragments based on naming convention.")

    num_fragments = len(clean_files_list)
    if num_fragments < 2:
        print("CRITICAL ERROR: Need at least 2 fragments to assemble.")
        sys.exit(1)

    clean_frags = {}
    full_frags = {}
    fragment_dowries = {}
    
    print("\n[1/4] Analyzing files and executing double verification of overlaps...")
    for i in range(1, num_fragments + 1):
        clean_file = clean_files_list[i-1]
        full_file = full_files_list[i-1]
        
        c_atoms, c_bonds = read_mol2(clean_file)
        f_atoms, f_bonds = read_mol2(full_file)
        
        clean_frags[i] = {'orig_file': clean_file, 'atoms': c_atoms, 'bonds': c_bonds}
        
        dowries = verify_and_extract_dowry(c_atoms, f_atoms, f_bonds, i)
        fragment_dowries[i] = dowries

    print("\n[2/4] Predicting conjugation bonds and confirming structure...")
    choice_indices = [0] * (num_fragments - 1)
    confirmed_links = []
    
    while True:
        inter_bonds_preview = []
        atom_offset = 0
        all_atoms_final = []
        all_bonds_final = []
        current_links = []
        
        for i in range(1, num_fragments + 1):
            frag = clean_frags[i]
            new_atoms = [a.copy() for a in frag['atoms']]
            
            new_bonds = [{'a1': b['a1'] + atom_offset, 'a2': b['a2'] + atom_offset, 'type': b['type']} for b in frag['bonds']]
            
            for na in new_atoms: na['global_id'] = na['id'] + atom_offset
            all_atoms_final.extend(new_atoms)
            all_bonds_final.extend(new_bonds)
            frag['global_offset'] = atom_offset
            atom_offset += len(frag['atoms'])

        print("\nIdentified inter-fragment junctions:")
        for i in range(1, num_fragments):
            f1_atoms = clean_frags[i]['atoms']
            f2_atoms = clean_frags[i+1]['atoms']
            offset1 = clean_frags[i]['global_offset']
            offset2 = clean_frags[i+1]['global_offset']
            
            candidates = []
            for at1 in f1_atoms:
                for at2 in f2_atoms:
                    if at1['type'].lower().startswith('c') and at2['type'].lower().startswith('c'):
                        d = calc_dist((at1['x'], at1['y'], at1['z']), (at2['x'], at2['y'], at2['z']))
                        if 1.3 <= d <= 1.8:
                            candidates.append((at1, at2, d))
                            
            candidates.sort(key=lambda x: abs(x[2] - 1.54))
            idx = min(choice_indices[i-1], len(candidates) - 1)
            best_at1, best_at2, d = candidates[idx]
            
            inter_bonds_preview.append({'a1': best_at1['id'] + offset1, 'a2': best_at2['id'] + offset2, 'type': '1'})
            
            current_links.append({
                'part_L': i, 'atom_L_id': best_at1['id'], 'atom_L_name': best_at1['name'],
                'part_R': i+1, 'atom_R_id': best_at2['id'], 'atom_R_name': best_at2['name']
            })
            print(f"  Link {i}: Part {i} ({best_at1['name']}) === Part {i+1} ({best_at2['name']}) | dist: {d:.3f} Å")

        all_bonds_final.extend(inter_bonds_preview)
        write_preview_mol2('preview_assembly.mol2', all_atoms_final, all_bonds_final)

        ans = input("\n>>> Geometry saved to 'preview_assembly.mol2'. Is the conjugation correct? (y/n): ").strip().lower()
        if ans in ['y', 'yes']:
            confirmed_links = current_links
            break
        else:
            wrong = input(f"Which Link number is incorrect? (1-{num_fragments-1}): ").strip()
            try:
                wrong_idx = int(wrong) - 1
                if 0 <= wrong_idx < len(choice_indices):
                    choice_indices[wrong_idx] += 1
                else:
                    print("Invalid link number.")
            except:
                choice_indices = [c+1 for c in choice_indices]

    print("\n[3/4] Executing symmetric charge equilibration algorithm...")
    for link in confirmed_links:
        pL, pR = link['part_L'], link['part_R']
        aL_id, aR_id = link['atom_L_id'], link['atom_R_id']
        
        q_dowry_L = fragment_dowries[pL].get(aL_id, 0.0)
        q_dowry_R = fragment_dowries[pR].get(aR_id, 0.0)
        
        joint_added_charge = q_dowry_L + q_dowry_R
        half_charge = joint_added_charge / 2.0
        
        print(f"\n  > Link {pL}-{pR} [{link['atom_L_name']} -- {link['atom_R_name']}]")
        print(f"    Left dowry  (Part {pL}): {q_dowry_L:+.6f}")
        print(f"    Right dowry (Part {pR}): {q_dowry_R:+.6f}")
        print(f"    Joint dowry pool: {joint_added_charge:+.6f}")
        print(f"    => Adding exactly half (+ {half_charge:+.6f}) to each carbon.")
        
        atom_L = next(a for a in clean_frags[pL]['atoms'] if a['id'] == aL_id)
        atom_R = next(a for a in clean_frags[pR]['atoms'] if a['id'] == aR_id)
        
        atom_L['charge'] += half_charge
        atom_R['charge'] += half_charge

    print("\n[4/4] Saving fragments and generating GROMACS-optimized TLEAP script...")
    total_system_charge = 0.0
    
    for i in range(1, num_fragments + 1):
        frag = clean_frags[i]
        # Витягуємо номер з назви для збереження вихідного формату
        frag_num = get_fragment_number(frag['orig_file'])
        if frag_num == 999: 
            frag_num = i # Якщо файли задані вручну і не мають формату lig_partX
            
        out_name = f"lig_part{frag_num}_symmetric.mol2"
        save_updated_mol2(frag['orig_file'], out_name, frag['atoms'])
        net_q = sum(a['charge'] for a in frag['atoms'])
        total_system_charge += net_q
        print(f"  Saved {out_name:<25} | Net Charge = {net_q:+.4f}")

    print("\n==================================================================")
    print(f" EXPECTED TOTAL SYSTEM CHARGE: {total_system_charge:+.6f}")
    print("==================================================================")
    
    with open('build_poly.in', 'w') as f:
        f.write("# TLEAP Assembly Script (Symmetric Charges)\n")
        f.write("source leaprc.gaff2\n\n")
        
        # Завантаження frcmod файлів
        for i in range(1, num_fragments + 1):
            frag_num = get_fragment_number(clean_files_list[i-1])
            if frag_num == 999: frag_num = i
            frcmod_file = f"lig_part{frag_num}_AC.frcmod"
            if os.path.exists(frcmod_file):
                f.write(f"loadamberparams {frcmod_file}\n")
        f.write("\n")
        
        # Завантаження mol2
        for i in range(1, num_fragments + 1):
            frag_num = get_fragment_number(clean_files_list[i-1])
            if frag_num == 999: frag_num = i
            f.write(f"P{i} = loadmol2 lig_part{frag_num}_symmetric.mol2\n")
            
        f.write("\nlig = combine { " + " ".join([f"P{i}" for i in range(1, num_fragments + 1)]) + " }\n")
        
        for link in confirmed_links:
            f.write(f"bond lig.{link['part_L']}.{link['atom_L_name']} lig.{link['part_R']}.{link['atom_R_name']}\n")
            
        f.write("\nsavepdb lig lig.pdb\n")
        f.write("saveamberparm lig lig.prmtop lig.inpcrd\nquit\n")

    print("\nSuccess! Now run: tleap -f build_poly.in")
    print("Then run: acpype -p lig.prmtop -x lig.inpcrd -o gmx -d")

if __name__ == "__main__":
    main()
