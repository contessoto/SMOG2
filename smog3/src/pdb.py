import sys
import os
from .utils import smog_quit, trim, what_am_i
from .template import RESIDUES

def convert_pdb_to_gro_ndx(input_pdb, output_gro, output_g96, output_gro_scm, output_ndx, user_provided_map, free_coor, box_buffer=1.0, center_system=False):
    if not os.path.exists(input_pdb):
        smog_quit(f"Can't find file {input_pdb} !! ")

    chain_hash = {}
    chain_counter = 1
    res_counter = 1
    res_num_curr = "null"
    atom_info = []

    x_min = x_max = y_min = y_max = z_min = z_max = None

    with open(input_pdb, 'r') as f:
        lines = f.readlines()

    atom_num = 0
    last_ter_line = -1
    line_num = 0

    for line in lines:
        line_num += 1
        if line.startswith("END"):
            break

        if line.startswith("TER"):
            if last_ter_line != line_num:
                chain_counter += 1
                res_counter += 1
            last_ter_line = line_num
            continue

        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            res_name = line[17:21].strip()

            # Check if residue and atom defined in templates (RESIDUES global from template.py)
            if res_name not in RESIDUES or atom_name not in RESIDUES[res_name]['atoms']:
                continue

            atom_num += 1

            if free_coor:
                coords_string = line[30:].strip()
                coords = coords_string.split()
                if len(coords) < 3:
                     smog_quit(f"Not enough coordinates found on line: {line}")
                try:
                    x, y, z = float(coords[0]), float(coords[1]), float(coords[2])
                except ValueError:
                    smog_quit(f"Unable to interpret coordinate in PDB file with -freecoor. Problematic line: {line}")
            else:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    smog_quit(f"Unable to interpret coordinate in PDB file. Problematic line: {line}")

            # Convert to nm
            x *= 0.1
            y *= 0.1
            z *= 0.1

            if x_min is None:
                x_min = x_max = x
                y_min = y_max = y
                z_min = z_max = z
            else:
                x_max = max(x_max, x)
                x_min = min(x_min, x)
                y_max = max(y_max, y)
                y_min = min(y_min, y)
                z_max = max(z_max, z)
                z_min = min(z_min, z)

            chain = chain_counter

            res_num_pdb = line[22:26].strip()
            if res_num_curr == "null": res_num_curr = res_num_pdb

            if res_num_pdb != res_num_curr:
                res_num_curr = res_num_pdb
                if last_ter_line != line_num: # Simplified check, logic matches perl roughly
                     res_counter += 1

            if chain not in chain_hash: chain_hash[chain] = []

            # Modulo arithmetic for large numbers
            atom_num_out = atom_num % 100000
            res_counter_out = res_counter % 100000

            chain_hash[chain].append(atom_num_out)

            # Formats
            # Gro: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
            gro_line = f"{res_counter_out:5d}{res_name:<5s}{atom_name:>5s}{atom_num_out:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
            # SCM: %5d%-5s%5s%5d %10.4f %10.4f %10.4f
            gro_scm_line = f"{res_counter_out:5d}{res_name:<5s}{atom_name:>5s}{atom_num_out:5d} {x:10.4f} {y:10.4f} {z:10.4f}\n"

            atom_info.append({
                'x': x, 'y': y, 'z': z,
                'line': gro_line,
                'scm_line': gro_scm_line,
                'res_counter': res_counter, # Use unique counter for internal logic
                'res_name': res_name,
                'atom_name': atom_name,
                'atom_num': atom_num # Use unique atom number for internal logic
            })

            last_ter_line = -1 # Reset TER flag

    if atom_num == 0:
        smog_quit(f"No atoms found in {input_pdb}.")

    # Box dimensions
    x_range = (x_max - x_min) + box_buffer * 2
    y_range = (y_max - y_min) + box_buffer * 2
    z_range = (z_max - z_min) + box_buffer * 2

    if center_system:
        print("Will shift coordinates so that the molecule is centered in the box\n")
        shift_x = x_min - box_buffer
        shift_y = y_min - box_buffer
        shift_z = z_min - box_buffer

        for atom in atom_info:
            atom['x'] -= shift_x
            atom['y'] -= shift_y
            atom['z'] -= shift_z
            # Reformat lines with new coords
            atom['line'] = f"{atom['res_counter']:5d}{atom['res_name']:<5s}{atom['atom_name']:>5s}{atom['atom_num']:5d}{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}\n"
            atom['scm_line'] = f"{atom['res_counter']:5d}{atom['res_name']:<5s}{atom['atom_name']:>5s}{atom['atom_num']:5d} {atom['x']:10.4f} {atom['y']:10.4f} {atom['z']:10.4f}\n"

    # Write Gro
    if output_gro:
        with open(output_gro, 'w') as f:
            f.write(f"Gro file for a structure based model, generated with SMOG v3\n")
            f.write(f"{len(atom_info)}\n")
            for atom in atom_info:
                f.write(atom['line'])
            f.write(f"{x_range:10.5f} {y_range:10.5f} {z_range:10.5f}\n")
        print(f"{os.path.abspath(output_gro)} written")

    # Write Gro SCM
    if output_gro_scm:
        with open(output_gro_scm, 'w') as f:
            f.write("Temp Gro file with PDB precision for SCM calculations.\n")
            f.write(f"{len(atom_info)}\n")
            for atom in atom_info:
                f.write(atom['scm_line'])
            f.write(f"{x_range:10.5f} {y_range:10.5f} {z_range:10.5f}\n")

    # Write NDX
    if output_ndx:
        with open(output_ndx, 'w') as f:
            for chain_id in sorted(chain_hash.keys()):
                f.write(f"[ {chain_id} ]\n")
                f.write("\n".join(map(str, chain_hash[chain_id])))
                f.write("\n\n")

    return atom_info
