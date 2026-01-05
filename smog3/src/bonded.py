import numpy as np
import sys
import re
from .config import VERSION
from .utils import norm, dot, cross, sub, add, mult, eval_sub
from .template import (
    RESIDUES, INTERACTIONS, BOND_FUNCTIONALS, DIHEDRAL_FUNCTIONALS,
    ANGLE_FUNCTIONALS, TERM_RATIOS, BOND_TYPES_USED, PAIR_TYPES_USED,
    NB_TYPES_PRESENT, INTERACTION_THRESHOLD, ENERGY_GROUPS
)

BOND_VECTORS = {}
FREE_ANGLES = {}
EXCLUDE_BONDED = {}

def norm(vec):
    return np.linalg.norm(vec)

def inner(vec1, vec2):
    return np.dot(vec1, vec2)

def crossp(vec1, vec2):
    return np.cross(vec1, vec2)

def evaluate_param(p, value):
    # p is string from template, value is the measured geometry (r, theta)
    if '?' in p:
        p = p.replace('?', str(value))

    # Try direct conversion
    try:
        return float(p)
    except:
        pass

    # Use safe eval with context
    # value might be needed as variable? usually '?' handled above.
    # We can pass typical math functions.
    val = eval_sub(p, {})
    if val is not None:
        return float(val)
    return p

def calculate_bonds(atom_info, index_handle):
    bond_cache = []

    # 1. Intra-residue bonds
    # Iterate over residues present in PDB (mapped by atom_info)
    # atom_info is a list of dicts. We need to group by residue.

    # Group atoms by residue (res_counter)
    residues_in_pdb = {}
    for atom in atom_info:
        res_id = atom['res_counter']
        if res_id not in residues_in_pdb:
            residues_in_pdb[res_id] = {'name': atom['res_name'], 'atoms': {}}
        residues_in_pdb[res_id]['atoms'][atom['atom_name']] = atom

    # 1. Intra-residue bonds
    for res_id in sorted(residues_in_pdb.keys()):
        res_name = residues_in_pdb[res_id]['name']
        res_atoms = residues_in_pdb[res_id]['atoms']

        if res_name not in RESIDUES: continue

        # Get bonds from template
        template_res = RESIDUES[res_name]

        for bond_pair, bond_data in template_res['bonds'].items():
            # bond_pair is "atomA-atomB"
            atom_a_name, atom_b_name = bond_pair.split('-')

            if atom_a_name in res_atoms and atom_b_name in res_atoms:
                atom_a = res_atoms[atom_a_name]
                atom_b = res_atoms[atom_b_name]

                vec_a = np.array([atom_a['x'], atom_a['y'], atom_a['z']])
                vec_b = np.array([atom_b['x'], atom_b['y'], atom_b['z']])

                vec = vec_a - vec_b
                dist = norm(vec) # nm

                # Store vector for angles/dihedrals
                idx_a = atom_a['atom_num']
                idx_b = atom_b['atom_num']
                BOND_VECTORS[f"{idx_a}-{idx_b}"] = vec / dist if dist > 0 else vec
                BOND_VECTORS[f"{idx_b}-{idx_a}"] = -vec / dist if dist > 0 else -vec

                # Determine function type
                # In Perl: logic matches bTypes to find interaction function
                b_type_a = template_res['atoms'][atom_a_name]['bType']
                b_type_b = template_res['atoms'][atom_b_name]['bType']

                func_str = get_bond_function(b_type_a, b_type_b)
                if func_str:
                    val = bond_output(func_str, idx_a, idx_b, dist, index_handle)
                    if val:
                        bond_cache.append({'i': idx_a, 'j': idx_b, 'v': val})

    # 2. Inter-residue bonds (connections)
    # SMOG 2 iterates over residues and checks connections defined in .bif
    # CONNECTIONS global from template.py
    # This requires identifying connected residues. In PDB, usually sequential residues in same chain.

    # We iterate residue i and i+1
    sorted_res_ids = sorted(residues_in_pdb.keys())
    for i in range(len(sorted_res_ids) - 1):
        res_id1 = sorted_res_ids[i]
        res_id2 = sorted_res_ids[i+1]

        # Check if they are in same chain?
        # atom_info list is ordered. If TER was processed, they might be separate.
        # But we don't have explicit chain info in residues_in_pdb easily unless we stored it.
        # convert_pdb_to_gro_ndx tracked chains. But didn't store in atom_info dicts explicitly beyond line context.
        # Wait, PDB format has chain ID but we ignored it in convert_pdb_to_gro_ndx logic (just chain_counter).
        # We can infer connectivity if atom indices are somewhat contiguous or based on distance?
        # SMOG 2 assumes connectivity between sequential residues in PDB unless TER.

        # We need to know if res_id1 and res_id2 are connected.
        # In convert_pdb_to_gro_ndx, we handled TER.
        # If we assume res_counter increments are continuous within a chain.
        # But convert_pdb_to_gro_ndx increments res_counter.

        # Let's assume sequential residues in the list are connected if they are chemically connected.
        # SMOG 2 checks CONNECTIONS defined in BIF.

        res1 = residues_in_pdb[res_id1]
        res2 = residues_in_pdb[res_id2]

        type1 = RESIDUES[res1['name']]['residueType']
        type2 = RESIDUES[res2['name']]['residueType']

        from .template import CONNECTIONS

        if type1 in CONNECTIONS and type2 in CONNECTIONS[type1]:
            connection = CONNECTIONS[type1][type2]
            # bond list
            for bond in connection.findall('bond'):
                atom_names = [a.text.strip() for a in bond.findall('atom')]
                if len(atom_names) == 2:
                    name_a = atom_names[0]
                    name_b = atom_names[1]

                    # Logic: Atoms in connections are usually defined as "C" (from first) and "N" (from second).
                    # But they could be any atoms.
                    # In Perl: uses indices relative to current residue or next?
                    # The BIF usually assumes "atom" refers to current residue, unless specified?
                    # Actually standard connections are like C of i and N of i+1.
                    # BIF format:
                    # <connection ...> <bond><atom>C</atom><atom>N</atom></bond> ...
                    # SMOG 2 logic matches names in res1 and res2?
                    # Actually, usually specific atoms are named.

                    # If atom name ends in +, it means next residue? Or - for prev?
                    # SMOG 2 BIF doesn't seem to use +/- for connections explicitly in name for standard peptide bond C-N?
                    # It just lists atom names.
                    # Wait, look at BIF example:
                    # <connection name="amino-amino" residueType1="amino" residueType2="amino">
                    #   <bond energyGroup="r_a"><atom>C</atom><atom>N</atom></bond>

                    # This implies atom C from residueType1 and atom N from residueType2.

                    if name_a in res1['atoms'] and name_b in res2['atoms']:
                        atom_a = res1['atoms'][name_a]
                        atom_b = res2['atoms'][name_b]

                        idx_a = atom_a['atom_num']
                        idx_b = atom_b['atom_num']

                        # Distance check?
                        vec_a = np.array([atom_a['x'], atom_a['y'], atom_a['z']])
                        vec_b = np.array([atom_b['x'], atom_b['y'], atom_b['z']])
                        vec = vec_a - vec_b
                        dist = norm(vec)

                        if dist > 0.25: # Arbitrary cutoff to avoid connecting across breaks/TER if logic failed?
                            # If dist is large, maybe not connected.
                            # SMOG 2 relies on user inputs (TER).
                            # If we assume sequential in atom_info list implies connection unless TER.
                            # We can check if `atom_info` indices are adjacent? No, residues.
                            pass

                        # We will process it.
                        BOND_VECTORS[f"{idx_a}-{idx_b}"] = vec / dist if dist > 0 else vec
                        BOND_VECTORS[f"{idx_b}-{idx_a}"] = -vec / dist if dist > 0 else -vec

                        # Get function
                        b_type_a = RESIDUES[res1['name']]['atoms'][name_a]['bType']
                        b_type_b = RESIDUES[res2['name']]['atoms'][name_b]['bType']

                        func_str = get_bond_function(b_type_a, b_type_b)
                        if func_str:
                            val = bond_output(func_str, idx_a, idx_b, dist, index_handle)
                            if val:
                                bond_cache.append({'i': idx_a, 'j': idx_b, 'v': val})

    return bond_cache

def get_bond_function(type_a, type_b):
    if type_a in INTERACTIONS['bonds'] and type_b in INTERACTIONS['bonds'][type_a]:
        return INTERACTIONS['bonds'][type_a][type_b]

    # Wildcards
    if type_a in INTERACTIONS['bonds'] and '*' in INTERACTIONS['bonds'][type_a]:
        return INTERACTIONS['bonds'][type_a]['*']

    if type_b in INTERACTIONS['bonds'] and '*' in INTERACTIONS['bonds'][type_b]:
        return INTERACTIONS['bonds'][type_b]['*']

    if '*' in INTERACTIONS['bonds'] and '*' in INTERACTIONS['bonds']['*']:
        return INTERACTIONS['bonds']['*']['*']

    return None

def bond_output(input_func, i, j, value, index_handle):
    # Parses func string like "bond_harmonic(?,10000)" and returns formatted line
    # Simplified parsing
    if not input_func: return None

    funcs = input_func.split('+')
    output_string = ""

    for f in funcs:
        f = f.strip()
        f_name, params = parse_func_str(f)

        # Get gromacs function type
        # Hardcoded map for now based on SMOGglobals
        f_type_map = {
            'bond_harmonic': 1,
            'bond_type6': 6
        }
        f_type = f_type_map.get(f_name, 1)

        # Evaluate params
        # param[0] is usually distance/angle
        # If param[0] contains ?, replace with value

        parsed_params = []
        for idx, p in enumerate(params):
            parsed_params.append(evaluate_param(p, value))

        # Format
        # i j func_type p1 p2 ...
        param_str = "\t".join([f"{p:12.9e}" if isinstance(p, (float, int)) else str(p) for p in parsed_params])
        output_string += f"{i}\t{j}\t{f_type}\t{param_str}\n"

    return output_string

def parse_func_str(func_str):
    # extract name(args)
    m = re.match(r'([^(]+)\((.*)\)', func_str)
    if m:
        name = m.group(1).strip()
        args = [a.strip() for a in m.group(2).split(',')]
        return name, args
    return func_str, []

def calculate_angles(atom_info, index_handle):
    adj_list = {}
    for bond_key in BOND_VECTORS.keys():
        idx_a, idx_b = map(int, bond_key.split('-'))
        if idx_a not in adj_list: adj_list[idx_a] = []
        adj_list[idx_a].append(idx_b)

    angle_cache = []

    # Map index to atom type
    # atom_info is list of dicts
    idx_to_info = {a['atom_num']: a for a in atom_info}

    for j, neighbors in adj_list.items():
        if len(neighbors) < 2: continue

        for idx1 in range(len(neighbors)):
            for idx2 in range(idx1 + 1, len(neighbors)):
                i = neighbors[idx1]
                k = neighbors[idx2]

                # Sort i and k to ensure canonical order?
                # SMOG 2 checks both orders i-j-k and k-j-i

                # Get bTypes
                if i not in idx_to_info or j not in idx_to_info or k not in idx_to_info: continue

                info_i = idx_to_info[i]
                info_j = idx_to_info[j]
                info_k = idx_to_info[k]

                res_i = RESIDUES[info_i['res_name']]
                res_j = RESIDUES[info_j['res_name']]
                res_k = RESIDUES[info_k['res_name']]

                type_i = res_i['atoms'][info_i['atom_name']]['bType']
                type_j = res_j['atoms'][info_j['atom_name']]['bType']
                type_k = res_k['atoms'][info_k['atom_name']]['bType']

                func_str = get_angle_function(type_i, type_j, type_k)
                if not func_str: continue

                # Calculate angle
                vec_ji = BOND_VECTORS[f"{j}-{i}"] # Vector from j to i
                vec_jk = BOND_VECTORS[f"{j}-{k}"] # Vector from j to k

                # Angle between ji and jk
                # cos_theta = dot(ji, jk) / (|ji| |jk|)
                # BOND_VECTORS stores unit vectors already?
                # "BOND_VECTORS[f"{idx_a}-{idx_b}"] = vec / dist if dist > 0 else vec"
                # Yes, normalized.

                cos_theta = inner(vec_ji, vec_jk)
                # Clamp to -1, 1
                cos_theta = max(-1.0, min(1.0, cos_theta))
                theta_rad = np.arccos(cos_theta)
                theta_deg = np.degrees(theta_rad)

                val = angle_output(func_str, i, j, k, theta_deg)
                if val:
                    angle_cache.append({'i': i, 'j': j, 'k': k, 'v': val})

    return angle_cache

def get_angle_function(type_a, type_b, type_c):
    key = f"{type_a}-{type_b}-{type_c}"
    if key in INTERACTIONS['angles']: return INTERACTIONS['angles'][key]
    key_rev = f"{type_c}-{type_b}-{type_a}"
    if key_rev in INTERACTIONS['angles']: return INTERACTIONS['angles'][key_rev]

    # Wildcards matching logic (simplified)
    # SMOG 2 checks permutations and wildcards.
    # Check A-B-*
    # Check *-B-C
    # Check *-B-*
    # Center atom type B must match or be *? No, center usually explicit.
    # Order matters.

    # Simple check for now
    for t_a in [type_a, '*']:
        for t_b in [type_b, '*']:
            for t_c in [type_c, '*']:
                k = f"{t_a}-{t_b}-{t_c}"
                if k in INTERACTIONS['angles']: return INTERACTIONS['angles'][k]
                k_rev = f"{t_c}-{t_b}-{t_a}"
                if k_rev in INTERACTIONS['angles']: return INTERACTIONS['angles'][k_rev]
    return None

def angle_output(input_func, i, j, k, value):
    if not input_func: return None
    funcs = input_func.split('+')
    output_string = ""
    for f in funcs:
        f = f.strip()
        f_name, params = parse_func_str(f)

        f_type_map = {
            'angle_harmonic': 1,
        }
        f_type = f_type_map.get(f_name, 1)

        parsed_params = []
        for p in params:
            parsed_params.append(evaluate_param(p, value))

        param_str = "\t".join([f"{p:12.9e}" if isinstance(p, (float, int)) else str(p) for p in parsed_params])
        output_string += f"{i}\t{j}\t{k}\t{f_type}\t{param_str}\n"
    return output_string

def calculate_dihedrals(atom_info, index_handle):
    adj_list = {}
    for bond_key in BOND_VECTORS.keys():
        idx_a, idx_b = map(int, bond_key.split('-'))
        if idx_a not in adj_list: adj_list[idx_a] = []
        adj_list[idx_a].append(idx_b)

    dihe_cache = []
    idx_to_info = {a['atom_num']: a for a in atom_info}

    # Traverse graph for 4 atoms i-j-k-l
    for j, neighbors_j in adj_list.items():
        for k in neighbors_j:
            if k <= j: continue # Avoid double counting bonds, process j-k once

            # j-k is central bond
            # Find i bonded to j (i != k)
            # Find l bonded to k (l != j)

            neighbors_i = [n for n in adj_list.get(j, []) if n != k]
            neighbors_l = [n for n in adj_list.get(k, []) if n != j]

            for i in neighbors_i:
                for l in neighbors_l:
                    # Dihedral i-j-k-l

                    if i not in idx_to_info or j not in idx_to_info or k not in idx_to_info or l not in idx_to_info: continue

                    info_i, info_j, info_k, info_l = idx_to_info[i], idx_to_info[j], idx_to_info[k], idx_to_info[l]

                    res_i = RESIDUES[info_i['res_name']]
                    res_j = RESIDUES[info_j['res_name']]
                    res_k = RESIDUES[info_k['res_name']]
                    res_l = RESIDUES[info_l['res_name']]

                    type_i = res_i['atoms'][info_i['atom_name']]['bType']
                    type_j = res_j['atoms'][info_j['atom_name']]['bType']
                    type_k = res_k['atoms'][info_k['atom_name']]['bType']
                    type_l = res_l['atoms'][info_l['atom_name']]['bType']

                    # Need energy group for the central bond j-k
                    # Energy group is defined in residue or connection
                    # Assuming j and k are in same residue or connected

                    eg = get_energy_group(info_j, info_k)

                    func_str = get_dihedral_function(type_i, type_j, type_k, type_l, eg)
                    if not func_str: continue

                    # Calculate angle
                    b1 = BOND_VECTORS[f"{j}-{i}"] # j->i (Wait, definition is usually b1=j-i, b2=k-j, b3=l-k?)
                    # Standard definition: b1 = r_j - r_i = - (r_i - r_j) = - vec_ji
                    # Actually standard: b1 = i->j, b2 = j->k, b3 = k->l

                    vec_ij = -BOND_VECTORS[f"{j}-{i}"]
                    vec_jk = BOND_VECTORS[f"{j}-{k}"] # same as k-j reversed? No, j->k
                    vec_kl = -BOND_VECTORS[f"{l}-{k}"] # k->l

                    # Using normals
                    n1 = np.cross(vec_ij, vec_jk)
                    n2 = np.cross(vec_jk, vec_kl)

                    norm_n1 = norm(n1)
                    norm_n2 = norm(n2)

                    if norm_n1 < 1e-6 or norm_n2 < 1e-6: continue # Linear

                    n1 /= norm_n1
                    n2 /= norm_n2

                    m1 = np.cross(n1, vec_jk / norm(vec_jk))
                    x = np.dot(n1, n2)
                    y = np.dot(m1, n2)

                    theta_rad = np.arctan2(y, x)
                    theta_deg = np.degrees(theta_rad) # -180 to 180

                    # SMOG convention for cosine dihedrals often includes a shift
                    # But output formatter handles it.

                    val = dihedral_output(func_str, i, j, k, l, theta_deg)
                    if val:
                        dihe_cache.append({'i': i, 'j': j, 'k': k, 'l': l, 'v': val})

    return dihe_cache

def get_energy_group(info_j, info_k):
    res_j_name = info_j['res_name']
    res_k_name = info_k['res_name']
    atom_j_name = info_j['atom_name']
    atom_k_name = info_k['atom_name']

    # Same residue
    if res_j_name == res_k_name:
        bond_key = f"{atom_j_name}-{atom_k_name}"
        if bond_key in RESIDUES[res_j_name]['energyGroups']:
            return RESIDUES[res_j_name]['energyGroups'][bond_key]
        bond_key_rev = f"{atom_k_name}-{atom_j_name}"
        if bond_key_rev in RESIDUES[res_j_name]['energyGroups']:
            return RESIDUES[res_j_name]['energyGroups'][bond_key_rev]
    else:
        # Connected residues
        type_j = RESIDUES[res_j_name]['residueType']
        type_k = RESIDUES[res_k_name]['residueType']

        from .template import CONNECTIONS
        if type_j in CONNECTIONS and type_k in CONNECTIONS[type_j]:
            # Assuming standard connection, get first bond's EG?
            # In SMOG2 it checks matching atom names?
            # Simplified: Use first bond EG in connection
            for bond in CONNECTIONS[type_j][type_k].findall('bond'):
                return bond.get('energyGroup')

    return None

def get_dihedral_function(type_a, type_b, type_c, type_d, eg):
    if not eg: return None
    if eg not in INTERACTIONS['dihedrals']: return None

    # Check permutations
    patterns = [
        [type_a, type_b, type_c, type_d],
        [type_d, type_c, type_b, type_a]
    ]

    for pat in patterns:
        ta, tb, tc, td = pat
        # Wildcard match (simplified)
        # Check explicit
        key = f"{ta}-{tb}-{tc}-{td}"
        if key in INTERACTIONS['dihedrals'][eg]: return INTERACTIONS['dihedrals'][eg][key]

        # Check wildcards
        # This is where SMOG 2 scoring comes in.
        # For now, check full wildcards *
        key_wild = r"*-\*-\*-\*" # SMOG2 uses *
        if '*-*-*-*' in INTERACTIONS['dihedrals'][eg]: return INTERACTIONS['dihedrals'][eg]['*-*-*-*']

        # Intermediate wildcards
        for k in INTERACTIONS['dihedrals'][eg].keys():
            # Check match
            parts = k.split('-')
            if len(parts) != 4: continue
            match = True
            for idx, part in enumerate(parts):
                if part != '*' and part != pat[idx]:
                    match = False
                    break
            if match: return INTERACTIONS['dihedrals'][eg][k]

    return None

def dihedral_output(input_func, i, j, k, l, value):
    if not input_func: return None
    funcs = input_func.split('+')
    output_string = ""
    for f in funcs:
        f = f.strip()
        f_name, params = parse_func_str(f)

        f_type_map = {
            'dihedral_cosine': 1,
            'dihedral_harmonic': 2,
            'dihedral_ncos': 1
        }
        f_type = f_type_map.get(f_name, 1) # Default to 1

        parsed_params = []
        # value is in degrees

        # For cosine: theta0, k, mult
        # For ncos: theta0, k, mult
        # For harmonic: theta0, k

        # Need to handle substitutions and unit conversions
        # SMOG 2: param[0] += 180 for ncos/cosine

        val_local = value

        for idx, p in enumerate(params):
            parsed_params.append(evaluate_param(p, val_local))

        if f_name in ['dihedral_cosine', 'dihedral_ncos']:
            # Adjust theta0 (param 0)
            # param[0] += 180
            parsed_params[0] += 180.0

        param_str = "\t".join([f"{p:12.9e}" if isinstance(p, (float, int)) else str(p) for p in parsed_params])
        output_string += f"{i}\t{j}\t{k}\t{l}\t{f_type}\t{param_str}\n"

    return output_string
