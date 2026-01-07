import numpy as np
import sys
import re
from .config import VERSION
from .utils import norm, dot, cross, sub, add, mult, eval_sub
from .template import (
    RESIDUES, INTERACTIONS, BOND_FUNCTIONALS, DIHEDRAL_FUNCTIONALS,
    ANGLE_FUNCTIONALS, TERM_RATIOS, BOND_TYPES_USED, PAIR_TYPES_USED,
    NB_TYPES_PRESENT, INTERACTION_THRESHOLD, ENERGY_GROUPS, COUNT_DIHEDRALS, NORMALIZE_VALS
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
    val = eval_sub(p, {})
    if val is not None:
        return float(val)
    return p

def calculate_bonds(atom_info, index_handle):
    bond_cache = []

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

        template_res = RESIDUES[res_name]

        for bond_pair, bond_data in template_res['bonds'].items():
            atom_a_name, atom_b_name = bond_pair.split('-')

            if atom_a_name in res_atoms and atom_b_name in res_atoms:
                atom_a = res_atoms[atom_a_name]
                atom_b = res_atoms[atom_b_name]

                vec_a = np.array([atom_a['x'], atom_a['y'], atom_a['z']])
                vec_b = np.array([atom_b['x'], atom_b['y'], atom_b['z']])

                vec = vec_a - vec_b
                dist = norm(vec) # nm

                idx_a = atom_a['atom_num']
                idx_b = atom_b['atom_num']
                BOND_VECTORS[f"{idx_a}-{idx_b}"] = vec / dist if dist > 0 else vec
                BOND_VECTORS[f"{idx_b}-{idx_a}"] = -vec / dist if dist > 0 else -vec

                b_type_a = template_res['atoms'][atom_a_name]['bType']
                b_type_b = template_res['atoms'][atom_b_name]['bType']

                func_str = get_bond_function(b_type_a, b_type_b)
                if func_str:
                    val = bond_output(func_str, idx_a, idx_b, dist, index_handle)
                    if val:
                        bond_cache.append({'i': idx_a, 'j': idx_b, 'v': val})

    # 2. Inter-residue bonds (connections)
    sorted_res_ids = sorted(residues_in_pdb.keys())
    for i in range(len(sorted_res_ids) - 1):
        res_id1 = sorted_res_ids[i]
        res_id2 = sorted_res_ids[i+1]

        res1 = residues_in_pdb[res_id1]
        res2 = residues_in_pdb[res_id2]

        type1 = RESIDUES[res1['name']]['residueType']
        type2 = RESIDUES[res2['name']]['residueType']

        from .template import CONNECTIONS

        if type1 in CONNECTIONS and type2 in CONNECTIONS[type1]:
            connection = CONNECTIONS[type1][type2]
            for bond in connection.findall('bond'):
                atom_names = [a.text.strip() for a in bond.findall('atom')]
                if len(atom_names) == 2:
                    name_a = atom_names[0]
                    name_b = atom_names[1]

                    if name_a in res1['atoms'] and name_b in res2['atoms']:
                        atom_a = res1['atoms'][name_a]
                        atom_b = res2['atoms'][name_b]

                        idx_a = atom_a['atom_num']
                        idx_b = atom_b['atom_num']

                        vec_a = np.array([atom_a['x'], atom_a['y'], atom_a['z']])
                        vec_b = np.array([atom_b['x'], atom_b['y'], atom_b['z']])
                        vec = vec_a - vec_b
                        dist = norm(vec)

                        if dist > 0.25: pass # Skip if too far?

                        BOND_VECTORS[f"{idx_a}-{idx_b}"] = vec / dist if dist > 0 else vec
                        BOND_VECTORS[f"{idx_b}-{idx_a}"] = -vec / dist if dist > 0 else -vec

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
    if type_a in INTERACTIONS['bonds'] and '*' in INTERACTIONS['bonds'][type_a]:
        return INTERACTIONS['bonds'][type_a]['*']
    if type_b in INTERACTIONS['bonds'] and '*' in INTERACTIONS['bonds'][type_b]:
        return INTERACTIONS['bonds'][type_b]['*']
    if '*' in INTERACTIONS['bonds'] and '*' in INTERACTIONS['bonds']['*']:
        return INTERACTIONS['bonds']['*']['*']
    return None

def bond_output(input_func, i, j, value, index_handle):
    if not input_func: return None
    funcs = input_func.split('+')
    output_string = ""
    for f in funcs:
        f = f.strip()
        f_name, params = parse_func_str(f)
        f_type_map = {'bond_harmonic': 1, 'bond_type6': 6}
        f_type = f_type_map.get(f_name, 1)
        parsed_params = []
        for idx, p in enumerate(params):
            parsed_params.append(evaluate_param(p, value))
        param_str = "\t".join([f"{p:12.9e}" if isinstance(p, (float, int)) else str(p) for p in parsed_params])
        output_string += f"{i}\t{j}\t{f_type}\t{param_str}\n"
    return output_string

def parse_func_str(func_str):
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
    idx_to_info = {a['atom_num']: a for a in atom_info}

    for j, neighbors in adj_list.items():
        if len(neighbors) < 2: continue
        for idx1 in range(len(neighbors)):
            for idx2 in range(idx1 + 1, len(neighbors)):
                i = neighbors[idx1]
                k = neighbors[idx2]
                if i not in idx_to_info or j not in idx_to_info or k not in idx_to_info: continue

                info_i, info_j, info_k = idx_to_info[i], idx_to_info[j], idx_to_info[k]
                res_i, res_j, res_k = RESIDUES[info_i['res_name']], RESIDUES[info_j['res_name']], RESIDUES[info_k['res_name']]
                type_i = res_i['atoms'][info_i['atom_name']]['bType']
                type_j = res_j['atoms'][info_j['atom_name']]['bType']
                type_k = res_k['atoms'][info_k['atom_name']]['bType']

                func_str = get_angle_function(type_i, type_j, type_k)
                if not func_str: continue

                vec_ji = BOND_VECTORS[f"{j}-{i}"]
                vec_jk = BOND_VECTORS[f"{j}-{k}"]
                cos_theta = inner(vec_ji, vec_jk)
                cos_theta = max(-1.0, min(1.0, cos_theta))
                theta_deg = np.degrees(np.arccos(cos_theta))

                val = angle_output(func_str, i, j, k, theta_deg)
                if val:
                    angle_cache.append({'i': i, 'j': j, 'k': k, 'v': val})
    return angle_cache

def get_angle_function(type_a, type_b, type_c):
    key = f"{type_a}-{type_b}-{type_c}"
    if key in INTERACTIONS['angles']: return INTERACTIONS['angles'][key]
    key_rev = f"{type_c}-{type_b}-{type_a}"
    if key_rev in INTERACTIONS['angles']: return INTERACTIONS['angles'][key_rev]
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
        f_type = {'angle_harmonic': 1}.get(f_name, 1)
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
    dihedral_counts = {} # Key: central bond "min-max"

    # Collect dihedrals first to count and normalize
    raw_dihedrals = []

    for j, neighbors_j in adj_list.items():
        for k in neighbors_j:
            if k <= j: continue

            neighbors_i = [n for n in adj_list.get(j, []) if n != k]
            neighbors_l = [n for n in adj_list.get(k, []) if n != j]

            bond_key = f"{j}-{k}"

            for i in neighbors_i:
                for l in neighbors_l:
                    if i not in idx_to_info or j not in idx_to_info or k not in idx_to_info or l not in idx_to_info: continue

                    info_i, info_j, info_k, info_l = idx_to_info[i], idx_to_info[j], idx_to_info[k], idx_to_info[l]
                    res_i, res_j, res_k, res_l = RESIDUES[info_i['res_name']], RESIDUES[info_j['res_name']], RESIDUES[info_k['res_name']], RESIDUES[info_l['res_name']]
                    type_i, type_j, type_k, type_l = res_i['atoms'][info_i['atom_name']]['bType'], res_j['atoms'][info_j['atom_name']]['bType'], res_k['atoms'][info_k['atom_name']]['bType'], res_l['atoms'][info_l['atom_name']]['bType']

                    eg = get_energy_group(info_j, info_k)
                    func_str = get_dihedral_function(type_i, type_j, type_k, type_l, eg)
                    if not func_str: continue

                    vec_ij = -BOND_VECTORS[f"{j}-{i}"]
                    vec_jk = BOND_VECTORS[f"{j}-{k}"]
                    vec_kl = -BOND_VECTORS[f"{l}-{k}"]

                    n1 = np.cross(vec_ij, vec_jk)
                    n2 = np.cross(vec_jk, vec_kl)
                    if norm(n1) < 1e-6 or norm(n2) < 1e-6: continue

                    n1 /= norm(n1)
                    n2 /= norm(n2)
                    m1 = np.cross(n1, vec_jk / norm(vec_jk))
                    theta_deg = np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))

                    raw_dihedrals.append({
                        'i': i, 'j': j, 'k': k, 'l': l,
                        'theta': theta_deg,
                        'func': func_str,
                        'eg': eg,
                        'res_type': res_j['residueType'] # Approx
                    })

                    if COUNT_DIHEDRALS:
                        if bond_key not in dihedral_counts: dihedral_counts[bond_key] = 0
                        dihedral_counts[bond_key] += 1

    # Process and Normalize
    for d in raw_dihedrals:
        # Determine scaling
        scaling = 1.0

        # Dihedral counting normalization
        bond_key = f"{d['j']}-{d['k']}"
        count = dihedral_counts.get(bond_key, 1)
        if count > 0 and COUNT_DIHEDRALS:
            scaling /= count

        # Intra-relative strength from SIF
        # TERM_RATIOS[res_type]['energyGroup'][eg_name]
        # Need res_type of the energy group context?
        # get_energy_group logic uses res_name.
        # We stored approx res_type.
        # But EG definition in SIF is per residue type.

        # We need the residue type that DEFINES the EG?
        # Or the residue type of the bond?
        # Usually it's the residue containing the bond.
        # If bond connects residues, it's a bit ambiguous which resType to use for SIF lookup if they differ.
        # But usually EG name is unique enough or consistent.

        eg_info = None
        if d['res_type'] in TERM_RATIOS and 'energyGroup' in TERM_RATIOS[d['res_type']]:
             if d['eg'] in TERM_RATIOS[d['res_type']]['energyGroup']:
                 eg_info = TERM_RATIOS[d['res_type']]['energyGroup'][d['eg']]

        if eg_info:
            if 'intraRelativeStrength' in eg_info and eg_info['intraRelativeStrength']:
                try:
                    scaling *= float(eg_info['intraRelativeStrength'])
                except:
                    pass # Variable?

            # Normalize flag in EG?
            # "normalize": 1 means apply normalization?
            # If COUNT_DIHEDRALS is set, we normalized by count.
            # SMOG 2: setRatios uses normalization if NORMALIZE_VALS is on.
            # And specific EG normalize flag.

            if NORMALIZE_VALS and eg_info.get('normalize', 0) == 1:
                 # Applied via scaling above?
                 # SMOG 2 applies 1/count if countDihedrals is on.
                 # And potentially other normalizations.
                 pass

        val = dihedral_output(d['func'], d['i'], d['j'], d['k'], d['l'], d['theta'], scaling)
        if val:
            dihe_cache.append({'i': d['i'], 'j': d['j'], 'k': d['k'], 'l': d['l'], 'v': val})

    # Add Impropers
    # Iterate RESIDUES/CONNECTIONS and find impropers
    # Impropers: Star topology centered at 'a'. (a-b, a-c, a-d)

    for j, neighbors in adj_list.items():
        if len(neighbors) < 3: continue
        # We need to find specific improper definitions.
        # Iterate over all permutations? Or check template?

        # Check if j is a center of an improper in template
        if j not in idx_to_info: continue
        info_j = idx_to_info[j]
        res_j = RESIDUES[info_j['res_name']]

        # Get impropers for this residue
        impropers = res_j.get('impropers', [])

        # Each improper in template is list of 4 atom names [center, b, c, d]?
        # Or [a, b, c, d] where a is center?
        # BIF: <improper><atom>CA</atom><atom>N</atom><atom>C</atom><atom>CB</atom></improper>
        # Usually first atom is center. Or central atom is bonded to other 3.
        # Logic in SMOG 2: checks connectivity.

        res_atoms = res_j['atoms'] # Map name -> data (including index)
        # But we need to map name -> atom_num in our system.
        # info_j['atom_num'] is j.
        # We need to find neighbors by NAME.

        # We need to check all impropers defined for this residue
        # AND connection impropers.

        # For current residue impropers:
        # Check if all 4 atoms exist in current residue (or connected?)
        # SMOG 2 iterates impropers defined in BIF and tries to map them to atoms.

        # Re-iterate residues to find impropers
        pass

    # Better approach for impropers:
    # Iterate all residues in system.
    # For each residue, iterate its impropers.
    # Map atom names to indices.
    # Check if they form improper (connected).
    # Calculate and output.

    # Need access to all atoms in residue.
    # residues_in_pdb structure from calculate_bonds would be useful here.
    # Rebuild it locally
    residues_in_pdb = {}
    for atom in atom_info:
        res_id = atom['res_counter']
        if res_id not in residues_in_pdb:
            residues_in_pdb[res_id] = {'name': atom['res_name'], 'atoms': {}}
        residues_in_pdb[res_id]['atoms'][atom['atom_name']] = atom

    for res_id in sorted(residues_in_pdb.keys()):
        res_name = residues_in_pdb[res_id]['name']
        res_atoms_map = residues_in_pdb[res_id]['atoms']

        template_res = RESIDUES.get(res_name)
        if not template_res: continue

        # Residue impropers
        if 'impropers' in template_res:
            for imp_atoms in template_res['impropers']:
                # imp_atoms is list of 4 names
                # Check if they exist

                # Handling + / - neighbors?
                # SMOG 2 usually defines impropers within residue or specific connections.
                # If names are simple, look in res_atoms_map.

                indices = []
                valid = True
                for name in imp_atoms:
                    target_res_id = res_id
                    target_name = name
                    if name.endswith('+'):
                        target_res_id += 1
                        target_name = name[:-1]
                    elif name.endswith('-'):
                        target_res_id -= 1
                        target_name = name[:-1]

                    if target_res_id in residues_in_pdb and target_name in residues_in_pdb[target_res_id]['atoms']:
                        indices.append(residues_in_pdb[target_res_id]['atoms'][target_name]['atom_num'])
                    else:
                        valid = False
                        break

                if valid:
                    # Calculate
                    # func? Improper usually harmonic.
                    # Defined in .b file (parse_bonds_xml).
                    # We need types.

                    # Get bTypes
                    types = []
                    for name in imp_atoms:
                        types.append(template_res['atoms'][name]['bType'])

                    func_str = get_improper_function(types)
                    if func_str:
                        # Calculate angle
                        # Definition: improper angle between planes (ijk) and (jkl)?
                        # Or standard dihedral definition for the 4 atoms?
                        # Gromacs improper is a dihedral.

                        i, j, k, l = indices
                        # Calc theta using coordinates directly to cover impropers where bonds might be implicit

                        a1 = idx_to_info[i]
                        a2 = idx_to_info[j]
                        a3 = idx_to_info[k]
                        a4 = idx_to_info[l]

                        v_ij = np.array([a2['x']-a1['x'], a2['y']-a1['y'], a2['z']-a1['z']])
                        v_jk = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
                        v_kl = np.array([a4['x']-a3['x'], a4['y']-a3['y'], a4['z']-a3['z']])

                        n1 = np.cross(v_ij, v_jk)
                        n2 = np.cross(v_jk, v_kl)

                        if norm(n1) > 1e-6 and norm(n2) > 1e-6:
                            n1 /= norm(n1)
                            n2 /= norm(n2)
                            m1 = np.cross(n1, v_jk / norm(v_jk))
                            theta = np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))

                            # Output
                            val = dihedral_output(func_str, i, j, k, l, theta, 1.0) # Scaling?
                            if val:
                                dihe_cache.append({'i': i, 'j': j, 'k': k, 'l': l, 'v': val})

    return dihe_cache

def get_energy_group(info_j, info_k):
    res_j_name = info_j['res_name']
    res_k_name = info_k['res_name']
    atom_j_name = info_j['atom_name']
    atom_k_name = info_k['atom_name']

    if res_j_name == res_k_name:
        bond_key = f"{atom_j_name}-{atom_k_name}"
        if bond_key in RESIDUES[res_j_name]['energyGroups']:
            return RESIDUES[res_j_name]['energyGroups'][bond_key]
        bond_key_rev = f"{atom_k_name}-{atom_j_name}"
        if bond_key_rev in RESIDUES[res_j_name]['energyGroups']:
            return RESIDUES[res_j_name]['energyGroups'][bond_key_rev]
    else:
        type_j = RESIDUES[res_j_name]['residueType']
        type_k = RESIDUES[res_k_name]['residueType']
        from .template import CONNECTIONS
        if type_j in CONNECTIONS and type_k in CONNECTIONS[type_j]:
            for bond in CONNECTIONS[type_j][type_k].findall('bond'):
                return bond.get('energyGroup')
    return None

def get_dihedral_function(type_a, type_b, type_c, type_d, eg):
    if not eg: return None
    if eg not in INTERACTIONS['dihedrals']: return None
    patterns = [[type_a, type_b, type_c, type_d], [type_d, type_c, type_b, type_a]]
    for pat in patterns:
        ta, tb, tc, td = pat
        key = f"{ta}-{tb}-{tc}-{td}"
        if key in INTERACTIONS['dihedrals'][eg]: return INTERACTIONS['dihedrals'][eg][key]
        if '*-*-*-*' in INTERACTIONS['dihedrals'][eg]: return INTERACTIONS['dihedrals'][eg]['*-*-*-*']
        for k in INTERACTIONS['dihedrals'][eg].keys():
            parts = k.split('-')
            if len(parts) != 4: continue
            match = True
            for idx, part in enumerate(parts):
                if part != '*' and part != pat[idx]:
                    match = False
                    break
            if match: return INTERACTIONS['dihedrals'][eg][k]
    return None

def get_improper_function(types):
    # types: list of 4 bTypes
    if 'impropers' not in INTERACTIONS: return None

    ta, tb, tc, td = types
    patterns = [[ta, tb, tc, td]] # Order specific for impropers? Usually.

    for pat in patterns:
        ta, tb, tc, td = pat
        key = f"{ta}-{tb}-{tc}-{td}"
        if key in INTERACTIONS['impropers']: return INTERACTIONS['impropers'][key]

        # Wildcards
        for k in INTERACTIONS['impropers'].keys():
            parts = k.split('-')
            if len(parts) != 4: continue
            match = True
            for idx, part in enumerate(parts):
                if part != '*' and part != pat[idx]:
                    match = False
                    break
            if match: return INTERACTIONS['impropers'][k]
    return None

def dihedral_output(input_func, i, j, k, l, value, scaling):
    if not input_func: return None
    funcs = input_func.split('+')
    output_string = ""
    for f in funcs:
        f = f.strip()
        f_name, params = parse_func_str(f)
        f_type_map = {'dihedral_cosine': 1, 'dihedral_harmonic': 2, 'dihedral_ncos': 1}
        f_type = f_type_map.get(f_name, 1)
        parsed_params = []
        val_local = value
        for idx, p in enumerate(params):
            parsed_params.append(evaluate_param(p, val_local))

        if f_name in ['dihedral_cosine', 'dihedral_ncos']:
            parsed_params[0] += 180.0

        # Apply scaling to force constant (usually param[1] or param[0] depending on function)
        # harmonic(theta0, k) -> k is idx 1
        # cosine(theta0, k, mult) -> k is idx 1
        # Check param definitions?
        # Assuming index 1 is force constant for these standard funcs
        if len(parsed_params) > 1:
            parsed_params[1] *= scaling

        param_str = "\t".join([f"{p:12.9e}" if isinstance(p, (float, int)) else str(p) for p in parsed_params])
        output_string += f"{i}\t{j}\t{k}\t{l}\t{f_type}\t{param_str}\n"
    return output_string
