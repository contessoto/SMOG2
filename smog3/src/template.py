import xml.etree.ElementTree as ET
import os
import sys
import copy
import re
from .utils import smog_quit, smog_note, check_comment, what_am_i, trim
from .config import NB_TYPES_PRESENT

# Global variables for parsed data
RESIDUES = {}
INTERACTIONS = {}
FUNCTIONS = {}
SETTINGS = {}
ENERGY_GROUPS = {}
CONTACT_SETTINGS = {}
INTERACTION_THRESHOLD = {}
ION_DEFS = {}
CONNECTIONS = {}
BOND_FUNCTIONALS = {}
DIHEDRAL_FUNCTIONALS = {}
ANGLE_FUNCTIONALS = {}
USED_FUNCTIONS = {}
COUNT_DIHEDRALS = 1
NORMALIZE_VALS = 0
TERM_RATIOS = {}
B_TYPES_PRESENT = {}
PAIR_TYPES_PRESENT = {}
PAIR_TYPES_USED = {}
BOND_TYPES_USED = {}
ATOM_NAME_HASH = {}
ATOM_NAME_HASH_IN_RES = {}
EG_IN_BIF = {}
EG_IN_SIF = {}

def set_input_file_names(bif, sif, bond, nbond, ion):
    global BIF_XML, SIF_XML, BOND_XML, NBOND_XML, ION_DEFS_FILE
    BIF_XML = bif
    SIF_XML = sif
    BOND_XML = bond
    NBOND_XML = nbond
    ION_DEFS_FILE = ion

def clear_bif_memory():
    global RESIDUES, FUNCTIONS, USED_FUNCTIONS, CONTACT_SETTINGS, TERM_RATIOS, INTERACTIONS
    global EXCLUDE_BONDED, BOND_FUNCTIONALS, DIHEDRAL_FUNCTIONALS, ANGLE_FUNCTIONALS, CONNECTIONS
    global DIHEDRAL_ADJ_LIST, NB_TYPES_PRESENT, B_TYPES_PRESENT, PAIR_TYPES_PRESENT, PAIR_TYPES_USED
    global BOND_TYPES_USED, ION_DEFS
    RESIDUES = {}
    FUNCTIONS = {}
    USED_FUNCTIONS = {}
    CONTACT_SETTINGS = {}
    TERM_RATIOS = {}
    INTERACTIONS = {}
    EXCLUDE_BONDED = {}
    BOND_FUNCTIONALS = {}
    DIHEDRAL_FUNCTIONALS = {}
    ANGLE_FUNCTIONALS = {}
    CONNECTIONS = {}
    DIHEDRAL_ADJ_LIST = {}
    NB_TYPES_PRESENT = {}
    B_TYPES_PRESENT = {}
    PAIR_TYPES_PRESENT = {}
    PAIR_TYPES_USED = {}
    BOND_TYPES_USED = {}
    ION_DEFS = {}

def parse_bif():
    if not os.path.exists(BIF_XML):
        smog_quit(f"BIF file not found: {BIF_XML}")
    try:
        tree = ET.parse(BIF_XML)
        root = tree.getroot()
        bif_residues(root)
        bif_connections(root)
    except Exception as e:
        smog_quit(f"Failed to parse {BIF_XML}: {e}")

def bif_residues(root):
    residues_elem = root.find('residues')
    if residues_elem is None: return

    for residue in residues_elem.findall('residue'):
        res_name = residue.get('name')

        # Atoms
        atoms_data = {}
        index = 0
        atoms_elem = residue.find('atoms')
        if atoms_elem is not None:
            seen = {}
            for atom in atoms_elem.findall('atom'):
                at_content = atom.text.strip() if atom.text else ""

                if at_content not in ATOM_NAME_HASH:
                    ATOM_NAME_HASH[at_content] = []
                    ATOM_NAME_HASH_IN_RES[at_content] = []
                ATOM_NAME_HASH[at_content].append(atom)
                ATOM_NAME_HASH_IN_RES[at_content].append(res_name)

                if at_content in seen:
                    smog_quit(f"Error in .bif. Duplicate declaration of atom {at_content} in residue {res_name}.")
                seen[at_content] = 1

                nb_type = atom.get('nbType')
                b_type = atom.get('bType')
                pair_type = atom.get('pairType')
                charge = atom.get('charge')

                if not re.match(r'^[a-zA-Z0-9_]+$', nb_type):
                    smog_quit(f"Only letters, numbers and _ can appear in nbType definitions. nbType \"{nb_type}\" found in residue {res_name}")
                if not re.match(r'^[a-zA-Z0-9_]+$', b_type):
                    smog_quit(f"Only letters, numbers and _ can appear in bType definitions. bType \"{b_type}\" found in residue {res_name}")
                if not re.match(r'^[a-zA-Z0-9_]+$', pair_type):
                    smog_quit(f"Only letters, numbers and _ can appear in pairType definitions. pairType \"{pair_type}\" found in residue {res_name}")

                NB_TYPES_PRESENT[nb_type] = 1
                B_TYPES_PRESENT[b_type] = 1
                PAIR_TYPES_PRESENT[pair_type] = 1

                atoms_data[at_content] = {
                    "index": index,
                    "nbType": nb_type,
                    "bType": b_type,
                    "pairType": pair_type,
                    "charge": float(charge) if charge else None
                }
                index += 1

        # Impropers
        impropers_data = []
        impropers_elem = residue.find('impropers')
        if impropers_elem is not None:
            seen_imp = {}
            for improper in impropers_elem.findall('improper'):
                atom_list = [a.text.strip() for a in improper.findall('atom')]
                if len(atom_list) != 4:
                    smog_quit(f"Declaration of residue {res_name} has an improper that does not have 4 atoms\n")

                atom_string = "-".join(atom_list)
                atom_string_rev = "-".join(atom_list[::-1])

                # Check for special chars in atom names
                for t in atom_list:
                    if re.search(r'[?^&!@#%()-]', t):
                        smog_quit(f"Special characters not permitted in \"connection\" atom names: {t} found.")

                # Verify atom uniqueness in improper (not implemented in original but maybe useful?)
                # Verify duplicates
                if atom_string in seen_imp or atom_string_rev in seen_imp:
                    smog_quit(f"Error in .bif. Duplicate declaration of improper dihedral {atom_string} for residue {res_name}.")
                seen_imp[atom_string] = 1
                seen_imp[atom_string_rev] = 1

                impropers_data.append(atom_list)

        # Bonds
        bonds_data = {}
        energy_groups = {}
        bonds_elem = residue.find('bonds')
        if bonds_elem is not None:
            for bond in bonds_elem.findall('bond'):
                atom_list = [a.text.strip() for a in bond.findall('atom')]
                if len(atom_list) != 2: continue

                atom_a, atom_b = atom_list
                bond_key = f"{atom_a}-{atom_b}"

                if bond_key in bonds_data:
                    smog_quit(f"Error in .bif. Duplicate declaration of bond between atoms {atom_a} and {atom_b} in residue {res_name}.")

                bonds_data[bond_key] = bond # Store xml element or processed data? Storing element for now
                eg = bond.get('energyGroup')
                energy_groups[bond_key] = eg
                energy_groups[f"{atom_b}-{atom_a}"] = eg
                EG_IN_BIF[eg] = 1

        atom_count = residue.get('atomCount')
        if atom_count is None: atom_count = -1
        else: atom_count = int(atom_count)

        total_charge = residue.get('totalcharge')
        if total_charge: total_charge = float(total_charge)

        RESIDUES[res_name] = {
            "residueType": residue.get('residueType'),
            "atoms": atoms_data,
            "impropers": impropers_data,
            "bonds": bonds_data,
            "energyGroups": energy_groups,
            "atomCount": atom_count,
            "totalcharge": total_charge,
            "connect": residue.get('connect'), # Not fully parsed in perl snippet provided but referenced
            "meta": residue.get('meta')
        }

def bif_connections(root):
    connections_elem = root.find('connections')
    if connections_elem is None: return

    for connection in connections_elem.findall('connection'):
        res_a = connection.get('residueType1')
        res_b = connection.get('residueType2')

        if res_a not in CONNECTIONS: CONNECTIONS[res_a] = {}

        if res_b in CONNECTIONS[res_a]:
            smog_quit(f"Multiple connections defined between residueTypes {res_a} and {res_b}")

        CONNECTIONS[res_a][res_b] = connection

        for bond in connection.findall('bond'):
            EG_IN_BIF[bond.get('energyGroup')] = 1

def parse_sif(os_ref=None):
    if not os.path.exists(SIF_XML):
        smog_quit(f"SIF file not found: {SIF_XML}")
    try:
        tree = ET.parse(SIF_XML)
        root = tree.getroot()
        # sifVersion(root) # To implement
        # OShashAddConstants(root, os_ref) # To implement
        # sifFunctions(root) # To implement
        sif_settings(root)
        # CustomNonBondedCheckAdjust(root, os_ref) # To implement
    except Exception as e:
        smog_quit(f"Failed to parse {SIF_XML}: {e}")

def sif_settings(root):
    global ENERGY_GROUPS, TERM_RATIOS, NORMALIZE_VALS, CONTACT_SETTINGS, INTERACTION_THRESHOLD, COUNT_DIHEDRALS

    settings_elem = root.find('settings')
    if settings_elem is None: return

    groups_elem = settings_elem.find('Groups')
    if groups_elem is not None:
        # Energy Groups
        eg_norm = 0
        for eg in groups_elem.findall('energyGroup'):
            name = eg.get('name')
            res_type = eg.get('residueType')
            intra_rel = eg.get('intraRelativeStrength')
            normalize = eg.get('normalize')

            if normalize: normalize = int(normalize)
            else: normalize = 0

            # Allow strings for intraRelativeStrength (variables)
            # if intra_rel: intra_rel = float(intra_rel)

            if res_type not in TERM_RATIOS: TERM_RATIOS[res_type] = {}
            if 'energyGroup' not in TERM_RATIOS[res_type]: TERM_RATIOS[res_type]['energyGroup'] = {}

            TERM_RATIOS[res_type]['energyGroup'][name] = {
                "normalize": normalize,
                "intraRelativeStrength": intra_rel
            }
            EG_IN_SIF[name] = 1

            if normalize == 1: eg_norm += 1

        # Contact Groups
        cg_norm = 0
        if 'contactGroup' not in TERM_RATIOS: TERM_RATIOS['contactGroup'] = {}
        for cg in groups_elem.findall('contactGroup'):
            name = cg.get('name')
            intra_rel = cg.get('intraRelativeStrength')
            normalize = cg.get('normalize')

            if normalize: normalize = int(normalize)
            else: normalize = 0

            # if intra_rel: intra_rel = float(intra_rel)

            TERM_RATIOS['contactGroup'][name] = {
                "normalize": normalize,
                "intraRelativeStrength": intra_rel
            }
            if normalize == 1: cg_norm += 1

        if eg_norm > 0: NORMALIZE_VALS = 1
        else: NORMALIZE_VALS = 0

        # Group Ratios
        group_ratios = groups_elem.find('groupRatios')
        if group_ratios is not None:
            contacts_ratio = group_ratios.get('contacts')
            dihedrals_ratio = group_ratios.get('dihedrals')

            # Allow string variables
            # contacts_ratio = float(contacts_ratio)
            # dihedrals_ratio = float(dihedrals_ratio)

            TERM_RATIOS["contactRelative"] = contacts_ratio
            TERM_RATIOS["energyRelative"] = dihedrals_ratio
            # TERM_RATIOS["interRelativeTotal"] = contacts_ratio + dihedrals_ratio # Only if floats

    # Contact Settings
    contacts_elem = settings_elem.find('Contacts')
    if contacts_elem is not None:
        CONTACT_SETTINGS = {
            "method": contacts_elem.get('method'),
            "shadowRadiusBonded": contacts_elem.get('shadowRadiusBonded'),
            "shadowRadius": contacts_elem.get('shadowRadius'),
            "contactDistance": contacts_elem.get('contactDistance'),
            "proteinDelta": contacts_elem.get('proteinDelta')
        }

        # Contact Scaling
        scaling_elem = contacts_elem.find('contactScaling')
        if scaling_elem is not None: # Note: structure in perl implies multiple contactScaling entries? No, keys %contactScaling implies it's a list or hash. XML::Simple behavior.
             # Assuming list of scalings if multiple
             pass # Implement scaling parsing

    # Thresholds
    for thresh in ['bondsThreshold', 'anglesThreshold', 'contactsThreshold', 'distanceThreshold']:
        elem = settings_elem.find(thresh)
        if elem is not None:
            INTERACTION_THRESHOLD[thresh.replace('Threshold', '')] = elem.attrib # Store attributes
            if thresh == 'distanceThreshold':
                 INTERACTION_THRESHOLD['distance'] = elem.attrib # Map correctly

    # Dihedral Normalization
    dn_elem = settings_elem.find('dihedralNormalization')
    if dn_elem is not None:
        dc = dn_elem.get('dihedralCounting')
        if dc: COUNT_DIHEDRALS = int(dc)
    else:
        COUNT_DIHEDRALS = 1

def parse_nbonds(os_ref=None):
    if not os.path.exists(NBOND_XML):
        smog_quit(f"NB file not found: {NBOND_XML}")
    try:
        tree = ET.parse(NBOND_XML)
        root = tree.getroot()
        get_nb_defaults(root)
        get_mol_info(root)
        non_bond_loop(root)
        contact_loop(root)
    except Exception as e:
        smog_quit(f"Failed to parse {NBOND_XML}: {e}")

def get_nb_defaults(root):
    defaults_elem = root.find('defaults')
    if defaults_elem is not None:
        INTERACTIONS['gen-pairs'] = 'yes' if defaults_elem.get('gen-pairs') == '1' else 'no'
        INTERACTIONS['nbfunc'] = int(defaults_elem.get('nbfunc', 1))
        INTERACTIONS['gmx-combination-rule'] = int(defaults_elem.get('gmx-combination-rule', 1))
        INTERACTIONS['fudgeLJ'] = float(defaults_elem.get('fudgeLJ', -1))
        INTERACTIONS['fudgeQQ'] = float(defaults_elem.get('fudgeQQ', -1))
    else:
        INTERACTIONS['gen-pairs'] = 'no'
        INTERACTIONS['nbfunc'] = 1
        INTERACTIONS['gmx-combination-rule'] = 1
        INTERACTIONS['fudgeLJ'] = -1
        INTERACTIONS['fudgeQQ'] = -1

def get_mol_info(root):
    mol_elem = root.find('moltype')
    if mol_elem is not None:
        INTERACTIONS['molname'] = mol_elem.get('molname')
        INTERACTIONS['nrexcl'] = int(mol_elem.get('nrexcl'))
    else:
        INTERACTIONS['molname'] = "Macromolecule"
        INTERACTIONS['nrexcl'] = 3

def non_bond_loop(root):
    if 'nonbonds' not in INTERACTIONS: INTERACTIONS['nonbonds'] = {}
    nonbond_elem = root.find('nonbond') # Check if this wraps list or if multiple nonbond elements exist at root level in XML::Simple structure vs ElementTree
    # In ElementTree, we iterate over children
    for nb in root.findall('nonbond'):
        nb_type = nb.find('nbType').text.strip()
        mass = float(nb.get('mass'))
        charge = float(nb.get('charge'))
        ptype = nb.get('ptype')

        func = {"mass": mass, "charge": charge, "ptype": ptype}

        if INTERACTIONS['gmx-combination-rule'] == 1:
            func['c6'] = float(nb.get('c6'))
            func['c12'] = float(nb.get('c12'))
        elif INTERACTIONS['gmx-combination-rule'] == 2:
            func['sigma'] = float(nb.get('sigma'))
            func['epsilon'] = float(nb.get('epsilon'))

        INTERACTIONS['nonbonds'][nb_type] = func

def contact_loop(root):
    if 'contacts' not in INTERACTIONS: INTERACTIONS['contacts'] = {'func': {}, 'contactGroup': {}}
    for contact in root.findall('contact'):
        pair_types = [pt.text.strip() for pt in contact.findall('pairType')]
        if len(pair_types) < 2:
            continue

        type_a = pair_types[0]
        type_b = pair_types[1]

        func = contact.get('func')
        cg = contact.get('contactGroup')

        if type_a not in INTERACTIONS['contacts']['func']: INTERACTIONS['contacts']['func'][type_a] = {}
        if type_b not in INTERACTIONS['contacts']['func']: INTERACTIONS['contacts']['func'][type_b] = {}

        INTERACTIONS['contacts']['func'][type_a][type_b] = func
        INTERACTIONS['contacts']['func'][type_b][type_a] = func

        if 'contactGroup' not in INTERACTIONS['contacts']: INTERACTIONS['contacts']['contactGroup'] = {}
        if type_a not in INTERACTIONS['contacts']['contactGroup']: INTERACTIONS['contacts']['contactGroup'][type_a] = {}
        if type_b not in INTERACTIONS['contacts']['contactGroup']: INTERACTIONS['contacts']['contactGroup'][type_b] = {}

        INTERACTIONS['contacts']['contactGroup'][type_a][type_b] = cg
        INTERACTIONS['contacts']['contactGroup'][type_b][type_a] = cg

def parse_bonds_xml():
    if not os.path.exists(BOND_XML):
        smog_quit(f"Bond file not found: {BOND_XML}")
    try:
        tree = ET.parse(BOND_XML)
        root = tree.getroot()
        b_bonds(root)
        b_angles(root)
        b_dihedrals(root)
        b_impropers(root)
    except Exception as e:
        smog_quit(f"Failed to parse {BOND_XML}: {e}")

def b_bonds(root):
    if 'bonds' not in INTERACTIONS: INTERACTIONS['bonds'] = {}
    bonds_container = root.find('bonds')
    if bonds_container is not None:
        for bond in bonds_container.findall('bond'):
            func = bond.get('func')
            b_types = [bt.text.strip() for bt in bond.findall('bType')]
            if len(b_types) < 2: continue

            type_a, type_b = b_types[0], b_types[1]
            BOND_TYPES_USED[type_a] = 1
            BOND_TYPES_USED[type_b] = 1

            if type_a not in INTERACTIONS['bonds']: INTERACTIONS['bonds'][type_a] = {}
            if type_b not in INTERACTIONS['bonds']: INTERACTIONS['bonds'][type_b] = {}

            INTERACTIONS['bonds'][type_a][type_b] = func
            INTERACTIONS['bonds'][type_b][type_a] = func

def b_angles(root):
    if 'angles' not in INTERACTIONS: INTERACTIONS['angles'] = {}
    angles_container = root.find('angles')
    if angles_container is not None:
        for angle in angles_container.findall('angle'):
            func = angle.get('func')
            b_types = [bt.text.strip() for bt in angle.findall('bType')]
            if len(b_types) < 3: continue

            type_a, type_b, type_c = b_types[0], b_types[1], b_types[2]
            for t in b_types: BOND_TYPES_USED[t] = 1

            key_string = f"{type_a}-{type_b}-{type_c}"
            INTERACTIONS['angles'][key_string] = func
            key_string_rev = f"{type_c}-{type_b}-{type_a}"
            INTERACTIONS['angles'][key_string_rev] = func

def b_dihedrals(root):
    if 'dihedrals' not in INTERACTIONS: INTERACTIONS['dihedrals'] = {}
    dihedrals_container = root.find('dihedrals')
    if dihedrals_container is not None:
        for dihedral in dihedrals_container.findall('dihedral'):
            func = dihedral.get('func')
            eg = dihedral.get('energyGroup')
            b_types = [bt.text.strip() for bt in dihedral.findall('bType')]
            if len(b_types) < 4: continue

            type_a, type_b, type_c, type_d = b_types[0], b_types[1], b_types[2], b_types[3]
            for t in b_types: BOND_TYPES_USED[t] = 1

            if eg not in INTERACTIONS['dihedrals']: INTERACTIONS['dihedrals'][eg] = {}

            key_string = f"{type_a}-{type_b}-{type_c}-{type_d}"
            INTERACTIONS['dihedrals'][eg][key_string] = func

            key_string_rev = f"{type_d}-{type_c}-{type_b}-{type_a}"
            INTERACTIONS['dihedrals'][eg][key_string_rev] = func

def b_impropers(root):
    if 'impropers' not in INTERACTIONS: INTERACTIONS['impropers'] = {}
    impropers_container = root.find('impropers')
    if impropers_container is not None:
        for improper in impropers_container.findall('improper'):
            func = improper.get('func')
            b_types = [bt.text.strip() for bt in improper.findall('bType')]
            if len(b_types) < 4: continue

            type_a, type_b, type_c, type_d = b_types[0], b_types[1], b_types[2], b_types[3]
            for t in b_types: BOND_TYPES_USED[t] = 1

            key_string = f"{type_a}-{type_b}-{type_c}-{type_d}"
            INTERACTIONS['impropers'][key_string] = func

            key_string_rev = f"{type_d}-{type_c}-{type_b}-{type_a}"
            INTERACTIONS['impropers'][key_string_rev] = func
