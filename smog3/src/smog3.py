import argparse
import sys
import os
import subprocess
from .utils import smog_quit, check_already_exists, check_suffix, smog_note
from .config import VERSION
from . import template as t
from . import pdb as p
from . import bonded as b
from . import topology as top
from . import opensmog
import shutil

def parse_args():
    parser = argparse.ArgumentParser(description=f"SMOG v{VERSION}")

    parser.add_argument('-i', dest='input_pdb', help='Input PDB file to generate Hamiltonian', default='molecule.pdb')
    parser.add_argument('-g', dest='gro_file', help='Output .gro file name')
    parser.add_argument('-o', dest='top_file', help='Output .top file name', default='smog.top')
    parser.add_argument('-s', dest='shadow_file', help='Output .contacts file name')
    parser.add_argument('-n', dest='ndx_file', help='Output .ndx file name')
    parser.add_argument('-c', dest='contact_file', help='Input contact map file')

    parser.add_argument('-t', dest='input_folder', help='Folder containing templates of molecular and interaction definitions')
    parser.add_argument('-tCG', '-tcg', dest='coarse_folder', help='Folder containing templates used for coarse graining')

    # Case insensitive flags (by adding aliases)
    parser.add_argument('-AA', '-aa', action='store_true', dest='AA', help='Use default All-atom model')
    parser.add_argument('-AAgaussian', '-aagaussian', action='store_true', dest='AAgaussian', help='Use default All-atom model with Gaussian contacts')
    parser.add_argument('-CA', '-ca', action='store_true', dest='CA', help='Use default Calpha protein model')
    parser.add_argument('-CAgaussian', '-cagaussian', action='store_true', dest='CAgaussian', help='Use default Calpha protein model with Gaussian contacts')

    parser.add_argument('-dname', dest='dname', default='smog', help='Default name to use for all output files')
    parser.add_argument('-backup', dest='backup', default='yes', help='Back up any pre-existing output files')

    parser.add_argument('-warn', dest='warn', type=int, default=0, help='Convert the first N fatal errors to warnings')
    parser.add_argument('-info', action='store_true', help='Show information about time and energy units')
    parser.add_argument('-v', action='store_true', help='Print version')
    parser.add_argument('-freecoor', action='store_true', help='Allow free-format PDB')
    parser.add_argument('-SCMorig', '-scmorig', action='store_true', dest='SCMorig', help='Directly save SCM contact map')
    parser.add_argument('-keep4SCM', '-keep4scm', action='store_true', dest='keep4SCM', help='Keep temporary files used for SCM')

    parser.add_argument('-OpenSMOG', '-opensmog', action='store_true', dest='OpenSMOG', help='Enable OpenSMOG mode')
    parser.add_argument('-OpenSMOGxml', '-opensmogxml', dest='opensmog_xml', help='Output OpenSMOG XML file name')

    args = parser.parse_args()

    if args.v:
        print(f"Version {VERSION}")
        sys.exit(0)

    return args

def main():
    print(f"******* ******* ******* ******* ******* SMOG v{VERSION} ******* ******* ******* ******* *******")
    args = parse_args()

    # Process args logic similar to smogv2
    if args.dname:
        if not args.top_file: args.top_file = args.dname
        if not args.gro_file: args.gro_file = args.dname
        if not args.shadow_file: args.shadow_file = args.dname
        if not args.ndx_file: args.ndx_file = args.dname

    args.top_file = check_suffix(args.top_file, ".top")
    if args.gro_file: args.gro_file = check_suffix(args.gro_file, ".gro")
    else: args.gro_file = args.dname + ".gro"

    args.shadow_file = check_suffix(args.shadow_file, ".contacts")
    args.ndx_file = check_suffix(args.ndx_file, ".ndx")

    if args.OpenSMOG:
        if not args.opensmog_xml: args.opensmog_xml = args.dname + ".xml"
        args.opensmog_xml = check_suffix(args.opensmog_xml, ".xml")

    gro_file_scm = args.gro_file + "4SCM.gro"
    top_file_scm = args.top_file + "4SCM.top"

    if args.backup == "yes":
        backup_files = [args.top_file, args.gro_file, args.shadow_file, args.ndx_file]
        if args.OpenSMOG: backup_files.append(args.opensmog_xml)
        for f in backup_files:
            if f: check_already_exists(f)

    # Template selection
    from .config import DEFAULT_TEMPLATE_PATH

    input_folder = args.input_folder

    # Map flags to subdirectories
    subdir = None
    if args.AA:
        subdir = "SBM_AA"
    elif args.AAgaussian:
        subdir = "SBM_AA+gaussian"
    elif args.CA:
        subdir = "SBM_calpha"
    elif args.CAgaussian:
        subdir = "SBM_calpha+gaussian"

    if subdir:
        # Check explicit path first (legacy support)
        if input_folder and os.path.isdir(os.path.join(input_folder, subdir)):
             input_folder = os.path.join(input_folder, subdir)
        else:
             # Use detected default path
             input_folder = os.path.join(DEFAULT_TEMPLATE_PATH, subdir)

    if not input_folder:
        smog_quit("No template folder specified. Use -t or a default flag like -AA.")

    if not os.path.exists(input_folder):
        # Fallback: check if user provided path exists relative to CWD
        if args.input_folder and os.path.exists(args.input_folder):
            input_folder = args.input_folder
        else:
            smog_quit(f"Template folder {input_folder} not found. (Search base: {DEFAULT_TEMPLATE_PATH})")

    print(f"Parsing templates from {input_folder}")

    # Locate files
    bif = None
    sif = None
    nb = None
    bond = None

    for f in os.listdir(input_folder):
        if f.endswith(".bif"): bif = os.path.join(input_folder, f)
        elif f.endswith(".sif"): sif = os.path.join(input_folder, f)
        elif f.endswith(".nb"): nb = os.path.join(input_folder, f)
        elif f.endswith(".b"): bond = os.path.join(input_folder, f)

    t.set_input_file_names(bif, sif, bond, nb, "")
    t.parse_bif()
    t.parse_sif()
    t.parse_nbonds()
    t.parse_bonds_xml()

    print(f"Reading {args.input_pdb}")
    atom_info = p.convert_pdb_to_gro_ndx(
        args.input_pdb,
        args.gro_file,
        None,
        gro_file_scm,
        args.ndx_file,
        args.contact_file is not None,
        args.freecoor
    )

    print("Calculating bonded interactions...")
    index_handle = {atom['atom_num']: atom for atom in atom_info} # Need mapping?
    # Logic in bonded expects just index handle for output formatting?
    # Actually bonded.py uses atom_info list mostly.

    bonds = b.calculate_bonds(atom_info, index_handle)
    angles = b.calculate_angles(atom_info, index_handle)
    dihedrals = b.calculate_dihedrals(atom_info, index_handle)

    pairs = []
    exclusions = []

    # Contact Map
    if args.contact_file:
        print(f"Using provided contact file: {args.contact_file}")
        # Need to parse and convert to pairs
        pairs, exclusions = parse_contacts(args.contact_file, atom_info, args.OpenSMOG)
    else:
        print("Generating contact map...")
        # Write simplified TOP for SCM?
        # SCM needs atoms and maybe some basic connectivity?
        # SMOG 2 generates a top file for SCM.
        # "printTop($topFile4SCM, 0)"

        # We need to write a temporary top file for SCM
        # SCM temp top should contain OpenSMOG terms? SMOG2 code includes Dihedrals in SCM temp top.
        # For now, pass all terms.
        top.write_topology(top_file_scm, atom_info, bonds, angles, dihedrals, [], [], t.INTERACTIONS, t.INTERACTIONS.get('molname', 'Macromolecule'), t.INTERACTIONS.get('nrexcl', 3))

        smog_path = os.environ.get("SMOG_PATH", ".")
        scm_jar = os.path.join(smog_path, "src/tools/SCM.jar")

        # "java $memoryMax -jar $ENV{SMOG_PATH}/src/tools/SCM.jar -g $groFile4SCM -freecoor -t $topFile -o $shadowFile -ch $ndxFile $SCMparams "

        # Construct SCM params from settings
        scm_params = []
        c_settings = t.CONTACT_SETTINGS
        if c_settings:
            method = c_settings.get('method', '')
            if 'shadow' in method:
                scm_params.extend(['-m', 'shadow'])
                scm_params.extend(['-c', c_settings.get('contactDistance', '4.0')])
                scm_params.extend(['-s', c_settings.get('shadowRadius', '1.0')])
                scm_params.extend(['-br', c_settings.get('shadowRadiusBonded', '0.0')])
            elif 'cutoff' in method:
                scm_params.extend(['-m', 'shadow']) # SCM uses shadow mode with 0 radius for cutoff? Or specific flag? SMOG2 code: "-m shadow -c $dist -s $radius -br $radiusBonded" where radius=0.
                scm_params.extend(['-c', c_settings.get('contactDistance', '4.0')])
                scm_params.extend(['-s', '0.0'])
                scm_params.extend(['-br', '0.0'])

            if c_settings.get('proteinDelta'):
                scm_params.extend(['-pd', c_settings['proteinDelta']])

        scm_params.extend(['--smog2output'])

        cmd = ['java', '-jar', scm_jar, '-g', gro_file_scm, '-freecoor', '-t', top_file_scm, '-o', args.shadow_file, '-ch', args.ndx_file] + scm_params

        print(" ".join(cmd))
        subprocess.check_call(cmd)

        if args.SCMorig:
            shutil.copy(args.shadow_file, args.shadow_file + ".ShadowOutput")

        # Parse contacts
        pairs, exclusions = parse_contacts(args.shadow_file, atom_info, args.OpenSMOG)

    # Write final topology
    top.write_topology(args.top_file, atom_info, bonds, angles, dihedrals, pairs, exclusions, t.INTERACTIONS, t.INTERACTIONS.get('molname', 'Macromolecule'), t.INTERACTIONS.get('nrexcl', 3), opensmog_enabled=args.OpenSMOG)

    if args.OpenSMOG:
        opensmog.write_opensmog_xml(args.opensmog_xml)

    print(f"\nYour Structure-based Model is ready!\nFiles generated:\n\t{args.top_file}")
    if args.OpenSMOG: print(f"\t{args.opensmog_xml}")
    print(f"\t{args.gro_file}\n\t{args.ndx_file}\n\t{args.shadow_file}")

def parse_contacts(contact_file, atom_info, opensmog_enabled=False):
    pairs = []
    exclusions = []

    if not os.path.exists(contact_file): return pairs, exclusions

    # Need access to atom info to get types
    # atom_info is list of dicts indexable by atom_num
    idx_to_info = {a['atom_num']: a for a in atom_info}

    with open(contact_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split()
            if len(parts) < 3: continue

            try:
                # SCM output with --smog2output is: chain1 serial1 chain2 serial2 [dist]
                # indices are 1-based in output but we need to check consistency with atom_info keys

                i = int(parts[1])
                j = int(parts[3])

                if len(parts) >= 5:
                    dist_ang = float(parts[4])
                    dist_nm = dist_ang * 0.1
                else:
                    # Calculate distance from structure
                    if i in idx_to_info and j in idx_to_info:
                        a1 = idx_to_info[i]
                        a2 = idx_to_info[j]
                        # atom_info coords are in nm (converted in pdb.py)
                        d2 = (a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2
                        dist_nm = d2**0.5
                    else:
                        continue

                if i not in idx_to_info or j not in idx_to_info: continue

                info_i = idx_to_info[i]
                info_j = idx_to_info[j]

                res_i = t.RESIDUES[info_i['res_name']]
                res_j = t.RESIDUES[info_j['res_name']]

                type_i = res_i['atoms'][info_i['atom_name']]['pairType']
                type_j = res_j['atoms'][info_j['atom_name']]['pairType']

                func, cg = get_contact_function(type_i, type_j)

                if not func:
                    continue

                f_type_map = {
                    'contact_1': 1,
                    'contact_2': 1,
                    'contact_gaussian': 6,
                    'bond_type6': 6
                }

                f_name, f_params = b.parse_func_str(func)
                f_type = f_type_map.get(f_name, 1)

                cg_info = t.TERM_RATIOS.get('contactGroup', {}).get(cg, {})
                relative_strength = cg_info.get('intraRelativeStrength')
                if relative_strength is None: relative_strength = 1.0
                epsilon = float(relative_strength)

                parsed_params = []

                # Define variables to ensure they exist for OpenSMOG logic
                m_exp = 0.0
                n_exp = 0.0

                if f_name == 'contact_1':
                    m_exp = float(f_params[0])
                    n_exp = float(f_params[1])

                    r0 = dist_nm

                    b_coef = -1.0 / (m_exp / n_exp - 1.0)
                    a_coef = (m_exp / n_exp) * b_coef

                    a_coef *= epsilon * (r0 ** n_exp)
                    b_coef *= epsilon * (r0 ** m_exp)

                    parsed_params = [b_coef, a_coef]

                elif f_name == 'contact_gaussian':
                    # contact_gaussian(eps, width, r0, ...)
                    parsed_params = []
                    for p in f_params:
                         if '?' in p:
                             p_val = p.replace('?', str(dist_nm))
                             try: p_val = float(eval(p_val))
                             except: pass
                             parsed_params.append(p_val)
                         else:
                             try: parsed_params.append(float(p))
                             except: parsed_params.append(p)

                else:
                    parsed_params = [dist_nm, epsilon]

                if opensmog_enabled:
                    # Add to OpenSMOG
                    # Determine function name for OS (e.g. contact_1-M-N)
                    os_func_name = f_name
                    if f_name == 'contact_1':
                        os_func_name = f"contact_1-{int(m_exp)}-{int(n_exp)}"
                        # Params: A, B. parsed_params is [B, A]. Swap for OpenSMOG
                        os_params = [parsed_params[1], parsed_params[0]]

                        # Add to OpenSMOG data
                        opensmog.add_interaction('contacts', os_func_name,
                            opensmog.OS_POTENTIALS['contact_1']['expression'].replace('N', str(int(n_exp))).replace('M', str(int(m_exp))),
                            opensmog.OS_POTENTIALS['contact_1']['parameters'],
                            {'i': i, 'j': j, 'A': os_params[0], 'B': os_params[1]},
                            exclusions=1
                        )
                    elif f_name == 'contact_gaussian':
                        # Params: A, r0, sigmaG, a
                        # parsed_params from loop above might be raw.
                        # Need specific mapping.
                        # Assuming template: contact_gaussian(eps_c, eps_nc, sigma, ?)
                        # parsed_params: [eps_c, eps_nc, sigma, r0]
                        # OpenSMOG expects: A, r0, sigmaG, a
                        # Mapping from Perl: A=p[0], r0=p[3], sigmaG=p[2], a=p[1]
                        if len(parsed_params) >= 4:
                            os_params = {
                                'A': parsed_params[0] * epsilon if isinstance(parsed_params[0], float) else parsed_params[0], # Scale eps?
                                'r0': parsed_params[3],
                                'sigmaG': parsed_params[2],
                                'a': parsed_params[1]
                            }
                            # Note: eps logic in Perl: p[0] *= epsilon.
                            opensmog.add_interaction('contacts', 'contact_gaussian',
                                opensmog.OS_POTENTIALS['contact_gaussian']['expression'],
                                opensmog.OS_POTENTIALS['contact_gaussian']['parameters'],
                                {'i': i, 'j': j, **os_params},
                                exclusions=1
                            )
                    # Add support for other functions if needed

                    # Do NOT add to pairs list for TOP file if OpenSMOG enabled for this term
                    pass
                else:
                    pair_line = f"{i}\t{j}\t{f_type}\t" + "\t".join([f"{p:12.9e}" for p in parsed_params]) + "\n"
                    pairs.append(pair_line)

                exclusions.append(f"{i}\t{j}\n")

            except Exception as e:
                pass

    return pairs, exclusions

def get_contact_function(type_a, type_b):
    if 'contacts' not in t.INTERACTIONS or 'func' not in t.INTERACTIONS['contacts']: return None, None
    funcs = t.INTERACTIONS['contacts']['func']
    groups = t.INTERACTIONS['contacts']['contactGroup']

    # Check explicit
    if type_a in funcs and type_b in funcs[type_a]:
        return funcs[type_a][type_b], groups[type_a][type_b]

    # Wildcards
    # A-*
    if type_a in funcs and '*' in funcs[type_a]:
        return funcs[type_a]['*'], groups[type_a]['*']
    if type_b in funcs and '*' in funcs[type_b]:
        return funcs[type_b]['*'], groups[type_b]['*']

    # *-*
    if '*' in funcs and '*' in funcs['*']:
        return funcs['*']['*'], groups['*']['*']

    return None, None

if __name__ == "__main__":
    main()
