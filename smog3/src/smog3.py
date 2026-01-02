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

def parse_args():
    parser = argparse.ArgumentParser(description=f"SMOG v{VERSION}")

    parser.add_argument('-i', dest='input_pdb', help='Input PDB file to generate Hamiltonian', default='molecule.pdb')
    parser.add_argument('-g', dest='gro_file', help='Output .gro file name')
    parser.add_argument('-o', dest='top_file', help='Output .top file name', default='smog.top')
    parser.add_argument('-s', dest='shadow_file', help='Output .contacts file name')
    parser.add_argument('-n', dest='ndx_file', help='Output .ndx file name')
    parser.add_argument('-c', dest='contact_file', help='Input contact map file')

    parser.add_argument('-t', dest='input_folder', help='Folder containing templates of molecular and interaction definitions')
    parser.add_argument('-tCG', dest='coarse_folder', help='Folder containing templates used for coarse graining')

    parser.add_argument('-AA', action='store_true', help='Use default All-atom model')
    parser.add_argument('-AAgaussian', action='store_true', help='Use default All-atom model with Gaussian contacts')
    parser.add_argument('-CA', action='store_true', help='Use default Calpha protein model')
    parser.add_argument('-CAgaussian', action='store_true', help='Use default Calpha protein model with Gaussian contacts')

    parser.add_argument('-dname', dest='dname', default='smog', help='Default name to use for all output files')
    parser.add_argument('-backup', dest='backup', default='yes', help='Back up any pre-existing output files')

    parser.add_argument('-warn', dest='warn', type=int, default=0, help='Convert the first N fatal errors to warnings')
    parser.add_argument('-info', action='store_true', help='Show information about time and energy units')
    parser.add_argument('-v', action='store_true', help='Print version')
    parser.add_argument('-freecoor', action='store_true', help='Allow free-format PDB')
    parser.add_argument('-SCMorig', action='store_true', help='Directly save SCM contact map')
    parser.add_argument('-keep4SCM', action='store_true', help='Keep temporary files used for SCM')

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

    gro_file_scm = args.gro_file + "4SCM.gro"
    top_file_scm = args.top_file + "4SCM.top"

    if args.backup == "yes":
        for f in [args.top_file, args.gro_file, args.shadow_file, args.ndx_file]:
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
        pass
    else:
        print("Generating contact map...")
        # Write simplified TOP for SCM?
        # SCM needs atoms and maybe some basic connectivity?
        # SMOG 2 generates a top file for SCM.
        # "printTop($topFile4SCM, 0)"

        # We need to write a temporary top file for SCM
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

        # Parse contacts
        pairs, exclusions = parse_contacts(args.shadow_file, atom_info)

    # Write final topology
    top.write_topology(args.top_file, atom_info, bonds, angles, dihedrals, pairs, exclusions, t.INTERACTIONS, t.INTERACTIONS.get('molname', 'Macromolecule'), t.INTERACTIONS.get('nrexcl', 3))

    print(f"\nYour Structure-based Model is ready!\nFiles generated:\n\t{args.top_file}\n\t{args.gro_file}\n\t{args.ndx_file}\n\t{args.shadow_file}")

def parse_contacts(contact_file, atom_info):
    pairs = []
    exclusions = []

    if not os.path.exists(contact_file): return pairs, exclusions

    # Need access to atom info to get types
    # atom_info is list of dicts indexable by atom_num
    idx_to_info = {a['atom_num']: a for a in atom_info}

    # SCM output format with --smog2output is "i j chaini chainj dist" or just "i j ..."?
    # The file we inspected earlier had "1 15 1 202" ?
    # Actually "1 15 1 202" looks like chain i, chain j, type, distance? No.
    # SMOG 2 SCM wrapper says: "print $contactFile "$i $j $chaini $chainj $dist\n";"
    # But wait, we saw "1 15 1 202" -> i=1, j=15?
    # The first line "1 15 1 202" likely: 1 15 1 202? No float?
    # Maybe "i j chaini chainj distance" ?
    # Let's assume standard SCM format if --smog2output is used.
    # In smogv2.pl: "parseCONTACT" reads it.

    with open(contact_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split()
            if len(parts) < 3: continue

            try:
                i, j = int(parts[0]), int(parts[1])
                # SCM.jar output: i j dist (if --smog2output used and format is simple)
                # But previously we saw "1 15 1 202"
                # If length is 4, maybe "i j chaini chainj" and no dist?
                # Or "i j chaini dist"?
                # SMOG 2 uses --smog2output.
                # In SCM code (ShadowMain.java):
                # if(ShadowSettings.SMOG2_OUTPUT_ON) ... System.out.println("\t"+s)
                # And s depends on ContactMap.java logic.
                # It writes "atom[i].getChain() ... atom[i].getPosition() ... dist?"
                # Actually earlier grep showed: output.write( ... getChain ... getResNum ... )
                # It seems complicated.
                # BUT `parseCONTACT` in smogv2 simply takes the LAST token as distance.
                # `my $dist=$tokens[$#tokens];`
                # So we will do the same.
                dist = float(parts[-1])

                # Note: if dist is > 100 or suspiciously large, check units.
                # SCM uses Angstroms?
                # SMOG 2 converts: `if($angToNano){ $dist *= 0.1; }`
                # $angToNano is 0.1.
                # So we assume Angstrom input and convert to nm.

                if i not in idx_to_info or j not in idx_to_info: continue

                info_i = idx_to_info[i]
                info_j = idx_to_info[j]

                res_i = t.RESIDUES[info_i['res_name']]
                res_j = t.RESIDUES[info_j['res_name']]

                type_i = res_i['atoms'][info_i['atom_name']]['pairType']
                type_j = res_j['atoms'][info_j['atom_name']]['pairType']

                func, cg = get_contact_function(type_i, type_j)

                if not func: continue

                dist_nm = dist * 0.1

                f_type_map = {
                    'contact_1': 1,
                    'contact_2': 1,
                    'contact_gaussian': 6,
                    'bond_type6': 6
                }

                f_name, f_params = b.parse_func_str(func)
                f_type = f_type_map.get(f_name, 1)

                cg_info = t.TERM_RATIOS.get('contactGroup', {}).get(cg, {})
                relative_strength = cg_info.get('intraRelativeStrength', 1.0)
                epsilon = float(relative_strength)

                parsed_params = []

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
                    # Check definition in SMOGglobals.pm: A, r0, sigmaG, a
                    # params: A (eps), a (?), r0?, sigmaG?
                    # "A*((1+a/(A*r^12))*(1-exp(-(r-r0)^2/(2*sigmaG^2)))-1)"
                    # Perl: ($A,$r0,$sigma_G,$a) = ($paramArr->[0],$paramArr->[3],$paramArr->[2],$paramArr->[1]);
                    # paramArr in Perl before shift: eps_c, eps_nc, sigma, r0.
                    # python f_params: same order as string.
                    # contact_gaussian(eps_c, eps_nc, sigma) usually r0 is implicit?
                    # If template has ?, it is replaced by r0.
                    # We need to map `f_params` to output.

                    # Assume simple Gaussian for now if template used it.
                    # Standard: A, r0, sigmaG

                    # Reusing contactParseGaussian logic from perl port would be best but complex here.
                    # Fallback to passing raw params with ? substitution.
                    pass
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

                pair_line = f"{i}\t{j}\t{f_type}\t" + "\t".join([f"{p:12.9e}" for p in parsed_params]) + "\n"
                pairs.append(pair_line)
                exclusions.append(f"{i}\t{j}\n")

            except Exception as e:
                # print(f"Error parsing contact line {line}: {e}")
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
