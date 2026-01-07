from .template import INTERACTIONS, RESIDUES
from .utils import smog_quit, smog_note
import xml.etree.ElementTree as ET
import xml.dom.minidom

# OpenSMOG data structure
# {
#   'contacts': { 'contacts_type': { name: { expression: ..., parameters: [], interaction: [ {i, j, ...} ] } } },
#   'dihedrals': { 'dihedrals_type': ... },
#   'angles': { ... },
#   'constants': { 'constant': { name: { value: ... } } },
#   'nonbond': { 'nonbond_bytype': { ... } }
# }

OPENSMOG_DATA = {}

def init_opensmog():
    global OPENSMOG_DATA
    OPENSMOG_DATA = {
        'contacts': {'contacts_type': {}},
        'dihedrals': {'dihedrals_type': {}},
        'angles': {'angles_type': {}},
        'externals': {'externals_type': {}},
        'constants': {'constant': {}},
        'nonbond': {'nonbond_bytype': {}}
    }

def add_interaction(type_key, name, expression, params, interaction_data, exclusions=None):
    # type_key: contacts, dihedrals, angles
    # name: function name (e.g. contact_gaussian)
    # expression: OpenMM expression
    # params: list of parameter names
    # interaction_data: dict of values {i, j, param1, ...}
    # exclusions: 1 to generate, 0 otherwise

    if type_key not in OPENSMOG_DATA:
        init_opensmog()

    group = OPENSMOG_DATA[type_key][f"{type_key}_type"]
    if name not in group:
        group[name] = {
            'expression': {'expr': expression},
            'parameter': params,
            'interaction': []
        }
        if exclusions is not None:
            group[name]['exclusions'] = {'generate': exclusions}

    group[name]['interaction'].append(interaction_data)

def write_opensmog_xml(filename):
    if not OPENSMOG_DATA:
        smog_note("No OpenSMOG data to write.")
        return

    root = ET.Element("OpenSMOGforces", OpenSMOGversion="1.3")

    # Process types in order?
    for type_key in ['contacts', 'dihedrals', 'angles', 'externals', 'constants', 'nonbond']:
        if type_key not in OPENSMOG_DATA: continue

        section_data = OPENSMOG_DATA[type_key]
        if not section_data: continue

        # Check if empty (e.g. contacts_type empty)
        if type_key == 'constants':
            # Handle constants
            const_section = section_data.get('constant', {})
            if not const_section: continue
            elem = ET.SubElement(root, type_key)
            for name, val in const_section.items():
                ET.SubElement(elem, "constant", name=name, value=str(val['value']))
        elif type_key == 'nonbond':
             # Handle nonbond
             # Currently only nonbond_bytype supported
             nb_section = section_data.get('nonbond_bytype', {})
             if not nb_section: continue
             elem = ET.SubElement(root, type_key)
             # ... implementation pending for nonbonds
             pass
        else:
            # Interactions
            subtype_key = f"{type_key}_type"
            if subtype_key not in section_data or not section_data[subtype_key]: continue

            elem = ET.SubElement(root, type_key)

            for func_name, func_data in section_data[subtype_key].items():
                if not func_data['interaction']: continue

                sub_elem = ET.SubElement(elem, subtype_key, name=func_name)
                ET.SubElement(sub_elem, "expression", expr=func_data['expression']['expr'])

                if 'exclusions' in func_data:
                    ET.SubElement(sub_elem, "exclusions", generate=str(func_data['exclusions']['generate']))

                for p in func_data['parameter']:
                    p_elem = ET.SubElement(sub_elem, "parameter")
                    p_elem.text = p

                for interact in func_data['interaction']:
                    # Format attributes
                    attr = {}
                    for k, v in interact.items():
                        if isinstance(v, int):
                            attr[k] = str(v)
                        else:
                            attr[k] = f"{v:.5e}" # Scientific notation
                    ET.SubElement(sub_elem, "interaction", **attr)

    tree = ET.ElementTree(root)
    # Pretty print
    xmlstr = xml.dom.minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")

    with open(filename, "w") as f:
        f.write(xmlstr)
    print(f"\t{filename}")

# Definitions of OpenSMOG potentials matching SMOG 2
# From OpenSMOG.pm / smogv2
OS_POTENTIALS = {
    'contact_1': {
        'expression': "A/r^N-B/r^M",
        'parameters': ["A", "B"] # N, M are substituted in name contact_1-M-N
    },
    'bond_type6': {
        'expression': "eps*0.5*(r-r0)^2",
        'parameters': ["eps", "r0"]
    },
    'contact_gaussian': {
        'expression': "A*((1+a/(A*r^12))*(1-exp(-(r-r0)^2/(2*sigmaG^2)))-1)",
        'parameters': ["A", "r0", "sigmaG", "a"]
    },
    'angle_harmonic': {
        'expression': "0.5*weight*(theta-theta0)^2",
        'parameters': ["theta0", "weight"],
        'defaults': {'theta0': 0, 'weight': 1}
    },
    'dihedral_cosine': {
        'expression': "weight*(1-cos(multiplicity*(theta-theta0)))",
        'parameters': ["theta0", "weight", "multiplicity"],
        'defaults': {'theta0': 0, 'weight': 1}
    },
    'dihedral_ncos': {
        'expression': "weight*(1-cos(multiplicity*theta-theta0))",
        'parameters': ["theta0", "weight", "multiplicity"],
        'defaults': {'theta0': 0, 'weight': 1}
    },
    'dihedral_harmonic': {
        'expression': "weight*(0.5*min(dtheta, 2*pi-dtheta)^2); dtheta = abs(theta-theta0); pi = 3.1415926535",
        'parameters': ["theta0", "weight"],
        'defaults': {'theta0': 0, 'weight': 1}
    }
}
