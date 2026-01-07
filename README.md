########################################################

              README for SMOG 2 / SMOG 3

SMOG 2 (and its Python successor SMOG 3) is a software package for
generating structure-based models, as provided by the
smog-server.org webserver (i.e. SMOG v1). These models
may be used to perform simulations using Gromacs, NAMD,
or OpenMM.

SMOG 3 is a direct port of SMOG 2 to Python 3, aiming for
improved portability, easier installation, and future extensibility,
while maintaining compatibility with SMOG 2 outputs and templates.

SMOG 2/3 is free software, distributed under the GNU
General Public License. See COPYING.

This distribution includes definitions for:
- C-alpha model (Clementi et al. JMB 2000)
- All-atom model (Whitford et al. Proteins 2009)
- Gaussian contact potentials

The software allows designing custom structure-based models without
modifying the source code. The `-OpenSMOG` option facilitates
usage with OpenMM.

For more information, consult the manual and tutorials at
smog-server.org.

Feedback: info@smog-server.org
Users Forum: https://mailman.rice.edu/mailman/listinfo/smog-users

      Jeffrey Noel, Jose' Onuchic and Paul Whitford
      (SMOG 3 port by Jules)

########################################################


########################################################

              INSTALLATION (SMOG 3)

SMOG 3 simplifies installation using standard Python tools.

Prerequisites:
- Python 3.6+
- Java (required for SCM contact map generation tool)
- pip

Installation from PyPI (if available):
> pip install smog3

Installation from source:
> git clone <repo_url>
> cd smog3
> pip install .

For development (editable install):
> pip install -e .

Verifying installation:
> smog3 -v
Version 2.7.0

########################################################


########################################################

              USAGE

The `smog3` command interface mirrors `smog2`.

Common Examples:

1. Default C-alpha model (proteins only):
   > smog3 -i protein.pdb -CA

2. Default All-atom model:
   > smog3 -i molecule.pdb -AA

3. All-atom model with Gaussian contacts:
   > smog3 -i molecule.pdb -AAgaussian

4. Using OpenSMOG (for OpenMM):
   > smog3 -i molecule.pdb -AA -OpenSMOG

   This generates `smog.xml` in addition to standard GROMACS files.

5. Custom Templates:
   > smog3 -i molecule.pdb -t /path/to/my/templates

Flags are case-insensitive (e.g., `-aa`, `-opensmog`).

Main Arguments:
  -i <file>       Input PDB file (default: molecule.pdb)
  -o <file>       Output topology file (default: smog.top)
  -g <file>       Output GRO file
  -n <file>       Output index file
  -t <dir>        Directory containing templates
  -dname <name>   Default base name for output files (default: smog)
  -warn <int>     Convert first N fatal errors to warnings

Models:
  -AA             All-Atom Model
  -CA             C-alpha Model
  -AAgaussian     All-Atom with Gaussian contacts
  -CAgaussian     C-alpha with Gaussian contacts

Advanced:
  -OpenSMOG       Enable OpenSMOG mode (generates XML)
  -contact_file   Use provided contact map instead of generating one
  -SCMorig        Save original SCM output

For a full list of options:
> smog3 --help

########################################################
