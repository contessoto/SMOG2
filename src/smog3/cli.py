from __future__ import annotations

import sys

from .compat import run_smog2, run_tool
import os
from .tablegen import main as tablegen_main
from .scale_energies_native import main as scale_energies_main
from .extract_native import main as extract_main
from .editgro_native import main as editgro_main
from .adjustpdb_native import main as adjustpdb_main
from .ions_native import main as ions_main
from .modifyxml_native import main as modifyxml_main
from .smog2_native import main as smog2_native_main
from .parity_direct import main as parity_direct_main


def smog3_main() -> int:
    return smog2_native_main(sys.argv[1:])


def smog2_main() -> int:
    rc = smog2_native_main(sys.argv[1:])
    if rc != 0 and os.environ.get("SMOG3_LEGACY_PERL_FALLBACK", "0") == "1":
        return run_smog2(sys.argv[1:])
    return rc


def smog_adjustpdb_main() -> int:
    return adjustpdb_main(sys.argv[1:])


def smog_editgro_main() -> int:
    return editgro_main(sys.argv[1:])


def smog_extract_main() -> int:
    return extract_main(sys.argv[1:])


def smog_ions_main() -> int:
    return ions_main(sys.argv[1:])


def smog_modifyxml_main() -> int:
    return modifyxml_main(sys.argv[1:])


def smog_scale_energies_main() -> int:
    return scale_energies_main(sys.argv[1:])


def smog_tablegen_main() -> int:
    return tablegen_main(sys.argv[1:])


def smog3_parity_direct_main() -> int:
    return parity_direct_main(sys.argv[1:])
