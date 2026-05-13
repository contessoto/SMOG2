"""Native OpenSMOG XML modification helper.

This small runtime tool applies the subset of ``smog_modifyXML`` edits used by
validation/tutorial workflows: removing force lines by group and scaling numeric
XML attributes for matching group pairs.  It edits XML text directly to preserve
SMOG2-compatible formatting and does not call Perl.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def main(argv: list[str]) -> int:
    """Apply group-based OpenSMOG XML removal or parameter scaling edits."""

    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-OpenSMOG", default="OpenSMOG.xml")
    p.add_argument("-OpenSMOGout", default="OpenSMOG.out.xml")
    p.add_argument("-modtype", default=None)
    p.add_argument("-modgrp1", default=None)
    p.add_argument("-modgrp2", default=None)
    p.add_argument("-modparam", default=None)
    p.add_argument("-modby", default=None)
    p.add_argument("-remove", action="store_true")
    p.add_argument("-help", "-?", action="store_true")
    ns, extra = p.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_modifyXML -OpenSMOG in.xml -OpenSMOGout out.xml [mod flags]")
        return 1

    txt = Path(ns.OpenSMOG).read_text()
    if ns.remove and ns.modgrp1:
        txt = re.sub(rf'.*group1="{re.escape(ns.modgrp1)}".*\n?', '', txt)
    elif ns.modparam and ns.modby is not None and ns.modgrp1:
        fac = float(ns.modby)
        pat = re.compile(rf'({re.escape(ns.modparam)}=")([0-9eE+\-.]+)(")')
        out_lines = []
        for ln in txt.splitlines(True):
            if f'group1="{ns.modgrp1}"' in ln and (ns.modgrp2 is None or f'group2="{ns.modgrp2}"' in ln):
                ln = pat.sub(lambda m: f'{m.group(1)}{float(m.group(2))*fac}{m.group(3)}', ln)
            out_lines.append(ln)
        txt = ''.join(out_lines)

    Path(ns.OpenSMOGout).write_text(txt)
    return 0
