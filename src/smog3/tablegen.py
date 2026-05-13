"""Native table generation for SMOG Debye-Huckel/nonbonded helper workflows.

The ``smog_tablegen``-compatible command writes GROMACS tabulated potentials
from analytic parameters used in SMOG tutorials.  It is a standalone Python
runtime helper and does not invoke Perl, SMOG2, or Docker.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


DR = 0.002


def _v(kappa: float, r: float) -> float:
    return math.exp(-kappa * r) / r


def _vprime(kappa: float, r: float) -> float:
    return -kappa * math.exp(-kappa * r) / r - math.exp(-kappa * r) / (r * r)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-N", type=int, default=6)
    p.add_argument("-M", type=int, default=12)
    p.add_argument("-ic", type=float, default=0.0)
    p.add_argument("-sd", type=float, default=1.0)
    p.add_argument("-sc", type=float, default=1.5)
    p.add_argument("-tl", type=float, default=5.0)
    p.add_argument("-table", default="table.xvg")
    p.add_argument("-unit", default=None)
    p.add_argument("-temp", type=float, default=None)
    p.add_argument("-help", "-?", action="store_true")
    return p


def main(argv: list[str]) -> int:
    """Generate a GROMACS ``table.xvg`` from SMOG-style table parameters."""

    parser = _build_parser()
    ns, extra = parser.parse_known_args(argv)
    if ns.help or extra:
        print("usage: smog_tablegen -N <int> -M <int> -ic <float> -sd <float> -sc <float> -tl <float> -table <name> -unit <kCal|kJ> -temp <float>")
        return 1

    n = ns.N
    m = ns.M
    conc = ns.ic
    r1 = ns.sd
    rc = ns.sc
    rtable = ns.tl
    outfilename = ns.table if ns.table.endswith(".xvg") else f"{ns.table}.xvg"

    temperature = 300.0 if ns.temp is None else ns.temp
    t_scale = temperature / 300.0

    unit_arg = ns.unit
    if unit_arg is None or "kcal" in unit_arg.lower():
        u_scale = 33.20 / 138.935485
        unit = "kCal"
    elif unit_arg.lower() == "kj":
        u_scale = 1.0
        unit = unit_arg
    else:
        raise SystemExit("Only kJ and kCal are acceptable values to use with -unit")

    prefactor = t_scale * u_scale

    if n >= m:
        raise SystemExit("For consistency with smog2, only M < N is supported")
    if rc <= r1 and rc != 0:
        raise SystemExit("switching distance must be shorter than truncated distance")
    if rc < 0 or r1 < 0:
        raise SystemExit("switching distance and truncated distance must be greater than, or equal to, zero")
    if rc > rtable:
        raise SystemExit("table must extend beyond the cutoff distance")

    if Path(outfilename).exists():
        raise SystemExit(f"File {outfilename} already exists")

    kappa = 3.2 * math.sqrt(conc)
    ntable = int(rtable / DR)

    if rc > 0:
        mcoef = -_v(kappa, rc)
        ncoef = (rc - r1) ** 2
        ocoef = (rc - r1) ** 3
        rcoef = -_vprime(kappa, rc)
        scoef = 2 * (rc - r1)
        tcoef = 3 * (rc - r1) ** 2
        acoef = ((mcoef / ocoef) - (rcoef / tcoef)) / ((ncoef / ocoef) - (scoef / tcoef))
        bcoef = rcoef / tcoef - acoef * scoef / tcoef
    else:
        acoef = 1.0
        bcoef = 1.0

    with open(outfilename, "w", encoding="utf-8") as f:
        f.write("# Force table for use with structure-based SMOG models\n")
        f.write("# Generated with smog_tablegen (version 2.7beta)\n")
        f.write("# Potential of the form U=K*q_i*q_j*exp(-kappa*r)/r-A/r^N+B/r^M\n")
        f.write(f"# kappa={kappa} nm^-1\n")
        f.write(f"# N={n}\n")
        f.write(f"# M={m}\n")
        f.write(f"# T={temperature}\n")
        f.write(f"# units: {unit}\n")
        f.write("# For electrostatics:\n")
        f.write(f"#    \tswitching distance={r1} nm\n")
        f.write(f"#\tcutoff distance={rc} nm\n")
        f.write(f"#\tion concentration={conc} M\n")

        for i in range(8):
            r = i * DR
            f.write(f"{r} 0 0 0.0 0.0 0.0 0.0\n")

        for i in range(8, ntable):
            r = i * DR
            r1v = -1 / r**n
            r2v = -n / r ** (n + 1)
            r3v = 1 / r**m
            r4v = m / r ** (m + 1)
            if r < r1 and rc != 0:
                ve = prefactor * _v(kappa, r)
                vpe = -prefactor * _vprime(kappa, r)
                f.write(f"{r} {ve} {vpe} {r1v} {r2v} {r3v} {r4v}\n")
            elif r1 <= r <= rc:
                ve = prefactor * (_v(kappa, r) + acoef * (r - r1) ** 2 + bcoef * (r - r1) ** 3)
                vpe = -prefactor * (_vprime(kappa, r) + acoef * 2 * (r - r1) + bcoef * 3 * (r - r1) ** 2)
                f.write(f"{r} {ve} {vpe} {r1v} {r2v} {r3v} {r4v}\n")
            else:
                f.write(f"{r} 0.0 0.0 {r1v} {r2v} {r3v} {r4v}\n")

        r = ntable * DR
        r1v = -1 / r**n
        r2v = -n / r ** (n + 1)
        r3v = 1 / r**m
        r4v = m / r ** (m + 1)
        f.write(f"{r} 0.0 0.0 {r1v} {r2v} {r3v} {r4v}\n")

    print(f"Parameter used to generate table file {outfilename}")
    print(f"\tM={m}\n\tN={n}\n\tmonovalent ion conc.={conc} M\n\tkappa={kappa}/nm\n\tbegin switching function={r1} nm\n\tend switching function={rc} nm\n\tLength of table={rtable} nm\n\tTemperature={temperature} K\n\tUnits: {unit}")
    print("\nTable file written successfully\n")
    return 0
