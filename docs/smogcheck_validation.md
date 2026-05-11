# SMOG-CHECK Validation

SMOG3 has been validated with the original SMOG-CHECK harness by placing a
temporary `smog2` wrapper first in `PATH`.  The legacy harness still prints
`SMOG2` because its messages are hardcoded, but each `smog2` invocation is
routed to Python SMOG3 code.

## How The Drop-In Wrapper Works

`scripts/run_smogcheck_dropin_smog3.sh` creates a run directory under
`smogcheck_dropin_runs/` and writes a temporary executable:

```text
smogcheck_dropin_runs/<run>/bin/smog2
```

That wrapper exports:

```bash
SMOG3_LEGACY_PERL_FALLBACK=0
PYTHONPATH=<repo>/src
```

and executes:

```bash
python3 -m smog3.smogcheck_dropin_smog2 "$@"
```

The Python module translates legacy SMOG-CHECK command lines into native SMOG3
arguments and then calls `smog3.smog2_native.main`.

## Validation Commands

```bash
PYTHONPATH=src pytest -q tests
bash scripts/run_selected_two_stage_parity.sh
bash scripts/run_smogcheck_dropin_smog3.sh 1 115
```

Expected final line:

```text
PASSED TESTS 1 to 115
```

## Wrapper Invocation Log

Each drop-in run writes:

```text
smogcheck_dropin_runs/<run>/smog3-wrapper-invocations.jsonl
```

Each JSON line records:

- current working directory
- process id
- original SMOG-CHECK arguments
- translated SMOG3 arguments

The verified full run recorded 138 SMOG3 invocations and zero SMOG3 Perl
fallback invocations.

## Perl Boundary

Perl is used by the original SMOG-CHECK harness.  Perl is not used by SMOG3
normal runtime.  Original SMOG2 is used only inside the official Docker image
when generating two-stage parity baselines.
