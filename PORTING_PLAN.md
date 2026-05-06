# SMOG2 -> smog3 Porting Plan

## Goals
- Preserve SMOG2 behavior closely enough that `SMOG-CHECK/smog-check` passes.
- Port Perl components incrementally to Python.
- Keep Java components (SCM/WHAM) via wrappers initially.
- Prepare for PyPI packaging and stable CLI compatibility.

## Current architecture summary
- Main executable logic: `src/smogv2` (Perl).
- Tool scripts: `src/tools/*` (Perl), with SCM/WHAM jars for Java functionality.
- Validation harness: `SMOG-CHECK/smog-check` -> `SMOG-CHECK/src/check.v2.pl` with
  large scenario matrix `SMOG-CHECK/share/settings/smog2.testlist`.

## Migration strategy
1. **Compatibility layer first**
   - Introduce Python package `smog3`.
   - Ship console scripts matching existing command names.
   - Delegate to current Perl scripts to preserve exact behavior.
2. **Incremental replacement**
   - Port lower-risk tools first (`tablegen`, `scale-energies`), then others.
   - Add output-parity tests comparing Perl vs Python on representative inputs.
3. **Core smog2 orchestration port**
   - Move parsing, topology generation, and writers into Python modules.
   - Keep SCM/WHAM calls via Java subprocess wrappers.
4. **Parity gating**
   - Run subset of smog-check continuously during migration.
   - Expand to full smog-check pass before retiring Perl entry points.

## CLI compatibility requirements
- Commands to preserve:
  - `smog2`
  - `smog_adjustPDB`
  - `smog_editgro`
  - `smog_extract`
  - `smog_ions`
  - `smog_modifyXML`
  - `smog_scale-energies`
  - `smog_tablegen`
- Keep default output names and key error/exit behavior.

## Packaging plan
- Use `pyproject.toml` + setuptools build backend.
- Source layout: `src/smog3/`.
- Console entry points in `[project.scripts]`.

## Testing plan
- Unit tests for wrapper/env/path behavior.
- Parity tests: Perl-vs-Python output comparison on selected cases.
- Integration tests that run smog-check subsets against smog3 commands.

## Acceptance criteria
- `smog-check` passes with smog3 CLI commands in place.
- Java-dependent features remain functional via wrappers.
- Documentation updated with installation and usage for smog3.


## Phase status (as of this branch)
- Phase 1 (packaging + compatibility wrappers): **completed**.
- Phase 2 (first native port: `smog_tablegen`): **completed**.
- Phase 3 (`smog-check` wrapper compatibility mode): **completed** via `smog3-install-compat`.
- Phase 4 (smog-check subset automation): **in progress** via `smog3-smogcheck-subset`; subset execution support added.
- Phase 5 (baseline-vs-candidate artifact parity reports): **in progress** via JSON digest report support in `smog3-smogcheck-subset --baseline-smog2 ...`.

- Phase 6 (native topology/index parsing foundations): **in progress** with `smog3.gmx` parser helpers and native `smog_scale-energies` and `smog_extract` now wired as default CLI paths in `smog3.cli`. `smog_extract` currently supports core extraction flow (`-f/-g/-n/-of/-og/-om`, `-group`, `-ndxorder`, `-nondx` with `-distfrom/-distval`) while OpenSMOG/restraints outputs remain pending.
