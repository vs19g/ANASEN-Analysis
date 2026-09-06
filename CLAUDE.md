# ANASEN Analysis — working notes for Claude

## Agent behavior
- ALWAYS output a brief, step-by-step plan before modifying files or running commands.
- NEVER push to remote. Commit only when explicitly asked.
- **Targeted edits only.** Never regenerate a whole file to make a local change.
- **Never read `MakeVertex.C`**
- **Never read `TrackRecon.C` end to end** (~4,250 lines; ~63k tokens each, and it stays in the prefix for the rest of the session). Grep for the symbol, then Read with offset/limit around the hit. If a question genuinely needs the whole file, say so first. Do not explore the codebase open-endedly.
- **Do not spawn subagents** for work that can be done inline — each starts cold and re-derives context already loaded. Worth it only for wide parallel searches where the conclusion is all that's needed.
- **You cannot verify your own work here.** No test suite, no CI, no data in the repo; validation is visual, via ROOT `TBrowser`, on a machine that has the data. Say explicitly when a change is unverified — never report a physics change as working because it compiled. The dangerous bugs are silent: a flipped sign in a kinematics term, an inverted spline direction, a wrong argument order, an off-by-one in a wire mask. None of them throw. The run completes, the plot renders, and the spectrum is wrong.
- **Review once per coherent change, not per edit.** Scale depth to blast radius: scratch macros, shell, and docs are shallow; `TrackRecon.C`, `MakeVertex.C`, `Armory/`, or any calibration fit is deep and should tolerate uncertain findings — silent failure is what's being hunted. Run `/security-review` only for file I/O or shell interpolation in the `run_*.sh` drivers; it is not a physics tool.

## Build & execution
- ROOT macros (`.C`) compile via ACLiC inside ROOT. Do not link into binaries.
- Standalone binaries (`EventBuilder`, `Mapper`): `cd Armory && make`
- Pre-compile for parallel runs: `root -q -l -b -e '.L TrackRecon.C++O'`
- Batch runs: `./run_17F.sh`, `./run_27Al.sh`, `./run_tr.sh`
- `run_tr.sh` gates its stages with `if [[ 1 -eq 0 ]]` / `if [[ 1 -eq 1 ]]`. Flip the literal to enable a block.
- Data is NOT in the repo. Mapped ROOT files live at `../ANASEN_analysis/data/${DATASET}_Data/${PREFIX}${run}_mapped.root`.

## Runtime configuration
All config is `getenv`-driven; the run scripts `export` it. Defaults are in `TrackRecon.C::Begin()`.

Physics / results-affecting:
`DATASET` `reactiondata` `CO2percent` `pressure_in_torr` `CATHODE_GAIN` `PC_ENERGY_CALIBRATION`
`source_vertex` `CUTLIST` `DITHER_SIGMA` `RNG_SEED` `BEAM_AXIS_X` `BEAM_AXIS_Y`
`timecut_low` `timecut_high` `DISABLE_BAD_ANODE_WIRES` `A1C1_LOWBAND_RFACTOR` `A1C1_CFRAC_SPLIT`
`A1C1_ANODEE_COFF` `A1C1_ANODEE_REF` `A1C1_MISSING_FMAX`
`A1C1_Z_SCALE_QQQ` `A1C1_Z_SCALE_SX3` `A1C1_Z_OFF_QQQ` `A1C1_Z_OFF_SX3`

Plumbing: `OUT_DIR` `RUN_NUMBER` `FLUSH_BARRIER` `MAX_RSS_MB` `MEMCHECK_STRIDE`

**`DEDX_SCALE` is read by `eloss_calculations/Eloss.py`, NOT by `TrackRecon.C`.**
It multiplies `catima.dedx` when the lookup tables are generated. Changing it invalidates every `.dat` table — re-run `Eloss.py` (the run scripts do this automatically before processing). Tables are keyed `eloss_calculations/<species>_lookup_<E>MeV_<P>torr_<CO2>pc.dat` (e.g. `alpha_lookup_50MeV_250torr_3pc.dat`).

## Footguns
- **`HistPlotter` keys `oMap` by histogram NAME only.** The folder argument does not namespace. Two `Fill` calls with the same name in different folders merge silently — no warning, no error, wrong plot. Always make names unique unless required.
- **Si–PC coincidence sits at negative dt.** The real band is roughly `-450 < t_si - t_pc < -200`, not near zero. Do not "correct" a timing gate toward zero.
- **RNG**: one shared `anasenRandom` (`TRandom3`, seed 4357, `RNG_SEED` to override) — do not introduce fresh `TRandom3(0)` instances; it breaks reproducibility.
- **Pass By Reference in Event Loops:** Always iterate via `const auto&` over cluster/hit vectors to prevent expensive struct copies in multi-million event runs.

## Invariants to check on any physics change
Units (MeV vs keV, cm vs mm); sign conventions on `z` and beam direction; eloss spline domain and extrapolation; energy conservation in the kinematics; dead-wire / neighbour masking; RNG seeding. When reviewing, include the `.dat` and `run_*.sh` parameter deltas — the most consequential recent changes have been constants, not code.

## Docs
`README.md` has the pipeline documentation and full file references.