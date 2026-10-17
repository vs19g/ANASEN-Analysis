# Cross-section pipeline (ANASEN active target)

This folder turns data into a cross section for any reaction in the ANASEN active target. There are three steps,
one script each, and one config file per study shared by all three:

```
root -l -b -q 'extract_yields.C+("examples/27Al_ap_runA.toml")'           # 1. ROOT -> yields_<det>.csv
~/myenv/bin/python -B geometric_efficiency.py examples/27Al_ap_runA.toml   # 2. detectors.py -> eps_<det>.csv  ~10 s
~/myenv/bin/python -B cross_section.py        examples/27Al_ap_runA.toml   # 3. -> xs_<det>.csv, xs.png      ~3 s
```

Run them in WSL from this folder (`~/ANASEN_analysis/cross_section`). Step 1 needs ROOT; steps 2 and 3 need
`~/myenv/bin/python`, not ROOT.

| file | what it does |
|---|---|
| `extract_yields.C` | **Step 1, the only step that reads the data.** It reads a TH1/TH2, or fills a histogram from a TTree with a cut, and writes a plain CSV histogram. |
| `geometric_efficiency.py` | Step 2: geometric efficiency per bin, a toy Monte Carlo of the reaction, with every hit rule taken from `detectors.py`. |
| `detectors.py` | **The detectors as they are now** (Run A): QQQ3, the four analysed SX3s, the counter, the beam tube and the PCB hole, plus the hit rules. To change the geometry, edit this file. |
| `cross_section.py` | Step 3: reads only CSVs (the yields and an efficiency) and computes σ(E_cm) or dσ/dΩ(θ_cm). |
| `common.py` | What the three scripts share: the config, the beam track E(z), the bins. |
| `examples/` | One commented config per study: `27Al_ap_runA.toml` and `paap_ta_runs_hist.toml`. |
| `examples/reference/` | Barnard 1964 (EXFOR A1113) for p + ⁴He. |
| `inputs/` | The beam track (`beam_tracks.csv`) and the TALYS states (`talys_27Al_ap.csv`). The p + ⁴He track is read from `../eloss_calculations/`. |
| `out/<study>/` | Each study's outputs: `yields_<det>.csv`, `eps_<det>.csv`, `xs_<det>.csv`, `xs.png`. |

## Starting a new study

1. Copy the closest example in `examples/` (excitation function: `27Al_ap_runA.toml`; angular distribution:
   `paap_ta_runs_hist.toml`) and give it its own `out_dir`.
2. Fill in `[extract]` for your histogram or tree, then run step 1 and look at `yields_<det>.csv`.
3. Fill in `[reaction]`, `[beam]` and `[bins]`, then run step 2. Check that `eps_<det>.csv` is non-zero where you
   expect events.
4. Fill in `[target]`, `[acceptance]` and any `[[reference]]`, then run step 3.

If the detectors have changed, edit `detectors.py`, not the configs.

![27Al(a,p)](out/27Al_ap/xs.png)

## The formula

| mode | result | formula |
|---|---|---|
| `excitation` | σ(E_cm), mb | Y / (N_b · n · dz · ε) |
| `angular` | dσ/dΩ(θ_cm) in one E_cm (or z) window, mb/sr | Y / (N_b · n · dz · ε · ΔΩ) |

- **Y:** the yield in the bin.
- **N_b:** the number of beam particles.
- **n:** target atoms per cm² per mm of path.
- **dz:** the length of beam path over which the beam's E_cm lies in the bin, read off the beam track E(z).
- **ε:** the efficiency of the bin.

## Step 1: `extract_yields.C`, data to CSV

`[extract] source` says where the events are:
- `"hist"`: a TH1 or TH2 in a ROOT file. A TH2 is projected on the axis that holds `variable`, or kept 2D with
  `angle_axis` for the angular mode.
- `"tree"`: a ROOT TTree (`file` may be a glob or a list of runs), filled into `bins` after `cut`.

Cuts use `and` / `or` (the macro turns them into `&&` / `||`) and ROOT's TTree::Draw syntax.

**Histogram paths** are `"folder/name"`.
- **They are case-sensitive.** Production files of different vintages spell the detector differently: the 27Al
  outputs use `_m27Alax+misc_qqq/...`, the proton runs' `Output_p/Al_protons.root` uses `elastic+misc_QQQ_a/...`.
  Put the spelling in `[det_vars]`:
  ```toml
  dets = ["qqq", "sx3"]                # the names detectors.py knows
  [det_vars.qqq]
  DET = "QQQ"
  [det_vars.sx3]
  DET = "SX3"
  [extract]
  hist = "elastic+misc_{DET}_a/elastic_VertexReconZ_a_{DET}"
  ```
- **A name that itself contains `/`** (e.g. `dE3_E_AnodeQQQR_A1C2_(vertex_fix_z/100)=0.000000`) is found too.
- **If the path doesn't exist,** the macro stops and names the closest histogram paths in the file.

The CSV it writes, `<out_dir>/yields_<det>.csv`, is the hand-off between the steps. `cross_section.py` reads nothing
else, so a yield made any other way works too as long as it is in this format:

```
1D:  <var>_lo, <var>_hi, counts, variance
2D:  <var>_lo, <var>_hi, theta_cm_lo, theta_cm_hi, counts, variance     (angular mode)
```

`<var>` is `z` (vertex z, mm), `e_beam` (MeV) or `ecm` (MeV). The bins may be finer than the analysis bins and need
not line up with them. Step 3 shares each one out by overlap.

It is run like this:

```
root -l -b -q 'extract_yields.C+("examples/27Al_ap_runA.toml")'
root -l -b -q 'extract_yields.C+("examples/27Al_ap_runA.toml", "out/elsewhere")'   # optional output folder
```

- **Trees:** it fills the histogram with `TTree::Draw` on a `TChain`, so `file` may be a glob or a list of runs.
- **Cuts:** write them with `and` / `or` (e.g. `"det == 1 and e_si > 0.5"`). The macro turns them into `&&` / `||`.
- **TOML:** the macro has a small reader that handles only what these configs use: strings, numbers, one-line lists
  and `{start, stop, step}` tables. It stops with an error on anything else.
- **Errors:** a missing file or tree, or a tree with no entries, stops it with exit code 1 rather than writing zeros.
- **One edge case:** an entry exactly on a tree histogram's last upper edge goes to ROOT's overflow.

## Step 2: `geometric_efficiency.py`, ε from the detector geometry

The toy runs once per bin:
- vertices uniform over the bin's path (angular mode: over the window);
- beam energy at each vertex from the nominal track;
- the residual state drawn from fixed weights, or from TALYS state cross sections at the bin's E_cm;
- isotropic CM emission (angular mode: inside the θ_cm bin);
- relativistic two-body kinematics and straight tracks.

The particle the silicon sees can be the ejectile or the recoil (`detect`). θ_cm is always the ejectile's, so with
`detect = "recoil"` the bins are still in the ejectile's CM angle. For p + ⁴He that is the proton's angle, Barnard's
convention, which TrackRecon's `elastic_VertexZ_vs_ThetaCM` histograms also use. The output has one column per
condition, each times the azimuthal coverage:

| column | the particle... |
|---|---|
| `eps_si` | reaches the analysed silicon (solid angle) |
| `eps_si_thr` | ...and has a vertex energy above the Si threshold. This is the ch5 "solid angle" convention: a kinematic cut, not geometry. |
| `eps_counter` | `eps_si`, and it crosses the guard, anode and cathode radii inside the counter |
| **`eps_geometric`** | `eps_counter`, and it clears the beam tube and the PCB hole. **This is the geometric efficiency.** |

Each column has a Monte Carlo error, `<column>_rel_err`. Gas energy loss, Si thresholds after the gas and the counter's
deposit threshold are **not** included, because they are not geometry. That is why the adopted `eps_tube` of
`../scratch/xs_27Al_lowE_2026-10-15` is lower than `eps_geometric`; to use it instead, see the commented
`[acceptance]` block in the example. The angular distribution is isotropic.

## Step 3: `cross_section.py`

`[acceptance] source` picks ε:
- `"geometric"`: step 2's file; `column` chooses which efficiency.
- `"csv"`: any file, matched on a bin edge or interpolated.
- `"value"`: one number for every bin.
- `"none"`: ε = 1, which gives σ_app.

`[[acceptance.extra]]` multiplies in further factors, such as an intrinsic efficiency. The output `xs_<det>.csv` keeps
every intermediate number: the bins, z ranges, dz, counts, σ_app, ε, σ and its errors (stat, ε, total), and the
references. The header lines record every input. The N_b and n uncertainties are a common scale, so the header gives
them as a percentage and they stay out of the error columns.

## The config (TOML)

Paths are relative to the config file, and `~` is the home folder. `{det}` is replaced by each entry of `dets`, and
so is every key of `[det_vars.<det>]`. Commented templates in `examples/`:
- `27Al_ap_runA.toml`: a TH2 from the merged Run A production file, all three steps.
- `paap_ta_runs_hist.toml`: p + ⁴He elastic **straight from a TrackRecon.C production file**, as for 27Al. It reads
  the TH2 `elastic_VertexZ_vs_ThetaCM_cath_exgate_a_<DET>` (vertex z × proton θ_cm). `protonAlphaElastic_core` fills
  that histogram since 2026-10-16 in every species folder, in four variants: all, `_cath`, `_exgate`, `_cath_exgate`.
  - Input: `~/ANASEN_analysis/Output_p/results_run022.root` (run 22 alone, run_tr.sh block 5 with CUTLIST),
    N_b = 1.35 pA × 672.65 s = 5.67 × 10⁹, vertex window z = −170 to −53 mm.
  - On 2026-10-16 the merged runs 18–22 agreed with a per-hit route to 2–4% in QQQ3 counts and about 3% in dσ/dΩ.

| section | used by | keys |
|---|---|---|
| top | all | `title`, `mode` (`excitation` / `angular`), `dets` (names in `detectors.py`: `qqq`, `sx3`), `out_dir` |
| `[extract]` histogram | 1 | `source = "hist"`, `file`, `hist`, `variable`, `axis` (the axis holding `variable`), `angle_axis` (keep it 2D), `other_range` |
| `[extract]` tree | 1 | `source = "tree"`, `file` (glob or list), `tree`, `variable`, `branch`, `bins`, `cut`, `weight_branch`, `angle_branch` + `angle_bins` |
| `[reaction]` | 2, 3 | `beam_mass_u`, `target_mass_u`; for step 2 also `ejectile_mass_u`, `recoil_mass_u`, `detect`, and `states` *or* `states_csv`; optional `ecm_per_ebeam` |
| `[beam]` | 2, 3 | track: `track_csv` + `track_e_col` (+ `nominal_e_col`, `track_z_col`), *or* `eloss_table` + `e_entrance_mev` + `z_entrance_mm`. N_b: `n_beam` *or* `charge_c` *or* `current_a` + `time_s`; `charge_state`, `n_beam_rel_err`. For Run A runs, `../runlog/beam_exposure.py` gives `n_beam` for any set of runs from the elog's cup readings |
| `[bins]` | 2, 3 | `ecm_edges` (a list or `{start, stop, step}`); angular: `theta_edges_deg` and `window_ecm` *or* `window_z_mm` |
| `[efficiency]` | 2 | `n_per_bin` (default 400 000), `seed` |
| `[target]` | 3 | `pressure_torr`, `temperature_k`, `mole_fraction`, `atoms_per_molecule`, `rel_err`; *or* `n_per_cm2_per_mm` |
| `[yield]` | 3 | `csv`: only to read yields from somewhere other than step 1's output |
| `[acceptance]` | 3 | `source` and its keys (see above), `[[acceptance.extra]]` |
| `[[reference]]` | 3 | `csv`, `x_col`, `y_col`, `label`, `scale`, `dets` |
| `[plot]` | 3 | `yscale`, `xlim`, `ylim`, `xlabel`, `ylabel` |

## Example: p + ⁴He elastic scattering, Ta run 22 (`examples/paap_ta_runs_hist.toml`)

```
root -l -b -q 'extract_yields.C+("examples/paap_ta_runs_hist.toml")'
~/myenv/bin/python -B geometric_efficiency.py examples/paap_ta_runs_hist.toml
~/myenv/bin/python -B cross_section.py        examples/paap_ta_runs_hist.toml
```

dσ/dΩ(θ_cm) from the α recoils (`detect = "recoil"`, bins in the proton's θ_cm, Barnard's convention), against
Barnard, Jones & Weil (1964). Histogram: `elastic_VertexZ_vs_ThetaCM_cath_exgate_a_<DET>` (cathode topologies,
graphical α cut, |E_x| < 1 MeV). The merged runs 18–22 gave QQQ3 0.44–0.54 × Barnard at 132–156°; SX3 is at 1–3%.
Not in `detectors.py`: the Ta holder rod (~28% of the azimuth downstream of each foil) and Si thresholds after the gas.

## Checks (run in the Windows repo, `/mnt/c/Users/sengi/ANASEN_analysis/cross_section`; the scripts are not in this copy)

- 27Al(α,p), Run A: yields, dz and nominal E_cm equal the earlier study's in every bin; ε_si_thr agrees with its
  independent toy within the MC errors.
- Synthetic events with a known σ, ε and dσ/dΩ are returned with unit-Gaussian pulls.
- `extract_yields.C` and the Python version gave byte-identical yields on every example.

**Not verified:** the physics inputs are whatever the config and `detectors.py` say. The checks show the arithmetic and
the bookkeeping are right, and that `detectors.py` reproduces the geometry of `anasen_detectability.py`. They cannot
show that the geometry, a track or an N_b is right. In particular:
- The dE/dx scale of a `.dat` table is not written in the file. The live 250 Torr / 3% tables were regenerated on
  2026-10-16 and match the 0.962 track to 3 keV.
- `detectors.py` gives QQQ3 0.967 of the azimuth, but TrackRecon drops QQQ pixels without a calibration, which leaves
  about 0.675. QQQ3 ε is therefore about 1.4× too high, and its cross sections about 1.4× too low.
- Where only a sliver of solid angle is open, ε is tiny and σ explodes. Example: QQQ3 at E_cm 5.4 MeV, ε = 1.5 × 10⁻⁵,
  behind the beam tube. The error bar shows it; no ε cut is applied.