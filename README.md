# ANASEN Analysis

Analysis code for the **Array for Nuclear Astrophysics and Structure with Exotic Nuclei (ANASEN)** detector at FSU. Processes raw FSUNSCL data through event building, channel mapping, calibration, and physics-level vertex reconstruction for transfer reaction experiments (e.g. ²⁷Al(α,p) and ¹⁷F(α,p)).

---

## Detector Overview

| Detector | Type | Role |
|----------|------|------|
| **SX3** | Double-sided silicon strip (barrel, ~88 mm radius) | Energy and Z-position of recoils/ejectiles |
| **QQQ** | Annular silicon (forward, z = 100 mm) | Energy and polar angle of forward-going particles |
| **PC** | Helical proportional counter | dE and azimuthal/Z position via 24 anode + 24 cathode wires |

The PC uses 24 twisted anode wires and 24 cathode wires. Wire geometry, crossover positions, pseudo-wire interpolation, and track-plane intersection are all handled by `Armory/ClassPW.h`.

---

## Full Analysis Chain

```
Raw .fsu files  (FSUNSCL digitizer output)
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│  1. EVENT BUILDING                                              │
│  Binary: EventBuilder  (Armory/EventBuilder.cpp)               │
│  Script: buildEvents.sh  or  ProcessRun.sh <run> <tw> 0        │
│  Input : *.fsu files                                           │
│  Output: Run_NNN_<timeWindow>ns.root  (raw, unmapped TTree)    │
│  Hits within a configurable time window are grouped as events. │
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│  3. PRE-ANALYSIS CHECKS  (optional)                            │
│  Macro: PreAnalysis.C                                          │
│  Plots raw rates and energy spectra per channel from unmapped  │
│  data — useful for identifying dead/noisy channels before      │
│  mapping.                                                      │
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│  2. CHANNEL MAPPING                                             │
│  Binary: Mapper  (Armory/Mapper.cpp)                           │
│  Script: ProcessRun.sh <run> <tw> 0  (calls Mapper internally) │
│  Config: mapping.h                                             │
│  Input : eventbuilt ROOT tree                                  │
│  Output: Run_NNN_mapped.root                                   │
│  Translates hardware (digitizer SN, channel) to logical        │
│  detector identity (SX3/QQQ/PC, strip/wire number).            │
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│  4. CALIBRATION                                                 │
│  ├── sx3cal/EXFit.C / EXFit2.C                                  │
│  │       Fit SX3 front-strip position vs back-strip energy      │
│  │       to extract front/back gain coefficients                │
│  ├── sx3cal/LRFit.C                                             │
│  │       Left-right ratio fit for SX3 position calibration      │
│  │       Output: sx3cal/{17F,27Al}/  (frontgains.dat,           │
│  │               backgains.dat, rightgains.dat per run set)     │
│  ├── GainMatchQQQ.C                                             │
│  │       QQQ ring/wedge gain matching                           │
│  │       Output: qqq_GainMatch.dat                              │
│  ├── Calibration.C                                              │
│  │        Final absolute energy calibration for all detectors   │
│  │        Output: qqq_Calib.dat                                 │
│  ├── PCGainMatch.C                                              │
│  │       PC anode and cathode gain matching                     │
│  │       Output: slope_intercept_*.dat                          │
│  ├── FitHistogramsWithTSpectrum_Sequential_Improved.C           │
│  │       Automated peak finding on PC wire spectra              │
│  │       done separately for anodes and cathodes then combined  │
│  │       (TSpectrum-based; writes centroids for gain match)     │
│  └──QQQ_Calcheck.C                                              │
│      Verifies QQQ calibration against known sources             │
└─────────────────────────────────────────────────────────────────┘
      │
      ▼
┌─────────────────────────────────────────────────────────────────┐
│  5.  VERTEX RECONSTRUCTION  (MakeVertex)                        │
│  Macro : MakeVertex.C  (ROOT TSelector)                         │
│  Batch : run_sx3.sh  (choose dataset, wireflip and offset)      │
│  Input : Run_NNN_mapped.root                                    │
│  Output: result_runNNN.root  (calibrated TTree                  │
│          + diagnostic histograms)                               │
│  Applies energy calibrations, reconstructs PC hit positions     │
│  (anode energy + cathode charge division), clusters wires,      │
│  and builds SX3/QQQ/PC event objects.                           │
│  3D track-reconstruction diagnostics using                      │
│  PC wire geometry and Si hit positions.                         │
└─────────────────────────────────────────────────────────────────┘
```

---

## Channel Mapping (`mapping.h` / `mapping_alpha.h`)

The Mapper binary reads the raw TTree (digitizer serial number + channel) and rewrites it with logical detector labels. `mapping.h` is used for standard beam runs; `mapping_alpha.h` for source/alpha calibration runs where the electronics configuration differs. The active mapping is selected by editing the `#include` at the top of `Armory/Mapper.cpp` before building.

`PrintMapping()` (defined in the mapping header) prints the full channel table at runtime for verification.

---

## MakeVertex — Vertex Reconstruction

`MakeVertex.C` is a ROOT `TSelector` that loops over the analyzed event tree and fills diagnostic and physics histograms. It loads energy-loss lookup tables from `eloss_calculations/` and uses `Armory/Kinematics.h` for two-body reaction kinematics.

### Diagnostic Section Toggles

Comment out any `#define` at the top of `MakeVertex.C` to skip that block at compile time:

```cpp
#define DIAG_WIREMULT   // anode/cathode cluster multiplicity plots
#define DIAG_1WIRE      // raw per-wire dE vs Si (no PC required)
#define DIAG_PC_SX3     // PC-SX3 coincidence analysis
#define DIAG_1A0C_SX3   // 1A0C single-wire vertex with SX3
#define DIAG_1A0C_QQQ   // 1A0C single-wire vertex with QQQ
#define DIAG_NA0C_SX3   // nA0C (n>=2) pseudo-wire vertex with SX3
#define DIAG_NA0C_QQQ   // nA0C (n>=2) pseudo-wire vertex with QQQ
#define DIAG_PC_QQQ     // PC-QQQ coincidence analysis
```

### Histogram Output Folders

| ROOT Folder | Content |
|-------------|---------|
| `wiremult` | Anode/cathode cluster size and multiplicity matrices |
| `1wire` | Raw single-wire dE vs SX3/QQQ energy and phi correlations |
| `hTiming` | PC-Si timing spectra |
| `hPCZSX3` | PC-Z reconstruction with SX3 coincidence |
| `hPCzQQQ` | PC-Z projection with QQQ coincidence |
| `1A0C`,`2A0C`, `3A0C`, … | n anodes (n≥1), 0 cathode: pseudo-wire vertex + excitation spectra |
| `ainterp_noc` | Pseudo-wire interpolation diagnostics (no cathode gate) |
| `phicut` | PC-QQQ with phi coincidence gate |

### Vertex Reconstruction Methods

**nA0C (n ≥ 1 anodes, `DIAG_NA0C_SX3/QQQ`):** All anode clusters are flattened and passed to `GetPseudoWire()` to produce an energy-weighted average wire. `getClosestWirePosAtWirePhi()` then finds the point on that pseudo-wire at the Si phi. Histograms land in a folder named `{n}A0C` (e.g. `2A0C`).

---

## Energy Loss (`eloss_calculations/`)

Lookup tables pre-computed with [pycatima](https://github.com/hrosiak/pycatima) for the isobutane/helium fill gas at 250 Torr. Tables cover protons, alphas, aluminum beam, fluorine, and oxygen. `Eloss.py` regenerates the tables; `make_eloss_table.C` is a ROOT-based alternative generator.

The `cm_to_MeV` / `MeV_to_cm` TGraph pairs loaded in `MakeVertex.C` provide fast range–energy conversion for correcting energy loss in the Si detectors.

---

## Geometry Utilities

| File | Purpose |
|------|---------|
| `grid_generate.py` | Generates `detector_geometry.dat` — calculates SX3, QQQ, and PC wire azimuthal angles from first principles given detector radii and strip positions |
| `shadowplay.py` | Exports ANASEN geometry (wires, SX3, QQQ) as a VTK file + `anasen_labels.csv` for visualisation in ParaView |
| `detector_geometry.dat` | Pre-generated geometry reference (SX3 strip angles, wire positions) |
| `anasen_labels.csv` | Label table for ParaView geometry display |

---

## FEM Simulations (`anasen_fem/`)

Finite-element simulations of the PC electric field using gmsh (meshing), Elmer (FEM solver), and ParaView (visualisation). See `anasen_fem/README.md` for full setup and installation instructions.

```bash
cd anasen_fem
python3 run.py    # mesh → FEM solve → field extraction → ParaView output
```

---

## Batch Processing Scripts

| Script | Purpose |
|--------|---------|
| `buildEvents.sh <run> <tw>` | Build events from raw .fsu for a single run |
| `ProcessRun.sh <run> <tw> <opt>` | Full single-run pipeline: event build (opt=0) or analysis (opt=1) |
| `BatchProcess.sh` | Parallel batch processing over a run range via GNU parallel |
| `process_mapped_run.sh <start> <end>` | Run `TrackRecon.C` in parallel over a range of mapped files |
| `run_sx3.sh/run_27Al.sh/run_17F.sh` | Run `VertexRecon.C` in parallel over a range of mapped files and hAdds them for analysis|

---

## Dependencies

- [ROOT](https://root.cern) ≥ 6.x (TSelector, TChain, TVector3, TF1, TGraph, TSpectrum)
- [GNU parallel](https://www.gnu.org/software/parallel/) — for batch scripts
- [pycatima](https://github.com/hrosiak/pycatima) — for regenerating energy-loss tables
- gmsh + Elmer + ParaView — for FEM simulations only

---

## Key Files Reference

| File | Purpose |
|------|---------|
| `mapping.h` / `mapping_alpha.h` | Hardware→logical channel mapping (beam / alpha-source configs) |
| `Armory/EventBuilder.cpp` | Time-window coincidence event builder (compiled binary) |
| `Armory/Mapper.cpp` | Channel remapper: SN+ch → detector ID (compiled binary) |
| `PreAnalysis.C / .h` | Raw rate and energy checks on mapped data |
| `sx3cal/EXFit.C`, `EXFit2.C` | SX3 front/back gain coefficient extraction |
| `sx3cal/LRFit.C` | SX3 left-right position calibration |
| `sx3cal/{17F,27Al}/` | Per-experiment gain coefficient files (frontgains, backgains, rightgains) |
| `GainMatchQQQ.C` | QQQ gain matching |
| `Calibration.C / .h` | Absolute energy calibration |
| `QQQ_Calcheck.C` | QQQ calibration verification |
| `PCGainMatch.C` | PC anode/cathode gain matching |
| `FitHistogramsWithTSpectrum_Sequential_Improved.C` | Automated PC anode peak fitting |
| `Analyzer.C / .h` | Event-level clustering and PC hit reconstruction, (older in need of rehaul)|
| `Analysis.C` | TChain wrapper to run Analyzer over a run range |
| `processRun.C` | Single-file wrapper to run Analyzer on one ROOT file |
| `MakeVertex.C / .h` | Vertex reconstruction and physics histogram production |
| `TrackRecon.C / .h` | Alternate to VertexRecon, to track changes through code, eventaully plan to clean up one for actual data analysis |
| `RunTimeSummary.C` | Run-by-run timing/rate summary across a run range |
| `Armory/ClassPW.h` | PC wire geometry: crossovers, pseudo-wire, track-plane intersection |
| `Armory/Kinematics.h/.cpp` | Two-body reaction kinematics (excitation energy, Q-value) |
| `Armory/HistPlotter.h` | Histogram management with named ROOT folders |
| `Armory/SX3Geom.h` | SX3 detector geometry |
| `Armory/PC_StepLadder_Correction.h` | PC Z-position slope correction |
| `eloss_calculations/Eloss.py` | Energy-loss table generator (pycatima) |
| `grid_generate.py` | Generates `detector_geometry.dat` from first principles |
| `shadowplay.py` | Exports ANASEN geometry to VTK + CSV for ParaView |
| `detector_geometry.dat` | Pre-generated detector geometry reference |
| `anasen_labels.csv` | Label table for ParaView geometry visualisation |
| `anasen_fem/` | FEM electric field simulations for the PC |