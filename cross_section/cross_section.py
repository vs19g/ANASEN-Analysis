"""Step 3: the cross section, from the yields CSV (step 1) and an efficiency (step 2, or any other).

    ~/myenv/bin/python -B cross_section.py examples/27Al_ap_runA.toml       # no ROOT needed; any python with pandas

    sigma(E_cm)        = Y / (N_b * n * dz * eps)                  mode = "excitation"
    dsigma/dOmega(th)  = Y / (N_b * n * dz * eps * dOmega)         mode = "angular"

    Y       yield in the bin, from <out_dir>/yields_<det>.csv (or [yield] csv); see extract_yields.py for the format
    N_b     number of beam particles                              [beam]
    n       target atoms per cm^2 per mm of beam path             [target]
    dz      length of beam path (mm) over which the beam's E_cm lies in the bin, from the beam track E(z)
    eps     efficiency of the bin                                 [acceptance]
    dOmega  solid angle of the theta_cm bin (angular mode)

Writes <out_dir>/xs_<det>.csv (every intermediate number; the header records every input) and <out_dir>/xs.png.

WHAT IT DOES, STEP BY STEP (see analyse_detector())
    1. beam track:     beam energy E at every z
    2. normalisation:  N_b and n
    3. bins:           E_cm bins -> their z ranges and dz  (or theta_cm bins inside one energy window)
    4. yields:         share the yield histogram's bins into those bins
    5. efficiency:     eps for each bin
    6. cross section:  the formula, with statistical and efficiency errors
"""
import argparse
import datetime
import re
from pathlib import Path

import numpy as np
import pandas as pd

import common

BOLTZMANN = 1.380649e-23             # J/K
PASCAL_PER_TORR = 133.322368
ELEMENTARY_CHARGE = 1.602176634e-19  # C
CM2_PER_MB = 1e-27
YIELD_VARIABLES = ("z", "e_beam", "ecm")


# ============================================================ normalisation

def number_of_beam_particles(beam):
    """N_b from: n_beam directly, or integrated charge, or current x time. Returns (N_b, relative error, text)."""
    charge_state = beam.get("charge_state", 1)
    if "n_beam" in beam:
        n_b, text = beam["n_beam"], "given"
    elif "charge_c" in beam:
        n_b = beam["charge_c"] / (charge_state * ELEMENTARY_CHARGE)
        text = f"{beam['charge_c']} C / ({charge_state} e)"
    elif "current_a" in beam:
        n_b = beam["current_a"] * beam["time_s"] / (charge_state * ELEMENTARY_CHARGE)
        text = f"{beam['current_a']} A x {beam['time_s']} s / ({charge_state} e)"
    else:
        raise KeyError("[beam] needs n_beam, charge_c, or current_a and time_s")
    return float(n_b), float(beam.get("n_beam_rel_err", 0.0)), text


def target_atoms_per_cm2_per_mm(target):
    """n = (gas molecules per cm^3) x (mole fraction of the target molecule) x (target atoms per molecule) x 0.1 cm/mm.
    E.g. He in 250 Torr He + 3% CO2 at 293.15 K: mole_fraction 0.97, atoms_per_molecule 1 -> 7.988e17.
    Returns (n, relative error, text)."""
    if "n_per_cm2_per_mm" in target:
        return float(target["n_per_cm2_per_mm"]), float(target.get("rel_err", 0.0)), "given"
    pressure_pa = target["pressure_torr"] * PASCAL_PER_TORR
    molecules_per_cm3 = pressure_pa / (BOLTZMANN * target["temperature_k"]) * 1e-6
    fraction = target.get("mole_fraction", 1.0)
    atoms = target.get("atoms_per_molecule", 1)
    n = molecules_per_cm3 * fraction * atoms * 0.1
    text = f"{target['pressure_torr']} Torr, {target['temperature_k']} K, mole fraction {fraction}, {atoms} atom(s)/molecule"
    return n, float(target.get("rel_err", 0.0)), text


# ============================================================ yields

def read_yields(file):
    """The yields CSV of extract_yields.py. Returns (table, variable), variable = z, e_beam or ecm."""
    table = pd.read_csv(file, comment="#")
    for variable in YIELD_VARIABLES:
        if f"{variable}_lo" in table.columns:
            return table, variable
    raise ValueError(f"{file}: no column z_lo, e_beam_lo or ecm_lo (see extract_yields.py for the format)")


def in_yield_units(variable, k, ecm_range, z_range):
    """Convert a bin, given both as an E_cm range and a z range, into the units of the yield's variable."""
    if variable == "ecm":
        return ecm_range
    if variable == "e_beam":
        return ecm_range[0] / k, ecm_range[1] / k
    return z_range


def fraction_inside(src_lo, src_hi, lo, hi):
    """Fraction of each yield bin [src_lo, src_hi) that lies inside [lo, hi). Used to share the yield's bins
    between the analysis bins when their edges do not line up."""
    return np.clip(np.minimum(src_hi, hi) - np.maximum(src_lo, lo), 0, None) / (src_hi - src_lo)


def yields_excitation(table, variable, t, k):
    """Yield and its error in each E_cm bin."""
    if "theta_cm_lo" in table.columns:          # a 2D yield: add up the theta bins first
        table = table.groupby([f"{variable}_lo", f"{variable}_hi"], as_index=False)[["counts", "variance"]].sum()
    src_lo, src_hi = table[f"{variable}_lo"].to_numpy(), table[f"{variable}_hi"].to_numpy()
    lo, hi = in_yield_units(variable, k, (t["ecm_lo_mev"], t["ecm_hi_mev"]), (t["z_up_mm"], t["z_down_mm"]))
    counts, errors = [], []
    for a, b in zip(lo, hi):
        f = fraction_inside(src_lo, src_hi, a, b)
        counts.append(np.sum(f * table["counts"]))
        errors.append(np.sqrt(np.sum(f ** 2 * table["variance"])))
    return np.array(counts), np.array(errors)


def yields_angular(table, variable, t, window, k):
    """Yield and its error in each theta_cm bin, for the part of the yield inside the energy (or z) window."""
    if "theta_cm_lo" not in table.columns:
        raise ValueError("angular mode needs a 2D yield (theta_cm_lo/hi columns): set angle_axis or angle_branch in [extract]")
    ecm_lo, ecm_hi, z_up, z_down = window
    w_lo, w_hi = in_yield_units(variable, k, (ecm_lo, ecm_hi), (z_up, z_down))
    in_window = fraction_inside(table[f"{variable}_lo"].to_numpy(), table[f"{variable}_hi"].to_numpy(), w_lo, w_hi)
    th_lo, th_hi = table["theta_cm_lo"].to_numpy(), table["theta_cm_hi"].to_numpy()
    counts, errors = [], []
    for a, b in zip(t["theta_lo_deg"], t["theta_hi_deg"]):
        f = in_window * fraction_inside(th_lo, th_hi, a, b)
        counts.append(np.sum(f * table["counts"]))
        errors.append(np.sqrt(np.sum(f ** 2 * table["variance"])))
    return np.array(counts), np.array(errors)


# ============================================================ efficiency

def efficiency(config, det, t, full_path):
    """eps and its relative error for every bin. [acceptance] source =
         "geometric": <out_dir>/eps_<det>.csv from geometric_efficiency.py; `column` picks the efficiency
                      (default eps_geometric), its Monte Carlo error comes with it
         "csv":       any file: csv, x_col, eps_col, match = "bin" (same bin edge: our column `key` = their x_col)
                      or "interp" (interpolate at our column `key`; NaN outside the file's range), rel_err_col or rel_err
         "value":     one number for every bin (value, rel_err)
         "none":      eps = 1, the apparent cross section (the default)
       [[acceptance.extra]] name, value, rel_err: further factors, multiplied in (e.g. an intrinsic efficiency)."""
    eff = config.get("acceptance", {})
    mode = config.get("mode", "excitation")
    source = eff.get("source", "none")
    n_bins = len(t)
    eps, rel_err = np.ones(n_bins), np.zeros(n_bins)
    text = "none (eps = 1: apparent cross section)"

    if source == "geometric":
        column = eff.get("column", "eps_geometric")
        file = common.output_file(config, full_path, "eps", det)
        eff = {"csv": str(file), "match": "bin", "eps_col": column, "rel_err_col": f"{column}_rel_err",
               "x_col": "ecm_lo_mev" if mode == "excitation" else "theta_lo_deg", "extra": eff.get("extra", [])}
        source = "csv"
    if source == "csv":
        table = pd.read_csv(full_path(eff["csv"]), comment="#")
        match = eff.get("match", "interp")
        default_key = {("excitation", "bin"): "ecm_lo_mev", ("excitation", "interp"): "ecm_nominal_mev",
                       ("angular", "bin"): "theta_lo_deg", ("angular", "interp"): "theta_c_deg"}[(mode, match)]
        key = eff.get("key", default_key)
        their_x = table[eff["x_col"]].to_numpy(float)
        their_eps = table[eff["eps_col"]].to_numpy(float)
        if "rel_err_col" in eff:
            their_err = table[eff["rel_err_col"]].to_numpy(float)
        else:
            their_err = np.full(len(table), eff.get("rel_err", 0.0))

        ours = t[key].to_numpy(float)
        if match == "bin":
            for i, x in enumerate(ours):
                same = np.flatnonzero(np.isclose(their_x, x, rtol=1e-6, atol=1e-9))
                eps[i] = their_eps[same[0]] if len(same) else np.nan
                rel_err[i] = their_err[same[0]] if len(same) else np.nan
        else:
            order = np.argsort(their_x)
            eps = np.interp(ours, their_x[order], their_eps[order], left=np.nan, right=np.nan)
            rel_err = np.interp(ours, their_x[order], their_err[order], left=np.nan, right=np.nan)
        text = f"{eff['csv']}, column {eff['eps_col']} ({match} on {eff['x_col']} = our {key})"
    elif source == "value":
        eps = np.full(n_bins, float(eff["value"]))
        rel_err = np.full(n_bins, float(eff.get("rel_err", 0.0)))
        text = f"constant {eff['value']}"
    elif source != "none":
        raise ValueError(f"[acceptance] source must be geometric, csv, value or none, not {source!r}")

    for extra in eff.get("extra", []):
        eps = eps * extra["value"]
        rel_err = np.hypot(rel_err, extra.get("rel_err", 0.0))
        text += f" x {extra.get('name', 'extra factor')} {extra['value']}"
    return eps, rel_err, text


# ============================================================ cross section

def cross_section(t, n_b, n, eps, eps_rel_err, name):
    """The formula. `name` is "sigma_mb" (excitation) or "dsdo_mb_sr" (angular). Adds columns to t:
         <name>_app, <name>_app_err      apparent cross section (eps = 1), statistical error
         <name>, <name>_err_stat         divided by eps
         <name>_err_eps, <name>_err      error from eps; total (stat and eps in quadrature)
       The N_b and n uncertainties are a common scale for all bins and are NOT in these columns."""
    exposure = n_b * n * t["dz_mm"] * t["d_omega_sr"] * CM2_PER_MB      # events per mb (per mb/sr)
    exposure = exposure.where(exposure > 0)                             # dz = 0 (past the stopping point) -> NaN
    eps = pd.Series(eps, index=t.index).where(lambda e: e > 0)          # eps = 0 -> NaN
    t[f"{name}_app"] = t["counts"] / exposure
    t[f"{name}_app_err"] = t["stat_err"] / exposure
    t["eps"] = eps
    t["eps_rel_err"] = eps_rel_err
    t[name] = t[f"{name}_app"] / eps
    t[f"{name}_err_stat"] = t[f"{name}_app_err"] / eps
    t[f"{name}_err_eps"] = t[name] * t["eps_rel_err"]
    t[f"{name}_err"] = np.hypot(t[f"{name}_err_stat"], t[f"{name}_err_eps"])
    return t


def add_references(t, references, det, x_col, full_path):
    """Reference curves ([[reference]]: csv, x_col, y_col, label, optional scale and dets) interpolated at our x
    (NaN outside their range). Returns the curves for the plot."""
    curves = []
    for ref in references or []:
        if "dets" in ref and det not in ref["dets"]:
            continue
        table = pd.read_csv(full_path(ref["csv"]), comment="#").sort_values(ref["x_col"])
        x = table[ref["x_col"]].to_numpy(float)
        y = ref.get("scale", 1.0) * table[ref["y_col"]].to_numpy(float)
        label = ref.get("label", Path(ref["csv"]).stem)
        column = "ref_" + re.sub(r"\W+", "_", label).strip("_").lower()
        t[column] = np.interp(t[x_col], x, y, left=np.nan, right=np.nan)
        curves.append((label, x, y))
    return curves


def analyse_detector(config, det, full_path):
    """Steps 1-6 for one detector. Returns the result table, its header lines, the x column, the cross-section
    column name and the reference curves."""
    mode = config.get("mode", "excitation")
    k = common.ecm_per_beam_energy(config["reaction"])

    # 1. beam track
    track = common.load_beam_track(config["beam"], full_path)

    # 2. normalisation
    n_b, n_b_rel_err, n_b_text = number_of_beam_particles(config["beam"])
    n, n_rel_err, n_text = target_atoms_per_cm2_per_mm(config["target"])

    # 3. bins and 4. yields
    yield_file = full_path(config["yield"]["csv"]) if "csv" in config.get("yield", {}) \
        else common.output_file(config, full_path, "yields", det)
    yields, variable = read_yields(yield_file)
    t, window = common.make_bins(config, track, k)
    if mode == "excitation":
        t["counts"], t["stat_err"] = yields_excitation(yields, variable, t, k)
        name, x_col, window_text = "sigma_mb", "ecm_nominal_mev", ""
    else:
        t["counts"], t["stat_err"] = yields_angular(yields, variable, t, window, k)
        name, x_col = "dsdo_mb_sr", "theta_c_deg"
        ecm_lo, ecm_hi, z_up, z_down = window
        window_text = f"; window E_cm {ecm_lo:.4g} to {ecm_hi:.4g} MeV, z {z_up:.2f} to {z_down:.2f} mm"

    # 5. efficiency and 6. cross section
    eps, eps_rel_err, eps_text = efficiency(config, det, t, full_path)
    t = cross_section(t, n_b, n, eps, eps_rel_err, name)
    curves = add_references(t, config.get("reference"), det, x_col, full_path)

    header = [
        f"{config.get('title', '')} {det.upper()}: {mode} mode{window_text}",
        f"yield: {yield_file} (binned in {variable})",
        f"beam track: {track.description}; E_cm = {k:.6g} x E_beam",
        f"N_b = {n_b:.6g} ({n_b_text}); n = {n:.6g} target atoms per cm2 per mm ({n_text})",
        f"eps: {eps_text}",
        f"common scale uncertainty from N_b and n: {100 * np.hypot(n_b_rel_err, n_rel_err):.3g} % (not in the error columns)",
        f"written by cross_section.py on {datetime.date.today().isoformat()}",
    ]
    return t, header, x_col, name, curves


def make_plot(results, config, out_file):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    opts = config.get("plot", {})
    excitation = config.get("mode", "excitation") == "excitation"
    fig, ax = plt.subplots(figsize=opts.get("figsize", (9, 6)))
    plotted_refs = set()
    for det, (t, x_col, name, curves) in results.items():
        ok = np.isfinite(t[name]) & (t[name] > 0)
        ax.errorbar(t.loc[ok, x_col], t.loc[ok, name], yerr=t.loc[ok, f"{name}_err"],
                    fmt="o", ms=4, capsize=2, label=det.upper() or "data")
        for label, x, y in curves:
            if label not in plotted_refs:
                ax.plot(x, y, "-", lw=1.2, label=label)
                plotted_refs.add(label)
    ax.set_xlabel(opts.get("xlabel", "E$_\\mathrm{cm}$ (MeV)" if excitation else "$\\theta_\\mathrm{cm}$ (deg)"))
    ax.set_ylabel(opts.get("ylabel", "$\\sigma$ (mb)" if excitation else "d$\\sigma$/d$\\Omega$ (mb/sr)"))
    ax.set_yscale(opts.get("yscale", "log"))
    if "xlim" in opts:
        ax.set_xlim(*opts["xlim"])
    if "ylim" in opts:
        ax.set_ylim(*opts["ylim"])
    ax.set_title(config.get("title", ""))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_file, dpi=110)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Step 3: cross section from the yields CSV and an efficiency.")
    parser.add_argument("config", help="TOML config file")
    args = parser.parse_args(argv)
    config, full_path = common.read_config(args.config)

    results = {}
    for det in config.get("dets", [""]):
        print(f"--- {det or 'data'}")
        det_config = common.config_for_detector(config, det)
        t, header, x_col, name, curves = analyse_detector(det_config, det, full_path)
        common.write_csv(t, header, common.output_file(config, full_path, "xs", det))
        results[det] = (t, x_col, name, curves)
        print("\n".join(header[:-1]))
        show = [x_col, "dz_mm" if name == "sigma_mb" else "d_omega_sr", "counts", f"{name}_app", "eps", name,
                f"{name}_err"] + [c for c in t.columns if c.startswith("ref_")]
        print(t[show].to_string(index=False, float_format=lambda v: f"{v:.4g}"), "\n")

    out_png = common.output_file(config, full_path, "xs", "").with_suffix(".png")
    make_plot(results, config, out_png)
    print(f"written to {out_png.parent}")
    return results


if __name__ == "__main__":
    main()
