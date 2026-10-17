"""Step 2: geometric efficiency of every bin, for the detectors as they are now (detectors.py).

    ~/myenv/bin/python -B geometric_efficiency.py examples/27Al_ap_runA.toml      # ~1 min; numpy and pandas only

A toy Monte Carlo per bin (same bins as the cross section):
  - vertices uniform over the bin's stretch of beam path (angular mode: over the window);
  - beam energy at each vertex from the nominal beam track;
  - the residual nucleus' state drawn from [reaction] states: fixed weights, or state cross sections
    from a table (e.g. TALYS) at the bin's nominal E_cm;
  - emission isotropic in the CM (angular mode: cos theta_cm uniform inside the theta bin, theta_cm of the ejectile);
  - relativistic two-body kinematics, straight tracks from the vertex.

Writes <out_dir>/eps_<det>.csv, one row per bin, with these efficiencies (each times the azimuthal coverage):
    eps_si          the particle reaches the analysed silicon                 (solid angle)
    eps_si_thr      eps_si, and its energy at the vertex is above the Si threshold   (the ch5 "solid angle" convention;
                    a kinematic cut, not geometry)
    eps_counter     eps_si, and it crosses the guard, anode and cathode radii inside the counter
    eps_geometric   eps_counter, and it clears the beam tube and the PCB hole    <- THE GEOMETRIC EFFICIENCY
and the Monte Carlo relative error of each (<name>_rel_err).
Not included, because they are not geometry: energy loss in the gas, the Si threshold after the gas, the counter's
deposit threshold, intrinsic efficiencies. They are why the adopted eps_tube of ../scratch/xs_27Al_lowE_2026-10-15 is lower.

CONFIG
    [reaction]    beam_mass_u, target_mass_u, ejectile_mass_u, recoil_mass_u   (atomic masses; Q from these)
                  detect = "ejectile" (default) or "recoil": the particle the silicon sees
                  states = [{ex_mev = 0.0, weight = 1.0}, ...]           fixed state weights, or
                  states_csv = file with ecm_mev and one column ex_<Ex in MeV> per state (state cross sections)
    [efficiency]  n_per_bin (default 400000), seed
    [beam], [bins], mode, dets as for the cross section.
"""
import argparse

import numpy as np
import pandas as pd

import common
import detectors

U = 931.49410242          # MeV per u
EFFICIENCIES = ("eps_si", "eps_si_thr", "eps_counter", "eps_geometric")


def two_body(t_beam, reaction, ex, cos_cm):
    """Relativistic two-body kinematics, beam on target at rest -> ejectile (at cos_cm in the CM) + recoil (state ex).
    Returns the lab kinetic energy (MeV) and polar angle (rad) of the detected particle, and whether the channel is
    open at this beam energy."""
    m1, m2 = reaction["beam_mass_u"] * U, reaction["target_mass_u"] * U
    m3, m4 = reaction["ejectile_mass_u"] * U, reaction["recoil_mass_u"] * U + ex
    s = (m1 + m2) ** 2 + 2 * m2 * t_beam                      # invariant mass squared
    sqrt_s = np.sqrt(s)
    is_open = sqrt_s > m3 + m4
    e_total = t_beam + m1 + m2
    beta = np.sqrt(t_beam ** 2 + 2 * t_beam * m1) / e_total   # of the CM in the lab
    gamma = e_total / sqrt_s
    p_cm = np.sqrt(np.clip((s - (m3 + m4) ** 2) * (s - (m3 - m4) ** 2), 0, None)) / (2 * sqrt_s)
    sin_cm = np.sqrt(1 - cos_cm ** 2)

    # the recoil goes the opposite way in the CM
    mass, sign = (m3, 1.0) if reaction.get("detect", "ejectile") == "ejectile" else (m4, -1.0)
    e_cm = np.sqrt(p_cm ** 2 + mass ** 2)
    p_z = gamma * (sign * p_cm * cos_cm + beta * e_cm)
    energy = gamma * (e_cm + beta * sign * p_cm * cos_cm)
    theta = np.arctan2(p_cm * sin_cm, p_z)
    return np.where(is_open, energy - mass, 0.0), theta, is_open


def q_value_mev(reaction):
    return (reaction["beam_mass_u"] + reaction["target_mass_u"]
            - reaction["ejectile_mass_u"] - reaction["recoil_mass_u"]) * U


def state_weights(reaction, ecm, states_table):
    """Excitation energies of the residual nucleus and their probabilities at this E_cm. States that are closed
    (Ex >= Q + E_cm) get weight 0; if none is open, the ground state alone."""
    if states_table is None:
        states = reaction.get("states", [{"ex_mev": 0.0, "weight": 1.0}])
        ex = np.array([s["ex_mev"] for s in states], float)
        w = np.array([s["weight"] for s in states], float)
    else:
        columns = [c for c in states_table.columns if c.startswith("ex_")]
        ex = np.array([float(c[3:]) for c in columns])
        # the table row nearest this E_cm, among rows with any cross section (below the first such row, its weights)
        rows = states_table[states_table[columns].sum(axis=1) > 0]
        row = rows.iloc[np.argmin(np.abs(rows["ecm_mev"].to_numpy() - max(ecm, rows["ecm_mev"].min())))]
        w = row[columns].to_numpy(float)
    w = np.where(ex < q_value_mev(reaction) + ecm, w, 0.0)
    if w.sum() <= 0:
        w = (ex == 0).astype(float)
    return ex, w / w.sum()


def efficiency_of_bin(det, reaction, track, z_range, cos_range, ecm, states_table, n, rng):
    """Toy for one bin: n reactions with vertices uniform in z_range and cos(theta_cm) uniform in cos_range."""
    z = rng.uniform(z_range[0], z_range[1], n)
    t_beam = track.energy_at(z, nominal=True)
    ex_states, weights = state_weights(reaction, ecm, states_table)
    ex = rng.choice(ex_states, size=n, p=weights)
    cos_cm = rng.uniform(cos_range[0], cos_range[1], n)
    energy, theta, is_open = two_body(t_beam, reaction, ex, cos_cm)
    is_open &= t_beam > 0

    on_si = is_open & detectors.hits_silicon(det, z, theta)
    passed = {"eps_si": on_si,
              "eps_si_thr": on_si & (energy > detectors.SI_THRESHOLD_MEV[det]),
              "eps_counter": on_si & detectors.crosses_counter(z, theta)}
    passed["eps_geometric"] = passed["eps_counter"] & detectors.clears_apertures(z, theta)

    phi = detectors.PHI_FRACTION[det]
    row = {"ground_state_weight": float(weights[ex_states == 0].sum()), "frac_open": float(is_open.mean())}
    for name, ok in passed.items():
        p = ok.mean()
        row[name] = phi * p
        row[f"{name}_rel_err"] = np.sqrt(p * (1 - p) / n) / p if p > 0 else np.nan     # binomial
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description="Step 2: geometric efficiency per bin (detectors.py).")
    parser.add_argument("config", help="TOML config file")
    args = parser.parse_args(argv)
    config, full_path = common.read_config(args.config)

    for det in config.get("dets", [""]):
        c = common.config_for_detector(config, det)
        reaction, settings = c["reaction"], c.get("efficiency", {})
        mode = c.get("mode", "excitation")
        k = common.ecm_per_beam_energy(reaction)
        track = common.load_beam_track(c["beam"], full_path)
        bins, window = common.make_bins(c, track, k)
        states_table = pd.read_csv(full_path(reaction["states_csv"]), comment="#") if "states_csv" in reaction else None
        rng = np.random.default_rng(settings.get("seed", 20261016))
        n = int(settings.get("n_per_bin", 400_000))

        rows = []
        for _, b in bins.iterrows():
            if b["dz_mm"] <= 0:                                  # past the beam's stopping point
                rows.append({name: np.nan for name in EFFICIENCIES})
                continue
            if mode == "excitation":
                ecm = b["ecm_nominal_mev"]
                cos_range = (-1.0, 1.0)
            else:
                ecm = k * track.energy_at(0.5 * (b["z_up_mm"] + b["z_down_mm"]), nominal=True)
                cos_range = (np.cos(np.radians(b["theta_hi_deg"])), np.cos(np.radians(b["theta_lo_deg"])))
            rows.append(efficiency_of_bin(det, reaction, track, (b["z_up_mm"], b["z_down_mm"]), cos_range,
                                          ecm, states_table, n, rng))
        table = pd.concat([bins.reset_index(drop=True), pd.DataFrame(rows)], axis=1)

        states_text = reaction.get("states_csv", reaction.get("states", "ground state only"))
        header = [f"{c.get('title', '')} {det.upper()}: geometric efficiency per bin ({mode} mode), detectors.py",
                  f"beam track: {track.description}; Q = {q_value_mev(reaction):.4f} MeV; detected: "
                  f"{reaction.get('detect', 'ejectile')}; states: {states_text}",
                  f"isotropic CM, straight tracks, {n} reactions per bin; azimuthal coverage {detectors.PHI_FRACTION[det]:.4f}",
                  "eps_geometric = silicon + counter + apertures; eps_si_thr adds the Si threshold at the vertex energy",
                  "written by geometric_efficiency.py"]
        out = common.output_file(config, full_path, "eps", det)
        common.write_csv(table, header, out)
        x_col = "ecm_nominal_mev" if mode == "excitation" else "theta_c_deg"
        print(f"--- {det}: -> {out}")
        print(table[[x_col, *EFFICIENCIES]].to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
