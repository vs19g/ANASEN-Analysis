"""Shared by extract_yields.py, geometric_efficiency.py and cross_section.py: the config, the beam track, the bins.

The cross section in this active target:
    sigma(E_cm)        = Y / (N_b * n * dz * eps)                  mode = "excitation"
    dsigma/dOmega(th)  = Y / (N_b * n * dz * eps * dOmega)         mode = "angular"
The beam loses energy along z, so a vertex at z means the reaction happened at the beam energy E(z) there.
An E_cm bin is the stretch of beam path dz over which the beam's E_cm lies in the bin.
"""
import os
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================ config

def read_config(config_file):
    """Read the TOML config. Returns the config and a function that turns paths in it into absolute paths
    (relative paths are relative to the config file's folder)."""
    config_file = Path(config_file).resolve()
    with open(config_file, "rb") as f:
        config = tomllib.load(f)

    def full_path(p):
        p = os.path.expanduser(str(p))           # "~/..." is the home folder, as in a shell
        if p.startswith("/") or Path(p).is_absolute():
            return p
        return str((config_file.parent / p).resolve())

    return config, full_path


def config_for_detector(config, det):
    """Copy of the config with "{det}" (and any key of [det_vars.<det>]) replaced in every string.
    This lets one config serve several detectors, e.g. hist = "..._{det}"."""
    replacements = {"det": det}
    replacements.update(config.get("det_vars", {}).get(det, {}))

    def replace(value):
        if isinstance(value, str):
            for key, new in replacements.items():
                value = value.replace("{" + key + "}", str(new))
            return value
        if isinstance(value, dict):
            return {k: replace(v) for k, v in value.items()}
        if isinstance(value, list):
            return [replace(v) for v in value]
        return value

    return replace(config)


def output_file(config, full_path, kind, det):
    """Where each step writes for one detector: <out_dir>/<kind>_<det>.csv, kind = yields, eps or xs."""
    out_dir = Path(full_path(config.get("out_dir", "out")))
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / (f"{kind}_{det}.csv" if det else f"{kind}.csv")


def make_edges(spec):
    """Bin edges from a list [e0, e1, ...] or from {start = .., stop = .., step = ..} (stop included)."""
    if isinstance(spec, dict):
        n_bins = int(round((spec["stop"] - spec["start"]) / spec["step"]))
        return spec["start"] + spec["step"] * np.arange(n_bins + 1)
    return np.array(spec, dtype=float)


def write_csv(table, header_lines, out_file):
    with open(out_file, "w", newline="\n") as f:
        for line in header_lines:
            f.write(f"# {line}\n")
        table.to_csv(f, index=False, float_format="%.6g")


# ============================================================ beam track

class BeamTrack:
    """Beam energy (MeV, lab) along z (mm). The energy falls with z and is 0 past the stopping point.

    `energy` is the track used to define the bins (the one the data were reconstructed with). `nominal` is an
    optional second track (e.g. another stopping-power scale): it gives each bin's nominal E_cm and the beam
    energy in the efficiency toy. It defaults to the same track."""

    def __init__(self, z, energy, nominal=None, description=""):
        self.z = np.asarray(z, dtype=float)
        self.energy = np.asarray(energy, dtype=float)
        self.nominal = self.energy if nominal is None else np.asarray(nominal, dtype=float)
        self.description = description

    def energy_at(self, z, nominal=False):
        return np.interp(z, self.z, self.nominal if nominal else self.energy)

    def z_where_energy_is(self, e):
        """The z where the beam has energy e.
        e <= 0 gives the stopping point; e above the entrance energy gives the entrance z."""
        if e <= 0:
            stopped = np.flatnonzero(self.energy <= 0)
            return float(self.z[stopped[0]]) if len(stopped) else float(self.z[-1])
        # np.interp needs increasing x: energy decreases along z, so interpolate in -energy
        return float(np.interp(-e, -self.energy, self.z))


def load_beam_track(beam, full_path):
    """Either a CSV with columns z and E (e.g. beam_tracks.csv), or an energy-loss table from
    eloss_calculations/ plus the beam energy at the entrance."""
    if "track_csv" in beam:
        table = pd.read_csv(full_path(beam["track_csv"]), comment="#")
        z = table[beam.get("track_z_col", "z_mm")]
        energy = table[beam["track_e_col"]]
        nominal = table[beam["nominal_e_col"]] if "nominal_e_col" in beam else None
        text = f"{beam['track_csv']}, column {beam['track_e_col']}"
        if nominal is not None:
            text += f" (nominal E_cm from {beam['nominal_e_col']})"
        return BeamTrack(z, energy, nominal, text)

    if "eloss_table" in beam:
        # Table columns: energy (MeV), thickness (mg/cm2), distance along the path (cm)
        table = np.loadtxt(full_path(beam["eloss_table"]), comments="#")
        table_energy, table_cm = table[:, 0], table[:, 2]
        e_entrance, z_entrance = float(beam["e_entrance_mev"]), float(beam["z_entrance_mm"])
        if e_entrance > table_energy[0]:
            raise ValueError(f"e_entrance_mev = {e_entrance} is above the table's start energy {table_energy[0]}")
        # Where in the table the beam starts, then walk along it in 0.1 mm steps
        start_cm = np.interp(-e_entrance, -table_energy, table_cm)
        z = np.arange(z_entrance, beam.get("z_max_mm", 400.0) + 0.05, 0.1)
        path_cm = start_cm + (z - z_entrance) / 10.0
        energy = np.interp(path_cm, table_cm, table_energy, right=0.0)
        return BeamTrack(z, energy, None, f"{beam['eloss_table']}, {e_entrance} MeV at z = {z_entrance} mm")

    raise KeyError("[beam] needs either track_csv or eloss_table")


def ecm_per_beam_energy(reaction):
    """E_cm = E_beam x m_target / (m_beam + m_target)   (non-relativistic; or give ecm_per_ebeam directly)."""
    if "ecm_per_ebeam" in reaction:
        return float(reaction["ecm_per_ebeam"])
    return reaction["target_mass_u"] / (reaction["beam_mass_u"] + reaction["target_mass_u"])


# ============================================================ bins

def energy_bins(bins, track, k):
    """E_cm bins and, for each, the stretch of beam path [z_up, z_down] where the beam's E_cm is inside it.
    The higher energy is upstream, so z_up belongs to the bin's top edge."""
    ecm_edges = make_edges(bins["ecm_edges"])
    t = pd.DataFrame({"ecm_lo_mev": ecm_edges[:-1], "ecm_hi_mev": ecm_edges[1:]})
    t["ecm_c_mev"] = 0.5 * (t["ecm_lo_mev"] + t["ecm_hi_mev"])
    t["z_up_mm"] = [track.z_where_energy_is(e / k) for e in t["ecm_hi_mev"]]
    t["z_c_mm"] = [track.z_where_energy_is(e / k) for e in t["ecm_c_mev"]]
    t["z_down_mm"] = [track.z_where_energy_is(e / k) for e in t["ecm_lo_mev"]]
    t["dz_mm"] = t["z_down_mm"] - t["z_up_mm"]
    t["ecm_nominal_mev"] = k * track.energy_at(t["z_c_mm"], nominal=True)
    # True if the beam enters the chamber already inside this bin (dz then starts at the entrance)
    t["partial_top"] = t["ecm_hi_mev"] / k > track.energy[0]
    t["d_omega_sr"] = 1.0                    # not used in excitation mode
    return t


def angle_bins(bins, track, k):
    """theta_cm bins inside one window, given in E_cm (window_ecm) or in vertex z (window_z_mm).
    Returns the bins and the window as (ecm_lo, ecm_hi, z_up, z_down)."""
    if "window_ecm" in bins:
        ecm_lo, ecm_hi = bins["window_ecm"]
        z_up, z_down = track.z_where_energy_is(ecm_hi / k), track.z_where_energy_is(ecm_lo / k)
    else:
        z_up, z_down = bins["window_z_mm"]
        ecm_lo, ecm_hi = k * track.energy_at(z_down), k * track.energy_at(z_up)
    theta_edges = make_edges(bins["theta_edges_deg"])
    t = pd.DataFrame({"theta_lo_deg": theta_edges[:-1], "theta_hi_deg": theta_edges[1:]})
    t["theta_c_deg"] = 0.5 * (t["theta_lo_deg"] + t["theta_hi_deg"])
    t["d_omega_sr"] = 2 * np.pi * (np.cos(np.radians(t["theta_lo_deg"])) - np.cos(np.radians(t["theta_hi_deg"])))
    t["z_up_mm"], t["z_down_mm"] = z_up, z_down
    t["dz_mm"] = z_down - z_up
    return t, (ecm_lo, ecm_hi, z_up, z_down)


def make_bins(config, track, k):
    """The bins of the config's mode. Returns (table, window or None)."""
    if config.get("mode", "excitation") == "excitation":
        return energy_bins(config["bins"], track, k), None
    return angle_bins(config["bins"], track, k)
