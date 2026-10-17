"""The ANASEN detectors as they are now (Run A, 27Al campaign), for the geometric efficiency.

geometric_efficiency.py reads everything about the detectors from here. Change a number here and every efficiency
computed afterwards follows it.

Sources:
  - Si and counter: anasen_detectability.py of ~/TALYS/al027_a_scripts_stoichgas_dedx0962_2026-10-08
    ("detector geometry as used by TrackRecon.C"; counter radii from ClassPW.h)
  - beam tube and PCB hole: the author, 2026-10-15 (../scratch/reaction_efficiency_tube_2026-10-15)
Geometry only: straight tracks from the vertex, no energy loss.
"""
import numpy as np

# ---------------------------------------------------------------- silicon
# QQQ3: the downstream annular detector
QQQ_Z_MM = 105.0                       # TrackRecon.C qqq_z
QQQ_R_MM = (50.0, 100.0)               # inner and outer active radius (16 rings: rho = 50 + 50/16 (ring + 0.5))
QQQ_PHI_FRACTION = 4 * 87.0 / 360.0    # 4 quadrants x 87 deg of wedges

# SX3: the barrel; only the four detectors TrackRecon.C builds events for, downstream half only
SX3_R_MM = 88.0
SX3_Z_MM = (-3.0, 72.0)                # z_local + 75/2 - 3
SX3_IDS = (1, 3, 7, 9)
SX3_HALF_ANGLE = np.arctan(4 * 40.30 / (8.0 * 88.0 * np.cos(np.radians(15.0))))   # outer strip edge, each side
SX3_PHI_FRACTION = len(SX3_IDS) * 2 * SX3_HALF_ANGLE / (2 * np.pi)

PHI_FRACTION = {"qqq": QQQ_PHI_FRACTION, "sx3": SX3_PHI_FRACTION}

# Si energy thresholds (MeV), TrackRecon.C reaction_ax_core si_ecut. Used only for the eps_si_thr column.
SI_THRESHOLD_MEV = {"qqq": 0.45, "sx3": 0.6}

# ---------------------------------------------------------------- proportional counter
COUNTER_HALF_LENGTH_MM = 174.3
COUNTER_RADII_MM = (32.0, 37.0, 42.0)  # guard, anode, cathode (ClassPW.h)

# ---------------------------------------------------------------- upstream apertures
# (name, radius mm, z mm of the face that binds). A track from an on-axis vertex only moves away from the axis,
# so an opening binds at its downstream end.
APERTURES = (("beam tube", 16.0, -189.0),
             ("PCB hole", 18.0, -174.3))


# ---------------------------------------------------------------- hit rules
# z_vertex in mm, theta = lab polar angle of the particle (radians). All return one True/False per track.

def hits_silicon(det, z_vertex, theta):
    """Does the straight track from z_vertex at angle theta reach the active silicon of `det` ("qqq" or "sx3")?
    Polar angle only: multiply by PHI_FRACTION[det] for the azimuthal coverage.
    A track that reaches the SX3 barrel's z range is counted as SX3, not QQQ."""
    s, c = np.sin(theta), np.cos(theta)
    with np.errstate(divide="ignore", invalid="ignore"):
        z_at_sx3 = z_vertex + SX3_R_MM * c / s                         # z where the track reaches the barrel radius
        on_sx3 = (z_at_sx3 >= SX3_Z_MM[0]) & (z_at_sx3 <= SX3_Z_MM[1])
        if det == "sx3":
            return on_sx3
        r_at_qqq = np.where(c > 0, (QQQ_Z_MM - z_vertex) * s / c, -1.0)  # radius where it reaches the QQQ plane
    if det == "qqq":
        return (QQQ_Z_MM > z_vertex) & (r_at_qqq >= QQQ_R_MM[0]) & (r_at_qqq <= QQQ_R_MM[1]) & ~on_sx3
    raise ValueError(f"unknown detector {det!r}: detectors.py defines 'qqq' and 'sx3'")


def crosses_counter(z_vertex, theta):
    """Does the track cross the guard, anode and cathode radii inside the counter's length?"""
    with np.errstate(divide="ignore", invalid="ignore"):
        cot = np.cos(theta) / np.sin(theta)
        ok = np.ones(np.shape(z_vertex), bool)
        for r in COUNTER_RADII_MM:
            ok &= np.abs(z_vertex + r * cot) <= COUNTER_HALF_LENGTH_MM
    return ok


def clears_apertures(z_vertex, theta):
    """Does the track pass the beam tube and the PCB hole? Only vertices upstream of an aperture can be blocked
    by it: blocked when the track is farther from the axis than the opening's radius at the opening's z."""
    s, c = np.sin(theta), np.cos(theta)
    tan_theta = s / np.where(c > 0, c, 1e-9)        # backward tracks from upstream never reach the detectors
    ok = np.ones(np.shape(z_vertex), bool)
    for _name, radius, z_face in APERTURES:
        ok &= ~((z_vertex < z_face) & ((z_face - z_vertex) * tan_theta > radius))
    return ok
