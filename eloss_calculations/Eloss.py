import pycatima as catima
import numpy as np

# --- 1. Constants ---
P_TORR = 250
TEMP_K = 293.15 
R = 8.3144
MEV2U = 1.0 / 931.494
P_CO2 = 4

# Gas Density Calculations
p_pa = P_TORR * 133.322
molar_density = p_pa / (R * TEMP_K)
m_he, m_c, m_o= 4.0026, 12.0000, 15.9949
m_mix_avg = ((1 - P_CO2 / 100) * m_he) + (P_CO2 / 100 * (m_c + 2*m_o))
rho_g_cm3 = (molar_density * m_mix_avg) / 1e6
print(f"Gas density at {P_TORR} Torr: {rho_g_cm3:.6e} g/cm^3")

# --- 2. Material & Step Setup ---
material_def = [(m_he, 2, (1 - P_CO2 / 100)), (m_c, 6, P_CO2 / 100), (m_o, 8, 2*P_CO2 / 100)]
gas_mix = catima.Material(material_def)
gas_mix.density(rho_g_cm3)

# Thickness step settings
step_mg_cm2 = 0.001           # 1 ug/cm2 steps as per your example -- kept fine for
                               # numerical accuracy of the dedx integration itself.
step_g_cm2 = step_mg_cm2 / 1000.0
max_steps = 1000000000             # Adjust based on how far you want to track

coarse_step_cm = 0.2   # row spacing over most of the track
fine_step_cm = 0.03     # row spacing near the Bragg peak
fine_zone_frac = 0.085   # fraction of the *total* range treated as "near the peak"

def generate_lookup(z, mass_u, e_start_mev, label):
    filename = f"{label}_lookup_{e_start_mev}MeV_{P_TORR}torr_{P_CO2}pc.dat"
    header = f"Energy(MeV) \tmg/cm2 \tcm\nStarting Energy: {e_start_mev} MeV"

    # Pass 1: integrate at full precision just to find the total range (needed
    # to know where the "last fine_zone_frac" of the track begins).
    projectile = catima.Projectile(mass_u, z)
    e_u = e_start_mev / mass_u
    total_thickness_g_cm2 = 0.0
    while e_u >= 0.0001:
        projectile.T(e_u)
        loss_mev = catima.dedx(projectile, gas_mix) * step_g_cm2
        e_u = (e_u * mass_u - loss_mev) / mass_u
        total_thickness_g_cm2 += step_g_cm2
    total_range_cm = total_thickness_g_cm2 / rho_g_cm3
    fine_zone_start_cm = total_range_cm * (1.0 - fine_zone_frac)

    projectile = catima.Projectile(mass_u, z)
    current_e_total = e_start_mev
    current_thickness_g_cm2 = 0.0
    next_checkpoint_cm = 0.0

    output = []
    last_dist_cm = None

    def append_row(e_total, thickness_g_cm2, dist_cm):
        if last_dist_cm is not None and dist_cm <= last_dist_cm:
            return False
        output.append([e_total, thickness_g_cm2 * 1000.0, dist_cm])
        return True

    for i in range(max_steps):
        dist_cm = current_thickness_g_cm2 / rho_g_cm3
        if dist_cm >= next_checkpoint_cm:
            if append_row(current_e_total, current_thickness_g_cm2, dist_cm):
                last_dist_cm = dist_cm
            step_cm = fine_step_cm if dist_cm >= fine_zone_start_cm else coarse_step_cm
            next_checkpoint_cm = dist_cm + step_cm

        e_u = current_e_total / mass_u
        if e_u < 0.0001:  # Stop at ATIMA limit
            append_row(current_e_total, current_thickness_g_cm2, dist_cm)
            break

        projectile.T(e_u)
        # dedx returns MeV / (g/cm2)
        loss_mev = catima.dedx(projectile, gas_mix) * step_g_cm2

        current_e_total = max(0.0, current_e_total - loss_mev)
        current_thickness_g_cm2 += step_g_cm2

    np.savetxt(filename, output, fmt='%.6f', delimiter='\t', header=header)
    print(f"Lookup table created: {filename} ({len(output)} rows, "
          f"range {total_range_cm:.2f} cm, fine zone below {fine_zone_start_cm:.2f} cm)")

# --- 3. Run ---
# Format: generate_lookup(Z, mass_u, E_start_MeV, label)
generate_lookup(1, 1.0078, 30, "proton")
generate_lookup(1, 2.01355, 30, "deutron")
generate_lookup(2, 4.0026, 50, "alpha")
generate_lookup(13,26.9815, 80, "aluminum")
generate_lookup(9,17.0021, 70, "fluorine")
generate_lookup(8,15.9949, 70, "oxygen")