import pycatima as catima
import numpy as np

# --- 1. Constants ---
P_TORR = 250
TEMP_K = 293.15 
R = 8.3144
MEV2U = 1.0 / 931.494

# Gas Density Calculations
p_pa = P_TORR * 133.322
molar_density = p_pa / (R * TEMP_K)
m_he, m_c, m_o= 4.0026, 12.0000, 15.9949
m_mix_avg = (0.96 * m_he) + (0.04 * (m_c + 2*m_o))
rho_g_cm3 = (molar_density * m_mix_avg) / 1e6
print(f"Gas density at {P_TORR} Torr: {rho_g_cm3:.6e} g/cm^3")

# --- 2. Material & Step Setup ---
material_def = [(m_he, 2, 0.96), (m_c, 6, 0.04), (m_o, 8, 0.08)]
gas_mix = catima.Material(material_def)
gas_mix.density(rho_g_cm3)

# Thickness step settings
step_mg_cm2 = 0.001           # 1 ug/cm2 steps as per your example
step_g_cm2 = step_mg_cm2 / 1000.0 
max_steps = 1000000000             # Adjust based on how far you want to track

def generate_lookup(z, mass_u, e_start_mev, label):
    filename = f"{label}_lookup_{e_start_mev}MeV.dat"
    projectile = catima.Projectile(mass_u, z)
    
    current_e_total = e_start_mev
    current_thickness_g_cm2 = 0.0
    
    output = []
    header = f"Energy(MeV) \tmg/cm2 \tcm\nStarting Energy: {e_start_mev} MeV"

    for i in range(max_steps):
        # 1. Record current state
        dist_cm = current_thickness_g_cm2 / rho_g_cm3
        output.append([current_e_total, current_thickness_g_cm2 * 1000.0, dist_cm])
        
        # 2. Calculate energy loss for the NEXT step
        e_u = current_e_total / mass_u
        if e_u < 0.0001: # Stop at ATIMA limit
            break
            
        projectile.T(e_u)
        # dedx returns MeV / (g/cm2)
        loss_mev = catima.dedx(projectile, gas_mix) * step_g_cm2
        
        # 3. Update values
        current_e_total -= loss_mev
        current_thickness_g_cm2 += step_g_cm2

    np.savetxt(filename, output, fmt='%.6f', delimiter='\t', header=header)
    print(f"Lookup table created: {filename}")

# --- 3. Run ---
# Format: generate_lookup(Z, mass_u, E_start_MeV, label)
generate_lookup(1, 1.0078, 20, "proton")
generate_lookup(2, 4.0026, 20, "alpha")
generate_lookup(13,26.9815, 80, "aluminum")
generate_lookup(9,17.0021, 70, "fluorine")
generate_lookup(8,15.9949, 70, "oxygen")