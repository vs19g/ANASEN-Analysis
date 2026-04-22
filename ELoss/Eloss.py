#VSitaraman 2024-06-20
# This script calculates the stopping power and range for protons and alpha particles in a gas mixture of 96% Helium and 4% CO2 at 250 Torr using the Catima library.
import numpy as np

# --- 1. Experimental & Physical Constants ---
P_TORR = 250
TEMP_K = 293.15 
R = 8.3144
MEV2U = 1.0 / 931.494  # Conversion factor for mass

# Pressure to Pascal and Molar Density
p_pa = P_TORR * 133.322
molar_density = p_pa / (R * TEMP_K)

# Nuclear Data (Masses in u)
m_he = 4.002603
m_c  = 12.00000
m_o  = 15.99491
m_co2 = m_c + (2 * m_o)

# Effective mass of mixture (96% He, 4% CO2)
m_mix_avg = (0.96 * m_he) + (0.04 * m_co2)
# Density in g/cm3
rho = (molar_density * m_mix_avg) / 1e6

# --- 2. Create the Material ---
# We treat 0.96 and 0.04 as the stoichiometry/abundance in the mixture
material_def = [
    (m_he, 2, 0.96), # Helium
    (m_c,  6, 0.04), # Carbon (from CO2)
    (m_o,  8, 0.08)  # Oxygen (2 * 0.04 from CO2)
]
gas_mix = catima.Material(material_def)
gas_mix.density(rho)

# --- 3. Generation Function ---
def save_loss_file(z, mass_u, label):
    filename = f"{label}_loss_250torr.dat"
    # Energy range from 0.05 MeV/u to 50 MeV/u
    energies_u = np.logspace(-1.3, 1.7, 250) 
    
    # Initialize Projectile: (Mass_u, Z)
    projectile = catima.Projectile(mass_u, z)
    
    header = (f"{label.capitalize()} Energy Loss | 96% He + 4% CO2 @ {P_TORR} Torr\n"
              f"Density: {rho:.6e} g/cm3\n"
              f"E_total(MeV) \tdE/dx(MeV/cm) \tRange(cm)")
    
    output_data = []
    for e_u in energies_u:
        # Set kinetic energy in MeV/u
        projectile.T(e_u)
        
        # Calculate stopping power (MeV / (g/cm2))
        # Multiply by density to get MeV/cm
        stop_power = catima.dedx(projectile, gas_mix)
        dedx_mev_cm = stop_power * rho
        
        # Calculate range (g/cm2) and convert to cm
        range_g_cm2 = catima.range(projectile, gas_mix)
        range_cm = range_g_cm2 / rho
        
        e_total = e_u * mass_u
        output_data.append([e_total, dedx_mev_cm, range_cm])
    
    np.savetxt(filename, output_data, fmt='%.6e', delimiter='\t', header=header)
    print(f"File created: {filename}")

# --- 4. Execute ---
# Proton: Z=1, Mass ~1.007 u
save_loss_file(1, 1.0078, "proton")

# Alpha: Z=2, Mass ~4.0026 u
save_loss_file(2, 4.0026, "alpha")

print(f"\nCompleted. Calculated gas density: {rho:.6e} g/cm3")