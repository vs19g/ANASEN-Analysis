import ROOT
import os
import sys

# 1. FIX: Manually load the Garfield library if it's not in the ROOT namespace
# Update this path to your actual installation location
garfield_lib_path = "/home/vsitaraman/garfieldpp/install/lib/libGarfield.so"

if os.path.exists(garfield_lib_path):
    ROOT.gSystem.Load(garfield_lib_path)
else:
    print(f"CRITICAL ERROR: {garfield_lib_path} not found.")
    sys.exit(1)

# Verify access
try:
    test_gas = ROOT.Garfield.MediumMagboltz()
except AttributeError:
    print("ERROR: Garfield shared library loaded, but 'Garfield' namespace not found.")
    print("Ensure Garfield was compiled with -DWITH_PYTHON=ON")
    sys.exit(1)

# --- 2. GAS SETUP (96% He, 4% CO2) ---
gas = ROOT.Garfield.MediumMagboltz()
gas_file = "He96_CO2_4_260Torr.gas"

if not os.path.exists(gas_file):
    print("Generating new Magboltz gas table (260 Torr)...")
    gas.SetComposition("he", 96.0, "co2", 4.0)
    gas.SetTemperature(293.15)
    gas.SetPressure(260.0) # Pressure in Torr
    
    # Grid must cover the field at the 18um wire surface (~100kV/cm)
    gas.SetFieldGrid(10., 150000., 20, True)
    
    # 10 iterations is usually plenty for a simple mix
    gas.GenerateGasTable(10)
    gas.WriteGasFile(gas_file)
else:
    print(f"Loading existing gas table: {gas_file}")
    gas.LoadGasFile(gas_file)

# --- 3. FIELD MAP SETUP ---
fm = ROOT.Garfield.ComponentElmer()

# Update these filenames to match your Elmer SIF output exactly
# Assuming ElmerGrid was run on 'wires2d' directory
fm.Initialise("wires2d/mesh.nodes", 
              "wires2d/mesh.elements", 
              "wires2d/mesh.boundary", 
              "wires2d/elfield_anasen.result", "mm")

# Set the medium (Body 13 from your Gmsh script)
fm.SetMedium(0, gas)

# --- 4. SENSOR AND DRIFT SETUP ---
sensor = ROOT.Garfield.Sensor()
sensor.AddComponent(fm)

# Heavy Ion Drift (RKF) - Best for the general track
drift = ROOT.Garfield.DriftLineRKF()
drift.SetSensor(sensor)

# Electron Avalanche (Microscopic) - Best for high-field gain
aval = ROOT.Garfield.AvalancheMicroscopic()
aval.SetSensor(sensor)

# --- 5. EXECUTION ---
# Starting position (e.g., near the IC wires at r=23mm or closer to Anodes)
x0, y0, z0, t0 = 35.0, 0.0, 0.0, 0.0

print(f"Simulating heavy ion drift from r={x0}...")
drift.DriftIon(x0, y0, z0, t0)

print(f"Simulating electron avalanche from r={x0}...")
# AvalancheElectron(x, y, z, t, energy, dx, dy, dz)
aval.AvalancheElectron(x0, y0, z0, t0, 0.1, 0.0, 0.0, 0.0)

print("Simulation complete.")