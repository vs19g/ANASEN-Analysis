import ROOT
import os
import sys

# 1. FIX: Manually load the Garfield library if it's not in the ROOT namespace
# Update this path to your actual installation location
garfield_lib_path = "/home/vs19g/garfieldpp/install/lib/libGarfield.so" #penguin path
# garfield_lib_path = "/home/vsitaraman/garfieldpp/install/lib/libGarfield.so" #laptop path

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
    sys.exit(1)

# --- 2. GAS SETUP (96% He, 4% CO2) ---
gas = ROOT.Garfield.MediumMagboltz()
gas_file = "He96_CO2_4_260Torr.gas"

if not os.path.exists(gas_file):
    print("Generating new Magboltz gas table (260 Torr)...")
    # --- Optimized Magboltz Settings ---
    gas.SetComposition("he", 96.0, "co2", 4.0)
    gas.SetTemperature(293.15)
    gas.SetPressure(260.0)

    # 1. Limit the Field Grid: 
    gas.SetFieldGrid(10., 80000., 20, True) 

    # 2. Reduce the precision slightly for the first run:
    gas.GenerateGasTable(8) 
    gas.WriteGasFile(gas_file)
else:
    print(f"Loading existing gas table: {gas_file}")
    gas.LoadGasFile(gas_file)

# --- 3. FIELD MAP SETUP ---
fm = ROOT.Garfield.ComponentElmer2d()


fm.Initialise("wires2d/mesh.header", 
              "wires2d/mesh.elements", 
              "wires2d/mesh.nodes", 
              "dielectrics.dat", 
              "wires2d/elstatics.result", 
              "mm")

# Set the medium (Body 13 from your Gmsh script)
fm.SetMedium(0, gas)

# --- 4. SENSOR AND DRIFT SETUP ---
sensor = ROOT.Garfield.Sensor()
sensor.AddComponent(fm)
sensor.SetArea(-50.0, -50.0, -5.0, 50.0, 50.0, 5.0) #hardcoding the sesnsor area to define a psuedo 3d geometry

# Heavy Ion Drift (RKF) - Best for the general track
drift = ROOT.Garfield.DriftLineRKF()
drift.SetSensor(sensor)

# Electron Avalanche (Microscopic) - Best for high-field gain
aval = ROOT.Garfield.AvalancheMicroscopic()
aval.SetSensor(sensor)

# --- 5. EXECUTION ---
# Starting position (e.g., near the IC wires at r=23mm or closer to Anodes)
x0, y0, z0, t0 = 3.50, 0.0, 0.0, 0.0

print(f"Simulating heavy ion drift from r={x0}...")
drift.DriftIon(x0, y0, z0, t0)

# Create a file to store the heavy ion track
with open("heavy_ion_track.csv", "w") as f:
    f.write("x,y,z,t\n")
    
    # After running drift.DriftIon(x0, y0, z0, t0):
    n_points = drift.GetNumberOfDriftLinePoints()
    for i in range(n_points):
        xi, yi, zi, ti = ROOT.double(0), ROOT.double(0), ROOT.double(0), ROOT.double(0)
        drift.GetDriftLinePoint(i, xi, yi, zi, ti)
        f.write(f"{xi},{yi},{zi},{ti}\n")

print(f"Simulating electron avalanche from r={x0}...")
# AvalancheElectron(x, y, z, t, energy, dx, dy, dz)
aval.AvalancheElectron(x0, y0, z0, t0, 0.1, 0.0, 0.0, 0.0)

with open("avalanche_endpoints.csv", "w") as f:
    f.write("x,y,z,t\n")
    
    # After aval.AvalancheElectron(...)
    n_endpoints = aval.GetNumberOfEndpoints()
    for i in range(n_endpoints):
        # Get start and end points of each electron in the avalanche
        x1, y1, z1, t1, e1 = ROOT.double(0), ROOT.double(0), ROOT.double(0), ROOT.double(0), ROOT.double(0)
        x2, y2, z2, t2, e2, status = ROOT.double(0), ROOT.double(0), ROOT.double(0), ROOT.double(0), ROOT.double(0), ROOT.int(0)
        
        aval.GetEndpoint(i, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status)
        # We save the endpoint (x2, y2, z2) where the electron was collected or attached
        f.write(f"{x2},{y2},{z2},{t2}\n")

print("Simulation complete.")