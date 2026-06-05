#!/home/vsitaraman/ParaView-6.1.0-RC1-MPI-Linux-Python3.12-x86_64/bin/pvbatch

import sys
import numpy as np
from paraview.simple import *
import paraview.servermanager as sm
from paraview.vtk.numpy_interface import dataset_adapter as dsa

# ==========================================
# COMMAND LINE ARGUMENTS
# ==========================================
if len(sys.argv) < 2:
 print("Usage: pvbatch 1dpotplot.py <z_locus in mm>")
 sys.exit(1)

z_loc = float(sys.argv[1])

INPUT_FILE = f"wires2d/elfield_anasen_t0001.vtu"
OUTPUT_CSV = f"anode_to_cathode_1d_{z_loc}.csv"
OUTPUT_IMAGE = f"1D_Variation_Plot_{z_loc}.png"

print(f"--- Starting 1D Potential Extraction for Z = {z_loc} mm ---")

# ==========================================
# 1. EXACT GEOMETRIC CALCULATION
# ==========================================
print("Calculating theoretical wire coordinates...")

wireShift = 4.0
k = 2 * np.pi / 24.0

offset_a1 = -10.0 * k
offset_c1 = -5.5 * k
offset_a2 = offset_a1 + (wireShift * k)
offset_c2 = offset_c1 - (wireShift * k)

t = (z_loc + 174.3) / (2 * 174.3)

# Anode 0 theoretical position
xa_loc = (37.0 * np.cos(offset_a1)) + t * (37.0 * np.cos(offset_a2) - 37.0 * np.cos(offset_a1))
ya_loc = (37.0 * np.sin(offset_a1)) + t * (37.0 * np.sin(offset_a2) - 37.0 * np.sin(offset_a1))

# All Cathodes theoretical positions
xc1_all = np.array([42.0 * np.cos(k * i + offset_c1) for i in range(24)])
yc1_all = np.array([42.0 * np.sin(k * i + offset_c1) for i in range(24)])
xc2_all = np.array([42.0 * np.cos(k * i + offset_c2) for i in range(24)])
yc2_all = np.array([42.0 * np.sin(k * i + offset_c2) for i in range(24)])

xc_loc_all = xc1_all + t * (xc2_all - xc1_all)
yc_loc_all = yc1_all + t * (yc2_all - yc1_all)

# Nearest Cathode theoretical position
distances = np.sqrt((xc_loc_all - xa_loc)**2 + (yc_loc_all - ya_loc)**2)
nearest_c_idx = np.argmin(distances)
xc_loc = xc_loc_all[nearest_c_idx]
yc_loc = yc_loc_all[nearest_c_idx]

# Theoretical targets in meters
math_anode = np.array([xa_loc / 1000.0, ya_loc / 1000.0, 0.0])
math_cathode = np.array([xc_loc / 1000.0, yc_loc / 1000.0, 0.0])
# ==========================================
# 2. MESH DATA SNAPPING (Electrode-Constrained)
# ==========================================
print(f"Loading VTU file: {INPUT_FILE} to snap to electrode nodes...")

reader = XMLUnstructuredGridReader(FileName=[INPUT_FILE])
UpdatePipeline()

# Fetch data to Python
raw_data = sm.Fetch(reader)
wrapped_data = dsa.WrapDataObject(raw_data)

pts = wrapped_data.Points
pot = wrapped_data.PointData['potential']

print("\nPotential range:")
print(" min =", np.min(pot))
print(" max =", np.max(pot))

# -------------------------------------------------
# Find actual electrode nodes
# -------------------------------------------------

anode_candidates = np.where(pot > 640.0)[0]

if len(anode_candidates) == 0:
 raise RuntimeError("No anode nodes found near 650 V")

cathode_candidates = np.where(np.abs(pot) < 1.0)[0]

if len(cathode_candidates) == 0:
 raise RuntimeError("No cathode nodes found near 0 V")

print(f"Found {len(anode_candidates)} anode nodes")
print(f"Found {len(cathode_candidates)} cathode nodes")

# -------------------------------------------------
# Snap theoretical anode location to nearest
# actual anode node
# -------------------------------------------------

anode_dist = np.linalg.norm(
 pts[anode_candidates] - math_anode,
 axis=1
)

anode_idx = anode_candidates[np.argmin(anode_dist)]

true_anode_pos = pts[anode_idx]
true_anode_pot = pot[anode_idx]

# -------------------------------------------------
# Snap theoretical cathode location to nearest
# actual cathode node
# -------------------------------------------------

cathode_dist = np.linalg.norm(
 pts[cathode_candidates] - math_cathode,
 axis=1
)

cathode_idx = cathode_candidates[np.argmin(cathode_dist)]

true_cathode_pos = pts[cathode_idx]
true_cathode_pot = pot[cathode_idx]

# -------------------------------------------------
# Diagnostics
# -------------------------------------------------

anode_snap_error = np.linalg.norm(
 true_anode_pos - math_anode
)

cathode_snap_error = np.linalg.norm(
 true_cathode_pos - math_cathode
)

print("\nSelected electrode nodes:")
print(f"Anode:")
print(f" position = {true_anode_pos}")
print(f" potential = {true_anode_pot:.3f} V")
print(f" distance = {anode_snap_error*1e6:.1f} um")

print(f"\nCathode:")
print(f" position = {true_cathode_pos}")
print(f" potential = {true_cathode_pot:.3f} V")
print(f" distance = {cathode_snap_error*1e6:.1f} um")
# ==========================================
# 3. PARAVIEW EXTRACTION & PLOTTING
# ==========================================
print("Extracting data along the snapped line...")
plot_over_line = PlotOverLine(Input=reader)

# Apply the SNAPPED coordinates (converted to standard Python lists)
plot_over_line.Point1 = true_anode_pos.tolist()
plot_over_line.Point2 = true_cathode_pos.tolist()
plot_over_line.Resolution = 1000 

SaveData(OUTPUT_CSV, proxy=plot_over_line)

print("Rendering plot...")
line_chart_view = CreateView('XYChartView')
line_chart_view.ViewSize = [1000, 600]
line_chart_view.ChartTitle = f'Potential & E-Field: Anode to Cathode (Z = {z_loc} mm)'
line_chart_view.LeftAxisTitle = 'Potential (V)'
line_chart_view.BottomAxisTitle = 'Distance from Anode (m)'

line_display = Show(plot_over_line, line_chart_view)

# Tell ParaView to draw BOTH variables
# line_display.SeriesVisibility = ['potential', 'electric field_Magnitude']
line_display.SeriesVisibility = ['potential']

# Style Potential (Red)
line_display.SeriesColor = ['potential', '1.0', '0.0', '0.0'] 
line_display.SeriesLineThickness = ['potential', '3']

# # Style Electric Field (Blue, Right Axis)
# line_display.SeriesColor = ['electric field_Magnitude', '0.0', '0.5', '1.0'] 
# line_display.SeriesLineThickness = ['electric field_Magnitude', '3']
# line_display.SeriesPlotCorner = ['electric field_Magnitude', '1']
# line_chart_view.RightAxisTitle = 'Electric Field (V/m)'

SaveScreenshot(OUTPUT_IMAGE, view=line_chart_view)
print(f"Saved plot screenshot to: {OUTPUT_IMAGE}")
print("--- Extraction Complete ---")