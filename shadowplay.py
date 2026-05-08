import math
import csv

# ==========================================
# CONFIGURATION
# ==========================================
VTK_FILE = "anasen_geometry.vtk"
CSV_FILE = "anasen_labels.csv"

Z_SOURCE = 53.0      # Source position
R_A = 37.0           # Anode radius (mm)
R_C = 43.0           # Cathode radius (mm)
R_SX3 = 88.0         # SX3 Barrel radius (mm)
Z_LEN = 348.6        # Total Z length (174.3 * 2 from ClassPW.h)
Z_HALF = Z_LEN / 2.0

points = []
lines = []
polygons = []
cell_colors = []
labels = [] 

def add_point(x_code, y_code, z_code):
    # TRANSFORM: Map the C++ Frame (+Y Down, +Z Into Page)
    # to the Visual ParaView Frame (+Y Up, +Z Out of Page)
    points.append((x_code, -y_code, -z_code))
    return len(points) - 1

def project_shadow(r_wire, phi_minusZ_deg, phi_plusZ_deg, color_code):
    shadow_pts = []
    steps = 50 
    for i in range(steps + 1):
        t = i / float(steps)
        # 1. Coordinate on the wire in CODE frame
        z_w = -Z_HALF + (t * Z_LEN)
        phi_w_deg = phi_minusZ_deg * (1 - t) + phi_plusZ_deg * t
        phi_w = math.radians(phi_w_deg)
        
        # 2. Ray from Source through the wire
        dx = r_wire * math.cos(phi_w)
        dy = r_wire * math.sin(phi_w)
        dz = z_w - Z_SOURCE
        
        # 3. Scale ray until it hits SX3 (R = 88)
        r_current = math.sqrt(dx**2 + dy**2)
        alpha = R_SX3 / r_current
        
        x_proj = alpha * dx
        y_proj = alpha * dy
        z_proj = Z_SOURCE + alpha * dz
        
        shadow_pts.append(add_point(x_proj, y_proj, z_proj))
        
    lines.append(shadow_pts)
    cell_colors.append(color_code)

# ==========================================
# 1. GENERATE PC WIRES & SHADOWS (LINES)
# ==========================================
k = 360.0 / 24.0

for i in range(24):
    # --- ANODES ---
    # From ClassPW.h: offset_a1 = -135 deg, offset_a2 = -90 deg
    phi_a_minusZ = -k * i - 135.0
    phi_a_plusZ  = -k * i - 90.0
    
    x_a1, y_a1 = R_A * math.cos(math.radians(phi_a_minusZ)), R_A * math.sin(math.radians(phi_a_minusZ))
    x_a2, y_a2 = R_A * math.cos(math.radians(phi_a_plusZ)),  R_A * math.sin(math.radians(phi_a_plusZ))
    
    pa1 = add_point(x_a1, y_a1, -Z_HALF)
    pa2 = add_point(x_a2, y_a2, Z_HALF)
    lines.append([pa1, pa2])
    cell_colors.append(1)  # Color 1 = Anodes
    project_shadow(R_A, phi_a_minusZ, phi_a_plusZ, 4)
    
    # Place label at +Z_visual (which maps from -Z_code in the transform)
    labels.append((x_a1 * 1.1, -y_a1 * 1.1, Z_HALF + 15, f"A{i}"))
    
    # --- CATHODES ---
    # From ClassPW.h: offset_c1 = -97.5 deg, offset_c2 = -142.5 deg
    phi_c_minusZ = k * i - 97.5
    phi_c_plusZ  = k * i - 142.5
    
    x_c1, y_c1 = R_C * math.cos(math.radians(phi_c_minusZ)), R_C * math.sin(math.radians(phi_c_minusZ))
    x_c2, y_c2 = R_C * math.cos(math.radians(phi_c_plusZ)),  R_C * math.sin(math.radians(phi_c_plusZ))
    
    pc1 = add_point(x_c1, y_c1, -Z_HALF)
    pc2 = add_point(x_c2, y_c2, Z_HALF)
    lines.append([pc1, pc2])
    cell_colors.append(2) # Color 2 = Cathodes
    project_shadow(R_C, phi_c_minusZ, phi_c_plusZ, 5)
    
    labels.append((x_c1 * 1.1, -y_c1 * 1.1, Z_HALF + 30, f"C{i}"))

# ==========================================
# 2. GENERATE SX3 BARREL (POLYGONS)
# ==========================================
# Calculate exact strip positions based on MakeVertex.C formulas
for det_id in range(12):
    for stripF in range(4):
        stripF_rev = 3 - stripF
        num = (2 * stripF_rev - 3) * 40.30
        den = 8.0 * R_SX3 * math.cos(math.radians(15.0))
        beta_n = 15.0 + math.degrees(math.atan2(num, den))
        
        # phi_n perfectly computes the C++ angle
        phi_code_deg = ((-det_id + 0.5) * 30.0 + beta_n) + 45.0
        phi_code = math.radians(phi_code_deg)
        
        d_phi = 10.0 / R_SX3  # ~6.5 degrees width per strip
        phi_L = phi_code - (d_phi / 2.0)
        phi_R = phi_code + (d_phi / 2.0)
        
        p1 = add_point(R_SX3 * math.cos(phi_L), R_SX3 * math.sin(phi_L), -Z_HALF)
        p2 = add_point(R_SX3 * math.cos(phi_R), R_SX3 * math.sin(phi_R), -Z_HALF)
        p3 = add_point(R_SX3 * math.cos(phi_R), R_SX3 * math.sin(phi_R), Z_HALF)
        p4 = add_point(R_SX3 * math.cos(phi_L), R_SX3 * math.sin(phi_L), Z_HALF)
        
        polygons.append([p1, p2, p3, p4])
        cell_colors.append(3) # SX3 Color
        
        # Place Labels correctly in visual space
        x_lbl, y_lbl = R_SX3 * 1.15 * math.cos(phi_code), R_SX3 * 1.15 * math.sin(phi_code)
        labels.append((x_lbl, -y_lbl, 0, f"D{det_id}_S{stripF}"))

# ==========================================
# 3. WRITE VTK AND CSV
# ==========================================
with open(VTK_FILE, "w") as f:
    f.write("# vtk DataFile Version 2.0\nANASEN Geometry\nASCII\nDATASET POLYDATA\n")
    f.write(f"POINTS {len(points)} float\n")
    for pt in points: f.write(f"{pt[0]:.4f} {pt[1]:.4f} {pt[2]:.4f}\n")
        
    line_data_size = sum(len(l) + 1 for l in lines)
    f.write(f"\nLINES {len(lines)} {line_data_size}\n")
    for l in lines: f.write(f"{len(l)} " + " ".join(map(str, l)) + "\n")
        
    poly_data_size = sum(len(p) + 1 for p in polygons)
    f.write(f"\nPOLYGONS {len(polygons)} {poly_data_size}\n")
    for p in polygons: f.write(f"{len(p)} " + " ".join(map(str, p)) + "\n")
        
    f.write(f"\nCELL_DATA {len(lines) + len(polygons)}\n")
    f.write("SCALARS EntityType int 1\nLOOKUP_TABLE default\n")
    for c in cell_colors: f.write(f"{c}\n")

with open(CSV_FILE, "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["X", "Y", "Z", "Label"])
    writer.writerows(labels)

print(f"Generated {VTK_FILE} and {CSV_FILE}.")