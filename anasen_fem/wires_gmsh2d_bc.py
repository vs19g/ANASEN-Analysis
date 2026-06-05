# This script generates a 2D Gmsh geometry for the ANASEN detector's wire planes, including an optional hot needle and IC wires. The geometry is designed for adaptive meshing, with finer mesh near the wires and coarser mesh in empty space. The script takes a z-locus as a command-line argument to determine the interpolation between the two wire planes.
import numpy as np
import gmsh, sys

# --- Configuration Flags ---
include_ic_wires = True
include_needle = True

gmsh.initialize()
gmsh.model.add("adaptive_mesh")
gmsh.option.setNumber("General.NumThreads", 10)
# gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfElements", 200000)
# gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfNodes", 200000)
# gmsh.option.setNumber("Mesh.Adapt.MaxIter",5)
# gmsh.option.setNumber("Mesh.MeshSizeMin", 5e-3)
# gmsh.option.setNumber("Mesh.MeshSizeMax", 10.0)
gmsh.option.setNumber("Geometry.Tolerance", 4e-2)
# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

lc = 0.04
# z_loc = -174.3

if len(sys.argv) < 2:
    print("Usage: python3 wires_gmsh2d_bc.py <z_locus in mm>")
    quit()

z_loc = float(sys.argv[1])

wireShift = 4.0
k = 2 * np.pi / 24.0
kg = k/2.0

# Plane 1 Offsets (-zmax/2)
# Anodes: -6*k (base) - 4*k (correction) = -10*k
offset_a1 = -6 * k - 4 * k
# Cathodes: -6*k (base) + 0.5*k (half-placement)
offset_c1 = -6 * k + (k / 2.0)
# Guard wires: aligned with cathodes
offset_g1 = offset_c1

# Plane 2 Offsets (+zmax/2) with Twist
offset_a2 = offset_a1 + (wireShift * k)
offset_c2 = offset_c1 - (wireShift * k)
offset_g2 = offset_c2

# 1 needle, 24 ic1, 24 ic2, 48 guard wires, 24 anodes, 24 cathodes

# needle at plane 1 at -zmax/2 no rotation
xarr_needle = np.array([0])
yarr_needle = np.array([0])

# ic1 wires, plane 1 at -zmax/2 no rotation
xarr_i11 = np.array([23 * np.cos(k * i) for i in range(24)])
yarr_i11 = np.array([23 * np.sin(k * i) for i in range(24)])

# ic1 wires, plane 1 at -zmax/2 no rotation
xarr_i21 = np.array([23 * np.cos(k * i + k / 2.0) for i in range(24)])
yarr_i21 = np.array([23 * np.sin(k * i + k / 2.0) for i in range(24)])

# --- Coordinate Arrays (Plane 1) ---
# Anodes: -k*i (Left-handed twist direction)
xarra_1 = np.array([37 * np.cos(-k * i + offset_a1) for i in range(24)])
yarra_1 = np.array([37 * np.sin(-k * i + offset_a1) for i in range(24)])

# Cathodes: +k*i (Right-handed twist direction)
xarrc_1 = np.array([42 * np.cos(k * i + offset_c1) for i in range(24)])
yarrc_1 = np.array([42 * np.sin(k * i + offset_c1) for i in range(24)])

# Guard Wires (48 wires, use kg spacing)
xarrg_1 = np.array([32 * np.cos(kg * i + offset_g1) for i in range(48)])
yarrg_1 = np.array([32 * np.sin(kg * i + offset_g1) for i in range(48)])

# needle at plane 2 at zmax/2
xarr_needle_2 = np.array([0])
yarr_needle_2 = np.array([0])

# #ic1 wires, plane 2 at zmax/2
xarr_i12 = np.array([23 * np.cos(k * i) for i in range(24)])
yarr_i12 = np.array([23 * np.sin(k * i) for i in range(24)])

# #ic2 wires, plane 2 at zmax/2
xarr_i22 = np.array([23 * np.cos(k * i + k / 2.0) for i in range(24)])
yarr_i22 = np.array([23 * np.sin(k * i + k / 2.0) for i in range(24)])

# --- Coordinate Arrays (Plane 2) ---
xarra_2 = np.array([37 * np.cos(-k * i + offset_a2) for i in range(24)])
yarra_2 = np.array([37 * np.sin(-k * i + offset_a2) for i in range(24)])

xarrc_2 = np.array([42 * np.cos(k * i + offset_c2) for i in range(24)])
yarrc_2 = np.array([42 * np.sin(k * i + offset_c2) for i in range(24)])

xarrg_2 = np.array([32 * np.cos(kg * i + offset_g2) for i in range(48)])
yarrg_2 = np.array([32 * np.sin(kg * i + offset_g2) for i in range(48)])

direction_needle_x = xarr_needle_2 - xarr_needle
direction_needle_y = yarr_needle_2 - yarr_needle

direction_ic1_x = xarr_i12 - xarr_i11
direction_ic1_y = yarr_i12 - yarr_i11

direction_ic2_x = xarr_i22 - xarr_i21
direction_ic2_y = yarr_i22 - yarr_i21

direction_guard_x = xarrg_2 - xarrg_1
direction_guard_y = yarrg_2 - yarrg_1

direction_anodes_x = xarra_2 - xarra_1
direction_anodes_y = yarra_2 - yarra_1

direction_cathodes_x = xarrc_2 - xarrc_1
direction_cathodes_y = yarrc_2 - yarrc_1

t = (z_loc + 174.3) / (2 * 174.3)  # z=-174.3 is 0, z=+174.3 is 1

xloc_needle = xarr_needle + t * direction_needle_x
yloc_needle = yarr_needle + t * direction_needle_y
xloc_i1 = xarr_i11 + t * direction_ic1_x
yloc_i1 = yarr_i11 + t * direction_ic1_y
xloc_i2 = xarr_i21 + t * direction_ic2_x
yloc_i2 = yarr_i21 + t * direction_ic2_y
xloc_g = xarrg_1 + t * direction_guard_x
yloc_g = yarrg_1 + t * direction_guard_y
xloc_a = xarra_1 + t * direction_anodes_x
yloc_a = yarra_1 + t * direction_anodes_y
xloc_c = xarrc_1 + t * direction_cathodes_x
yloc_c = yarrc_1 + t * direction_cathodes_y

# wire_radius_a = 0.018 #mm
# wire_radius_c = 0.0762 #mm
# wire_radius_g = 0.0762 #mm
wire_radius = 0.254  # mm
needle = []
ic1_wires = []
ic2_wires = []
guard_wires = []
anode_wires = []
cathode_wires = []
iw1_tags = [(3, i) for i in range(24)]
iw2_tags = [(3, i + 24) for i in range(24)]
gw_tags = [(3, i + 48) for i in range(48)]
aw_tags = [(3, i) for i in range(24)]
cw_tags = [(3, i + 24) for i in range(24)]

# for i,[xa,ya,xc,yc] in enumerate(zip(xarra_1,yarra_1,xarrc_1,yarrc_1)):
# create Hot Needle (1 total)
for i, (xn, yn) in enumerate(zip(xloc_needle, yloc_needle)):
    if include_needle:
        ndisk = gmsh.model.occ.addDisk(xn, yn, 0, wire_radius, wire_radius)
        needle.append(ndisk)

# create Guard Wires (48 total)
for xg, yg in zip(xloc_g, yloc_g):
    gdisk = gmsh.model.occ.addDisk(xg, yg, 0, wire_radius, wire_radius)
    guard_wires.append(gdisk)

# create Cathode Wires (24 total)
for xc, yc in zip(xloc_c, yloc_c):
    cdisk = gmsh.model.occ.addDisk(xc, yc, 0, wire_radius, wire_radius)
    cathode_wires.append(cdisk)

# create IC Anode and Cathode Wires (24 total each)
for i, (xa, ya) in enumerate(zip(xloc_a, yloc_a)):
    adisk = gmsh.model.occ.addDisk(xa, ya, 0, wire_radius, wire_radius)
    anode_wires.append(adisk)

    # Place IC wires only if flag is True
    if include_ic_wires:
        i1disk = gmsh.model.occ.addDisk(
            xloc_i1[i], yloc_i1[i], 0, wire_radius, wire_radius
        )
        i2disk = gmsh.model.occ.addDisk(
            xloc_i2[i], yloc_i2[i], 0, wire_radius, wire_radius
        )
        ic1_wires.append(i1disk)
        ic2_wires.append(i2disk)

anasen_barrel = gmsh.model.occ.addDisk(0, 0, 0, 500, 500)
# gmsh.model.occ.synchronize()
# gmsh.model.mesh.embed(1,anode_wires+cathode_wires,2,anasen_barrel)

gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
gmsh.model.occ.synchronize()


# --- Surface Extraction ---
def get_surfs(disks):
    surfs = []
    for d in disks:
        surfs += [
            s[1] for s in gmsh.model.getBoundary([(2, d)], oriented=False) if s[0] == 1
        ]
    return surfs


needle_surfs = get_surfs(needle) if include_needle else []
gwire_surfs = get_surfs(guard_wires)
awire_surfs = get_surfs(anode_wires)
cwire_surfs = get_surfs(cathode_wires)
i1wire_surfs = get_surfs(ic1_wires) if include_ic_wires else []
i2wire_surfs = get_surfs(ic2_wires) if include_ic_wires else []


all_active_wire_surfs = (
    needle_surfs + gwire_surfs + awire_surfs + cwire_surfs + i1wire_surfs + i2wire_surfs
)
gmsh.model.mesh.embed(1, all_active_wire_surfs, 2, anasen_barrel)

f1 = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f1, "CurvesList", all_active_wire_surfs)

f2 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f2, "InField", f1)
gmsh.model.mesh.field.setNumber(f2, "SizeMin", 0.05)  # Fine mesh near wires
gmsh.model.mesh.field.setNumber(f2, "SizeMax", 5.0)  # Large mesh in empty space
gmsh.model.mesh.field.setNumber(f2, "DistMin", 0.5)  # Apply SizeMin within 1mm
gmsh.model.mesh.field.setNumber(f2, "DistMax", 15.0)  # Transition to SizeMax by 20mm

gmsh.model.mesh.field.setAsBackgroundMesh(f2)


# --- Physical Groups ---
# Needle
if include_needle:
    gmsh.model.addPhysicalGroup(1, needle_surfs, tag=1, name="hot_needle")

# IC Wires
if include_ic_wires:
    gmsh.model.addPhysicalGroup(1, i1wire_surfs, tag=2, name="ic_wire_1")
    gmsh.model.addPhysicalGroup(1, i2wire_surfs, tag=3, name="ic_wire_2")

# Proportional Counter Wires
gmsh.model.addPhysicalGroup(1, gwire_surfs, tag=10, name="guard_wires")
gmsh.model.addPhysicalGroup(1, awire_surfs, tag=20, name="anode_wires")
gmsh.model.addPhysicalGroup(1, cwire_surfs, tag=30, name="cathode_wires")

# Gas Volume (2D)
gmsh.model.addPhysicalGroup(2, [anasen_barrel], tag=13, name="gas")

gmsh.option.setNumber("Mesh.Algorithm", 6)

gmsh.model.mesh.generate(dim=2)
# gmsh.model.mesh.refine()
# gmsh.model.mesh.refine()
gmsh.model.mesh.setOrder(2)
gmsh.write("wires2d.msh")
# gmsh.fltk.run()
gmsh.finalize()
