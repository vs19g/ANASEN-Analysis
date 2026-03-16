import numpy as np
import gmsh,sys

gmsh.initialize()
gmsh.model.add("adaptive_mesh")
gmsh.option.setNumber('General.NumThreads', 4)
#gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfElements", 200000)
#gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfNodes", 200000)
#gmsh.option.setNumber("Mesh.Adapt.MaxIter",5)
#gmsh.option.setNumber("Mesh.MeshSizeMin", 5e-3)
#gmsh.option.setNumber("Mesh.MeshSizeMax", 10.0)
gmsh.option.setNumber("Geometry.Tolerance", 4e-2)
#gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

lc = 0.04
#z_loc = -178.3

if len(sys.argv) < 2:
    print("Usage: python3 wires_gmsh2d_bc.py <z_locus in mm>")
    quit()

z_loc = float(sys.argv[1])

k=(2*np.pi/24.)

#anodes, plane 1 at -zmax/2
k=-2*np.pi/24.
offset = 6*k + 3*k #-pi/2
xarra_1 = np.array([37*np.cos(k*i+offset) for i in np.arange(0,24)])
yarra_1 = np.array([37*np.sin(k*i+offset) for i in np.arange(0,24)])

#cathodes, plane 1 at -zmax/2
kc=2*np.pi/24.
offsetc = -4*kc + 2*kc -  np.pi/24 #-pi/4
xarrc_1 = np.array([42*np.cos(kc*i+offsetc) for i in np.arange(0,24)])
yarrc_1 = np.array([42*np.sin(kc*i+offsetc) for i in np.arange(0,24)])

#anodes, plane 2 at +zmax/2
offset = offset-3*k
xarra_2 = np.array([37*np.cos(k*i+offset) for i in np.arange(0,24)])
yarra_2 = np.array([37*np.sin(k*i+offset) for i in np.arange(0,24)])

#cathodes, plane2 at +zmax/2
offsetc = offsetc-3*kc
xarrc_2 = np.array([42*np.cos(kc*i+offsetc) for i in np.arange(0,24)])
yarrc_2 = np.array([42*np.sin(kc*i+offsetc) for i in np.arange(0,24)])

direction_anodes_x = xarra_2 - xarra_1
direction_anodes_y = yarra_2 - yarra_1

direction_cathodes_x = xarrc_2 - xarrc_1
direction_cathodes_y = yarrc_2 - yarrc_1

t = (z_loc+178.3)/(2*178.3) #z=-178.3 is 0, z=+178.3 is 1
xloc_a = xarra_1 + t*direction_anodes_x
yloc_a = yarra_1 + t*direction_anodes_y
xloc_c = xarrc_1 + t*direction_cathodes_x
yloc_c = yarrc_1 + t*direction_cathodes_y

wire_radius_a = 0.018 #mm
wire_radius_c = 0.0762 #mm
anode_wires = []
cathode_wires = []
aw_tags = [(3,i) for i in range(24)]
cw_tags = [(3,i+24) for i in range(24)]

#for i,[xa,ya,xc,yc] in enumerate(zip(xarra_1,yarra_1,xarrc_1,yarrc_1)):
for i,[xa,ya,xc,yc] in enumerate(zip(xloc_a,yloc_a,xloc_c,yloc_c)):
    print(i,xa,ya,-178.3,xc,yc,-178.3)
    adisk = gmsh.model.occ.addDisk(xa,ya,0,wire_radius_a,wire_radius_a)
    cdisk = gmsh.model.occ.addDisk(xc,yc,0,wire_radius_c,wire_radius_c)
    anode_wires.append(adisk)
    cathode_wires.append(cdisk)

anasen_barrel = gmsh.model.occ.addDisk(0,0,0,500,500)
#gmsh.model.occ.synchronize()
#gmsh.model.mesh.embed(1,anode_wires+cathode_wires,2,anasen_barrel)

gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
gmsh.model.occ.synchronize()

awire_surfs = []
for w in anode_wires:
    awire_surfs += [s[1] for s in gmsh.model.getBoundary([(2,w)], oriented=False) if s[0] == 1]

cwire_surfs = []
for w in cathode_wires:
    cwire_surfs += [s[1] for s in gmsh.model.getBoundary([(2,w)], oriented=False) if s[0] == 1]
gmsh.model.mesh.embed(1,cwire_surfs+awire_surfs,2,anasen_barrel)

for s in gmsh.model.getBoundary([(2,w)],oriented=False):
    if s[0] == 1:
        anasen_bdry=s[1]


f1 = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f1,"CurvesList",cwire_surfs+awire_surfs)

f2 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f2,"InField",f1)
gmsh.model.mesh.field.setNumber(f2,"SizeMin",0.1)
gmsh.model.mesh.field.setNumber(f2,"SizeMax",5)
gmsh.model.mesh.field.setNumber(f2,"DistMin",1)
gmsh.model.mesh.field.setNumber(f2,"DistMax",20)

gmsh.model.mesh.field.setAsBackgroundMesh(f2)

gmsh.model.addPhysicalGroup(1, awire_surfs, tag=10)
gmsh.model.setPhysicalName(1,10,"anode_wires")

gmsh.model.addPhysicalGroup(1, cwire_surfs, tag=20)
gmsh.model.setPhysicalName(1,20,"cathode_wires")

#gmsh.model.addPhysicalGroup(1, [anasen_bdry], tag=30)
#gmsh.model.setPhysicalName(1,30,"barrel_boundary")

gmsh.model.addPhysicalGroup(2,[anasen_barrel],tag=13)
gmsh.model.setPhysicalName(1,13,"gas")

gmsh.option.setNumber("Mesh.Algorithm", 6)

gmsh.model.mesh.generate(dim=2)
gmsh.model.mesh.refine()
gmsh.write("wires2d.msh")
#gmsh.fltk.run()
gmsh.finalize()
