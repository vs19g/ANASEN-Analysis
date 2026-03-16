import numpy as np
import gmsh

gmsh.initialize()
gmsh.model.add("adaptive_mesh")
gmsh.option.setNumber('General.NumThreads', 4)
#gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfElements", 200000)
#gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfNodes", 200000)
#gmsh.option.setNumber("Mesh.Adapt.MaxIter",5)
#gmsh.option.setNumber("Mesh.MeshSizeMin", 5e-3)
#gmsh.option.setNumber("Mesh.MeshSizeMax", 10.0)
#gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)

lc = 0.01

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

pa1 = []
pa2 = []
pc1 = []
pc2 = []

wire_radius = 0.254 #mm
anode_wires = []
cathode_wires = []
aw_tags = [(3,i) for i in range(24)]
cw_tags = [(3,i+24) for i in range(24)]

for i,[xa,ya,xc,yc,xa2,ya2,xc2,yc2] in enumerate(zip(xarra_1,yarra_1,xarrc_1,yarrc_1,xarra_2,yarra_2,xarrc_2,yarrc_2)):
    print(i,xa,ya,-178.3,xc,yc,-178.3,xa2,ya2,178.3,xc2,yc2,178.3)
    anode_wires.append(gmsh.model.occ.addCylinder(xa,ya,-178.3,(xa2-xa),(ya2-ya),178.3*2,wire_radius,i)) #x,y,z of first face center, dx,dy,dz of the axis, then the wire radius
    cathode_wires.append(gmsh.model.occ.addCylinder(xc,yc,-178.3,(xc2-xc),(yc2-yc),178.3*2,wire_radius,i+24)) #cathode tags 24-47, anode 0-23


anasen_barrel = gmsh.model.occ.addCylinder(0,0,-500,0,0,500+605,300,1234) #tag 1234
#anasen_barrel = gmsh.model.occ.addCylinder(0,0,-500,0,0,500+605,300,1234) #tag 1234

gmsh.model.occ.synchronize()

all_wires = aw_tags+cw_tags
gmsh.model.occ.fragment([(3,1234)],all_wires)
gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()
gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
gmsh.option.setNumber("Geometry.OCCFixDegenerated", 1)
gmsh.option.setNumber("Geometry.OCCFixSmallEdges", 1)
gmsh.option.setNumber("Geometry.OCCFixSmallFaces", 1)

wire_surfs = []
for w in anode_wires + cathode_wires:
    wire_surfs += [s[1] for s in gmsh.model.getBoundary([(3,w)], oriented=False) if s[0] == 2]
#'''
f1 = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f1, "FacesList", wire_surfs) # Example curves
f2 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f2, "InField", f1)
gmsh.model.mesh.field.setNumber(f2, "SizeMin", 0.05)
gmsh.model.mesh.field.setNumber(f2, "SizeMax", 5.)
gmsh.model.mesh.field.setNumber(f2, "DistMin", 1.)
gmsh.model.mesh.field.setNumber(f2, "DistMax", 20.)
gmsh.model.mesh.field.setAsBackgroundMesh(f2)
#'''
gmsh.option.setNumber("Mesh.Algorithm", 2)
gmsh.option.setNumber("Mesh.Algorithm3D", 1) # For 3D meshes

gmsh.model.mesh.generate(dim=3)
#gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()
gmsh.fltk.run()
gmsh.finalize()
