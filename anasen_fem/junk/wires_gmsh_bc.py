import numpy as np
import gmsh

gmsh.initialize()
#gmsh.model.add("adaptive_mesh")
gmsh.option.setNumber('General.NumThreads', 4)
#gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfElements", 200000)
#gmsh.option.setNumber("Mesh.Adapt.MaxNumberOfNodes", 200000)
#gmsh.option.setNumber("Mesh.Adapt.MaxIter",5)
#gmsh.option.setNumber("Mesh.MeshSizeMin", 5e-3)
#gmsh.option.setNumber("Mesh.MeshSizeMax", 10.0)
gmsh.option.setNumber("Geometry.Tolerance", 1e-2)
#gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)

lc = 0.04
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

wire_radius = 0.254 #mm
anode_wires = []
cathode_wires = []
aw_tags = [(3,i) for i in range(24)]
cw_tags = [(3,i+24) for i in range(24)]

for i,[xa,ya,xc,yc,xa2,ya2,xc2,yc2] in enumerate(zip(xarra_1,yarra_1,xarrc_1,yarrc_1,xarra_2,yarra_2,xarrc_2,yarrc_2)):
    print(i,xa,ya,-178.3,xc,yc,-178.3,xa2,ya2,178.3,xc2,yc2,178.3)
    pa1 = gmsh.model.occ.addPoint(xa,ya,-178.3,lc)
    pa2 = gmsh.model.occ.addPoint(xa2,ya2,178.3,lc)
    pc1 = gmsh.model.occ.addPoint(xc,yc,-178.3,lc)
    pc2 = gmsh.model.occ.addPoint(xc2,yc2,178.3,lc)
    linea = gmsh.model.occ.addLine(pa1,pa2)
    linec = gmsh.model.occ.addLine(pc1,pc2)
    anode_wires.append(linea)
    cathode_wires.append(linec)

#anasen_barrel = gmsh.model.occ.addCylinder(0,0,-500,0,0,500+605,300,1234) #tag 1234
anasen_barrel = gmsh.model.occ.addCylinder(0,0,-200, 0,0,400, 300) #tag 1234

gmsh.model.occ.synchronize()
gmsh.model.mesh.embed(1,anode_wires+cathode_wires,3,anasen_barrel)

f1 = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(f1,"CurvesList",anode_wires+cathode_wires)

f2 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(f2,"InField",f1)
gmsh.model.mesh.field.setNumber(f2,"SizeMin",0.08)
gmsh.model.mesh.field.setNumber(f2,"SizeMax",5)
gmsh.model.mesh.field.setNumber(f2,"DistMin",1)
gmsh.model.mesh.field.setNumber(f2,"DistMax",20)

gmsh.model.mesh.field.setAsBackgroundMesh(f2)

gmsh.model.addPhysicalGroup(1, anode_wires, tag=10)
gmsh.model.setPhysicalName(1,10,"anode_wires")

gmsh.model.addPhysicalGroup(1, cathode_wires, tag=20)
gmsh.model.setPhysicalName(1,20,"cathode_wires")

gmsh.option.setNumber("Mesh.Algorithm",6)
gmsh.option.setNumber("Mesh.Algorithm3D", 10) # For 3D meshes

gmsh.model.mesh.generate(dim=3)
gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()
#gmsh.model.mesh.refine()
gmsh.fltk.run()
gmsh.finalize()
