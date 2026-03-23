import os

#val=-178.3
val=17.83
count=11    
while val<178.3+0.1:
    print(val)
    os.system("python3 wires_gmsh2d_bc.py "+str(val))
    os.system("ElmerGrid 14 2 wires2d.msh")
    os.system("ElmerSolver wires2d.sif")
    os.system("./paraview_plotter.py")
    os.system("cp contour_output.png contour_output_z_%02d_%1.4f.png"%(count,val))
    val=val+17.83
    count = count + 1

