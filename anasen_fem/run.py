import code
import os

val=-178.3
# val=17.83
count=0    
while val<17.83+0.1:
    print(val)
    os.system("python3 wires_gmsh2d_bc.py "+str(val))
    os.system("cp wires2d/mesh_files/wires2d.msh wires2d/mesh_files/%02d_%1.4f.msh"%(count,val))
    os.system("ElmerGrid 14 2 wires2d/mesh_files/wires2d_%02d_%1.4f.msh"%(count,val))
    os.system("cp wires2d.sif wires2d/sif_files/wires2d_%02d_%1.4f.sif"%(count,val))
    os.system("ElmerSolver wires2d/sif_files/wires2d_%02d_%1.4f.sif"%(count,val))
    os.system("./paraview_plotter.py")
    os.system("cp contour_output.png png/contour_output_z_%02d_%1.4f.png"%(count,val))
    # os.system("cp field_output.png field_ouput_z_%02d_%1.4f.png"%(count,val))

    val=val+17.83
    count = count + 1
    break

