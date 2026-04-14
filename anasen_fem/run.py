import code
import os

# val=-178.3
val=89.15
count=11   
while val<178.3+0.1:
    print(val)
    os.system("python3 wires_gmsh2d_bc.py "+str(val))
    os.system("ElmerGrid 14 2 wires2d.msh")
    os.system("ElmerSolver wires2d.sif")
    os.system("./paraview_plotter.py")
    os.system("cp wires2d.msh wires2d/mesh_files/wires2d%02d_%1.4f.msh"%(count,val))
    os.system("cp wires2d.sif wires2d/sif_files/wires2d_%02d_%1.4f.sif"%(count,val))
    os.system("cp wires2d/elfield_anasen_t0001.vtu wires2d/vtu_files/elfield_anasen_%02d_%1.4f.vtu"%(count,val))
    os.system("cp contour_output.png png/Contour_output_z_%02d_%1.4f.png"%(count,val))
    os.system("cp Field_output.png png/Field_ouput_z_%02d_%1.4f.png"%(count,val))
    val=val+17.83
    count = count + 1
    break

# os.system("tar -cvzf wiress2d/mesh.tar.gz wires2d/mesh_files")
# os.system("rm -rf wires2d/mesh_files/*")
# os.system("tar -cvzf wires2d/sif.tar.gz wires2d/sif_files")
# os.system("rm -rf wires2d/sif_files/*")
# os.system("tar -cvzf wires2d/vtu.tar.gz wires2d/vtu_files")
# os.system("rm -rf wires2d/vtu_files/*")