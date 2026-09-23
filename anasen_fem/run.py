import os
import sys

# Cathode wire (0..23) that wires2d_weight.sif drives to 1 V. Passed to the
# mesher, which puts that one wire in its own physical group, tag 40.
SELECTED_CATHODE = 1

# 0 = everything in one pass
# 1 = mesh and the physical field only
# 2 = weighting field and dot product only, reusing stage 1's mesh
STAGE = 0

# gmsh comes from apt and lives in the system python3; meshio is pip-installed
# into the venv, which cannot see apt packages. Neither interpreter has both,
# so name them explicitly rather than relying on whichever python3 is on PATH.
MESH_PY = "/usr/bin/python3"
DOT_PY = os.path.expanduser("~/myenv/bin/python3")


def run(cmd):
    """os.system that stops on failure.

    Worth the four lines: a solve that quietly did nothing once produced a
    field of zeros, and every later step ran happily on it.
    """
    print(cmd)
    if os.system(cmd) != 0:
        sys.exit("failed: " + cmd)


# z-loci to run, in mm.
#
# The wires twist: over z = 0 -> 174.3 the anode lattice turns 60 deg against
# the cathode lattice, which is exactly 4 cathode pitches of 15 deg. So the
# field pattern goes through the same cycle four times across that range, while
# the wire radii drift slowly and monotonically (anodes 32.0 -> 37.0 mm).
#
# A uniform 17.43 step samples the fast cycle 2.75 times per period, which
# aliases it -- consecutive slices land at relative azimuths 7.5, 14.1, 5.7,
# 12.2, 3.5 deg, jumping around rather than sweeping. Sample the two variations
# separately instead:
#   - 8 slices across one pitch, which walks the azimuth in even ~2 deg steps
#     at essentially fixed radius
#   - 5 slices one whole pitch apart, which holds the azimuth near-fixed and
#     walks the radius over its full range
# 12 z-loci in total, against 11 for the uniform sweep.
PITCH_Z = 174.3 / 4.0

Z_VALUES = sorted(set(
    [round(i * PITCH_Z / 8.0, 4) for i in range(8)]
    + [round(i * PITCH_Z, 4) for i in range(5)]
))

# Running only two loci for now; delete this line for the full sweep.
# 19.794 is the crossover: the relative azimuth wraps through zero there, so an
# anode sits directly radially outside a cathode, and the closest anode-cathode
# approach over the whole range (4.339 mm) is at 19.707. It is deliberately off
# the grid above, whose nearest point is 21.7875.
Z_VALUES = [0.0, 19.794]

# Stage 2 reuses whatever single mesh is sitting in wires2d/, so it cannot span
# several z. Without this it would happily solve the weighting field on one z's
# mesh and archive the result under another z's name -- and dotproduct.py could
# not catch it, since the node sets would agree perfectly.
if STAGE == 2 and len(Z_VALUES) > 1:
    sys.exit(
        "STAGE 2 applies to the one mesh in wires2d/, which is whichever z\n"
        "stage 1 ran last. Set Z_VALUES to that single z."
    )

for count, val in enumerate(Z_VALUES, start=10):
    print(val)
    os.system("mkdir -p wires2d/mesh_files wires2d/sif_files wires2d/vtu_files png")

    # Each stage archives its own output before the next one starts, so a
    # failure in stage 2 does not throw away what stage 1 already produced.
    if STAGE in (0, 1):
        run("%s wires_gmsh2d_bc.py %s %d" % (MESH_PY, val, SELECTED_CATHODE))
        run("ElmerGrid 14 2 wires2d.msh -2d")
        run("ElmerSolver wires2d.sif")
        run("./paraview_plotter.py")

        os.system("cp wires2d.msh wires2d/mesh_files/wires2d%02d_%1.4f.msh"%(count,val))
        os.system("cp wires2d.sif wires2d/sif_files/wires2d_%02d_%1.4f.sif"%(count,val))
        os.system("cp wires2d/elfield_anasen_t0001.vtu wires2d/vtu_files/elfield_anasen_%02d_%1.4f.vtu"%(count,val))
        os.system("cp contour_output.png png/Contour_output_z_%02d_%1.4f.png"%(count,val))
        os.system("cp contour_quarter_output.png png/Contour_output_z_%02d_%1.4f_quarter.png"%(count,val))
        os.system("cp Field_output.png png/Field_ouput_z_%02d_%1.4f.png"%(count,val))
        os.system("cp Streamlines_quarter_output.png png/Streamlines_output_z_%02d_%1.4f_quarter.png"%(count,val))

    # Stage 2 does not re-mesh: it reuses wires2d/ as stage 1 left it, which is
    # what makes the node-by-node dot product valid.
    if STAGE in (0, 2):
        run("ElmerSolver wires2d_weight.sif")
        run("%s dotproduct.py" % DOT_PY)
        run("./paraview_dotproduct.py")

        os.system("cp wires2d_weight.sif wires2d/sif_files/wires2d_weight_%02d_%1.4f.sif"%(count,val))
        os.system("cp wires2d/elfield_weight_t0001.vtu wires2d/vtu_files/elfield_weight_%02d_%1.4f.vtu"%(count,val))
        os.system("cp wires2d/dotproduct.vtu wires2d/vtu_files/dotproduct_%02d_%1.4f.vtu"%(count,val))
        os.system("cp DotProduct_output.png png/DotProduct_output_z_%02d_%1.4f.png"%(count,val))
        os.system("cp DotProduct_quarter_output.png png/DotProduct_output_z_%02d_%1.4f_quarter.png"%(count,val))

    # os.system("python3 garfield_sim.py")

    # (a bare "break" here runs only the first z, as the original did)

# os.system("tar -cvzf wiress2d/mesh.tar.gz wires2d/mesh_files")
# os.system("rm -rf wires2d/mesh_files/*")
# os.system("tar -cvzf wires2d/sif.tar.gz wires2d/sif_files")
# os.system("rm -rf wires2d/sif_files/*")
# os.system("tar -cvzf wires2d/vtu.tar.gz wires2d/vtu_files")
# os.system("rm -rf wires2d/vtu_files/*")
