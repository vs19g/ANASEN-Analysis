#!/home/vsitaraman/ParaView-6.1.0-RC1-MPI-Linux-Python3.12-x86_64/bin/pvbatch
# #######!/home/vsitaraman/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch
import numpy as np
import sys
from paraview.simple import *

# Optional argument so an archived VTU can be re-plotted without re-solving:
#   ./paraview_plotter.py wires2d/vtu_files/elfield_anasen_10_0.0000.vtu
INPUT_VTU = sys.argv[1] if len(sys.argv) > 1 else "wires2d/elfield_anasen_t0001.vtu"

# --- contour settings ---
# The whole-PC view keeps the original sparse, evenly spaced levels; at that
# scale more lines just fill in solid. The quadrant gets the denser, biased
# set, since that is where the structure between wires is actually readable.
POT_MIN, POT_MAX = 0.0, 660.0
N_CONTOURS_FULL, CONTOUR_BIAS_FULL = 40, 1.0
N_CONTOURS_QUARTER, CONTOUR_BIAS_QUARTER = 80, 1.0

# Line widths. Thin lines wash out their own colour, so the quadrant views draw
# heavier than the full-PC one; tubes keep them smooth at 2000 px.
LINE_WIDTH_FULL = 3.0
LINE_WIDTH_QUARTER = 2.5
STREAM_LINE_WIDTH = 2.5

# Arrows are glyphed onto the contour points, so doubling the quadrant contour
# density doubled the arrow count too. This offsets that.
GLYPH_STRIDE = 96

# --- quadrant view ---
# Centre sits at radius 35.4 mm, near the outer edge of the wire band, which
# pushes the wires to one side and leaves the far corner empty. Pulled in to
# 31.8 mm so the band sits more centrally.
ZOOM_CENTER = [-0.0225, 0.0225, 0.0]
ZOOM_SCALE = 0.022   # half-height of the view in metres; larger = more in frame


def contour_levels(n, bias):
    """Isosurface values.

    phi ~ ln(r) around a wire, so evenly spaced levels crowd onto the wire
    surfaces and thin out in the gas. Raising a symmetric parameter to `bias`
    pulls levels towards mid-range, where the saddles between wires are,
    without adding more of them. 1.0 is plain linear spacing; much above 1.5
    and the middle levels start landing on top of each other.
    """
    u = np.linspace(-1.0, 1.0, n)
    u = np.sign(u) * np.abs(u) ** bias
    mid = 0.5 * (POT_MIN + POT_MAX)
    return list(mid + 0.5 * (POT_MAX - POT_MIN) * u)


reader = XMLUnstructuredGridReader(FileName=[INPUT_VTU])

contour_filter = Contour(Input=reader,ContourBy = 'potential')
contour_filter.Isosurfaces = contour_levels(N_CONTOURS_FULL, CONTOUR_BIAS_FULL)

renderView = GetActiveViewOrCreate('RenderView')
renderView.ViewSize = [2000,2000]
renderView.OrientationAxesVisibility = 0 # Hide axis
renderView.UseColorPaletteForBackground=0
renderView.Background = [0.1, 0.1, 0.1] # Set background to dark gray (RGB 0-1)

renderView.MultiSamples = 8  # 0 disables it, 4-8 is usually sufficient

ResetCamera()

contour_display = Show(contour_filter, renderView)
contour_display.LineWidth = LINE_WIDTH_FULL
contour_display.RenderLinesAsTubes = 1    # Makes lines look smoother at high res
#colorbar
contour_display_potentialLUT = GetColorTransferFunction('potential', contour_display, separate=True)
contour_display_potentialLUT.ApplyPreset('Cool to Warm', True)
contour_display.SetScalarBarVisibility(renderView, True)

#axesGrid = renderView.AxesGridrfcxgdtv
#axesGrid.Visibility = 1
#axesGrid.XTitle = "x (mm)"
#axesGrid.YTitle = "y (mm)"

# 1. Get the active view
view = GetActiveView()

# 2. Define your desired coordinate ranges (x_min, x_max, y_min, y_max, z_min, z_max)
x_min, x_max = -0.05, 0.05
y_min, y_max = -0.05, 0.05
z_min, z_max = -0.05, 0.05

# 3. Calculate Center, Position, and Parallel Scale
center = [(x_min + x_max) / 2.0, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0]
# Position the camera far away along Z to look at the center
position = [center[0], center[1], 1.0] 
# Parallel scale defines how much of the scene is visible. 
# It is usually half the height of the viewed area.
view.CameraParallelScale = max((x_max - x_min), (y_max - y_min))/1.6

# 4. Apply settings
view.CenterOfRotation = center
view.CameraPosition = position
view.CameraFocalPoint = center
view.CameraViewUp = [0.0, 1.0, 0.0] # Y-axis is up

# 5. Enable Parallel Projection (optional, often better for exact mapping)
view.CameraParallelProjection = 1

#ResetCamera()
Render()

SaveScreenshot("contour_output.png")

# Everything from here on is a quadrant view: denser levels, heavier lines.
contour_filter.Isosurfaces = contour_levels(N_CONTOURS_QUARTER, CONTOUR_BIAS_QUARTER)
contour_display.LineWidth = LINE_WIDTH_QUARTER
contour_display.RenderLinesAsTubes = 1

# 1. Get the active view
view = GetActiveView()

# 1. Set the Focal Point to the middle of the quadrant in metres
zoom_center = ZOOM_CENTER

# 2. Tighten the Parallel Scale 
view.CameraParallelScale = ZOOM_SCALE

# 3. Position the Camera (0.5m away is fine)
view.CameraPosition = [zoom_center[0], zoom_center[1], 0.5]
view.CameraFocalPoint = zoom_center
view.CameraViewUp = [0.0, 1.0, 0.0]

# Equipotentials over one quadrant, saved before the arrows go on so this and
# Field_output.png below share a camera and overlay exactly.
Render()
SaveScreenshot("contour_quarter_output.png")

# pot_threshold = Threshold(Input=reader)
# pot_threshold.Scalars = ['POINTS', 'potential']
# pot_threshold.ThresholdMethod = 'Above Upper Threshold'
# pot_threshold.UpperThreshold = 100.0  

# --- 2. Create the Glyph Filter (The Arrows) ---
# IMPORTANT: Use 'pot_threshold' as the Input, not the 'reader'
glyph = Glyph(Input=contour_filter, GlyphType='Arrow') #
# glyph = Glyph(Input=reader, GlyphType='Arrow') #this uses all field line snot just the ones from the equipotential lines shown

# Orientation Array: Use the 'electric field' vector from Elmer
glyph.OrientationArray = ['POINTS', 'electric field']
glyph.ScaleArray = ['POINTS', 'No scale array']
glyph.ScaleFactor = 0.001  

glyph.GlyphMode = 'Every Nth Point'
glyph.Stride = GLYPH_STRIDE

# --- 3. Display the Glyphs ---
glyph_display = Show(glyph, renderView)

# Set the representation to Surface so we see the full arrow colors
glyph_display.Representation = 'Surface'


# This is the critical line: Color the arrows by the 'potential' scalar
ColorBy(glyph_display, ('POINTS', 'potential'))
glyph_display.LookupTable = contour_display_potentialLUT
contour_display_potentialLUT.RescaleTransferFunction(POT_MIN, POT_MAX)

# Optional: Disable the scalar bar for the arrows to avoid cluttering 
# the existing 'potential' scalar bar.
glyph_display.SetScalarBarVisibility(renderView, False)

# --- 4. Final Render ---
Render()
SaveScreenshot("Field_output.png")


# --- Streamlines over the same quadrant ---
# Field lines rather than arrows: with no B field the drift velocity is
# parallel to E, so these are drift paths up to diffusion, i.e. a map of which
# wire collects charge from where. The glyphs above are fixed length, so they
# show direction only. Quadrant only; over the whole PC it is an unreadable mat.
# Equipotentials the streamlines start from. One level near the anodes seeds
# lines only in a ring around each anode, so with the lines truncated short the
# rest of the frame comes out empty. Spanning the range instead puts seeds
# right across the gas: every equipotential is a closed curve winding between
# the wires, so together they cover the volume.
SEED_POTENTIALS = [660.0, 600.0, 550.0, 500.0, 450.0, 400.0,350.0, 300.0,250.0, 250.0, 200.0, 150.0, 100.0, 50.0]

# Keep 1 seed point in N along the seed contours. Measured point counts per
# level, for working out the total:
#
#     600 V   4970      300 V  29971
#     500 V   8861      200 V  51558
#     400 V  16331      100 V  51645
#                        50 V  27259
#
# 250 V and 150 V are not measured; from the trend they are of order 40000 and
# 50000, putting the nine levels near 280000 points. Stride 3 is then ~94000
# seeds, of which roughly an eighth land in the rendered quadrant. That is
# three times the previous seed count, but each line is now 4 mm instead of
# 20 mm, so the integration work is comparable.
SEED_STRIDE = 3
# Arc length each line is allowed to run, in metres. Field lines terminate on
# electrodes, so any line long enough to reach one will converge with every
# other line reaching the same wire -- that convergence is the field, not a
# plotting artefact, and the only way to avoid drawing it is to stop the lines
# short. The anode-cathode gap runs 4.3-6.2 mm depending on z, so 3 mm keeps
# every line clear of the cathodes. Raise it towards 0.005 for longer lines,
# and they will start meeting again at the wires.
# Short lines, and the gaps between them closed by seeding more planes rather
# than by making each line longer. 4 mm is just under the 4.3 mm closest
# anode-cathode approach, so no line quite reaches a wire and they stay
# separate instead of converging.
STREAM_MAX_LENGTH = 0.0005

# StreamTracer counts steps as well as length, and its step sizes are in
# cell-length units rather than metres: 0.2 of a cell by default, and cells are
# 0.05 mm at the wires. 2000 steps therefore buys about 20 mm at best and less
# wherever the mesh is fine, so the default caps the lines before
# STREAM_MAX_LENGTH ever does. Raise it well past what the length needs.
MAX_STEPS = 20000
STREAM_DIRECTION = 'BOTH'   # 'BOTH' also draws the half that runs into the anode
STREAM_RANGE = [-200.0, 650.0]   # full electrode span, needle to anode

# Clear the quadrant view of everything else so only the field lines show.
Hide(glyph, renderView)
Hide(contour_filter, renderView)

seed_contour = Contour(Input=reader, ContourBy=['POINTS', 'potential'])
seed_contour.Isosurfaces = SEED_POTENTIALS

# StreamTracer's own Point Cloud and Line seeds are 3D and this mesh is a
# single plane at z = 0, so they would mostly land outside it. Seeding from a
# contour of the mesh keeps every seed in the plane by construction.
seed_points = MaskPoints(Input=seed_contour)
seed_points.OnRatio = SEED_STRIDE
seed_points.RandomSampling = 0

stream = StreamTracerWithCustomSource(Input=reader, SeedSource=seed_points)
stream.Vectors = ['POINTS', 'electric field']
# FORWARD follows E, which runs from the anode outwards, so each line leaves
# its seed and heads for a cathode. With BOTH, the backward half of every line
# also climbs to the anode the seed ring encircles, and all of them meet there
# -- that is the joining. BOTH is still the honest full field line if you want
# it; this only changes what gets drawn.
stream.IntegrationDirection = STREAM_DIRECTION
stream.MaximumStreamlineLength = STREAM_MAX_LENGTH
stream.MaximumSteps = MAX_STEPS

# Keep integration in the z = 0 plane. Newer ParaView only, so tolerate it
# being absent.
try:
    stream.SurfaceStreamlines = 1
except Exception as exc:
    print("note: SurfaceStreamlines unavailable (%s)" % exc)

stream_display = Show(stream, renderView)
stream_display.Representation = 'Surface'
stream_display.LineWidth = STREAM_LINE_WIDTH
stream_display.RenderLinesAsTubes = 1

ColorBy(stream_display, ('POINTS', 'potential'))
streamLUT = GetColorTransferFunction('potential', stream_display, separate=True)
streamLUT.ApplyPreset('Cool to Warm', True)
streamLUT.RescaleTransferFunction(STREAM_RANGE[0], STREAM_RANGE[1])
stream_display.LookupTable = streamLUT
stream_display.SetScalarBarVisibility(renderView, True)

renderView.Background = [0.0, 0.0, 0.0]

Render()
SaveScreenshot("Streamlines_quarter_output.png")
