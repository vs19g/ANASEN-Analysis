#!/home/vs19g/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch
import numpy as np
import sys
from paraview.simple import *

reader = XMLUnstructuredGridReader(FileName=["wires2d/elfield_anasen_t0001.vtu"])

contour_filter = Contour(Input=reader,ContourBy = 'potential')
contour_filter.Isosurfaces = [i for i in np.arange(0,660,660/40.)]

renderView = GetActiveViewOrCreate('RenderView')
renderView.ViewSize = [2000,2000]
renderView.OrientationAxesVisibility = 0 # Hide axis
renderView.UseColorPaletteForBackground=0
renderView.Background = [0.1, 0.1, 0.1] # Set background to dark gray (RGB 0-1)

renderView.MultiSamples = 8  # 0 disables it, 4-8 is usually sufficient

ResetCamera()

contour_display = Show(contour_filter, renderView)
contour_display.LineWidth = 3.0          # Increase this for thicker lines
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

#make glyps for field lines
contour_display.LineWidth = 1.0          # Increase this for thicker lines
contour_display.RenderLinesAsTubes = 0    # Makes lines look smoother at high res

# 1. Get the active view
view = GetActiveView()

# 1. Set the Focal Point to the middle of the quadrant in metres
zoom_center = [-0.025, 0.025, 0.0] 

# 2. Tighten the Parallel Scale 
view.CameraParallelScale = 0.015 

# 3. Position the Camera (0.5m away is fine)
view.CameraPosition = [zoom_center[0], zoom_center[1], 0.5]
view.CameraFocalPoint = zoom_center
view.CameraViewUp = [0.0, 1.0, 0.0]
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
glyph.Stride = 24

# --- 3. Display the Glyphs ---
glyph_display = Show(glyph, renderView)

# Set the representation to Surface so we see the full arrow colors
glyph_display.Representation = 'Surface'


# This is the critical line: Color the arrows by the 'potential' scalar
ColorBy(glyph_display, ('POINTS', 'potential'))
glyph_display.LookupTable = contour_display_potentialLUT
contour_display_potentialLUT.RescaleTransferFunction(0.0, 660.0)

# Optional: Disable the scalar bar for the arrows to avoid cluttering 
# the existing 'potential' scalar bar.
glyph_display.SetScalarBarVisibility(renderView, False)

# --- 4. Final Render ---
Render()
SaveScreenshot("Field_output.png")
