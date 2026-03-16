#!~/ParaView-6.1.0-RC1-MPI-Linux-Python3.12-x86_64/bin/pvbatch
import numpy as np
import sys
from paraview.simple import *

reader = XMLUnstructuredGridReader(FileName=["wires2d/elfield_anasen_t0001.vtu"])

contour_filter = Contour(Input=reader,ContourBy = 'potential')
contour_filter.Isosurfaces = [i for i in np.arange(0,660,650/24.)]

renderView = GetActiveViewOrCreate('RenderView')
renderView.ViewSize = [800,800]
renderView.OrientationAxesVisibility = 0 # Hide axis
renderView.UseColorPaletteForBackground=0
renderView.Background = [0.1, 0.1, 0.1] # Set background to dark gray (RGB 0-1)

renderView.MultiSamples = 8  # 0 disables it, 4-8 is usually sufficient

ResetCamera()

contour_display = Show(contour_filter, renderView)

#colorbar
contour_display_potentialLUT = GetColorTransferFunction('potential', contour_display, separate=True)
contour_display_potentialLUT.ApplyPreset('Cool to Warm', True)
contour_display.SetScalarBarVisibility(renderView, True)

#axesGrid = renderView.AxesGrid
#axesGrid.Visibility = 1
#axesGrid.XTitle = "x (mm)"
#axesGrid.YTitle = "y (mm)"

# 1. Get the active view
view = GetActiveView()

# 2. Define your desired coordinate ranges (x_min, x_max, y_min, y_max, z_min, z_max)
# Example: Look at a box from -10 to 10 in all dimensions
x_min, x_max = -50.0, 50.0
y_min, y_max = -50.0, 50.0
z_min, z_max = -50.0, 50.0

# 3. Calculate Center, Position, and Parallel Scale
center = [(x_min + x_max) / 2.0, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0]
# Position the camera far away along Z to look at the center
position = [center[0], center[1], z_min - 30.0] 
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
