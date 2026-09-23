#!/home/vsitaraman/ParaView-6.1.0-RC1-MPI-Linux-Python3.12-x86_64/bin/pvbatch
# Point the shebang above at your own pvbatch/pvpython, as for paraview_plotter.py.
#
# Renders the E . E_w map written by dotproduct.py.
#
# By Shockley-Ramo the dot product describes the current induced on ONE
# electrode: the cathode that wires2d_weight.sif held at 1 V. Everywhere else
# the number is just the tail of that electrode's weighting field, so the view
# is centred on that wire and nothing else. Its position is read out of the
# weighting solve rather than recomputed from the geometry -- the driven
# electrode is wherever the weighting potential reaches 1.
#
# Usage (from anasen_fem/):
#   ./paraview_dotproduct.py [dotproduct.vtu] [output.png] [weight.vtu]
#
# Re-plotting an archived z needs its archived weighting solve as well, since
# that is where the driven cathode's position is read from:
#   ./paraview_dotproduct.py wires2d/vtu_files/dotproduct_10_0.0000.vtu \
#       png/DotProduct_output_z_10_0.0000.png \
#       wires2d/vtu_files/elfield_weight_10_0.0000.vtu

import sys

import numpy as np
from paraview.simple import *

INPUT_VTU = sys.argv[1] if len(sys.argv) > 1 else "wires2d/dotproduct.vtu"
OUTPUT_IMAGE = sys.argv[2] if len(sys.argv) > 2 else "DotProduct_output.png"

# Weighting solve, used only to locate the driven cathode. Must be the one that
# belongs to this dot product -- the wire moves with z.
WEIGHT_VTU = sys.argv[3] if len(sys.argv) > 3 else "wires2d/elfield_weight_t0001.vtu"

# Half-height of the view in metres. Cathodes sit every 15 degrees at a radius
# of about 36.4 mm, so the neighbours either side are 9.5 mm away; 14 mm takes
# them both in with a margin.
# Image size. The structure round the driven cathode is much wider than it is
# tall, so a square frame spends most of its area on empty gas.
VIEW_SIZE = [2000, 1200]

# CameraParallelScale, which is the half-HEIGHT of the view in metres. The
# width follows from the aspect: half-height * (VIEW_SIZE[0] / VIEW_SIZE[1]).
# So widening the frame alone does not crop anything, it just shows more
# sideways; the scale has to come down with it. 0.0084 at 2000x1200 keeps the
# same 28 mm width as the old square 0.014 and trims the height to 16.8 mm.
ZOOM_SCALE = 0.0084

# The driven cathode swings round the axis as z changes, so with a fixed "up"
# every z comes out rotated by a different amount and the shots cannot be laid
# side by side. Aligning up with the outward radial direction through the
# cathode puts the detector axis below and the barrel wall above in every
# frame, whatever z. Set False for plain screen-aligned axes.
ALIGN_RADIAL = True

# Shift of the view centre along the outward radial, in metres. Up in the frame
# is the outward radial, so a NEGATIVE value walks the camera in towards the
# detector axis and the cathode rides higher in the picture. Positive drops it.
# Purely framing; it does not move the wire.
CENTER_SHIFT = -0.0036

# Second, wider shot over the same quadrant paraview_plotter.py frames, so the
# dot product can be laid beside the field and contour views of that z. Centre
# is at x < 0, y > 0; negate the y for the lower-left quadrant instead.

# Colour range for DotProduct. Both ends negative: E and E_w point in broadly
# opposite senses over most of the drift volume.
# Most of the painted patch lies between -2e7 and -2e6. Against a range
# reaching -3.1e8 that is the top few percent of the scale, so it all came out
# at the white end of Greens. Narrowing the range spreads those values over the
# green; anything stronger simply clamps to the darkest end, which only affects
# the last fraction of a millimetre at the wire surface.
DOT_MIN, DOT_MAX = -8.0e7, -3.0e6

# Opacity of the painted map. Anything at least as strong as OPAQUE_BELOW is
# drawn solid; from there it ramps down to transparent at FADE_TO, so the grey
# shows through with a soft round edge rather than the ragged one a hard cut
# leaves.
#
# The plateau matters more than either endpoint. Ramping all the way from
# DOT_MIN (-3.1e8) meant the opacity at -5e7 was only about 0.13, so the whole
# patch read as washed out however the fade point was set.
OPAQUE_BELOW = -2.0e7
FADE_TO = -2.0e6

# Coarse cut, only to keep the far field out of the render. Must be nearer zero
# than FADE_TO so its hard edge lands where the surface is already fully
# transparent and cannot be seen.
PAINT_CUT = -5.0e5

# "signed-log" maps sign and magnitude together through a signed logarithm and
# a diverging colour map; "linear" is the earlier one-sided Greens scale, which
# the constants above still drive.
COLOR_MODE = "signed-log"
# X0 sets where the mapping turns from linear to logarithmic, and so how much
# of the weak far field survives at all: at 1e6 the data map to -2.51 .. 1.62,
# against -3.51 .. 2.61 at 1e5. Raise it to paint less.
SYMLOG_X0 = 1.0e6
SYMLOG_LIMIT = 2.5       # symmetric colour limit; match it to the mapped range
SYMLOG_OPAQUE = 0.6      # |value| at which the map is fully opaque
SYMLOG_FADE = 0.3        # |value| below which it is transparent

# Shape of the two opacity ramps, as opposed to where they sit.
#
# A piecewise function's midpoint describes the segment running RIGHTWARD from
# each node, and the two ramps run in opposite senses -- falling on the
# negative side, rising on the positive. So a single midpoint below 0.5 applied
# to both makes the negative lobe fade out sooner while the positive lobe goes
# opaque sooner, rendering the two polarities at different effective
# thresholds. Mirroring it keeps them symmetric: below 0.5 both sides go opaque
# sooner, so more is painted, equally on each.
#
# Sharpness is the independent dial for how abrupt the edge is: 0 smooth, 1 a
# hard step. 0.5 is already close to a step.
SYMLOG_MIDPOINT = 0.3
SYMLOG_SHARPNESS = 0.0
DIVERGING_PRESET = "PRGn"   # purple / white / green; zero sits at white

BACKGROUND = [0.35, 0.34, 0.33]   # grey

# Equipotential lines over the whole frame, painted gas or not. Off by default:
# the map is the subject here, and the lines are already in the contour plots.
SHOW_CONTOURS = False
CONTOUR_BY = "potential"
CONTOUR_LEVELS = list(np.arange(0, 660, 660 / 40.0))


def threshold_between(source, array, low, high):
    """Threshold filter, tolerating the pre-5.10 property names."""
    t = Threshold(Input=source, Scalars=["POINTS", array])
    try:
        t.LowerThreshold = low
        t.UpperThreshold = high
        t.ThresholdMethod = "Between"
    except Exception:
        t.ThresholdRange = [low, high]
    return t


# --- locate the driven cathode: the weighting potential is 1 on it ---
weight_reader = XMLUnstructuredGridReader(FileName=[WEIGHT_VTU])
driven = threshold_between(weight_reader, "potential", 0.99, 1.01)
driven.UpdatePipeline()

info = driven.GetDataInformation()
if info.GetNumberOfPoints() == 0:
    sys.exit(
        "no nodes with weighting potential ~1 in %s -- the weighting solve did\n"
        "not drive any electrode, so there is nothing to centre on." % WEIGHT_VTU
    )
b = info.GetBounds()
CATHODE = [0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.0]
print("driven cathode at x = %.5f m, y = %.5f m" % (CATHODE[0], CATHODE[1]))

radius = (CATHODE[0] ** 2 + CATHODE[1] ** 2) ** 0.5
outward = [CATHODE[0] / radius, CATHODE[1] / radius, 0.0] if radius > 0 else [0.0, 1.0, 0.0]

# Frame centre, offset from the wire along that radial.
CENTER = [
    CATHODE[0] + CENTER_SHIFT * outward[0],
    CATHODE[1] + CENTER_SHIFT * outward[1],
    0.0,
]

Delete(driven)
Delete(weight_reader)

# --- the dot product map ---
reader = XMLUnstructuredGridReader(FileName=[INPUT_VTU])

renderView = GetActiveViewOrCreate("RenderView")
renderView.ViewSize = VIEW_SIZE
renderView.OrientationAxesVisibility = 1
renderView.UseColorPaletteForBackground = 0
renderView.Background = BACKGROUND
renderView.MultiSamples = 8

if COLOR_MODE == "signed-log":
    # E . E_w spans about 3.5 decades negative and 2.6 positive, and 38% of the
    # nodes are positive -- charge drifting there induces current of the
    # opposite polarity on the cathode. A linear scale over one sign shows
    # neither the decades nor the sign, so map through a signed log:
    #
    #     sign(d) * log10(1 + |d| / SYMLOG_X0)
    #
    # sign() is not in the calculator namespace, hence d/(|d|+1); the +1 is
    # negligible beside the 1e5 scale and the expression tends to 0 anyway
    # where d does.
    calc = PythonCalculator(Input=reader)
    calc.ArrayName = "SignedLogDot"
    calc.ArrayAssociation = "Point Data"
    calc.Expression = (
        "(DotProduct/(abs(DotProduct)+1.0))*log10(1.0+abs(DotProduct)/%g)"
        % SYMLOG_X0
    )
    painted = calc
    array = "SignedLogDot"
    colour_range = (-SYMLOG_LIMIT, SYMLOG_LIMIT)   # symmetric, so zero is neutral
    preset = DIVERGING_PRESET
    # Transparent either side of zero, opaque past SYMLOG_OPAQUE, both signs.
    _m, _s = SYMLOG_MIDPOINT, SYMLOG_SHARPNESS
    opacity_points = [
        -SYMLOG_LIMIT, 1.0, 0.5, _s,
        -SYMLOG_OPAQUE, 1.0, 1.0 - _m, _s,   # mirrored: this ramp falls
        -SYMLOG_FADE, 0.0, 0.5, _s,
        SYMLOG_FADE, 0.0, _m, _s,            # this one rises
        SYMLOG_OPAQUE, 1.0, 0.5, _s,
        SYMLOG_LIMIT, 1.0, 0.5, _s,
    ]
    bar_title = "sign(E.Ew) log10(1+|E.Ew|/%g)" % SYMLOG_X0
else:
    painted = threshold_between(reader, "DotProduct", -1.0e30, PAINT_CUT)
    array = "DotProduct"
    colour_range = (DOT_MIN, DOT_MAX)
    preset = "Greens"
    opacity_points = [
        DOT_MIN, 1.0, 0.5, 0.0,
        OPAQUE_BELOW, 1.0, 0.5, 0.0,
        FADE_TO, 0.0, 0.5, 0.0,
    ]
    bar_title = "E . E_w  (V/m^2)"

surface_display = Show(painted, renderView)
surface_display.Representation = "Surface"
ColorBy(surface_display, ("POINTS", array))

dotLUT = GetColorTransferFunction(array)
dotLUT.ApplyPreset(preset, True)
dotLUT.RescaleTransferFunction(colour_range[0], colour_range[1])

# Points go (value, opacity, midpoint, sharpness) in ascending value.
dotPWF = GetOpacityTransferFunction(array)
dotPWF.Points = opacity_points
dotLUT.ScalarOpacityFunction = dotPWF
dotLUT.EnableOpacityMapping = 1

surface_display.SetScalarBarVisibility(renderView, True)
try:
    GetScalarBar(dotLUT, renderView).Title = bar_title
except Exception as exc:
    print("note: could not set the scalar bar title (%s)" % exc)

# --- contour lines over the whole field, including the unpainted gas ---
if SHOW_CONTOURS:
    contour_filter = Contour(Input=reader, ContourBy=["POINTS", CONTOUR_BY])
    contour_filter.Isosurfaces = CONTOUR_LEVELS

    contour_display = Show(contour_filter, renderView)
    contour_display.LineWidth = 1.0
    contour_display.DiffuseColor = [1.0, 1.0, 1.0]
    # ColorBy(rep, None) is rejected by ParaView 6.1 -- it forwards the
    # association as the string 'NONE' and raises -- so set the property.
    contour_display.ColorArrayName = [None, ""]
    contour_display.SetScalarBarVisibility(renderView, False)

# --- camera: parallel projection, centred on the driven cathode ---
# ParaView resets the camera on a view's FIRST render, throwing away anything
# set beforehand, so burn that render here and set the camera after it.
# Measured: scale 0.008 set before the first Render comes back as 0.021; set
# again afterwards it holds.
Render()

view = GetActiveView()
view.CameraParallelScale = ZOOM_SCALE
view.CenterOfRotation = CENTER
view.CameraPosition = [CENTER[0], CENTER[1], 1.0]
view.CameraFocalPoint = CENTER
if ALIGN_RADIAL:
    up = outward
    print("view rotated %.1f deg so the outward radial is up"
          % np.degrees(np.arctan2(up[0], up[1])))
else:
    up = [0.0, 1.0, 0.0]
view.CameraViewUp = up
view.CameraParallelProjection = 1

Render()
SaveScreenshot(OUTPUT_IMAGE)
print("Saved %s" % OUTPUT_IMAGE)

