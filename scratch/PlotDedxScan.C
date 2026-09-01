/***************************************************
 *
 * PlotDedxScan.C
 *
 * Builds the comma-separated file/label lists make_pretty.C's
 * multi-file overlay function expects, straight from your
 * DEDX_SCALE scan -- so you can superimpose every scale's histogram
 * for a given gate (a1c2fix, a1c1c2, etc.) in one call, without
 * typing out each file path and label by hand.
 *
 * Only load THIS file -- it already #includes make_pretty.C.
 *
 * Usage:
 *
 *   .L PlotDedxScan.C
 **   std::vector<double> scales = {0.70, 0.75, 0.80, 0.85, 0.87, 0.90, 0.92, 0.95, 1.00, 1.05, 1.10, 1.15};
 *   PlotDedxScan(scales,
 *       "_m27Alax+misc_sx3_p/m27Alax_Ex_from_p_a1c1c2_sx3");
 *
 *   // or for the a1c2fix gate:
 *   PlotDedxScan(scales,
 *       "_m27Alax+misc_sx3_p/m27Alax_Ex_from_p_a1c2fix_sx3");
 *
 ***************************************************/

#ifndef PlotDedxScan_C
#define PlotDedxScan_C

#include "scratch/make_prettyplots.C"
#include <vector>

// ============================================================
// scales       : the DEDX_SCALE values from your bash scan, e.g.
//                {0.70, 0.75, ..., 1.20} -- MUST match your folder
//                names exactly (same decimal formatting), since each
//                path is built as folderPrefix + scale + "/" + fileName
// histPath     : histogram path inside each file (same for every
//                scale), e.g. "_m27Alax+misc_sx3_p/m27Alax_Ex_from_p_a1c1c2_sx3"
// folderPrefix : e.g. "Output_27Al_"
// fileName     : ROOT file name inside each folder, e.g. "output_27Al.root"
// xlabel,ylabel: axis titles for the overlay plot
// xMin,xMax    : x-axis zoom range -- set these to the region where your
//                peaks actually live (e.g. -3, 9) to cut out the flat,
//                empty tails and make the overlaid curves easier to read.
//                Leave at -9999 for the full auto range.
// yMin,yMax    : optional y-axis range; leave at -9999 for auto
// canvasW,canvasH : output image size in pixels. Defaults here (3000x2100)
//                are larger than make_pretty.C's own defaults (2100x1575)
//                since a busy multi-curve overlay benefits from the extra
//                resolution -- especially once zoomed in with xMin/xMax.
// maxOverlay    : if scales.size() exceeds this, evenly subsample down to
//                maxOverlay curves (always keeping the first and last)
//                rather than plotting everything -- past ~5-6 overlapping
//                histograms in the same narrow x-range, more curves stops
//                helping and just makes the plot harder to read regardless
//                of resolution. Set to 0 to disable subsampling entirely.
// ============================================================
void PlotDedxScan(
    std::vector<double> scales,
    TString              histPath,
    TString              folderPrefix = "Output_27Al_",
    TString               fileName    = "output_27Al.root",
    TString              xlabel       = "Excitation Energy [MeV]",
    TString              ylabel       = "Counts",
    double               xMin         = -4.0,
    double               xMax         = 8.0,
    double               yMin         = -9999.0,
    double               yMax         = -9999.0,
    int                  canvasW      = 3000,
    int                  canvasH      = 2100,
    int                  maxOverlay   = 6
){
  if (maxOverlay > 0 && (int)scales.size() > maxOverlay) {
    std::vector<double> subset;
    int n = (int) scales.size();
    for (int k = 0; k < maxOverlay; k++) {
      int idx = (maxOverlay == 1) ? 0 : (int) std::round(k * (n - 1) / double(maxOverlay - 1));
      subset.push_back(scales[idx]);
    }
    printf("PlotDedxScan: %d scales given, subsampling to %d for readability: ", n, maxOverlay);
    for (double s : subset) printf("%.2f ", s);
    printf("\n(pass maxOverlay=0, or maxOverlay >= %d, to plot all of them)\n", n);
    scales = subset;
  }
  if (scales.empty()) {
    printf("PlotDedxScan: no scale values given.\n");
    return;
  }

  TString filesCSV, labelsCSV;
  for (size_t i = 0; i < scales.size(); i++) {
    TString folder = Form("%s%.2f", folderPrefix.Data(), scales[i]);
    TString path   = folder + "/" + fileName;

    if (i > 0) { filesCSV += ","; labelsCSV += ","; }
    filesCSV  += path;
    labelsCSV += Form("%.2f", scales[i]);
  }

  printf("Overlaying %d scale value(s) for histogram: %s\n", (int) scales.size(), histPath.Data());

  // matches make_pretty.C's multi-file overlay overload exactly:
  // (filesCSV, labelsCSV, histName, xAxisLabel, yAxisLabel, yMin, yMax, xMin, xMax, canvasW, canvasH)
  make_prettyplots(filesCSV, labelsCSV, histPath, xlabel, ylabel, yMin, yMax, xMin, xMax, canvasW, canvasH);
}

#endif