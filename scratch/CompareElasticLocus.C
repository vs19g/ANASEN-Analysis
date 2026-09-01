/***************************************************
 *
 * CompareElasticLocus.C
 *
 * DEDX_SCALE calibration check that sidesteps the excitation-energy
 * peak-fitting degeneracy entirely (see ScanDedxScalePosition.C's
 * pinning/boundary problems for why that approach got compromised).
 *
 * The idea: overlay a THEORETICAL kinematic curve on top of a
 * MEASURED locus.
 *
 *   MEASURED  : m27Alax_Ef_vs_theta_p_sx3 -- proton angle (purely
 *               geometric) vs Efix (built from the UNSCALED proton
 *               tables, per Eloss.py's mass_u>=10 condition). Neither
 *               ingredient depends on DEDX_SCALE, so this locus
 *               should be essentially identical across every scale
 *               folder -- this macro cross-checks that directly.
 *
 *   THEORETICAL: for each trial DEDX_SCALE, pull that scale's own
 *               reconstructed beam_energy_at_vertex from
 *               m27Alax_BeamEnergy_vs_VertexZ_sx3 (this DOES depend
 *               on DEDX_SCALE, via the aluminum/beam table), build a
 *               Kinematics object with it, and root-find the Ex=0
 *               (elastic/ground-state) locus across angle using
 *               predictElasticEnergy() -- copied verbatim from
 *               TrackRecon.C so the physics matches exactly.
 *
 * The DEDX_SCALE whose theoretical curve best tracks the measured
 * (DEDX_SCALE-independent) locus is your best calibration -- a
 * direct kinematic comparison, no spectral fitting involved.
 *
 * Only load THIS file. Requires Armory/Kinematics.h to be reachable
 * from wherever you compile this -- adjust the #include path below
 * if your directory layout differs from TrackRecon.C's.
 *
 * Usage:
 *
 *   .L CompareElasticLocus.C+
 *
 *  std::vector<double> scales = {0.70, 0.75, 0.80, 0.85, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.95, 1.00, 1.05, 1.10, 1.15};
 *
 *  ScoreElasticLocusEdge(scales, -1, "Output_27Al_", "output_27Al.root","_m27Alax+misc_sx3_p/m27Alax_Ef_vs_theta_p_sx3", "_m27Alax+misc_sx3/m27Alax_BeamEnergy_vs_VertexZ_sx3", 15, 50, 5);
 *
 *   CompareElasticLocus(scales);
 *
 ***************************************************/

#ifndef CompareElasticLocus_C
#define CompareElasticLocus_C

#include "Armory/Kinematics.h"

#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <cmath>
#include <vector>

// ---- mass constants, copied from TrackRecon.C to match exactly ----
static const double mass_27Al = 26.981538;
static const double mass_4He  = 4.002603254;
static const double mass_1H   = 1.007825032;
static const double mass_30Si = 29.973770;

// ---- predictElasticEnergy, copied verbatim from TrackRecon.C ----
// Root-finds the ejectile kinetic energy t3 at a given angle such
// that Kinematics::getExc(t3, angle) == 0 (the elastic/ground-state
// locus). Returns -1 if no single root exists in [t3_lo, t3_hi]
// (kinematically forbidden angle, or an ambiguous multi-valued locus).
inline double predictElasticEnergy(Kinematics &kin, double angle3_deg, double t3_lo = 0.001, double t3_hi = 60.0, int iters = 60)
{
  const int N = 200;
  double dt = (t3_hi - t3_lo) / N;
  int n_sign_changes = 0;
  double seg_lo = t3_lo, seg_hi = t3_hi;
  double prev = kin.getExc(t3_lo, angle3_deg);
  for (int k = 1; k <= N; ++k)
  {
    double t = t3_lo + k * dt;
    double cur = kin.getExc(t, angle3_deg);
    if (std::isfinite(prev) && std::isfinite(cur) && prev * cur < 0.0)
    {
      ++n_sign_changes;
      seg_lo = t - dt;
      seg_hi = t;
    }
    if (std::isfinite(cur))
      prev = cur;
  }
  if (n_sign_changes == 0)
    return -1.0; // no root in range (e.g. kinematically forbidden angle)
  if (n_sign_changes > 1)
    return -1.0; // ambiguous (multi-valued) locus -> reject
  double f_lo = kin.getExc(seg_lo, angle3_deg);
  for (int i = 0; i < iters; ++i)
  {
    double t3_mid = 0.5 * (seg_lo + seg_hi);
    double f_mid = kin.getExc(t3_mid, angle3_deg);
    if (!std::isfinite(f_mid))
      return -1.0;
    if (f_mid * f_lo <= 0.0)
      seg_hi = t3_mid;
    else
    {
      seg_lo = t3_mid;
      f_lo = f_mid;
    }
  }
  return 0.5 * (seg_lo + seg_hi);
}

// ---- shared helper: pull a scale's representative beam_energy_at_vertex ----
// Returns -1 on any failure (file/histogram missing, invalid profile bin).
inline double GetBeamEnergyAtVertex(TString folder, TString fileName, TString beamHist, double &repZ_out) {
  TFile *f = TFile::Open(folder + "/" + fileName, "READ");
  if (!f || f->IsZombie()) return -1;
  TH2 *hBeam = (TH2*) f->Get(beamHist);
  if (!hBeam) { f->ls(); return -1; }
  TProfile *prof = hBeam->ProfileX(Form("beamProf_%s", folder.Data()));
  double repZ = prof->GetMean();
  int zBin = prof->FindBin(repZ);
  double beamE = prof->GetBinContent(zBin);
  repZ_out = repZ;
  return beamE;
}

// ============================================================
// scales         : DEDX_SCALE values to overlay theoretical curves for
// referenceScale : which scale's folder to pull the MEASURED locus
//                  from -- shouldn't matter which, since the locus is
//                  independent of DEDX_SCALE (this is cross-checked
//                  automatically against a second folder, see below)
// folderPrefix, fileName : same convention as the other scan macros
// lociHist       : the measured (theta, Ef) 2D histogram
// beamHist       : the (vertex_z, beam_energy_at_vertex) 2D histogram
//                  used to read off each scale's beam_energy_at_vertex
// thetaMin,thetaMax : angle range (degrees) to draw theoretical curves over
// ============================================================
void CompareElasticLocus(
    std::vector<double> scales,
    double               referenceScale = -1,   // -1 = use scales[0]
    TString              folderPrefix   = "Output_27Al_",
    TString              fileName       = "output_27Al.root",
    TString              lociHist       = "_m27Alax+misc_sx3_p/m27Alax_Ef_vs_theta_p_sx3",
    TString              beamHist       = "_m27Alax+misc_sx3/m27Alax_BeamEnergy_vs_VertexZ_sx3",
    double               thetaMin       = 20,
    double               thetaMax       = 160
){
  if (scales.empty()) { printf("CompareElasticLocus: no scale values given.\n"); return; }
  if (referenceScale < 0) referenceScale = scales[0];

  // ---------- 1) measured locus from the reference folder ----------
  TString refFolder = Form("%s%.2f", folderPrefix.Data(), referenceScale);
  TFile *fref = TFile::Open(refFolder + "/" + fileName, "READ");
  if (!fref || fref->IsZombie()) {
    printf("ERROR: could not open reference folder file: %s\n", (refFolder + "/" + fileName).Data());
    return;
  }
  TH2 *hLocus = (TH2*) fref->Get(lociHist);
  if (!hLocus) {
    printf("ERROR: '%s' not found in reference folder\n", lociHist.Data());
    fref->ls();
    return;
  }
  hLocus->SetDirectory(0);

  // ---------- cross-check: does the locus actually look the same in
  //            a different scale's folder? (validates the core
  //            assumption this whole approach rests on) ----------
  if (scales.size() > 1) {
    double otherScale = (scales[0] == referenceScale && scales.size() > 1) ? scales[1] : scales[0];
    TString otherFolder = Form("%s%.2f", folderPrefix.Data(), otherScale);
    TFile *fother = TFile::Open(otherFolder + "/" + fileName, "READ");
    if (fother && !fother->IsZombie()) {
      TH2 *hOther = (TH2*) fother->Get(lociHist);
      if (hOther) {
        double n1 = hLocus->GetEntries(), n2 = hOther->GetEntries();
        double m1x = hLocus->GetMean(1), m2x = hOther->GetMean(1);
        double m1y = hLocus->GetMean(2), m2y = hOther->GetMean(2);
        printf("Cross-check: locus entries/means at scale %.2f vs %.2f:\n", referenceScale, otherScale);
        printf("  entries: %.0f vs %.0f\n", n1, n2);
        printf("  <theta>: %.3f vs %.3f\n", m1x, m2x);
        printf("  <Ef>   : %.3f vs %.3f\n", m1y, m2y);
        if (n1 > 0 && std::abs(n1 - n2) / n1 > 0.05)
          printf("  NOTE: entry counts differ by >5%% -- something upstream of this\n"
                 "  histogram (a cut, a gate) may depend on DEDX_SCALE after all;\n"
                 "  worth investigating before trusting the overlay below.\n");
      }
      fother->Close();
    }
  }

  // ---------- 2) draw the measured locus ----------
  TCanvas *c = new TCanvas("cElasticLocus", "Measured Ef vs theta with theoretical DEDX_SCALE curves", 1000, 700);
  hLocus->SetStats(0);
  hLocus->Draw("colz");

  int colors[] = {kRed, kOrange+7, kSpring+4, kGreen+2, kCyan+2, kAzure+1, kBlue, kViolet, kMagenta+1, kPink+1, kGray+2, kBlack};
  int nColors = 12;

  TLegend *leg = new TLegend(0.15, 0.60, 0.4, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetHeader("DEDX_SCALE");

  // ---------- 3) for each scale: get beam_energy_at_vertex, build the
  //            theoretical curve, overlay it ----------
  printf("\n%-10s %10s %16s\n", "scale", "rep. Z", "beamE@vertex");
  int idx = 0;
  for (double scale : scales) {
    TString folder = Form("%s%.2f", folderPrefix.Data(), scale);
    TFile *f = TFile::Open(folder + "/" + fileName, "READ");
    if (!f || f->IsZombie()) {
      printf("%-10.2f  ERROR: could not open folder\n", scale);
      continue;
    }
    TH2 *hBeam = (TH2*) f->Get(beamHist);
    if (!hBeam) {
      printf("%-10.2f  ERROR: '%s' not found\n", scale, beamHist.Data());
      f->ls();
      continue;
    }

    TProfile *prof = hBeam->ProfileX(Form("beamProf_%.2f", scale));
    double repZ = prof->GetMean(); // entries-weighted mean vertex_z
    int zBin = prof->FindBin(repZ);
    double beamE = prof->GetBinContent(zBin);
    f->Close();

    if (beamE <= 0) {
      printf("%-10.2f  ERROR: invalid beam energy at representative Z=%.1f\n", scale, repZ);
      continue;
    }
    printf("%-10.2f %10.1f %16.4f\n", scale, repZ, beamE);

    Kinematics kin(mass_27Al, mass_4He, mass_1H, mass_30Si, beamE / mass_27Al);

    TGraph *g = new TGraph();
    for (double th = thetaMin; th <= thetaMax; th += 1.0) {
      double Ef = predictElasticEnergy(kin, th);
      if (Ef > 0) g->SetPoint(g->GetN(), th, Ef);
    }
    if (g->GetN() < 2) {
      printf("%-10.2f  WARNING: theoretical curve has <2 valid points over [%.0f,%.0f] deg\n",
             scale, thetaMin, thetaMax);
    }
    g->SetLineColor(colors[idx % nColors]);
    g->SetLineWidth(2);
    g->Draw("L SAME");
    leg->AddEntry(g, Form("%.2f", scale), "l");
    idx++;
  }
  leg->Draw();

  printf("\nThe DEDX_SCALE whose colored curve best tracks the underlying measured\n"
         "(colz) density band is your best calibration candidate -- this comparison\n"
         "doesn't involve fitting the excitation-energy spectrum at all, so it's\n"
         "immune to the mean/sigma degeneracies we ran into scanning peak positions.\n\n");
}

// ============================================================
// ScoreElasticLocusEdge
//
// Quantifies what CompareElasticLocus's plot asks you to eyeball: in
// narrow angle slices across the DISCRIMINATING region (where
// theoretical curves for different DEDX_SCALE actually separate --
// typically low angle; check your CompareElasticLocus plot to see
// where curves diverge vs collapse together before trusting the
// default range here), extract the measured edge (a high percentile
// of Ef, since the Ex=0 locus is a boundary/edge feature, not the
// bulk of the statistics) and compare it against each scale's
// theoretical curve at the same angles. Reports summed squared
// residual vs DEDX_SCALE -- the minimum is your best candidate,
// as an actual number instead of a judgment call.
//
// thetaMin,thetaMax,thetaStep : angle slices to score at -- restrict
//                to wherever CompareElasticLocus showed real
//                separation between curves
// edgePercentile : the measured "edge" in each angle slice is defined
//                as the Ef below which this fraction of that slice's
//                counts lie (0.98 default -- near the top of the
//                distribution without being thrown off by single
//                stray high-Ef outlier bins)
// minEntriesPerSlice : angle slices with fewer total counts than this
//                are skipped (too little data to define an edge)
// ============================================================
void ScoreElasticLocusEdge(
    std::vector<double> scales,
    double               referenceScale     = -1,
    TString              folderPrefix       = "Output_27Al_",
    TString              fileName           = "output_27Al.root",
    TString              lociHist           = "m27Alax_Ef_vs_theta_p_sx3",
    TString              beamHist           = "m27Alax_BeamEnergy_vs_VertexZ_sx3",
    double               thetaMin           = 15,
    double               thetaMax           = 50,
    double               thetaStep          = 5,
    double               edgePercentile     = 0.98,
    int                  minEntriesPerSlice = 50
){
  if (scales.empty()) { printf("ScoreElasticLocusEdge: no scale values given.\n"); return; }
  if (referenceScale < 0) referenceScale = scales[0];

  TString refFolder = Form("%s%.2f", folderPrefix.Data(), referenceScale);
  TFile *fref = TFile::Open(refFolder + "/" + fileName, "READ");
  if (!fref || fref->IsZombie()) {
    printf("ERROR: could not open reference folder file: %s\n", (refFolder + "/" + fileName).Data());
    return;
  }
  TH2 *hLocus = (TH2*) fref->Get(lociHist);
  if (!hLocus) {
    printf("ERROR: '%s' not found in reference folder\n", lociHist.Data());
    fref->ls();
    return;
  }
  hLocus->SetDirectory(0);

  // ---------- 1) extract the measured edge in each angle slice ----------
  std::vector<double> thetaSlices, measuredEdge;
  TAxis *xax = hLocus->GetXaxis();
  for (double th = thetaMin; th <= thetaMax; th += thetaStep) {
    int b1 = xax->FindBin(th - thetaStep/2.0);
    int b2 = xax->FindBin(th + thetaStep/2.0);
    TH1D *slice = hLocus->ProjectionY(Form("slice_%.1f", th), b1, b2);

    double total = slice->Integral();
    if (total < minEntriesPerSlice) {
      printf("theta=%.1f: skipped (only %.0f entries, need >= %d)\n", th, total, minEntriesPerSlice);
      delete slice;
      continue;
    }

    double cum = 0, edge = -1;
    int nb = slice->GetNbinsX();
    for (int b = 1; b <= nb; b++) {
      cum += slice->GetBinContent(b);
      if (cum / total >= edgePercentile) { edge = slice->GetBinCenter(b); break; }
    }
    delete slice;

    if (edge <= 0) continue;
    thetaSlices.push_back(th);
    measuredEdge.push_back(edge);
    printf("theta=%.1f: measured edge (p%.0f) = %.3f MeV\n", th, edgePercentile*100, edge);
  }

  if (thetaSlices.empty()) {
    printf("ERROR: no usable angle slices -- widen [thetaMin,thetaMax], lower minEntriesPerSlice,\n"
           "or lower edgePercentile.\n");
    return;
  }

  // ---------- 2) score each scale against the measured edge ----------
  TGraph *gScore = new TGraph();
  gScore->SetTitle("Sum-squared edge residual vs DEDX_SCALE;DEDX_SCALE;#Sigma(theory - measured edge)^{2} [MeV^{2}]");
  gScore->SetMarkerStyle(20);

  printf("\n%-10s %16s\n", "scale", "sum sq. resid.");
  double bestScale = -1, bestScore = 1e18;
  for (double scale : scales) {
    TString folder = Form("%s%.2f", folderPrefix.Data(), scale);
    double repZ;
    double beamE = GetBeamEnergyAtVertex(folder, fileName, beamHist, repZ);
    if (beamE <= 0) {
      printf("%-10.2f  ERROR: invalid/missing beam energy\n", scale);
      continue;
    }

    Kinematics kin(mass_27Al, mass_4He, mass_1H, mass_30Si, beamE / mass_27Al);

    double sumSq = 0;
    int nUsed = 0;
    for (size_t i = 0; i < thetaSlices.size(); i++) {
      double theory = predictElasticEnergy(kin, thetaSlices[i]);
      if (theory <= 0) continue; // kinematically forbidden / ambiguous at this angle
      double resid = theory - measuredEdge[i];
      sumSq += resid * resid;
      nUsed++;
    }

    if (nUsed == 0) {
      printf("%-10.2f  ERROR: theoretical curve invalid at every scored angle\n", scale);
      continue;
    }

    printf("%-10.2f %16.4f  (%d/%d angles used)\n", scale, sumSq, nUsed, (int)thetaSlices.size());
    gScore->SetPoint(gScore->GetN(), scale, sumSq);
    if (sumSq < bestScore) { bestScore = sumSq; bestScale = scale; }
  }

  TCanvas *c = new TCanvas("cScoreElasticLocusEdge", "DEDX_SCALE via elastic-locus edge matching", 900, 650);
  gScore->Draw("APL");

  if (bestScale > 0) {
    printf("\nBest-matching DEDX_SCALE among tested values (minimum summed residual): %.2f\n", bestScale);

    // ---------- parabolic refinement near the minimum ----------
    // Fit only the points close to the discrete minimum, rather than
    // the whole curve, since the score isn't parabolic far from the
    // minimum (see the steep rise at the edges of your scan) -- a
    // global parabola fit would be pulled around by those points.
    double fitWindow = 3 * (scales.size() > 1 ? std::abs(scales[1]-scales[0]) : 0.1);
    // widen the window a bit so at least a handful of points are
    // typically included even with uneven scale spacing
    fitWindow = std::max(fitWindow, 0.08);

    TF1 *parab = new TF1("parab", "pol2", bestScale - fitWindow, bestScale + fitWindow);
    TFitResultPtr fr = gScore->Fit(parab, "RSQ"); // R: restrict to window, S: get result, Q: quiet
    parab->SetLineColor(kRed);
    parab->SetLineStyle(2);
    parab->Draw("SAME");

    double a = parab->GetParameter(2), b = parab->GetParameter(1);
    if (a > 0) { // sanity: should open upward near a true minimum
      double vertexScale = -b / (2*a);
      printf("Parabolic refinement (fit window: scale in [%.3f, %.3f]):\n", bestScale-fitWindow, bestScale+fitWindow);
      printf("  analytic minimum at DEDX_SCALE = %.4f\n\n", vertexScale);
    } else {
      printf("Parabolic fit near the minimum did not open upward (a=%.3g) -- the\n"
             "points there may be too flat/noisy for a reliable sub-grid estimate;\n"
             "trust the discrete best value (%.2f) instead.\n\n", a, bestScale);
    }
  }
}

#endif