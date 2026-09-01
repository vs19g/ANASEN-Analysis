/***************************************************
 *
 * ScanDedxGasCalib.C
 *
 * Direct first-principles DEDX_SCALE calibration check, using
 * TrackRecon.C's "BeamEnergy_ETrack_vs_EKin" plot -- NOT
 * dEgasPred_vs_dEgasCalib (an earlier dead end: that histogram is
 * built entirely from the UNSCALED proton tables, so it's blind to
 * DEDX_SCALE and would look identical at every scale).
 *
 * From TrackRecon.C:
 *   double snapped_level = snapToNearestLevel(Ex, levels_30Si_MeV, level_residual);
 *   double ebeam_kin_MeV = invertBeamEnergyMeV(m_beam, mass_4He, m3, m4,
 *                              Efix, theta * 180/M_PI, snapped_level);
 *   Fill2D(..., beam_energy_at_vertex, ebeam_kin_MeV, ...);
 *
 * x = beam_energy_at_vertex : built from the ALUMINUM table --
 *                              DEPENDS on DEDX_SCALE
 * y = ebeam_kin_MeV         : the beam energy at the vertex REQUIRED
 *                              for this event's measured (Efix, theta)
 *                              to exactly match a known literature
 *                              level, solved purely from Efix
 *                              (unscaled proton tables) and theta
 *                              (geometry) -- INDEPENDENT of DEDX_SCALE
 *
 * If DEDX_SCALE is correct, x should equal y for every event -- a
 * clean y=x diagonal. This covers FIVE literature levels at once
 * (levels_30Si_MeV = {0.0, 2.235, 3.498, 6.550, 6.870}), not just the
 * ground state the way CompareElasticLocus.C did, giving more
 * statistics and letting you check whether the y=x offset is
 * constant (a flat DEDX_SCALE correction is the right model) or
 * drifts with beam energy (it isn't -- see the QQQ/SX3 discussion).
 *
 * Usage:
 *
 *   .L ScanDedxScale.C+
 *
 *   std::vector<double> scales = {0.70,0.75,0.80,0.85,0.87,0.88,0.89,0.90,0.91,0.92,0.95,1.00,1.05,1.10,1.15};
 *
 *   ScanDedxScale(scales);                                    // SX3
 *   ScanDedxScale(scales, "Output_27Al_", "output_27Al.root", "_m27Alax+misc_qqq_p/m27Alax_BeamEnergy_ETrack_vs_EKin_p_a2c0_qqq");          // QQQ
 *
 ***************************************************/

#ifndef ScanDedxScale_C
#define ScanDedxScale_C

#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <vector>

// ============================================================
// scales       : the DEDX_SCALE values from your bash scan -- MUST
//                match your folder names exactly (same decimal
//                formatting)
// folderPrefix : e.g. "Output_27Al_"
// fileName     : ROOT file name inside each folder, e.g. "output_27Al.root"
// histName     : the BeamEnergy_ETrack_vs_EKin histogram. Swap "sx3"
//                for "qqq" (and adjust the gate tag if needed) to run
//                the QQQ side with the exact same macro.
// fitMin,fitMax: x-range (beam_energy_at_vertex, in MeV) to fit the
//                profile over. Histogram range is 0 to beamE0*1.5
//                (~84 MeV) with 400 bins, but real events cluster in
//                a much narrower band -- check your printed
//                "beamE@vertex" values from earlier scans (roughly
//                6-27 MeV across DEDX_SCALE 0.70-1.15) before
//                trusting these defaults; widen/narrow as needed.
// minEntriesPerBin : profile bins with fewer raw entries than this
//                are dropped from the line fit (noisy tails can
//                otherwise pull the slope around)
// ============================================================
void ScanDedxScale(
    std::vector<double> scales,
    TString              folderPrefix     = "Output_27Al_",
    TString              fileName         = "output_27Al.root",
    TString              histName         = "_m27Alax+misc_sx3_p/m27Alax_BeamEnergy_ETrack_vs_EKin_p_a2c0_sx3",
    double               fitMin           = 5.0,
    double               fitMax           = 25.0,
    int                  minEntriesPerBin = 20
){
  int nScales = (int) scales.size();
  if (nScales == 0) { printf("ScanDedxGasCalib: no scale values given.\n"); return; }

  TGraphErrors *gSlope = new TGraphErrors();
  gSlope->SetTitle("Slope (ebeam_kin/beam_energy_at_vertex) vs DEDX_SCALE;DEDX_SCALE;slope");
  gSlope->SetMarkerStyle(20);

  TGraphErrors *gIntercept = new TGraphErrors();
  gIntercept->SetTitle("Intercept vs DEDX_SCALE;DEDX_SCALE;intercept [MeV]");
  gIntercept->SetMarkerStyle(20);
  gIntercept->SetMarkerColor(kRed+1);
  gIntercept->SetLineColor(kRed+1);

  printf("\n%-8s %12s %12s %12s %10s\n", "scale", "slope", "+/-", "intercept", "chi2/ndf");

  for (int s = 0; s < nScales; s++) {
    TString folder = Form("%s%.2f", folderPrefix.Data(), scales[s]);
    TString path   = folder + "/" + fileName;

    TFile *f = TFile::Open(path, "READ");
    if (!f || f->IsZombie()) {
      printf("%-8.2f  ERROR: could not open %s\n", scales[s], path.Data());
      continue;
    }

    TH2 *h2 = (TH2*) f->Get(histName);
    if (!h2) {
      printf("%-8.2f  ERROR: histogram '%s' not found\n", scales[s], histName.Data());
      f->ls();
      continue;
    }

    // Profile: mean ebeam_kin_MeV (y) in bins of beam_energy_at_vertex (x)
    TProfile *profile = h2->ProfileX(Form("profile_%.2f", scales[s]));

    // Build a filtered graph: only bins inside [fitMin,fitMax] with at
    // least minEntriesPerBin raw entries, so sparsely populated bins
    // (noisy tails) don't pull the line fit around.
    TGraphErrors *gFit = new TGraphErrors();
    int nb = profile->GetNbinsX();
    for (int b = 1; b <= nb; b++) {
      double xc = profile->GetBinCenter(b);
      if (xc < fitMin || xc > fitMax) continue;
      if (profile->GetBinEntries(b) < minEntriesPerBin) continue;
      int gi = gFit->GetN();
      gFit->SetPoint(gi, xc, profile->GetBinContent(b));
      gFit->SetPointError(gi, 0, profile->GetBinError(b));
    }

    if (gFit->GetN() < 2) {
      printf("%-8.2f  ERROR: fewer than 2 usable bins after entry-count filtering "
             "(try lowering minEntriesPerBin or widening [fitMin,fitMax])\n", scales[s]);
      continue;
    }

    TF1 *line = new TF1(Form("line_%.2f", scales[s]), "pol1", fitMin, fitMax);
    gFit->Fit(line, "RQ"); // Q: quiet, don't spam per-scale fit output

    double slope    = line->GetParameter(1);
    double slopeErr = line->GetParError(1);
    double icept    = line->GetParameter(0);
    double icptErr  = line->GetParError(0);
    int ndf = line->GetNDF();
    double chi2 = line->GetChisquare();

    printf("%-8.2f %12.4f %12.4f %12.4f %10.3f\n",
           scales[s], slope, slopeErr, icept, ndf > 0 ? chi2/ndf : -1.0);

    int n = gSlope->GetN();
    gSlope->SetPoint(n, scales[s], slope);
    gSlope->SetPointError(n, 0, slopeErr);
    gIntercept->SetPoint(n, scales[s], icept);
    gIntercept->SetPointError(n, 0, icptErr);
  }

  TCanvas *c = new TCanvas("ScanDedxScale", "DEDX_SCALE via BeamEnergy_ETrack_vs_EKin", 1000, 700);
  c->Divide(1, 2);

  c->cd(1);
  gSlope->Draw("APL");
  TF1 *targetSlope = new TF1("targetSlope", "1", scales.front(), scales.back());
  targetSlope->SetLineColor(kGray+1);
  targetSlope->SetLineStyle(2);
  targetSlope->Draw("SAME");

  c->cd(2);
  gIntercept->Draw("APL");
  TF1 *targetIcept = new TF1("targetIcept", "0", scales.front(), scales.back());
  targetIcept->SetLineColor(kGray+1);
  targetIcept->SetLineStyle(2);
  targetIcept->Draw("SAME");

  printf("\nLook for where the slope curve (top) crosses the dashed slope=1\n"
         "line AND the intercept curve (bottom) crosses dashed intercept=0 --\n"
         "ideally near the same DEDX_SCALE. If they cross at different scales,\n"
         "a single multiplicative DEDX_SCALE may not fully capture the real\n"
         "correction (e.g. an additive offset might also be needed).\n\n");
}

#endif