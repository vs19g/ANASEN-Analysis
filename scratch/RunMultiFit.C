/***************************************************
 *
 * RunMultiFit.C
 *
 * Two fitting functions for multi-Gaussian excitation-energy spectra:
 *
 *   RunMultiFit(...)          -- wrapper around AutoFit.C's fitNGaussPol():
 *     1) opens a ROOT file living anywhere on disk
 *     2) pulls out a histogram, including ones nested inside
 *        TDirectories within the file (e.g. "subdir/histname")
 *     3) auto-generates the AutoFit-style parameter text file
 *        from a plain list of peak positions
 *     4) runs the n-Gauss + polynomial-background fit
 *     5) prints peak-amplitude correlations (via PrintCorrelations)
 *
 *   RunMultiFitWithTail(...)  -- same idea, but builds its own TF1
 *     directly (bypassing fitNGaussPol) to add an extra background
 *     term -- an exponential or Fermi/sigmoid tail -- for modeling a
 *     fusion-evaporation-like continuum on the low-Ex side, which
 *     fitNGaussPol's gaus(i)+pol(deg)-only parameter file has no way
 *     to express.
 *
 * Only load THIS file -- it already #includes AutoFit.C,
 * so don't separately .L AutoFit.C or you'll double-load it.
 *
 * Usage (from the ROOT prompt, with both files in the same dir):
 *
 *   .L RunMultiFit.C+
 *
 *   std::vector<double> peaks = {0, 2.2, 3.4, 4.8, 5.6, 6.55};
 *
 *   RunMultiFit(
 *       "/full/path/to/your/folder/yourfile.root",
 *       "_m27Alax+misc_sx3_p/m27Alax_EX_from_p_a2c0_sx3",
 *       peaks
 *   );
 *
 *   // or, with a fusion-evaporation tail on the low-Ex background:
 *   RunMultiFitWithTail(
 *       "/full/path/to/your/folder/yourfile.root",
 *       "_m27Alax+misc_sx3_p/m27Alax_EX_from_p_a2c0_sx3",
 *       peaks
 *   );
 *
 ***************************************************/

#ifndef RunMultiFit_C
#define RunMultiFit_C

// AutoFit.C uses TH1F*, ifstream, TCanvas, gROOT, gStyle, TLatex, and TList
// throughout but only gets forward declarations (or nothing at all) for
// several of them via its own includes -- fine for cling's interpreter,
// not enough for ACLiC's real compiler. Pull in full definitions here,
// before AutoFit.C is included, so they're available everywhere in this
// translation unit.
#include <fstream>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TList.h>

#include "AutoFit.C"

#include <TFile.h>
#include <TVirtualFitter.h>
#include <fstream>
#include <vector>

// ============================================================
// PrintCorrelations: prints the correlation coefficient between every
// pair of peak amplitudes from the most recently run fit. Deliberately
// uses TVirtualFitter rather than the TFitResultPtr returned by
// h->Fit(f,"S") -- on some ROOT/cling builds, holding a TFitResultPtr
// at the interactive prompt crashes with a
// "cling::runtime::internal::LifetimeHandler" linker error. Reading
// straight from TVirtualFitter's covariance matrix gets the same
// numbers without touching that code path.
//
// f            : the fit TF1 (e.g. from h->GetFunction("fit"))
// nPeaks       : number of Gaussian peaks in the fit (not counting the
//                background parameters)
// warnThreshold: |correlation| at or above this gets flagged -- pairs
//                this correlated shouldn't be trusted individually;
//                treat their combined area as the meaningful number
// stride       : parameters per peak before the shared/background block
//                -- 3 (default) for independent amp/mean/sigma per peak,
//                2 for RunMultiFitWithTail's sharedSigma=true mode
//                (amp/mean per peak, sigma shared separately)
// ============================================================
void PrintCorrelations(TF1 *f, int nPeaks, double warnThreshold = 0.8, int stride = 3) {
  if (!f) {
    printf("PrintCorrelations: no fit function given (was the fit successful?).\n");
    return;
  }

  TVirtualFitter *fitter = TVirtualFitter::GetFitter();
  if (!fitter) {
    printf("PrintCorrelations: no active fitter found -- run the fit first.\n");
    return;
  }

  printf("\n================ Peak-amplitude correlations ================\n");
  printf("%-8s %-8s %10s\n", "peakA", "peakB", "corr");

  std::vector<TString> flagged;
  for (int i = 0; i < nPeaks; i++) {
    for (int j = i + 1; j < nPeaks; j++) {
      int pi = stride * i, pj = stride * j; // amplitude is the first of each peak's params
      double covij = fitter->GetCovarianceMatrixElement(pi, pj);
      double si = f->GetParError(pi);
      double sj = f->GetParError(pj);
      double corr = (si > 0 && sj > 0) ? covij / (si * sj) : 0;

      bool isHigh = TMath::Abs(corr) >= warnThreshold;
      printf("%-8d %-8d %10.3f%s\n", i, j, corr, isHigh ? "   <-- high" : "");
      if (isHigh) flagged.push_back(Form("peak %d <-> peak %d (corr = %.3f)", i, j, corr));
    }
  }

  if (!flagged.empty()) {
    printf("\n  %d pair(s) at |corr| >= %.2f -- their individual areas are not\n"
           "  reliably separable; treat the combined area as the trustworthy number:\n",
           (int) flagged.size(), warnThreshold);
    for (auto &s : flagged) printf("    %s\n", s.Data());
  }
  printf("===============================================================\n\n");
}

// ============================================================
// fileName      : full path to the .root file (the "folder" part
//                 of the path is just the directory it lives in --
//                 TFile::Open handles that directly, see below)
// histPath      : name of the histogram inside the file. If the
//                 histogram lives inside a TDirectory, include the
//                 directory the same way you would in TBrowser,
//                 e.g. "dirName/histName" (works for nested dirs too:
//                 "dir1/dir2/histName")
// peakEnergies  : your peak positions, e.g. {0,2.2,3.4,4.8,5.6,6.55}
// sigmaGuess    : default initial sigma (same units as x-axis) used
//                 for every peak unless overridden per-peak below
// meanWindow    : half-width allowed around each mean during the fit
//                 (used only for peaks not covered by meanWindowOverride).
//                   0    -> auto: half the distance to the nearest
//                           neighboring peak on each side (tighter
//                           windows automatically appear for closely
//                           spaced peaks, e.g. two peaks 0.8 apart
//                           each get +/-0.4 max)
//                   >0   -> fixed +/- meanWindow
// degPol        : degree of the polynomial background (0=flat,
//                 1=linear, 2=quadratic, ...)
// xMin, xMax    : fit range; leave both at 0 to use the full
//                 histogram range
// fixMeans      : true locks every mean at its input value instead
//                 of letting it float within its window
// fixSigmas     : true locks every sigma at its guess instead of
//                 letting it float up to that value
// paraFile      : where the auto-generated parameter file is written
//                 (plain text, human-readable/editable afterward)
// sigmaOverride : optional per-peak sigma guesses; if given, overrides
//                 sigmaGuess for peaks with a matching index
// printCorr     : if true (default), automatically prints the peak-
//                 amplitude correlation matrix after fitting, flagging
//                 any pair whose areas aren't reliably separable
// corrWarnThreshold : |correlation| at/above this gets flagged in the
//                 printout (see PrintCorrelations above)
// ============================================================
void RunMultiFit(
    TString              fileName,
    TString              histPath,
    std::vector<double>  peakEnergies,
    double               sigmaGuess        = 0.15,
    double               meanWindow        = 0,
    int                  degPol            = 1,
    double               xMin              = 0,
    double               xMax              = 0,
    bool                 fixMeans          = false,
    bool                 fixSigmas         = false,
    TString              paraFile          = "AutoFit_para_auto.txt",
    std::vector<double>  sigmaOverride     = {},
    bool                 printCorr         = true,
    double               corrWarnThreshold = 0.8,
    std::vector<double>  meanWindowOverride = {}
){
  // ---------- 1) open the file ----------
  // TFile::Open takes the full path, folder included, e.g.
  // "/home/user/data/run074/analysis.root" -- no separate step
  // needed to "open the folder" first.
  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    printf("ERROR: could not open file: %s\n", fileName.Data());
    return;
  }

  // ---------- 2) get the histogram ----------
  TH1F *h = (TH1F*) f->Get(histPath);
  if (!h) {
    printf("ERROR: histogram not found at path: %s\n", histPath.Data());
    printf("---- top-level contents of the file, for reference ----\n");
    f->ls();
    return;
  }
  h->SetDirectory(0); // detach from file so it's safe even after f closes

  // ---------- 3) build the AutoFit parameter file ----------
  int n = (int) peakEnergies.size();
  if (n == 0) {
    printf("ERROR: no peak energies given.\n");
    return;
  }

  std::ofstream out(paraFile.Data());
  out << "# energy    lowE       highE      eFlag   sigma    sFlag   height\n";

  for (int i = 0; i < n; i++) {
    double e = peakEnergies[i];

    double lo, hi;
    double thisWindow = (i < (int)meanWindowOverride.size()) ? meanWindowOverride[i] : meanWindow;
    if (thisWindow > 0) {
      lo = e - thisWindow;
      hi = e + thisWindow;
    } else {
      // half the gap to each neighbor; edge peaks reuse their
      // only neighbor's gap on both sides
      double dLeft  = (i == 0)      ? (peakEnergies[i+1] - e)     : (e - peakEnergies[i-1]);
      double dRight = (i == n - 1)  ? (e - peakEnergies[i-1])     : (peakEnergies[i+1] - e);
      lo = e - dLeft  / 2.0;
      hi = e + dRight / 2.0;
    }

    double sig = (i < (int)sigmaOverride.size()) ? sigmaOverride[i] : sigmaGuess;

    double guess = h->GetBinContent(h->FindBin(e));
    if (guess <= 0) guess = h->GetMaximum() * 0.05; // avoid a zero starting amplitude

    out << e   << "  "
        << lo  << "  "
        << hi  << "  "
        << (fixMeans  ? 1 : 0) << "  "
        << sig << "  "
        << (fixSigmas ? 1 : 0) << "  "
        << guess << "\n";
  }
  out.close();

  printf("Wrote %d-peak parameter file: %s\n", n, paraFile.Data());
  for (int i = 0; i < n; i++) {
    printf("  peak %d : E = %.4f\n", i, peakEnergies[i]);
  }

  // ---------- 4) run the fit ----------
  fitNGaussPol(h, degPol, paraFile, xMin, xMax);

  // ---------- 5) report peak-amplitude correlations ----------
  // fitNGaussPol names its internal TF1 "fit" and fits directly on h,
  // so ROOT attaches it to h's list of functions -- retrieve it from
  // there rather than needing fitNGaussPol to return anything.
  if (printCorr) {
    TF1 *fit = h->GetFunction("fit");
    PrintCorrelations(fit, n, corrWarnThreshold);
  }
}

// ============================================================
// RunMultiFitWithTail
//
// Fits N Gaussians + a polynomial background + an extra background
// term representing a fusion-evaporation-like continuum tail on the
// low-Ex side of the spectrum.
//
// AutoFit.C's fitNGaussPol has no way to express this extra term --
// its parameter file only knows gaus(i) + pol(deg) -- so this builds
// and fits its own TF1 directly instead of routing through it, then
// reuses PrintCorrelations above (which works with any TF1*, not
// tied to fitNGaussPol internals).
//
// tailShape:
//   "fermi" (default) : A / (1 + exp((x - x0)/d))
//       Smooth sigmoid: ~flat at A for low x, turns off toward 0
//       above x0. BOUNDED for all x -- numerically safe regardless
//       of fit range. x0 is interpretable as roughly where the
//       continuum "turns off"; d is the turn-off width.
//   "exp"              : A * exp(-x / tau)
//       Simpler, more directly matches "evaporation spectrum"
//       language, but can blow up for x well below 0 if tau is
//       small (e.g. tau=0.1 at x=-3 gives exp(30)). Only use this
//       if your fit range doesn't extend far below the peaks, or if
//       you've confirmed tau stays large enough to behave.
//
// fileName, histPath   : same as RunMultiFit
// peakEnergies         : peak means (fixed by default -- see fixMeans)
// sigmaGuess           : default initial sigma used for every peak,
//                         unless overridden per-peak by sigmaOverride
// degPol               : background polynomial degree
// xMin, xMax           : fit range (required here -- no "full range"
//                         default, since the tail shape depends on
//                         where the range actually starts)
// tailAmpGuess          : initial tail amplitude; -1 = auto (0.5x hist max)
// tailScaleGuess        : tau (exp) or d/width (fermi)
// tailX0Guess           : only used for "fermi" -- initial turn-off location
// fixMeans              : true (default) locks peak means at peakEnergies,
//                         consistent with how you've been running RunMultiFit
// fixSigmas             : true locks every sigma at its guess/override
//                         instead of letting it float (default false)
// sigmaOverride         : optional per-peak sigma guesses, e.g.
//                         {1.2, 1.0, 0.3, 0.3, 0.4, 0.4, 0.4} -- overrides
//                         sigmaGuess for peaks with a matching index. Use
//                         this the same way you've been using it in
//                         RunMultiFit to give wide peaks (like your 0 and
//                         2.2 MeV ones) more room without loosening
//                         everyone else's ceiling too. IGNORED if
//                         sharedSigma=true (see below).
// sharedSigma           : if true, every peak shares ONE sigma parameter
//                         instead of getting its own -- tests the
//                         hypothesis that peak width is genuinely
//                         constant across the spectrum, rather than
//                         assuming it. sigmaOverride is ignored in this
//                         mode (there's only one sigma to set, from
//                         sigmaGuess); fixSigmas still works to lock
//                         that single shared value. Compare chi2/ndf and
//                         the residuals against a free-sigma run: if
//                         shared-sigma fits comparably well, your
//                         instinct that width shouldn't vary was right;
//                         if it fits noticeably worse -- especially
//                         around peaks you know are blended multi-level
//                         clusters (4.8, 5.6, 6.55) -- that's evidence
//                         the width differences are physically real,
//                         not fit noise.
// ============================================================
void RunMultiFitWithTail(
    TString              fileName,
    TString              histPath,
    std::vector<double>  peakEnergies,
    double               sigmaGuess     = 0.5,
    int                  degPol         = 1,
    double               xMin           = -3,
    double               xMax           = 9,
    TString              tailShape      = "fermi",
    double               tailAmpGuess   = -1,
    double               tailScaleGuess = 1.0,
    double               tailX0Guess    = 0.0,
    bool                 fixMeans       = true,
    bool                 fixSigmas      = false,
    std::vector<double>  sigmaOverride  = {},
    bool                 sharedSigma    = false
){
  // ---------- open file + histogram ----------
  TFile *f = TFile::Open(fileName, "READ");
  if (!f || f->IsZombie()) {
    printf("ERROR: could not open file: %s\n", fileName.Data());
    return;
  }

  TH1 *h = (TH1*) f->Get(histPath);
  if (!h) {
    printf("ERROR: histogram not found at path: %s\n", histPath.Data());
    printf("---- top-level contents of the file, for reference ----\n");
    f->ls();
    return;
  }
  h->SetDirectory(0);

  int nPeaks = (int) peakEnergies.size();
  if (nPeaks == 0) { printf("ERROR: no peak energies given.\n"); return; }
  if (tailShape != "fermi" && tailShape != "exp") {
    printf("ERROR: tailShape must be \"fermi\" or \"exp\", got \"%s\"\n", tailShape.Data());
    return;
  }

  // ---------- index layout ----------
  // Normal mode: each peak gets 3 params [amp, mean, sigma] at 3*i.
  // Shared-sigma mode: each peak gets 2 params [amp, mean] at 2*i, and
  // ALL peaks reference the same single sigma parameter index, placed
  // right after the last peak's amp/mean pair.
  int stride = sharedSigma ? 2 : 3;
  auto ampIdx  = [&](int i){ return stride*i; };
  auto meanIdx = [&](int i){ return stride*i + 1; };
  int sharedSigmaIdx = stride*nPeaks; // only meaningful if sharedSigma
  auto sigIdx  = [&](int i){ return sharedSigma ? sharedSigmaIdx : (stride*i + 2); };

  if (sharedSigma && !sigmaOverride.empty())
    printf("NOTE: sharedSigma=true -- sigmaOverride is ignored; using sigmaGuess (%.3f) for the single shared sigma.\n", sigmaGuess);

  // ---------- build the formula: gaussians + poly + tail ----------
  TString formula;
  for (int i = 0; i < nPeaks; i++) {
    if (i > 0) formula += "+";
    formula += Form("[%d]*exp(-0.5*((x-[%d])/[%d])*((x-[%d])/[%d]))",
                     ampIdx(i), meanIdx(i), sigIdx(i), meanIdx(i), sigIdx(i));
  }
  int polOffset = sharedSigma ? (stride*nPeaks + 1) : (stride*nPeaks);
  formula += Form("+pol%d(%d)", degPol, polOffset);

  int tailOffset = polOffset + (degPol + 1);
  int nTailPar = (tailShape == "exp") ? 2 : 3;
  if (tailShape == "exp")
    formula += Form("+[%d]*exp(-x/[%d])", tailOffset, tailOffset+1);
  else
    formula += Form("+[%d]/(1+exp((x-[%d])/[%d]))", tailOffset, tailOffset+1, tailOffset+2);

  TF1 *fitT = new TF1("fitWithTail", formula, xMin, xMax);

  // ---------- initial guesses: peaks ----------
  for (int i = 0; i < nPeaks; i++) {
    double e = peakEnergies[i];
    double guess = h->GetBinContent(h->FindBin(e));
    if (guess <= 0) guess = h->GetMaximum() * 0.3;

    fitT->SetParameter(ampIdx(i), guess);
    fitT->SetParLimits(ampIdx(i), 0, h->GetMaximum() * 2);

    fitT->SetParameter(meanIdx(i), e);
    if (fixMeans) fitT->FixParameter(meanIdx(i), e);
    else          fitT->SetParLimits(meanIdx(i), e - 0.2, e + 0.2);

    if (!sharedSigma) {
      double sig = (i < (int)sigmaOverride.size()) ? sigmaOverride[i] : sigmaGuess;
      fitT->SetParameter(sigIdx(i), sig);
      if (fixSigmas) fitT->FixParameter(sigIdx(i), sig);
      else           fitT->SetParLimits(sigIdx(i), 0.02, 3.0);
    }
  }

  if (sharedSigma) {
    fitT->SetParameter(sharedSigmaIdx, sigmaGuess);
    if (fixSigmas) fitT->FixParameter(sharedSigmaIdx, sigmaGuess);
    else           fitT->SetParLimits(sharedSigmaIdx, 0.02, 3.0);
  }

  // ---------- initial guesses: polynomial ----------
  for (int i = 0; i <= degPol; i++)
    fitT->SetParameter(polOffset + i, (i == 0) ? 1.0 : 0.0);

  // ---------- initial guesses: tail ----------
  double ampGuess = (tailAmpGuess > 0) ? tailAmpGuess : h->GetMaximum() * 0.5;
  fitT->SetParameter(tailOffset, ampGuess);
  fitT->SetParLimits(tailOffset, 0, h->GetMaximum() * 3);

  if (tailShape == "exp") {
    fitT->SetParameter(tailOffset+1, tailScaleGuess);
    fitT->SetParLimits(tailOffset+1, 0.3, 10);   // kept away from very small tau -- see header caution
  } else {
    fitT->SetParameter(tailOffset+1, tailX0Guess);
    fitT->SetParLimits(tailOffset+1, xMin, xMax);
    fitT->SetParameter(tailOffset+2, tailScaleGuess);
    fitT->SetParLimits(tailOffset+2, 0.02, 5);
  }

  // ---------- fit ----------
  // "R0": R restricts to the given range; 0 suppresses TH1::Fit's default
  // behavior of auto-drawing the fit function (in red) onto whatever pad
  // is currently active -- without this, that auto-draw can land on a
  // leftover pad from a previous call (e.g. if you're comparing a
  // free-sigma run against a sharedSigma run back-to-back) and show up
  // as a stray extra curve wherever gPad happened to be pointing.
  h->Fit(fitT, "R0");

  // unique run tag so re-running with different settings (e.g. comparing
  // sharedSigma true/false) gets its own canvas instead of colliding by
  // name with a previous run's
  TString runTag = Form("%s_%s", tailShape.Data(), sharedSigma ? "shared" : "free");

  // ---------- draw: two-pad canvas (fit + residual), matching the
  //            fitNGaussPol look you've been using throughout ----------
  TCanvas *c = new TCanvas("cFitWithTail_" + runTag, "Fit with fusion-evap tail: " + runTag, 900, 900);
  c->Divide(1, 2);

  // ---- top pad: histogram + total fit + individual peaks + background ----
  c->cd(1);
  gPad->SetPad(0, 0.3, 1, 1);
  h->SetStats(0);
  h->Draw();

  int colors[] = {kRed, kGreen+2, kBlue, kMagenta+1, kOrange+7, kCyan+2, kViolet, kSpring+4, kPink+1, kAzure+1};
  int nColors = 10;

  std::vector<TF1*> peakCurves(nPeaks);
  for (int i = 0; i < nPeaks; i++) {
    peakCurves[i] = new TF1(Form("peak%d_%s", i, runTag.Data()), "gaus", xMin, xMax);
    peakCurves[i]->SetParameter(0, fitT->GetParameter(ampIdx(i)));
    peakCurves[i]->SetParameter(1, fitT->GetParameter(meanIdx(i)));
    peakCurves[i]->SetParameter(2, fitT->GetParameter(sigIdx(i)));
    peakCurves[i]->SetLineColor(colors[i % nColors]);
    peakCurves[i]->SetLineWidth(2);
    peakCurves[i]->Draw("SAME");
  }

  // combined background (polynomial + tail), so its shape under the peaks is visible
  TString bgFormula = Form("pol%d(0)", degPol);
  int bgTailOffset = degPol + 1;
  if (tailShape == "exp")
    bgFormula += Form("+[%d]*exp(-x/[%d])", bgTailOffset, bgTailOffset+1);
  else
    bgFormula += Form("+[%d]/(1+exp((x-[%d])/[%d]))", bgTailOffset, bgTailOffset+1, bgTailOffset+2);
  TF1 *bgCurve = new TF1("bgCurve_" + runTag, bgFormula, xMin, xMax);
  for (int i = 0; i <= degPol; i++) bgCurve->SetParameter(i, fitT->GetParameter(polOffset+i));
  for (int i = 0; i < nTailPar; i++) bgCurve->SetParameter(bgTailOffset+i, fitT->GetParameter(tailOffset+i));
  bgCurve->SetLineColor(kGray+2);
  bgCurve->SetLineStyle(2);
  bgCurve->SetLineWidth(2);
  bgCurve->Draw("SAME");

  fitT->SetLineColor(kBlack);
  fitT->SetLineWidth(3);
  fitT->Draw("SAME");

  double chi2 = fitT->GetChisquare();
  int ndf = fitT->GetNDF();
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.06);
  latex.DrawLatex(0.15, 0.85, Form("#bar{#chi}^{2} : %.3f", ndf > 0 ? chi2/ndf : chi2));

  // ---- bottom pad: residual (Hist - fit), same convention as fitNGaussPol ----
  c->cd(2);
  gPad->SetPad(0, 0, 1, 0.3);
  TH1 *hRes = (TH1*) h->Clone("hRes_" + runTag);
  hRes->Add(fitT, -1);
  hRes->SetTitle("Residual");
  hRes->GetYaxis()->SetTitle("Hist - fit");
  hRes->SetStats(0);
  hRes->Draw();

  // ---------- report ----------
  printf("\n---- tail parameters (%s) ----\n", tailShape.Data());
  if (tailShape == "exp") {
    printf("  amplitude = %.4f +/- %.4f\n", fitT->GetParameter(tailOffset),   fitT->GetParError(tailOffset));
    printf("  tau       = %.4f +/- %.4f\n", fitT->GetParameter(tailOffset+1), fitT->GetParError(tailOffset+1));
  } else {
    printf("  amplitude     = %.4f +/- %.4f\n", fitT->GetParameter(tailOffset),   fitT->GetParError(tailOffset));
    printf("  x0 (turn-off) = %.4f +/- %.4f\n", fitT->GetParameter(tailOffset+1), fitT->GetParError(tailOffset+1));
    printf("  d  (width)    = %.4f +/- %.4f\n", fitT->GetParameter(tailOffset+2), fitT->GetParError(tailOffset+2));
  }
  printf("chi2/ndf = %.3f / %d = %.3f\n\n",
         fitT->GetChisquare(), fitT->GetNDF(), fitT->GetChisquare()/fitT->GetNDF());

  if (sharedSigma) {
    printf("shared sigma (all peaks) = %.4f +/- %.4f\n\n",
           fitT->GetParameter(sharedSigmaIdx), fitT->GetParError(sharedSigmaIdx));
  }

  PrintCorrelations(fitT, nPeaks, 0.8, stride);
}

#endif