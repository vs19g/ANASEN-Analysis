#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TF1.h>
#include <TAxis.h>
#include <TFile.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TLine.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>

// Optional dataset_filter (e.g., "17F" or "27Al") prevents mixing different
// gas pressure/temperature environments which causes gain smearing.
// Leave blank ("") to pool all files.
//
// known_dead_wires: wire indices (0-23 anode, 24-47 cathode) known to be
// non-functional regardless of what their calibration data looks like. Whatever
// slope their own points would produce -- least-squares or robust -- isn't a real
// calibration and is excluded here rather than trusted. Update this list as more
// dead channels are identified; currently anode 9 and 12.
void fit_pc_energy_calibration(const std::string& dataset_filter = "",
                                const std::vector<int>& known_dead_wires = {9, 12})
{
  std::vector<bool> isDead(48, false);
  for (int w : known_dead_wires)
    if (w >= 0 && w < 48)
      isDead[w] = true;

  std::vector<std::pair<double, double>> pts[48]; // [wire] -> (ADC, dE_gas MeV)

  TSystemDirectory dir("pc_calib_raw", "pc_calib_raw");
  TList *files = dir.GetListOfFiles();
  if (!files)
  {
    std::cerr << "fit_pc_energy_calibration: pc_calib_raw/ not found or empty -- "
              << "run TrackRecon.C with PC energy calibration enabled first." << std::endl;
    return;
  }

  int nFiles = 0;
  long long nOverflowCut = 0;
  long long nOverflowCutPerWire[48] = {0};
  TIter next(files);
  TSystemFile *f;
  
  // Read all points from all files into the single pooled array
  while ((f = (TSystemFile *)next()))
  {
    TString name = f->GetName();
    if (f->IsDirectory() || !name.BeginsWith("points_") || !name.EndsWith(".dat"))
      continue;
      
    // Apply dataset safeguard if requested
    if (!dataset_filter.empty() && !name.Contains(dataset_filter.c_str()))
      continue;

    std::ifstream infile(std::string("pc_calib_raw/") + name.Data());
    if (!infile.is_open())
      continue;
      
    int wire;
    double adc, dE_gas;
    while (infile >> wire >> adc >> dE_gas)
    {
      if (wire >= 0 && wire < 48)
      {
        if (adc < 64000.0) {
          pts[wire].push_back({adc, dE_gas});
        } else {
          nOverflowCut++;
          nOverflowCutPerWire[wire]++;
        }
      }
    }
    ++nFiles;
  }
  
  std::cout << "fit_pc_energy_calibration: read " << nFiles 
            << " run file(s) from pc_calib_raw/ (Filter: '" << dataset_filter << "')" << std::endl;
  std::cout << "fit_pc_energy_calibration: cut " << nOverflowCut 
            << " points due to ADC >= 64k overflow." << std::endl;
  std::cout << "fit_pc_energy_calibration: per-wire overflow breakdown (wire: cut / kept):" << std::endl;
  for (int wire = 0; wire < 48; ++wire)
  {
    long long kept = static_cast<long long>(pts[wire].size());
    if (nOverflowCutPerWire[wire] == 0 && kept == 0)
      continue;
    double fracCut = (nOverflowCutPerWire[wire] + kept > 0)
                          ? 100.0 * nOverflowCutPerWire[wire] / (nOverflowCutPerWire[wire] + kept)
                          : 0.0;
    std::cout << "  " << (wire < 24 ? "anode " : "cathode ") << (wire < 24 ? wire : wire - 24)
              << ": " << nOverflowCutPerWire[wire] << " cut / " << kept << " kept ("
              << fracCut << "% cut)" << std::endl;
  }

  // --- Adaptive lower-bound cut (anode only) -----------------------------
  // Raw ADC per wire can be genuinely bimodal (pileup, accidental coincidences,
  // wrong-topology leakage) with a lower population that would otherwise bias
  // both the least-squares slope and the median-based peak-matching fallback.
  // Rather than one fixed ADC cutoff -- which can't be right for every wire at
  // once when gains differ (the valley between the two populations sits at a
  // genuinely different absolute ADC per wire) -- this finds each wire's own
  // dominant peak from a coarse histogram of its own points, preferring the
  // highest-ADC prominent peak over a lower one, and cuts at the valley just
  // below it. Points below that valley are dropped before any fit/median runs.
  const double kPeakProminenceFrac = 0.25; // a local max must reach this fraction
                                            // of the tallest bin to count as a peak
  const int kFloorHistBins = 60;
  auto findAdaptiveFloor = [&](const std::vector<double> &adcVals) -> double
  {
    if (adcVals.size() < 10)
      return 0.0; // too few points to say anything about shape -- don't cut
    double maxADC = *std::max_element(adcVals.begin(), adcVals.end());
    if (maxADC <= 0.0)
      return 0.0;
    std::vector<int> hist(kFloorHistBins, 0);
    double binW = maxADC / kFloorHistBins;
    for (double v : adcVals)
    {
      int b = std::min(kFloorHistBins - 1, static_cast<int>(v / binW));
      if (b >= 0) hist[b]++;
    }
    int tallest = *std::max_element(hist.begin(), hist.end());
    if (tallest <= 0)
      return 0.0;

    // Scan from the highest-ADC bin down; the first local maximum that's
    // prominent enough is taken as "the" peak.
    int peakBin = -1;
    for (int b = kFloorHistBins - 1; b >= 1; --b)
    {
      if (hist[b] < kPeakProminenceFrac * tallest)
        continue;
      bool isLocalMax = hist[b] >= hist[b - 1] && (b == kFloorHistBins - 1 || hist[b] >= hist[b + 1]);
      if (isLocalMax)
      {
        peakBin = b;
        break;
      }
    }
    if (peakBin <= 0)
      return 0.0; // no clear peak below the top edge -- nothing to cut against

    // Walk down from the peak to the valley: the floor is the bin where the
    // count stops falling and starts rising again (the start of a lower
    // population), or the histogram runs out.
    int valleyBin = peakBin;
    for (int b = peakBin - 1; b >= 0; --b)
    {
      if (hist[b] > hist[valleyBin])
        break;
      valleyBin = b;
    }
    return valleyBin * binW;
  };

  std::vector<std::pair<double, double>> ptsFit[48]; // math uses this; pts[] stays raw for display
  for (int wire = 0; wire < 48; ++wire)
    ptsFit[wire] = pts[wire];

  std::vector<double> floorADC(48, 0.0);
  for (int wire = 0; wire < 24; ++wire) // anode only
  {
    std::vector<double> adcVals;
    adcVals.reserve(pts[wire].size());
    for (const auto &p : pts[wire])
      adcVals.push_back(p.first);
    double floor = findAdaptiveFloor(adcVals);
    if (floor <= 0.0)
      continue;
    floorADC[wire] = floor;
    size_t before = pts[wire].size();
    std::vector<std::pair<double, double>> kept;
    kept.reserve(pts[wire].size());
    for (const auto &p : pts[wire])
      if (p.first >= floor)
        kept.push_back(p);
    ptsFit[wire] = std::move(kept);
    std::cout << "fit_pc_energy_calibration: anode wire " << wire << " adaptive floor = "
              << floor << " ADC -- kept " << ptsFit[wire].size() << "/" << before << " point(s)" << std::endl;
  }

  // --- Pass 1: per-wire stats and the standard least-squares fit ---
  // Split out from plotting so we can look at ALL wires' fitted slopes before
  // deciding whether any individual one needs the robust fallback below.
  std::vector<double> slope_lsq(48, 1.0), maxX_arr(48, 0.0), maxY_arr(48, 0.0), n_arr(48, 0.0);
  std::vector<bool> ok_arr(48, false);

  for (int wire = 0; wire < 48; ++wire)
  {
    double n = static_cast<double>(ptsFit[wire].size());
    n_arr[wire] = n;
    bool ok = (n >= 2);

    double maxX = 0.0, maxY = 0.0; // display bounds: from the raw (unfiltered) data
    for (const auto &p : pts[wire])
    {
      if (p.first > maxX) maxX = p.first;
      if (p.second > maxY) maxY = p.second;
    }
    maxX_arr[wire] = maxX;
    maxY_arr[wire] = maxY;

    double sxx = 0, sxy = 0; // fit sums: from the filtered data
    for (const auto &p : ptsFit[wire])
    {
      sxx += p.first * p.first;
      sxy += p.first * p.second;
    }

    if (ok)
    {
      if (std::isfinite(sxx) && std::abs(sxx) > 1e-12)
        slope_lsq[wire] = sxy / sxx;
      else
        ok = false;
    }
    ok_arr[wire] = ok;
  }

  auto medianOf = [](std::vector<double> v) -> double
  {
    if (v.empty())
      return 0.0;
    std::sort(v.begin(), v.end());
    size_t n = v.size();
    return (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
  };

  // Median least-squares slope on each side (anode vs cathode gain structures differ,
  // so they're compared separately, not pooled) -- the reference point for flagging
  // any single wire's fit as an outlier.
  std::vector<double> anodeGoodSlopes, cathodeGoodSlopes;
  for (int wire = 0; wire < 48; ++wire)
  {
    if (!ok_arr[wire] || isDead[wire])
      continue;
    (wire < 24 ? anodeGoodSlopes : cathodeGoodSlopes).push_back(slope_lsq[wire]);
  }
  double medianAnodeSlope = medianOf(anodeGoodSlopes);
  double medianCathodeSlope = medianOf(cathodeGoodSlopes);

  // A least-squares slope more than this factor away from its side's median (or a
  // wire with too few points to fit at all) gets a peak-matched slope instead: that
  // wire's own ADC peak lined up against the A1C2 consensus dE_gas peak pooled from
  // every trusted wire on its side (computed below). Flagged wires are marked
  // method=2 in the output file and "PEAK-MATCHED (a1c2 consensus)" in their
  // diagnostic plot -- check those plots rather than trusting this blindly, since a
  // genuinely multi-modal raw distribution needs a different fix, not a different average.
  const double kOutlierRatioHi = 1.75;
  const double kOutlierRatioLo = 1.0 / kOutlierRatioHi;

  std::cout << "fit_pc_energy_calibration: median anode slope = " << medianAnodeSlope
            << ", median cathode slope = " << medianCathodeSlope << std::endl;

  // A1C2 consensus dE_gas peak: the target quantity (Ee - Ex from source-to-SX3
  // geometry) doesn't depend on which wire recorded the ADC -- it's the same
  // physical prediction for every event. Pooling it across every wire that DID get
  // a trusted direct fit ("calibrated using the a1c2 method") gives a far more
  // precise reference than any single low-statistics wire's own handful of points
  // could, which is what the wires below actually get matched against.
  std::vector<bool> isOutlier(48, false);
  for (int wire = 0; wire < 48; ++wire)
  {
    if (!ok_arr[wire] || isDead[wire])
      continue;
    double medianSide = (wire < 24) ? medianAnodeSlope : medianCathodeSlope;
    double ratio = (medianSide > 0.0) ? (slope_lsq[wire] / medianSide) : 1.0;
    isOutlier[wire] = (medianSide > 0.0) && (ratio > kOutlierRatioHi || ratio < kOutlierRatioLo);
  }
  std::vector<double> anodeConsensusPts, cathodeConsensusPts;
  int nTrustedAnode = 0, nTrustedCathode = 0;
  for (int wire = 0; wire < 48; ++wire)
  {
    if (!ok_arr[wire] || isDead[wire] || isOutlier[wire])
      continue;
    auto &bucket = (wire < 24) ? anodeConsensusPts : cathodeConsensusPts;
    for (const auto &p : ptsFit[wire])
      bucket.push_back(p.second);
    (wire < 24 ? nTrustedAnode : nTrustedCathode)++;
  }
  double consensusAnodeDEgas = medianOf(anodeConsensusPts);
  double consensusCathodeDEgas = medianOf(cathodeConsensusPts);
  std::cout << "fit_pc_energy_calibration: A1C2 consensus anode dE_gas peak = "
            << consensusAnodeDEgas << " MeV from " << anodeConsensusPts.size()
            << " point(s) across " << nTrustedAnode << " trusted anode wire(s)" << std::endl;
  std::cout << "fit_pc_energy_calibration: consensus cathode dE_gas peak = "
            << consensusCathodeDEgas << " MeV from " << cathodeConsensusPts.size()
            << " point(s) across " << nTrustedCathode << " trusted cathode wire(s)" << std::endl;

  // --- Setup ROOT Diagnostic Graphics ---
  gROOT->SetBatch(kTRUE); // Run silently without popping up windows
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird); 
  
  // Create output directory for individual PNGs
  gSystem->mkdir("pc_calib_plots", kTRUE);
  
  TFile *fOut = new TFile("pc_calib_diagnostics.root", "RECREATE");
  
  TCanvas *cAnodes = new TCanvas("cAnodes", "Anode Calibrations", 1800, 1200);
  cAnodes->Divide(6, 4, 0.01, 0.01);
  
  TCanvas *cCathodes = new TCanvas("cCathodes", "Cathode Calibrations", 1800, 1200);
  cCathodes->Divide(6, 4, 0.01, 0.01);

  std::string outFilename = dataset_filter.empty()
                                ? "pc_energy_calibration.dat"
                                : "pc_energy_calibration_" + dataset_filter + ".dat";
  std::ofstream outfile(outFilename);
  outfile << std::scientific << std::setprecision(6);
  
  // --- Pass 2: finalize each wire's slope (lsq / robust fallback / identity / known-dead), plot, write ---
  for (int wire = 0; wire < 48; ++wire)
  {
    double n = n_arr[wire];
    bool ok = ok_arr[wire];
    double slope = 1.0, intercept = 0.0;
    int method = 0; // 0 = identity/no data, 1 = least-squares, 2 = peak-matched (a1c2 consensus), 3 = known dead

    if (isDead[wire])
    {
      method = 3;
      std::cerr << "fit_pc_energy_calibration: wire " << wire
                << " is in known_dead_wires -- writing identity regardless of what its "
                << n << " raw point(s) would otherwise fit to." << std::endl;
    }
    else if (ok && !isOutlier[wire])
    {
      slope = slope_lsq[wire];
      method = 1;
    }
    else if (n >= 1)
    {
      // Either the least-squares slope is a >1.75x outlier vs its side's median, or
      // there weren't even enough points (n<2) to attempt a direct fit at all. Either
      // way: take this wire's own ADC peak (median of whatever points it has) and line
      // it up with the A1C2-consensus dE_gas peak from the trusted wires on its side,
      // rather than trusting a noisy few-point regression or falling back to identity.
      std::vector<double> xs;
      xs.reserve(ptsFit[wire].size());
      for (const auto &p : ptsFit[wire])
        xs.push_back(p.first);
      double medX = medianOf(xs);
      double targetY = (wire < 24) ? consensusAnodeDEgas : consensusCathodeDEgas;

      if (medX > 1e-9 && targetY > 0.0)
      {
        slope = targetY / medX;
        method = 2;
        std::cerr << "fit_pc_energy_calibration: wire " << wire << " "
                  << (ok ? "least-squares slope (" + std::to_string(slope_lsq[wire]) + ") is an outlier vs its side's median"
                         : "has too few points (" + std::to_string(static_cast<int>(n)) + ") for a direct fit")
                  << " -- using peak-matched slope (own ADC median=" << medX
                  << " vs consensus dE_gas=" << targetY << ") = " << slope
                  << " instead. Check pc_calib_plots/wire_" << Form("%02d", wire)
                  << ".png to confirm this makes sense." << std::endl;
      }
      else if (ok)
      {
        // No usable consensus target (e.g. this side has no trusted wires at all) --
        // fall back to whatever least-squares gave, even though it was flagged.
        slope = slope_lsq[wire];
        method = 1;
      }
      // else: n>=1 but medX<=0 and no lsq fit either -- falls through to identity below.
    }

    if (method == 0 && !isDead[wire]) {
      std::cerr << "fit_pc_energy_calibration: wire " << wire << " has " << n
                << " point(s) and no usable consensus target -- writing identity (slope=1, intercept=0)" << std::endl;
    }

    outfile << wire << " " << slope << " " << intercept << " " << method << "\n";

    // 2. Generate and Fill 2D Density Histogram
    double maxX = maxX_arr[wire], maxY = maxY_arr[wire];
    if (maxX <= 0) maxX = 64000.0;
    if (maxY <= 0) maxY = 10.0;

    TString wName = Form("%s %02d", wire < 24 ? "Anode" : "Cathode", wire < 24 ? wire : wire - 24);
    TH2D *h2 = new TH2D(Form("h2_wire_%d", wire), 
                        wName + "; ADC; dE_{gas} (MeV)", 
                        150, 0, maxX * 1.05, 
                        150, 0, maxY * 1.05);

    for (const auto &p : pts[wire]) {
      h2->Fill(p.first, p.second);
    }

    // 3. Setup Fit Line and Stats Box
    // Known-dead wires never get a fit line drawn, even if they have raw points
    // (wire 9 does) -- whatever structure is in that data isn't trusted as real
    // calibration, but the raw scatter is still worth keeping for reference/crosstalk
    // diagnosis.
    TF1 *fitLine = nullptr;
    TPaveText *pt = nullptr;
    TLine *floorLine = nullptr;
    bool drawFit = (n > 0) && !isDead[wire];

    if (drawFit) {
        fitLine = new TF1(Form("fit_%d", wire), "[0]*x", 0, maxX * 1.05); // Formula is strictly y = m*x
        fitLine->SetParameter(0, slope);
        fitLine->SetLineColor(method == 2 ? kMagenta : kRed);
        fitLine->SetLineWidth(2);

        pt = new TPaveText(0.15, 0.72, 0.55, 0.88, "NDC");
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->AddText(Form("N = %.0f", n));
        pt->AddText(Form("m = %.2e", slope));
        pt->AddText("b = 0 (Fixed)");
        if (method == 2)
          pt->AddText("PEAK-MATCHED (a1c2 consensus)");
        if (floorADC[wire] > 0.0)
          pt->AddText(Form("floor cut @ %.0f ADC", floorADC[wire]));
    }
    // Adaptive floor marker: drawn whenever a cut was found, even for wires
    // that ended up on the identity/dead path, so it's visible on every plot
    // where it was computed, not just the ones that used it in a fit.
    if (wire < 24 && floorADC[wire] > 0.0)
    {
      floorLine = new TLine(floorADC[wire], 0, floorADC[wire], maxY * 1.05);
      floorLine->SetLineColor(kGreen + 2);
      floorLine->SetLineStyle(2);
      floorLine->SetLineWidth(2);
    }

    // ------------------------------------------------------------
    // 4. Save High-Res Individual PNG
    // ------------------------------------------------------------
    TCanvas cTemp("cTemp", "cTemp", 800, 600);
    cTemp.SetGridx();
    cTemp.SetGridy();
    if (drawFit) {
        h2->Draw("COLZ");
        fitLine->Draw("SAME");
        pt->Draw();
        if (floorLine) floorLine->Draw("SAME");
    } else if (isDead[wire]) {
        h2->SetTitle(wName + (n > 0 ? " (KNOWN DEAD -- excluded from fit)" : " (KNOWN DEAD, NO DATA)"));
        h2->Draw(n > 0 ? "COLZ" : "");
        if (floorLine) floorLine->Draw("SAME");
    } else {
        h2->SetTitle(wName + " (DEAD/NO DATA)");
        h2->Draw();
        if (floorLine) floorLine->Draw("SAME");
    }
    cTemp.SaveAs(Form("pc_calib_plots/wire_%02d.png", wire));

    // ------------------------------------------------------------
    // 5. Draw onto Global PDF Canvas
    // ------------------------------------------------------------
    TVirtualPad* pad = (wire < 24) ? cAnodes->cd(wire + 1) : cCathodes->cd(wire - 24 + 1);
    pad->SetGridx();
    pad->SetGridy();

    if (drawFit) 
    {
      h2->DrawClone("COLZ");
      fitLine->DrawClone("SAME");
      pt->DrawClone();
      if (floorLine) floorLine->DrawClone("SAME");
    }
    else if (n > 0)
    {
      h2->DrawClone("COLZ"); // known-dead wire with data: show the raw scatter, no fit
      if (floorLine) floorLine->DrawClone("SAME");
    }
    else 
    {
      h2->DrawClone(); 
    }
    
    // Write histogram to ROOT file
    fOut->cd();
    h2->Write();
  }
  
  outfile.close();
  
  fOut->Close();

  std::cout << "fit_pc_energy_calibration: wrote " << outFilename << std::endl;
  std::cout << "fit_pc_energy_calibration: individual high-res PNGs saved to pc_calib_plots/" << std::endl;
}