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
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

// Optional dataset_filter (e.g., "17F" or "27Al") prevents mixing different
// gas pressure/temperature environments which causes gain smearing.
// Leave blank ("") to pool all files.
void fit_pc_energy_calibration(const std::string& dataset_filter = "")
{
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
        }
      }
    }
    ++nFiles;
  }
  
  std::cout << "fit_pc_energy_calibration: read " << nFiles 
            << " run file(s) from pc_calib_raw/ (Filter: '" << dataset_filter << "')" << std::endl;
  std::cout << "fit_pc_energy_calibration: cut " << nOverflowCut 
            << " points due to ADC >= 64k overflow." << std::endl;

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

  std::ofstream outfile("pc_energy_calibration.dat");
  outfile << std::scientific << std::setprecision(6);
  
  // --- Fit and Plot each wire ---
  for (int wire = 0; wire < 48; ++wire)
  {
    double slope = 1.0, intercept = 0.0;
    double n = static_cast<double>(pts[wire].size());
    bool ok = (n >= 2);
    
    // 1. Calculate Standard Linear Least Squares Fit (y = mx) where intercept is forced to 0
    double maxX = 0.0, maxY = 0.0;
    double sxx = 0, sxy = 0; // Only need sum(x^2) and sum(x*y) for fixed-0 intercept

    for (const auto &p : pts[wire])
    {
      if (p.first > maxX) maxX = p.first;
      if (p.second > maxY) maxY = p.second;
      sxx += p.first * p.first;
      sxy += p.first * p.second;
    }

    if (ok)
    {
      if (std::isfinite(sxx) && std::abs(sxx) > 1e-12)
      {
        slope = sxy / sxx;
        intercept = 0.0; // Forced mathematically
      }
      else
      {
        ok = false;
      }
    }

    if (!ok) {
      std::cerr << "fit_pc_energy_calibration: wire " << wire << " has too few points (" << n
                << ") to fit -- writing identity (slope=1, intercept=0)" << std::endl;
    }

    outfile << wire << " " << slope << " " << intercept << "\n";

    // 2. Generate and Fill 2D Density Histogram
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
    TF1 *fitLine = nullptr;
    TPaveText *pt = nullptr;
    
    if (n > 0) {
        fitLine = new TF1(Form("fit_%d", wire), "[0]*x", 0, maxX * 1.05); // Formula is strictly y = m*x
        fitLine->SetParameter(0, slope);
        fitLine->SetLineColor(kRed);
        fitLine->SetLineWidth(2);

        pt = new TPaveText(0.15, 0.75, 0.55, 0.88, "NDC");
        pt->SetFillColor(kWhite);
        pt->SetBorderSize(1);
        pt->AddText(Form("N = %.0f", n));
        pt->AddText(Form("m = %.2e", slope));
        pt->AddText("b = 0 (Fixed)");
    }

    // ------------------------------------------------------------
    // 4. Save High-Res Individual PNG
    // ------------------------------------------------------------
    TCanvas cTemp("cTemp", "cTemp", 800, 600);
    cTemp.SetGridx();
    cTemp.SetGridy();
    if (n > 0) {
        h2->Draw("COLZ");
        fitLine->Draw("SAME");
        pt->Draw();
    } else {
        h2->SetTitle(wName + " (DEAD/NO DATA)");
        h2->Draw();
    }
    cTemp.SaveAs(Form("pc_calib_plots/wire_%02d.png", wire));

    // ------------------------------------------------------------
    // 5. Draw onto Global PDF Canvas
    // ------------------------------------------------------------
    TVirtualPad* pad = (wire < 24) ? cAnodes->cd(wire + 1) : cCathodes->cd(wire - 24 + 1);
    pad->SetGridx();
    pad->SetGridy();

    if (n > 0) 
    {
      h2->DrawClone("COLZ");
      fitLine->DrawClone("SAME");
      pt->DrawClone();
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

  std::cout << "fit_pc_energy_calibration: wrote pc_energy_calibration.dat" << std::endl;
  std::cout << "fit_pc_energy_calibration: individual high-res PNGs saved to pc_calib_plots/" << std::endl;
  std::cout << "fit_pc_energy_calibration: global contact sheets saved to pc_calib_anodes.pdf and pc_calib_cathodes.pdf" << std::endl;
}