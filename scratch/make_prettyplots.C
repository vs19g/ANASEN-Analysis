// =============================================================================
// make_pretty.C
//
// Feed it a ROOT file and a histogram name, get a publication-quality PNG.
//
// Usage:
//   root -l -b -q 'make_pretty.C("myfile.root", "histName")'
//   root -l -b -q 'make_pretty.C("myfile.root", "histName", "x label", "y label")'
//
// Supports TH1F, TH1D, TH2F, TH2D — type is detected automatically.
// Output PNG is written to the same directory as the input file.
// =============================================================================

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TROOT.h"

#include <iostream>
#include <string>

// ---------------------------------------------------------------------------
// Style — called once before anything is drawn
// ---------------------------------------------------------------------------
void SetStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");

    gStyle->SetLabelSize(0.045, "xyz");
    gStyle->SetTitleSize(0.050, "xyz");
    gStyle->SetTitleOffset(1.15, "x");
    gStyle->SetTitleOffset(1.20, "y");

    gStyle->SetTickLength(0.025, "xy");
    gStyle->SetNdivisions(510, "xy");
    TGaxis::SetMaxDigits(4);

    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);

    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadTopMargin(0.06);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // kBird: perceptually uniform, grayscale-safe, colorblind-friendly
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(255);
}

void make_prettyplots(const char *rootFile,
                      const char *histName,
                      const char *xlabel = "",
                      const char *ylabel = "",
                      double minZ = -9999.0)
{

    SetStyle();

    // --- Open file ----------------------------------------------------------
    TFile *f = TFile::Open(rootFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "ERROR: Cannot open " << rootFile << "\n";
        return;
    }

    TObject *obj = f->Get(histName);
    if (!obj)
    {
        std::cerr << "ERROR: '" << histName << "' not found in " << rootFile << "\n";
        f->Close();
        return;
    }

    TObject *clone = obj->Clone(Form("%s_clone", histName));
    if (!clone)
    {
        std::cerr << "ERROR: Clone failed for '" << histName << "'\n";
        f->Close();
        return;
    }
    if (clone->InheritsFrom(TH1::Class()))
        ((TH1 *)clone)->SetDirectory(0); // detach — survives file close
    f->Close();

    // --- Detect dimension ---------------------------------------------------
    bool is2D = clone->InheritsFrom(TH2::Class());

    // --- Canvas: 2100 px wide = 7 in at 300 DPI -----------------------------
    int canvasW = 2100;
    int canvasH = is2D ? 1800 : 1575; // square-ish for 2D, 4:3 for 1D

    TCanvas c("c", "", 0, 0, canvasW, canvasH);
    c.cd();

    // Widen right margin for the colz palette bar
    if (is2D)
        gPad->SetRightMargin(0.13);

    if (is2D)
    {
        TH2 *h = (TH2 *)clone;
        h->SetStats(0);

        // --- Apply Minimum Z if provided ---
        if (minZ != -9999.0)
        {
            h->SetMinimum(minZ);
        }

        h->GetXaxis()->SetTitle(xlabel);
        h->GetYaxis()->SetTitle(ylabel);
        h->GetXaxis()->SetTitleOffset(1);
        h->GetYaxis()->SetTitleOffset(1.4);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->Draw("colz");
    }
    else
    {
        TH1 *h = (TH1 *)clone;
        h->SetStats(0);
        h->SetLineColor(kAzure - 5);
        h->SetLineWidth(3);
        h->GetXaxis()->SetTitle(xlabel);
        h->GetYaxis()->SetTitle(ylabel);
        h->GetXaxis()->SetTitleOffset(1);
        h->GetYaxis()->SetTitleOffset(1.4);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->Draw("hist");
    }

    // --- Save ---------------------------------------------------------------
    // Build output path: same directory as input file, named after the histogram
    std::string inPath(rootFile);
    std::string dir = inPath.substr(0, inPath.find_last_of("/\\"));
    if (dir == inPath)
        dir = "."; // no directory component — use cwd

    std::string outPath = dir + "/" + std::string(histName) + ".png";
    // Replace any "/" inside histName (e.g. "folder/hist") with "_"
    for (char &ch : outPath)
        if (ch == '/')
            ch = '_';

    c.Modified();
    c.Update();
    c.SaveAs(outPath.c_str());

    std::cout << "Saved: " << outPath << "\n";
}