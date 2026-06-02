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
void SetStyle() {
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

// ---------------------------------------------------------------------------
// make_pretty()
// ---------------------------------------------------------------------------
void make_pretty(const char* rootFile,
                 const char* histName,
                 const char* xlabel = "",
                 const char* ylabel = "") {

    SetStyle();

    // --- Open file ----------------------------------------------------------
    TFile* f = TFile::Open(rootFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << rootFile << "\n";
        return;
    }

    TObject* obj = f->Get(histName);
    if (!obj) {
        std::cerr << "ERROR: '" << histName << "' not found in " << rootFile << "\n";
        f->Close();
        return;
    }

    // --- Clone and detach BEFORE closing the file ---------------------------
    // ROOT ties every histogram to the TDirectory of the file it came from.
    // When f->Close() is called, ROOT deletes all objects in that directory —
    // including the clone — unless SetDirectory(0) is called first to detach
    // it. Drawing a deleted object causes the segfault you are seeing.
    TObject* clone = obj->Clone(Form("%s_clone", histName));
    if (!clone) {
        std::cerr << "ERROR: Clone failed for '" << histName << "'\n";
        f->Close();
        return;
    }
    if (clone->InheritsFrom(TH1::Class()))
        ((TH1*)clone)->SetDirectory(0);   // detach — survives file close
    f->Close();

    // --- Detect dimension ---------------------------------------------------
    bool is2D = clone->InheritsFrom(TH2::Class());

    // --- Canvas: 2100 px wide = 7 in at 300 DPI -----------------------------
    int canvasW = 2100;
    int canvasH = is2D ? 1800 : 1575;   // square-ish for 2D, 4:3 for 1D

    TCanvas c("c", "", 0, 0, canvasW, canvasH);
    c.cd();

    // Widen right margin for the colz palette bar
    if (is2D) gPad->SetRightMargin(0.13);

    // --- Draw ---------------------------------------------------------------
    if (is2D) {
        TH2* h = (TH2*)clone;
        h->SetStats(0);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->Draw("colz");
    } else {
        TH1* h = (TH1*)clone;
        h->SetStats(0);
        h->SetLineColor(kBlue+1);
        h->SetFillColorAlpha(kBlue+1, 0.25);
        h->SetFillStyle(1001);
        h->GetXaxis()->SetTitle("");
        h->GetYaxis()->SetTitle("");
        h->Draw("hist");
    }

    // --- Axis labels via TLatex (full typographic control) ------------------
    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.050);

    // X: y=0.06 sits just below the tick numbers inside the bottom margin (0.13).
    // Y: x=0.08 sits just left of the tick numbers inside the left margin (0.14).
    // Increase either value to push the label further from the axis.
    if (strlen(xlabel) > 0)
        tex.DrawLatex(0.46, 0.05, xlabel);

    if (strlen(ylabel) > 0) {
        tex.SetTextAngle(90);
        tex.DrawLatex(0.08, 0.40, ylabel);
        tex.SetTextAngle(0);
    }

    // --- Save ---------------------------------------------------------------
    // Build output path: same directory as input file, named after the histogram
    std::string inPath(rootFile);
    std::string dir = inPath.substr(0, inPath.find_last_of("/\\"));
    if (dir == inPath) dir = ".";   // no directory component — use cwd

    std::string outPath = dir + "/" + std::string(histName) + ".png";
    // Replace any "/" inside histName (e.g. "folder/hist") with "_"
    for (char& ch : outPath)
        if (ch == '/') ch = '_';

    c.Modified();
    c.Update();
    c.SaveAs(outPath.c_str());

    std::cout << "Saved: " << outPath << "\n";
}