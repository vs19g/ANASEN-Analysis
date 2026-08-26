// =============================================================================
// make_pretty.C
//
// Feed it a ROOT file and a histogram name, get a publication-quality PNG.
//
// Usage (Single File):
//   root -l -b -q 'make_pretty.C("myfile.root", "histName")'
//
// Usage (Multi-File Overlay):
//   root -l -b -q 'make_pretty.C("file1.root,file2.root", "Run 1,Run 2", "histName", "X-Axis", "Y-Axis")'
//
// =============================================================================

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <iostream>
#include <string>
#include <vector>

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

// ---------------------------------------------------------------------------
// Single Plot Function
// ---------------------------------------------------------------------------
void make_prettyplots(const char *rootFile,
                      const char *histName,
                      const char *xlabel = "",
                      const char *ylabel = "",
                      double minZ = -9999.0,
                      double minY = -9999.0,
                      double minX = -9999.0,
                      double maxX = -9999.0,
                      double maxY = -9999.0)
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
            h->GetXaxis()->SetRangeUser(minX, h->GetXaxis()->GetXmax());
            h->GetYaxis()->SetRangeUser(minY, h->GetYaxis()->GetXmax());
            h->SetMinimum(minZ);
        }
        if (minY != -9999.0)
            h->GetYaxis()->SetRangeUser(minY, h->GetYaxis()->GetXmax());
        if (minX != -9999.0)
            h->GetXaxis()->SetRangeUser(minX, h->GetXaxis()->GetXmax());

        if (maxX != -9999.0 && maxY != -9999.0)
        {
            h->GetXaxis()->SetRangeUser(minX, maxX);
            h->GetYaxis()->SetRangeUser(minY, maxY);
        }

        h->GetXaxis()->SetTitle(xlabel);
        h->GetYaxis()->SetTitle(ylabel);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.4);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->Draw("colz");
    }
    else
    {
        TH1 *h = (TH1 *)clone;
        if (minX != -9999.0)
            h->GetXaxis()->SetRangeUser(minX, h->GetXaxis()->GetXmax());
        if (maxX != -9999.0)
            h->GetXaxis()->SetRangeUser(minX, maxX);

        h->SetStats(0);
        h->SetLineColor(kAzure - 5);
        h->SetLineWidth(3);
        h->GetXaxis()->SetTitle(xlabel);
        h->GetYaxis()->SetTitle(ylabel);
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.4);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->Draw("hist");
    }

    // --- Save ---------------------------------------------------------------
    std::string inPath(rootFile);
    std::string dir = inPath.substr(0, inPath.find_last_of("/\\"));
    if (dir == inPath)
        dir = ".";

    std::string safeHistName = histName;
    for (char &ch : safeHistName)
        if (ch == '/')
            ch = '_';

    std::string outPath = dir + "/" + safeHistName + ".png";
    c.SaveAs(outPath.c_str());
    std::cout << "Saved: " << outPath << "\n";
}

// ---------------------------------------------------------------------------
// Multi-File Overlay Function
// ---------------------------------------------------------------------------
// NOTE: canvasW/canvasH added as trailing optional args (default = old
// hardcoded values, so every existing call site keeps working unchanged).
// Pass larger values for a bigger/higher-resolution output image.
void make_prettyplots(TString filesCSV, TString labelsCSV, TString histName, TString xAxisLabel="", TString yAxisLabel="", double yMin=-9999, double yMax=-9999, double xMin=-9999, double xMax=-9999, int canvasW=2100, int canvasH=1575)
{
    SetStyle();

    // Parse the comma-separated lists natively using ROOT's TString
    TObjArray* fileArr = filesCSV.Tokenize(",");
    TObjArray* labelArr = labelsCSV.Tokenize(",");

    if (fileArr->GetEntriesFast() != labelArr->GetEntriesFast())
    {
        std::cerr << "ERROR: Number of files does not match number of labels!\n";
        return;
    }

    // Array of distinct colors to cycle through
    int colors[] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7, kTeal - 1, kPink + 2};
    int nColors = 8;

    std::vector<TH1 *> hists;
    double globalMaxY = -999999.0;

    // Load all histograms and determine the maximum peak
    for (int i = 0; i < fileArr->GetEntriesFast(); i++)
    {
        TString currentFile = ((TObjString*)fileArr->At(i))->GetString();
        TString currentLabel = ((TObjString*)labelArr->At(i))->GetString();

        TFile *f = TFile::Open(currentFile, "READ");
        if (!f || f->IsZombie())
            continue;

        TH1 *h = (TH1 *)f->Get(histName);
        if (!h)
        {
            std::cerr << "WARNING: '" << histName << "' not found in " << currentFile << ". Skipping...\n";
            f->Close();
            continue;
        }

        TH1 *clone = (TH1 *)h->Clone(Form("h_%d", i));
        clone->SetDirectory(0); // Detach from file
        f->Close();

        // Apply X-range BEFORE checking the max Y-height, otherwise peaks
        // outside the viewing range might scale the Y-axis unnecessarily!
        if (xMin != -9999.0 && xMax != -9999.0)
            clone->GetXaxis()->SetRangeUser(xMin, xMax);
        else if (xMin != -9999.0)
            clone->GetXaxis()->SetRangeUser(xMin, clone->GetXaxis()->GetXmax());

        if (clone->GetMaximum() > globalMaxY)
        {
            globalMaxY = clone->GetMaximum();
        }

        hists.push_back(clone);
    }

    if (hists.empty())
    {
        std::cerr << "ERROR: No valid histograms found to plot.\n";
        return;
    }

    // --- Draw ---
    TCanvas c("c", "", 0, 0, canvasW, canvasH);
    c.cd();

    // Create Legend (Positioned top right)
    TLegend *leg = new TLegend(0.65, 0.70, 0.92, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // Transparent background

    for (size_t i = 0; i < hists.size(); i++)
    {
        hists[i]->SetStats(0);
        hists[i]->SetLineColor(colors[i % nColors]);
        hists[i]->SetLineWidth(3);

        if (i == 0)
        {
            // First histogram dictates the frame layout
            if (yMin != -9999.0 && yMax != -9999.0) {
                hists[i]->GetYaxis()->SetRangeUser(yMin, yMax);
            } else {
                hists[i]->SetMaximum(globalMaxY * 1.15); // Add 15% headroom automatically
            }

            if (xAxisLabel != "") hists[i]->GetXaxis()->SetTitle(xAxisLabel);
            if (yAxisLabel != "") hists[i]->GetYaxis()->SetTitle(yAxisLabel);
            
            hists[i]->GetXaxis()->SetTitleOffset(1.1);
            hists[i]->GetYaxis()->SetTitleOffset(1.4);
            hists[i]->GetXaxis()->CenterTitle(true);
            hists[i]->GetYaxis()->CenterTitle(true);
            hists[i]->Draw("hist");
        }
        else
        {
            // Subsequent histograms stack on top
            hists[i]->Draw("hist same");
        }

        TString currentLabel = ((TObjString*)labelArr->At(i))->GetString();
        leg->AddEntry(hists[i], currentLabel.Data(), "l");
    }

    leg->Draw();

    // --- Save ---
    TString safeHistName = histName;
    safeHistName.ReplaceAll("/", "_");

    TString outPath = safeHistName + "_overlay.png";
    c.SaveAs(outPath.Data());
    std::cout << "Saved: " << outPath.Data() << "\n";
}