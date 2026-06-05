#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <vector>


void scan_rf_timing(int startRun = 350, int endRun = 370) {
    TCanvas *c = new TCanvas("c1", "RF-MCP Timing Images", 0, 0, 1600, 800);
    c->Divide(2, 1);

    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+1, 
                    kOrange+7, kSpring+4, kViolet+2, kAzure+2, kTeal-1, kPink+2};
    int nColors = 12;

    std::vector<TFile*> files;
    int colorIdx = 0;

    for (int i = startRun; i <= endRun; i++) {
        TString filename = Form("../17F_output/results_run%d.root", i);
        TFile* f = TFile::Open(filename, "READ");
        
        if (!f || f->IsZombie()) {
            std::cout << "Skipping missing file: " << filename << std::endl;
            continue;
        }
        files.push_back(f); 

        int color = colors[colorIdx % nColors];
        colorIdx++;

        // --- Pad 1: Inner Ring  ---
        c->cd(1);
        c->GetPad(1)->SetGrid(1, 1);
        TH1F *h1 = (TH1F*)(f->Get("misc/dt_rf_mcp_qqq_innerring1"));
        if (h1) {
            h1->SetDirectory(0);
            h1->Rebin(2);
            h1->SetTitle(Form("Run %d: Inner Ring Timing", i));
            h1->SetLineColor(color);
            h1->SetLineWidth(3);
            h1->Draw("hist");
            
            c->Update();
            TLine *line1 = new TLine(105.189, gPad->GetUymin(), 105.189, gPad->GetUymax());
            line1->SetLineColor(kRed);
            line1->SetLineStyle(2); // Dashed
            line1->SetLineWidth(2);
            line1->Draw();
        }

        // --- Pad 2: Outer Ring  ---
        c->cd(2);
        c->GetPad(2)->SetGrid(1, 1);
        TH1F *h0 = (TH1F*)(f->Get("misc/dt_rf_mcp_qqq_innerring0"));
        if (h0) {
            h0->SetDirectory(0); 
            h0->Rebin(2); 
            h0->SetTitle(Form("Run %d: Outer Ring Timing", i));
            h0->SetLineColor(color);
            h0->SetLineWidth(3);
            h0->Draw("hist");
            
            c->Update();
            TLine *line0 = new TLine(105.189, gPad->GetUymin(), 105.189, gPad->GetUymax());
            line0->SetLineColor(kRed);
            line0->SetLineStyle(2); // Dashed
            line0->SetLineWidth(2);
            line0->Draw();
        }

        c->cd(1); c->Modified(); c->Update();
        c->cd(2); c->Modified(); c->Update();

        TString outName = Form("rf_timing_run%03d.png", i);
        c->SaveAs(outName);
        std::cout << "Saved: " << outName << std::endl;
    }

    std::cout << "\n=== IMAGE GENERATION COMPLETE ===" << std::endl;

    for (auto file : files) {
        if (file) file->Close();
    }
}