#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <vector>

void scan_rf_timing(int startRun = 350, int endRun = 370) {
    TCanvas *c = new TCanvas("c1", "RF-MCP Timing Animation", 0, 0, 1600, 800);
    c->Divide(2, 1);

    // Array of distinct colors to cycle through
    int colors[] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+1, 
                    kOrange+7, kSpring+4, kViolet+2, kAzure+2, kTeal-1, kPink+2};
    int nColors = 12;

    std::vector<TFile*> files;
    
    // Define the output file name
    const char* gifFile = "rf_timing_animation.gif";
    
    // Safety check: ROOT appends to GIFs. We must delete any old version 
    // before starting, otherwise it will just add new frames to the old file!
    gSystem->Unlink(gifFile); 

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

        // --- Pad 1: Inner Ring 1 ---
        c->cd(1);
        c->GetPad(1)->SetGrid(1, 1);
        TH1F *h1 = (TH1F*)(f->Get("misc/dt_rf_mcp_qqq_innerring1"));
        if (h1) {
            h1->SetDirectory(0);
            h1->Rebin(2);
            h1->SetTitle(Form("Run %d: Inner Ring 1 Timing", i));
            h1->SetLineColor(color);
            h1->SetLineWidth(3);
            h1->Draw("hist");
        }

        // --- Pad 2: Inner Ring 0 ---
        c->cd(2);
        c->GetPad(2)->SetGrid(1, 1);
        TH1F *h0 = (TH1F*)(f->Get("misc/dt_rf_mcp_qqq_innerring0"));
        if (h0) {
            h0->SetDirectory(0); 
            h0->Rebin(2); 
            h0->SetTitle(Form("Run %d: Inner Ring 0 Timing", i));
            h0->SetLineColor(color);
            h0->SetLineWidth(3);
            h0->Draw("hist");
        }

        // Update the canvas to reflect the new drawn histograms
        c->cd(1); c->Modified(); c->Update();
        c->cd(2); c->Modified(); c->Update();

        std::cout << "Adding Run " << i << " to GIF..." << std::endl;
        
        // --- ADD FRAME TO GIF ---
        // The "+50" tells ROOT to append this frame and wait 50 centiseconds (0.5 seconds)
        c->Print(Form("%s+50", gifFile));
    }

    // --- SEAL THE GIF ---
    // The "++" tells ROOT we are done adding frames and it should finalize the file
    c->Print(Form("%s++", gifFile));

    std::cout << "\n=== GIF GENERATION COMPLETE ===" << std::endl;
    std::cout << "Saved as: " << gifFile << std::endl;

    // Safely clean up memory
    for (auto file : files) {
        if (file) file->Close();
    }
}