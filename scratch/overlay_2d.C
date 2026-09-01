#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TF1.h"
#include <iostream>
#include <vector>

void overlay_2d(TString rootFile, TString histsCSV, TString labelsCSV, TString xAxisLabel, TString yAxisLabel) {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    TObjArray* histArr = histsCSV.Tokenize(",");
    TObjArray* labelArr = labelsCSV.Tokenize(",");
    
    // Standard ANASEN color sequence
    int colors[] = {kBlack, kRed+1, kAzure+2, kGreen+2, kMagenta+1};
    std::vector<TH2*> hists;
    
    TFile *f = TFile::Open(rootFile, "READ");
    if (!f || f->IsZombie()) return;

    for (int i = 0; i < histArr->GetEntriesFast(); i++) {
        TString hName = ((TObjString*)histArr->At(i))->GetString();
        TH2 *h = (TH2 *)f->Get(hName);
        if (h) {
            TH2 *clone = (TH2 *)h->Clone(Form("h_%d", i));
            clone->SetDirectory(0); // Detach from file
            hists.push_back(clone);
        } else {
            std::cerr << "Warning: Could not find " << hName << "\n";
        }
    }
    f->Close();

    if (hists.empty()) return;

    TCanvas *c = new TCanvas("c", "", 1800, 1800);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);
    
    // Turn on the X and Y grid lines
    c->SetGridx(1);
    c->SetGridy(1);

    TLegend *leg = new TLegend(0.65, 0.75, 0.9, 0.9);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001); // Solid background so grid doesn't bleed through
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.03);

    for (size_t i = 0; i < hists.size(); i++) {
        hists[i]->SetMarkerColor(colors[i % 5]);
        hists[i]->SetMarkerStyle(6); // Small dot for scatter
        
        if (i == 0) {
            hists[i]->GetXaxis()->SetTitle(xAxisLabel);
            hists[i]->GetYaxis()->SetTitle(yAxisLabel);
            hists[i]->GetXaxis()->CenterTitle();
            hists[i]->GetYaxis()->CenterTitle();
            hists[i]->Draw("scat");
            
            // Draw the y = x diagonal dashed line
            TF1 *diag = new TF1("diag", "x", -1000, 1000); 
            diag->SetLineStyle(2); 
            diag->SetLineColor(kGray+2);
            diag->Draw("same");
        } else {
            hists[i]->Draw("scat same");
        }
        
        TString label = ((TObjString*)labelArr->At(i))->GetString();
        leg->AddEntry(hists[i], label.Data(), "p");
    }
    leg->Draw();
    
    c->SaveAs("kinematic_states_overlay.png");
    std::cout << "Saved: kinematic_states_overlay.png\n";
}