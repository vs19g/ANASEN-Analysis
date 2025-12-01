#include "GainMatchQQQ.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <fstream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "Armory/HistPlotter.h"
#include "TVector3.h"

void makeplots()
{
    TCanvas *c1 = new TCanvas("c1", "c1", 0, 0, 1024, 768);
    c1->Divide(2, 2);

    TH2F *h1, *h2, *h3, *h4 = nullptr;

    std::string type = "Cal_cut";
    TFile *inFile = TFile::Open("Cal_checkQQQ.root");
    TF1 f1("x","x",0 ,80);


    for (int ring = 0; ring < 16; ring++)
    {
        for (int wedge = 0; wedge < 16; wedge++)
        {

            c1->cd(1);
            h1 = (TH2F *)inFile->Get(Form("E%s/h%s_qqq%d_ring%d_wedge%d", type.c_str(), type.c_str(), 0, ring, wedge));
            gPad->SetGrid();
            if (h1)
            {
                h1->SetTitle(Form("QQQ%d Ring %d Wedge %d %s E", 0, ring, wedge, type.c_str()));
                h1->Draw("COLZ");
                f1.DrawClone("same"); 
                
                gPad->Modified();
                gPad->Update();
            }
            
            c1->cd(2);
            h2 = (TH2F *)inFile->Get(Form("E%s/h%s_qqq%d_ring%d_wedge%d", type.c_str(), type.c_str(), 1, ring, wedge));
            gPad->SetGrid();
            if (h2)
            {
                std::cout<<h2<<std::endl;
                h2->SetTitle(Form("QQQ%d Ring %d Wedge %d %s E", 1, ring, wedge, type.c_str()));
                h2->Draw("COLZ");
                f1.DrawClone("same");
                gPad->Modified();
                gPad->Update();
            }
            
            c1->cd(3);
            h3 = (TH2F *)inFile->Get(Form("E%s/h%s_qqq%d_ring%d_wedge%d", type.c_str(), type.c_str(), 2, ring, wedge));
            gPad->SetGrid();
            if (h3)
            {
                h3->SetTitle(Form("QQQ%d Ring %d Wedge %d %s E", 2, ring, wedge, type.c_str()));
                h3->Draw("COLZ");
                f1.DrawClone("same");
                gPad->Modified();
                gPad->Update();
            }
            
            c1->cd(4);
            h4 = (TH2F *)inFile->Get(Form("E%s/h%s_qqq%d_ring%d_wedge%d", type.c_str(), type.c_str(), 3, ring, wedge));
            gPad->SetGrid();
            if (h4)
            {
                h4->SetTitle(Form("QQQ%d Ring %d Wedge %d %s E", 3, ring, wedge, type.c_str()));
                h4->Draw("COLZ");
                f1.DrawClone("same");
                gPad->Modified();
                gPad->Update();
            }
            while (gPad->WaitPrimitive());
        }
    }

    inFile->Close();
}