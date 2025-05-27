#define GainMatchSX3_cxx

#include "GainMatchSX3.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <fstream>
#include <utility>
#include <algorithm>
#include <TProfile.h>
#include "Armory/ClassSX3.h"

#include "TVector3.h"

TH2F *hSX3FvsB;
TH2F *hSX3FvsB_g;
TH2F *hsx3IndexVE;
TH2F *hsx3IndexVE_g;
int padID = 0;

SX3 sx3_contr;
TCutG *cut;
std::map<std::tuple<int, int, int>, std::vector<std::tuple<double, double, double>>> dataPoints;

void GainMatchSX3::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hSX3FvsB_g = new TH2F("hSX3FvsB_g", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hsx3IndexVE_g = new TH2F("hsx3IndexVE_g", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);

    sx3_contr.ConstructGeo();

    // Load the TCutG object
    TFile *cutFile = TFile::Open("sx3cut.root");
    if (!cutFile || cutFile->IsZombie())
    {
        std::cerr << "Error: Could not open sx3cut.root" << std::endl;
        return;
    }
    cut = dynamic_cast<TCutG *>(cutFile->Get("sx3cut"));
    if (!cut)
    {
        std::cerr << "Error: Could not find TCutG named 'sx3cut' in sx3cut.root" << std::endl;
        return;
    }
    cut->SetName("sx3cut"); // Ensure the cut has the correct name
}

Bool_t GainMatchSX3::Process(Long64_t entry)
{

    b_sx3Multi->GetEntry(entry);
    b_sx3ID->GetEntry(entry);
    b_sx3Ch->GetEntry(entry);
    b_sx3E->GetEntry(entry);
    b_sx3T->GetEntry(entry);
    b_qqqMulti->GetEntry(entry);
    b_qqqID->GetEntry(entry);
    b_qqqCh->GetEntry(entry);
    b_qqqE->GetEntry(entry);
    b_qqqT->GetEntry(entry);
    b_pcMulti->GetEntry(entry);
    b_pcID->GetEntry(entry);
    b_pcCh->GetEntry(entry);
    b_pcE->GetEntry(entry);
    b_pcT->GetEntry(entry);

    sx3.CalIndex();
    qqq.CalIndex();
    pc.CalIndex();

    std::vector<std::pair<int, int>> ID;
    for (int i = 0; i < sx3.multi; i++)
    {
        if (sx3.e[i] > 100)
        {
            ID.push_back(std::pair<int, int>(sx3.id[i], i));
            hsx3IndexVE->Fill(sx3.index[i], sx3.e[i]);
        }
    }

    if (ID.size() > 0)
    {
        std::sort(ID.begin(), ID.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                  { return a.first < b.first; });

        std::vector<std::pair<int, int>> sx3ID;
        sx3ID.push_back(ID[0]);
        bool found = false;

        for (size_t i = 1; i < ID.size(); i++)
        {
            if (ID[i].first == sx3ID.back().first)
            {
                sx3ID.push_back(ID[i]);
                if (sx3ID.size() >= 3)
                {
                    found = true;
                }
            }
            else
            {
                if (!found)
                {
                    sx3ID.clear();
                    sx3ID.push_back(ID[i]);
                }
            }
        }

        if (found)
        {
            int sx3ChUp = -1, sx3ChDn = -1, sx3ChBk = -1;
            float sx3EUp = 0.0, sx3EDn = 0.0, sx3EBk = 0.0;

            for (size_t i = 0; i < sx3ID.size(); i++)
            {
                int index = sx3ID[i].second;
                if (sx3.ch[index] < 8)
                {
                    if (sx3.ch[index] % 2 == 0)
                    {
                        sx3ChDn = sx3.ch[index] / 2;
                        sx3EDn = sx3.e[index];
                    }
                    else
                    {
                        sx3ChUp = sx3.ch[index] / 2;
                        sx3EUp = sx3.e[index];
                    }
                }
                else 
                {
                    sx3ChBk = sx3.ch[index] - 8;
                    // if (sx3ChBk == 2)
                    //     printf("Found back channel Det %d Back %d \n", sx3.id[index], sx3ChBk);
                    sx3EBk = sx3.e[index];
                }
            }

            // Fill the histogram for the front vs back
            hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);

            for (int i = 0; i < sx3.multi; i++)
            {
                if (sx3.id[i] == 3)
                {
                    // Fill the histogram for the front vs back with gain correction
                    hSX3FvsB_g->Fill(sx3EUp + sx3EDn, sx3EBk);
                    // Fill the index vs energy histogram
                    hsx3IndexVE_g->Fill(sx3.index[i], sx3.e[i]);
                // }
                // {
                    TString histName = Form("hSX3FVB_id%d_U%d_D%d_B%d", sx3.id[i], sx3ChUp, sx3ChDn, sx3ChBk);
                    TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                    if (!hist2d)
                    {
                        hist2d = new TH2F(histName, Form("hSX3FVB_id%d_U%d_D%d_B%d", sx3.id[i], sx3ChUp, sx3ChDn, sx3ChBk), 400, 0, 16000, 400, 0, 16000);
                    }

                    // if (sx3ChBk == 2)
                    //     printf("Found back channel Det %d Back %d \n", sx3.id[i], sx3ChBk);
                    // hsx3IndexVE_g->Fill(sx3.index[i], sx3.e[i]);
                    // hSX3FvsB_g->Fill(sx3EUp + sx3EDn, sx3EBk);

                    hist2d->Fill(sx3EUp + sx3EDn, sx3EBk);
                    
                    // if (cut && cut->IsInside(sx3EUp + sx3EDn, sx3EBk))
                    // if (sx3.id[i] < 24 && sx3ChUp < 4 && sx3ChBk < 4 && std::isfinite(sx3EUp) && std::isfinite(sx3EDn) && std::isfinite(sx3EBk))
                    {
                        // Accumulate data for gain matching
                        dataPoints[{sx3.id[i], sx3ChUp, sx3ChBk}].emplace_back(sx3EBk, sx3EUp, sx3EDn);
                    }
                }
            }
        }
    }

    return kTRUE;
}

void GainMatchSX3::Terminate()
{
    const int MAX_DET = 24;
    const int MAX_UP = 4;
    const int MAX_DOWN = 4;
    const int MAX_BK = 4;

    double gainArray[MAX_DET][MAX_UP][MAX_BK] = {{{0}}};
    bool gainValid[MAX_DET][MAX_UP][MAX_BK] = {{{false}}};
    std::map<int, TH2F *> updn2DHistos;
    std::map<int, double> upCorrFactor;

    std::ofstream outFile("sx3_GainMatch.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Gain fit using up+dn vs bk
    for (const auto &kv : dataPoints)

    {
        auto [id, ud, bk] = kv.first;
        const auto &pts = kv.second;
        if (pts.size() < 5)
            continue;

        std::vector<double> bkE, udE;
        for (const auto &pr : pts)
        {
            double eUp, eDn, eBk;
            std::tie(eBk, eUp, eDn) = pr;
            bkE.push_back(eBk);
            udE.push_back(eUp + eDn);
        }

        TGraph g(bkE.size(), bkE.data(), udE.data());
        TF1 f("f", "[0]*x", 0, 16000);
        g.Fit(&f, "QNR");
        gainArray[id][ud][bk] = f.GetParameter(0);
        gainValid[id][ud][bk] = true;
    }

    // Output results
    for (int id = 0; id < MAX_DET; ++id)
    {
        for (int bk = 0; bk < MAX_BK; ++bk)
        {
            for (int ud = 0; ud < MAX_UP; ++ud)
            {
                if (gainValid[id][ud][bk])
                {
                    outFile << id << " " << bk << " " << ud << " " << gainArray[id][ud][bk] << std::endl;
                    printf("Gain match Det%d Up+Dn%d Back%d → %.4f \n", id, ud, bk, gainArray[id][ud][bk]);
                }
            }
        }
    }

    for (int bk = 0; bk < MAX_BK; ++bk)
    {
        TString name = Form("hUpDnVsBk_%d", bk);
        TString title = Form("Up/Bk vs Dn/Bk for Back %d;Dn/Bk;Up/Bk", bk);
        updn2DHistos[bk] = new TH2F(name, title, 400, 0, 1, 400, 0, 1);
    }

    outFile.close();
    std::cout << "Gain matching complete." << std::endl;

    // === Create histograms ===
    TH2F *hFVB = new TH2F("hFVB", "Corrected Up+Dn vs Corrected Back;Corrected Back E;Up+Dn E",
                          400, 0, 16000, 400, 0, 16000);
    TH2F *hAsym = new TH2F("hAsym", "Up vs Dn dvide corrected back;Up/Back E;Dn/Back E",
                           400, 0.0, 1.0, 400, 0.0, 1.0);

    // Fill histograms
    for (const auto &kv : dataPoints)
    {
        auto [id, ud, bk] = kv.first;
        if (!gainValid[id][ud][bk])
            continue;
        double gain = gainArray[id][ud][bk];

        for (const auto &pr : kv.second)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pr;
            
            double updn = eUp + eDn;
            if (updn == 0|| eBk == 0)
                continue;

            double asym = (eUp - eDn) / updn;
            double correctedBack = eBk * gain;

            hFVB->Fill(correctedBack, updn);
            hAsym->Fill(eUp / correctedBack, eDn / correctedBack);
            updn2DHistos[bk]->Fill(eUp / correctedBack, eDn / correctedBack);
            // hAsym->()
        }
    }

    for (auto &[bk, h2] : updn2DHistos)
    {
        // Project along diagonal
        TProfile *prof = h2->ProfileY(Form("prof_bk%d", bk), 1, h2->GetNbinsX(), "s");
        TF1 *fitLine = new TF1("fitLine", "[0]*x", 0, 1);
        fitLine->SetParameters(-1, 0.5); // initial guess: slope -1
        prof->Fit(fitLine, "QNR");

        double slope = fitLine->GetParameter(0);
        if (slope == 0)
        {
            slope = 1e-6; // prevent div-by-zero
        }

        // double corr = slope / -1.0;
        upCorrFactor[bk] = -1.0 / slope;

        printf("Back %d: Fit slope = %.4f → Up correction factor = %.4f\n", bk, slope, upCorrFactor[bk]);
    }

    TH2F *hFVB_Corr = new TH2F("hFVB_Corr", "Corrected Up+Dn vs Back;Back E;Corrected Up+Dn E",
                               400, 0, 16000, 400, 0, 16000);
    TH2F *hAsym_Corr = new TH2F("hAsym_Corr", "Corrected Up/Bk vs Dn/Bk;Dn/Bk;Up/Bk",
                                400, 0, 1.5, 400, 0, 1.5);

    for (const auto &kv : dataPoints)
    {
        auto [id, ud, bk] = kv.first;
        double factor = upCorrFactor.count(bk) ? upCorrFactor[bk] : 1.0;

        for (const auto &pr : kv.second)
        {
            double correctedBack, eUp, eDn;
            std::tie(correctedBack, eUp, eDn) = pr;

            double eUpCorr = eUp * factor;
            double eDnCorr = eDn;
            double eSumCorr = eUpCorr + eDnCorr;

            hFVB_Corr->Fill(correctedBack, eSumCorr);
            hAsym_Corr->Fill(eDnCorr / correctedBack, eUpCorr / correctedBack);
        }
    }
    // Optional: save histograms to a file
    // TFile *outHist = new TFile("sx3_gainmatch_hists.root", "RECREATE");
    // hFVB->Write();
    // hAsym->Write();
    // outHist->Close();
}
