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
#include "TGraphErrors.h"
#include "TMultiDimFit.h"

#include "TVector3.h"

TH2F *hSX3FvsB;
TH2F *hSX3FvsB_g;
TH2F *hsx3IndexVE;
TH2F *hsx3IndexVE_g;
TH2F *hSX3;
TH2F *hsx3Coin;

int padID = 0;

SX3 sx3_contr;
TCutG *cut;
TCutG *cut1;
std::map<std::tuple<int, int, int, int>, std::vector<std::tuple<double, double, double>>> dataPoints;


void GainMatchSX3::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hSX3FvsB_g = new TH2F("hSX3FvsB_g", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hsx3IndexVE_g = new TH2F("hsx3IndexVE_g", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hSX3 = new TH2F("hSX3", "SX3 Front v Back; Fronts; Backs", 8, 0, 8, 4, 0, 4);

    hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);

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

        // Load the TCutG object
    TFile *cutFile1 = TFile::Open("UvD.root");
    if (!cutFile1 || cutFile1->IsZombie())
    {
        std::cerr << "Error: Could not open UvD.root" << std::endl;
        return;
    }
    cut1 = dynamic_cast<TCutG *>(cutFile1->Get("UvD"));
    if (!cut1)
    {
        std::cerr << "Error: Could not find TCutG named 'UvD' in UvD.root" << std::endl;
        return;
    }
    cut1->SetName("UvD");
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

        for (int j = i + 1; j < sx3.multi; j++)
        {
            if (sx3.id[i] == 3)
                hsx3Coin->Fill(sx3.index[i], sx3.index[j]);
        }
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

        // start with the first entry in the sorted array: channels that belong to the same detector are together in sequenmce
        std::vector<std::pair<int, int>> sx3ID;
        sx3ID.push_back(ID[0]);
        bool found = false;

        for (size_t i = 1; i < ID.size(); i++)
        { // Check if id of i belongs to the same detector and then add it to the detector ID vector
            if (ID[i].first == sx3ID.back().first)
            { // count the nunmber of hits that belong to the same detector
                sx3ID.push_back(ID[i]);

                if (sx3ID.size() >= 3)
                {
                    found = true;
                }
            }
            else
            { // the next event does not belong to the same detector, abandon the first event and continue with the next one
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
                // Check the channel number and assign it to the appropriate channel type
                if (sx3.ch[index] < 8)
                {
                    if (sx3.ch[index] % 2 == 0)
                    {
                        sx3ChDn = sx3.ch[index];
                        sx3EDn = sx3.e[index];
                    }
                    else
                    {
                        sx3ChUp = sx3.ch[index];
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
            // If we have a valid front and back channel, fill the histograms
            hSX3->Fill(sx3ChDn, sx3ChBk);
            hSX3->Fill(sx3ChUp, sx3ChBk);

            // Fill the histogram for the front vs back
            hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);

            for (int i = 0; i < sx3.multi; i++)
            {
                if (sx3.id[i] == 3 && sx3.e[i] > 100)
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

                    if (cut && cut->IsInside(sx3EUp + sx3EDn, sx3EBk))
                    // if (sx3.id[i] < 24 && sx3ChUp < 4 && sx3ChBk < 4 && std::isfinite(sx3EUp) && std::isfinite(sx3EDn) && std::isfinite(sx3EBk))
                    {
                        // Accumulate data for gain matching
                        dataPoints[{sx3.id[i], sx3ChBk, sx3ChUp, sx3ChDn}].emplace_back(sx3EBk, sx3EUp, sx3EDn);
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

    double gainArray[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
    bool gainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

    std::ofstream outFile("sx3_MultiDimFit_results.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // === Loop over all (id, bk, up, dn) combinations ===
    for (const auto &kv : dataPoints)
    {
        auto [id, bk, u, d] = kv.first;
        const auto &pts = kv.second;

        if (pts.size() < 5)
            continue;

        int N = pts.size();
        double *x0 = new double[N];
        double *x1 = new double[N];
        double *x2 = new double[N];
        double *y  = new double[N];
        int mPowers[]   = {1,1,1};

        // Fill arrays
        for (int i = 0; i < N; ++i)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pts[i];
            x0[i] = eBk;
            x1[i] = eUp;
            x2[i] = eDn;
            y[i]  = eUp + eDn;  // Target is front sum
        }

        // Build MultiDim Fit
        TMultiDimFit *mdf = new TMultiDimFit(3, TMultiDimFit::kMonomials, "v");
        mdf->SetMaxPowers(mPowers); // Up to quadratic terms
        mdf->SetMaxTerms(3); // Limit number of terms kept

        // Add points
        for (int i = 0; i < N; ++i)
        {
            double vars[3] = {x0[i], x1[i], x2[i]};
            mdf->AddRow(vars, y[i]);
        }

        mdf->MakeHistograms();
        mdf->FindParameterization();
        mdf->Print("a"); // Print coefficients
        mdf->
        // Save result string
        TString formula;
        mdf->GetFunctions(formula);
        outFile << id << " " << bk << " " << u << " " << d << " " << formula.Data() << std::endl;
        printf("Det %d Bk %d Up %d Dn %d — MultiDimFit formula: %s\n", id, bk, u, d, formula.Data());

        gainValid[id][bk][u][d] = true; // Save as "valid" so we can use it

        // Clean up
        delete[] x0;
        delete[] x1;
        delete[] x2;
        delete[] y;
        delete mdf;
    }

    outFile.close();
    std::cout << "MultiDim fits complete.\n";

    // === Example histogram after correction ===
    TH2F *hFVB = new TH2F("hFVB", "Corrected Up+Dn vs Back;Back E;Corrected Up+Dn E",
                          400, 0, 16000, 400, 0, 16000);

    for (const auto &kv : dataPoints)
    {
        auto [id, bk, u, d] = kv.first;
        if (!gainValid[id][bk][u][d])
            continue;

        // Recreate the fitted model to evaluate
        TMultiDimFit *mdf = new TMultiDimFit(3, TMultiDimFit::kMonomials, "v");
        mdf->SetMaxPowers(2);
        mdf->SetMaxTerms(10);

        // Re-fill points to refit (if needed — or you can serialize coefficients instead)
        const auto &pts = kv.second;
        for (const auto &pr : pts)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pr;
            double vars[3] = {eBk, eUp, eDn};
            double y = eUp + eDn;
            mdf->AddRow(vars, y);
        }

        mdf->FindParameterization();

        // Fill histogram with corrected "front" from model
        for (const auto &pr : kv.second)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pr;
            double vars[3] = {eBk, eUp, eDn};
            double correctedFront = mdf->Eval(vars);
            if (eBk == 0 || correctedFront == 0)
                continue;
            hFVB->Fill(eBk, correctedFront);
        }

        delete mdf;
    }

    // Save histogram if needed
    // TFile *outHist = new TFile("sx3_multidimfit_hists.root", "RECREATE");
    // hFVB->Write();
    // outHist->Close();
}
