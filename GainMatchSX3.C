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
#include <TGraphErrors.h>
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
std::map<std::tuple<int, int, int, int>, int> comboCounts;

const int MAX_DET = 24;
const int MAX_UP = 4;
const int MAX_DOWN = 4;
const int MAX_BK = 4;

double frontGain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool frontGainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

// ==== Configuration Flags ====
const bool interactiveMode = false; // If true: show canvas + wait for user
const bool verboseFit = true;       // If true: print fit summary and chi²
const bool drawCanvases = true;     // If false: canvases won't be drawn at all

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
    bool cut1Loaded = (cut1 != nullptr);
    cut1 = dynamic_cast<TCutG *>(cutFile1->Get("UvD"));
    if (!cut1)
    {
        std::cerr << "Error: Could not find TCutG named 'UvD' in UvD.root" << std::endl;
        return;
    }
    cut1->SetName("UvD");
    std::string filename = "sx3_GainMatchfront.txt";

    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error opening " << filename << "!" << std::endl;
        return;
    }

    int id, bk, u, d;
    double gain;
    while (infile >> id >> bk >> u >> d >> gain)
    {
        frontGain[id][bk][u][d] = gain;
        frontGainValid[id][bk][u][d] = true;
    }
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

        // for (int j = i + 1; j < sx3.multi; j++)
        // {
        //     if (sx3.id[i] == 3)
        //         hsx3Coin->Fill(sx3.index[i], sx3.index[j]);
        // }
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
                if (sx3.e[i] > 100)
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
            }
            for (int i = 0; i < sx3.multi; i++)
            {
                auto key = std::make_tuple(sx3.id[i], sx3ChBk, sx3ChUp, sx3ChDn);
                comboCounts[key]++;
                // If we have a valid front and back channel, fill the histograms
                hSX3->Fill(sx3ChDn+4, sx3ChBk);
                hSX3->Fill(sx3ChUp, sx3ChBk);
    
                // Fill the histogram for the front vs back
                hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);
            }

            
            for (int i = 0; i < sx3.multi; i++)
            {
                // if (sx3.id[i] == 4)
                {
                    auto key = std::make_tuple(sx3.id[i], sx3ChBk, sx3ChUp, sx3ChDn);

                    // Only continue if this combo has enough entries
                    if (comboCounts[key] < 100 || sx3EBk < 100 || sx3EUp < 100 || sx3EDn < 100)
                        continue;
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

                    // if (cut && cut->IsInside(sx3EUp + sx3EDn, sx3EBk))// && cut1 && cut1->IsInside(sx3EUp / sx3EBk, sx3EDn / sx3EBk))
                    {
                        // Accumulate data for gain matching
                        // if (frontGainValid[sx3.id[i]][sx3ChBk][sx3ChUp][sx3ChDn])
                        // {
                        //     sx3EUp *= frontGain[sx3.id[i]][sx3ChBk][sx3ChUp][sx3ChDn];
                        // }
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
    double fbgain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
    bool fbgainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

    // std::map<int, TH2F *> updn2DHistos;
    std::map<int, double> upCorrFactor;

    // === Gain matching ===

    std::ofstream outFile("sx3_GainMatchback.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Gain fit using up+dn vs bk
    for (const auto &kv : dataPoints)

    {
        // kv.first is a tuple of (id, up, bk)
        // kv.second is a vector of tuples (bkE, upE, dnE)
        auto [id, bk, u, d] = kv.first;
        const auto &pts = kv.second;
        // Check if we have enough points for fitting
        if (pts.size() < 5)
            continue;

        std::vector<double> bkE, udE;

        for (const auto &pr : pts)
        {
            double eUp, eDn, eBk;
            std::tie(eBk, eUp, eDn) = pr;
            if ((eBk < 100) || (eUp < 100) || (eDn < 100))
                continue; // Skip if any energy is less than 100
            bkE.push_back(eBk);
            udE.push_back(eUp + eDn);
        }

        // Fill the TGraph with bkE and udE
        // TGraph g(bkE.size(), bkE.data(), udE.data());
        // Fit the graph to a linear function
        if (bkE.size() < 5)
            continue; // Ensure we have enough points for fitting

        // TF1 f("f", "[0]*x", 0, 16000);
        // g.Fit(&f, "NR");

        // if (TMath::Abs(f.GetParameter(0) - 1) > 3.0)
        //     continue;

        const double fixedError = 10.0; // in ADC channels

        std::vector<double> xVals, yVals, exVals, eyVals;

        // Build data with fixed error
        for (size_t i = 0; i < udE.size(); ++i)
        {
            double x = udE[i]; // front energy
            double y = bkE[i]; // back energy

            xVals.push_back(x);
            yVals.push_back(y);
            exVals.push_back(fixedError); // error in front energy
            eyVals.push_back(fixedError); // error in back energy
        }

        // Build TGraphErrors with errors
        TGraphErrors g(xVals.size(), xVals.data(), yVals.data(), exVals.data(), eyVals.data());

        TF1 f("f", "[0]*x", 0, 16000);
        f.SetParameter(0, 1.0); // Initial guess

        if (drawCanvases)
        {
            TCanvas *c = new TCanvas(Form("c_%d_%d_%d_%d", id, bk, u, d), "Fit", 800, 600);
            g.SetTitle(Form("Detector %d: U%d D%d B%d", id, u, d, bk));
            g.SetMarkerStyle(20);
            g.SetMarkerColor(kBlue);
            g.Draw("AP");

            g.Fit(&f, interactiveMode ? "Q" : "QNR"); // 'R' avoids refit, 'N' skips drawing

            if (verboseFit)
            {
                double chi2 = f.GetChisquare();
                int ndf = f.GetNDF();
                double reducedChi2 = (ndf != 0) ? chi2 / ndf : -1;

                std::cout << Form("Det%d U%d D%d B%d → Gain: %.4f | χ²/ndf = %.2f/%d = %.2f",
                                  id, u, d, bk, f.GetParameter(0), chi2, ndf, reducedChi2)
                          << std::endl;
            }

            if (interactiveMode)
            {
                c->Update();
                gPad->WaitPrimitive();
            }
            else
            {
                c->Close(); // Optionally avoid clutter in batch
            }
        }
        else
        {
            g.Fit(&f, "QNR");
        }

        gainArray[id][bk][u][d] = f.GetParameter(0);
        gainValid[id][bk][u][d] = true;
        // }

        // // Output results
        // for (int id = 0; id < MAX_DET; ++id)
        // {
        //     for (int bk = 0; bk < MAX_BK; ++bk)
        //     {
        //         for (int u = 0; u < MAX_UP; ++u)
        //         {
        //             for (int d = 0; d < MAX_DOWN; ++d)
        //             {
        //                 // Check if the gain is valid for this detector, back, up, and down
        //                 if (gainValid[id][bk][u][d])
        //                 {
        if (TMath::Abs(gainArray[id][u][d][bk] - 1) < 0.3)
        {
            printf("Gain match Det%d Up%dDn%d Backs%d → %.4f \n", id, u, d, bk, gainArray[id][u][d][bk]);
            outFile << id << " " << bk << " " << u << " " << d << " " << gainArray[id][u][d][bk] << std::endl;
        }
        else if (gainArray[id][u][d][bk] != 0)
        {
            std::cerr << "Warning: Gain value out of range for Det " << id << " Up " << u << " Dn " << d << " Back " << bk << ": "
                      << gainArray[id][u][d][bk] << std::endl;
        }
    }
    //             }
    //         }
    //     }
    // }
    // }

    // for (int bk = 0; bk < MAX_BK; ++bk)
    // {
    //     TString name = Form("hUpDnVsBk_%d", bk);
    //     TString title = Form("Up/Bk vs Dn/Bk for Back %d;Dn/Bk;Up/Bk", bk);
    //     updn2DHistos[bk] = new TH2F(name, title, 400, 0, 1, 400, 0, 1);
    // }

    outFile.close();
    std::cout << "Gain matching complete." << std::endl;

    // === Create histograms ===
    TH2F *hFVB = new TH2F("hFVB", "Corrected Up+Dn vs Corrected Back;Corrected Back E;Up+Dn E",
                          600, 0, 16000, 600, 0, 16000);
    TH2F *hAsym = new TH2F("hAsym", "Up vs Dn dvide corrected back;Up/Back E;Dn/Back E",
                           400, 0.0, 1.0, 400, 0.0, 1.0);

    // Fill histograms
    for (const auto &kv : dataPoints)
    {
        auto [id, u, d, bk] = kv.first;
        if (!gainValid[id][u][d][bk])
            continue;
        double gain = gainArray[id][u][d][bk];

        // Prepare vectors to hold the points for TGraph
        std::vector<double> xVals;
        std::vector<double> yVals;

        for (const auto &pr : kv.second)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pr;

            double updn = eUp + eDn;
            if (updn == 0 || eBk == 0)
                continue;

            double asym = (eUp - eDn) / updn;
            double correctedBack = eBk * gain;

            hFVB->Fill(correctedBack, updn);
            hAsym->Fill(eUp / correctedBack, eDn / correctedBack);
        }
    }
}