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
const bool drawCanvases = false;     // If false: canvases won't be drawn at all

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

    // std::string filename = "sx3_GainMatchfront.txt";

    // std::ifstream infile(filename);
    // if (!infile.is_open())
    // {
    //     std::cerr << "Error opening " << filename << "!" << std::endl;
    //     return;
    // }

    // int id, bk, u, d;
    // double gain;
    // while (infile >> id >> bk >> u >> d >> gain)
    // {
    //     frontGain[id][bk][u][d] = gain;
    //     frontGainValid[id][bk][u][d] = true;
    // }
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

            // Build the correlated set once
            for (size_t i = 0; i < sx3ID.size(); i++)
            {
                if (sx3.e[i] > 100)
                {
                    int index = sx3ID[i].second;
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
                        sx3EBk = sx3.e[index];
                    }
                }
            }

            // Only if we found all three channels do we proceed
            if (sx3ChUp >= 0 && sx3ChDn >= 0 && sx3ChBk >= 0)
            {
                // Fill once per correlated set
                hSX3->Fill(sx3ChDn + 4, sx3ChBk);
                hSX3->Fill(sx3ChUp, sx3ChBk);
                hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);

                // Pick detector ID from one of the correlated hits (all same detector)
                int detID = sx3ID[0].first;

                TString histName = Form("hSX3FVB_id%d_U%d_D%d_B%d",
                                        detID, sx3ChUp, sx3ChDn, sx3ChBk);
                TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                if (!hist2d)
                {
                    hist2d = new TH2F(histName, histName,
                                      400, 0, 16000, 400, 0, 16000);
                }

                if (sx3EBk > 100 || sx3EUp > 100 || sx3EDn > 100)
                {
                    hSX3FvsB_g->Fill(sx3EUp + sx3EDn, sx3EBk);

                    // Use the correlated triplet directly
                    dataPoints[{detID, sx3ChBk, sx3ChUp, sx3ChDn}]
                        .emplace_back(sx3EBk, sx3EUp, sx3EDn);
                }

                hist2d->Fill(sx3EUp + sx3EDn, sx3EBk);
            }
        }
    }

    return kTRUE;
}

const double GAIN_ACCEPTANCE_THRESHOLD = 0.3;
void GainMatchSX3::Terminate()
{
    double backSlope[MAX_DET][MAX_BK] = {{0}};
    bool backSlopeValid[MAX_DET][MAX_BK] = {{false}};

    std::ofstream outFile("sx3_BackGains.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening sx3_BackGains.txt for writing!" << std::endl;
        return;
    }

    // === Gain fit: (Up+Dn) vs Back, grouped by [id][bk] ===
    for (int id = 0; id < MAX_DET; id++)
    {
        for (int bk = 0; bk < MAX_BK; bk++)
        {
            std::vector<double> bkE, udE;

            // Collect all (Up+Dn, Back) for this id,bk
            for (const auto &kv : dataPoints)
            {
                auto [cid, cbk, u, d] = kv.first;
                if (cid != id || cbk != bk)
                    continue;

                for (const auto &pr : kv.second)
                {
                    double eBk, eUp, eDn;
                    std::tie(eBk, eUp, eDn) = pr;
                    if ((eBk < 100) || (eUp < 100) || (eDn < 100))
                        continue;

                    bkE.push_back(eBk);
                    udE.push_back(eUp + eDn);
                }
            }

            if (bkE.size() < 5)
                continue; // not enough statistics

            // Build graph with errors
            const double fixedError = 0.0;                     // ADC channels
            std::vector<double> exVals(udE.size(), 0.0);        // no x error
            std::vector<double> eyVals(udE.size(), fixedError); // constant y error

            TGraphErrors g(udE.size(), udE.data(), bkE.data(),
                           exVals.data(), eyVals.data());

            TF1 f("f", "[0]*x", 0, 16000);
            // f.SetParameter(0, 1.0); // initial slope

            if (drawCanvases)
            {
                TCanvas *c = new TCanvas(Form("c_%d_%d", id, bk), "Back Fit", 800, 600);
                g.SetTitle(Form("Detector %d Back %d: (Up+Dn) vs Back", id, bk));
                g.SetMarkerStyle(20);
                g.SetMarkerColor(kBlue);
                g.Draw("AP");

                g.Fit(&f, interactiveMode ? "Q" : "QNR");

                if (verboseFit)
                {
                    double chi2 = f.GetChisquare();
                    int ndf = f.GetNDF();
                    double reducedChi2 = (ndf != 0) ? chi2 / ndf : -1;

                    std::cout << Form("Det%d Back%d → Slope: %.4f | χ²/ndf = %.2f/%d = %.2f",
                                      id, bk, f.GetParameter(0), chi2, ndf, reducedChi2)
                              << std::endl;
                }

                if (interactiveMode)
                {
                    c->Update();
                    gPad->WaitPrimitive();
                }
                else
                {
                    c->Close();
                }
            }
            else
            {
                g.Fit(&f, "QNR");
            }

            double slope = 1/f.GetParameter(0);
            if (std::abs(slope - 1.0) < 0.3) // sanity check
            {
                backSlope[id][bk] = slope;
                backSlopeValid[id][bk] = true;
                outFile << id << " " << bk << " " << slope << "\n";
                printf("Back slope Det%d Bk%d → %.4f\n", id, bk, slope);
            }
            else
            {
                std::cerr << "Warning: Bad slope for Det" << id << " Bk" << bk
                          << " slope=" << slope << std::endl;
            }
        }
    }

    outFile.close();
    std::cout << "Back gain matching complete." << std::endl;

    // === Create histograms ===
    TH2F *hFVB = new TH2F("hFVB", "Corrected Up+Dn vs Corrected Back;Up+Dn E;Corrected Back E",
                          600, 0, 16000, 600, 0, 16000);
    TH2F *hAsym = new TH2F("hAsym", "Up vs Dn divide corrected back;Up/Back E;Dn/Back E",
                           400, 0.0, 1.0, 400, 0.0, 1.0);

    // Fill histograms using corrected back energies
    for (const auto &kv : dataPoints)
    {
        auto [id, bk, u, d] = kv.first;
        if (!backSlopeValid[id][bk])
            continue;

        double slope = backSlope[id][bk];

        for (const auto &pr : kv.second)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pr;

            double updn = eUp + eDn;
            if (updn == 0 || eBk == 0)
                continue;

            double correctedBack = eBk * slope;
            double asym = (eUp - eDn) / updn;

            hFVB->Fill(updn,correctedBack );
            hAsym->Fill(eUp / correctedBack, eDn / correctedBack);
        }
    }
}
