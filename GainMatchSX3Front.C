#define GainMatchSX3Front_cxx

#include "GainMatchSX3Front.h"
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

// Gain arrays

const int MAX_DET = 24;
const int MAX_UP = 4;
const int MAX_DOWN = 4;
const int MAX_BK = 4;
double backGain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool backGainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};
double frontGain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool frontGainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

void GainMatchSX3Front::Begin(TTree * /*tree*/)
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
    bool cutLoaded = (cut != nullptr);
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
    std::string filename = "sx3_GainMatchback.txt";

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
        backGain[id][bk][u][d] = gain;
        backGainValid[id][bk][u][d] = true;
    }

    infile.close();
    std::cout << "Loaded back gains from " << filename << std::endl;
    SX3 sx3_contr;
}

Bool_t GainMatchSX3Front::Process(Long64_t entry)
{

    b_sx3Multi->GetEntry(entry);
    b_sx3ID->GetEntry(entry);
    b_sx3Ch->GetEntry(entry);
    b_sx3E->GetEntry(entry);
    b_sx3T->GetEntry(entry);

    sx3.CalIndex();

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
            // If we have a valid front and back channel, fill the histograms
            hSX3->Fill(sx3ChDn + 4, sx3ChBk);
            hSX3->Fill(sx3ChUp, sx3ChBk);

            // Fill the histogram for the front vs back
            hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);

            for (int i = 0; i < sx3.multi; i++)
            {
                if (sx3.id[i] == 3 && sx3.e[i] > 100)
                {
                    // back gain correction

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

                    if (cut && cut->IsInside(sx3EUp + sx3EDn, sx3EBk) &&
                        cut1 && cut1->IsInside(sx3EUp / sx3EBk, sx3EDn / sx3EBk))
                    {

                        if (backGainValid[sx3.id[i]][sx3ChBk][sx3ChUp][sx3ChDn])
                        {
                            sx3EBk *= backGain[sx3.id[i]][sx3ChBk][sx3ChUp][sx3ChDn];
                        }
                        // Accumulate data for gain matching
                        dataPoints[{sx3.id[i], sx3ChBk, sx3ChUp, sx3ChDn}].emplace_back(sx3EBk, sx3EUp, sx3EDn);
                    }
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

void GainMatchSX3Front::Terminate()
{

    std::map<std::tuple<int, int, int, int>, TVectorD> fitCoefficients;

    // === Gain matching ===

    std::ofstream outFile("sx3_GainMatchfront.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    TH2F *hUvD = new TH2F("hUvD", " UvD; Up/CorrBack; Down/CorrBack", 600, 0, 1, 600, 0, 1);

    for (const auto &kv : dataPoints)
    {
        auto [id, bk, u, d] = kv.first;
        const auto &pts = kv.second;

        if (pts.size() < 5)
            continue;

        std::vector<double> udE, corrBkE;

        for (const auto &pr : pts)
        {
            double eBkCorr, eUp, eDn;
            std::tie(eBkCorr, eUp, eDn) = pr;
            udE.push_back(eUp + eDn);
            corrBkE.push_back(eBkCorr);
            hUvD->Fill(eUp / eBkCorr, eDn / eBkCorr);
        }

        TGraph g(udE.size(), udE.data(), corrBkE.data());
        TF1 f("f", "[0]*x", 0, 40000);
        g.Fit(&f, "QNR");

        frontGain[id][bk][u][d] = f.GetParameter(0);
        frontGainValid[id][bk][u][d] = true;

        outFile << id << " " << bk << " " << u << " " << d << " " << frontGain[id][bk][u][d] << std::endl;
        printf("Front gain Det%d Back%d Up%dDn%d → %.4f\n", id, bk, u, d, frontGain[id][bk][u][d]);
    }

    outFile.close();
    std::cout << "Gain matching complete." << std::endl;

    // === Stage 3: Create corrected histogram ===
    TH2F *hCorrectedFvB = new TH2F("hCorrectedFvB", "Corrected;Corrected Front Sum;Corrected Back", 800, 0, 16000, 800, 0, 16000);
    TH2F *hCorrectedUvD = new TH2F("hCorrectedUvD", "Corrected UvD; UvD Up; UvD Down", 600, 0, 1, 600, 0, 1);

    for (const auto &kv : dataPoints)
    {

        auto [id, bk, u, d] = kv.first;
        double front = frontGain[id][bk][u][d];

        for (const auto &pr : kv.second)
        {
            double eBk, eUp, eDn;
            std::tie(eBk, eUp, eDn) = pr;
            double corrUp = eUp * front;
            double corrDn = eDn * front;

            hCorrectedFvB->Fill(corrUp + corrDn, eBk);
            hCorrectedUvD->Fill(corrUp / eBk, corrDn / eBk);
        }
    }

    // === Final canvas ===
    gStyle->SetOptStat(1110);
    TCanvas *c1 = new TCanvas("c1", "Gain Correction Results", 1200, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    hSX3FvsB_g->SetTitle("Before Correction (Gated)");
    hSX3FvsB_g->GetXaxis()->SetTitle("Measured Front Sum (E_Up + E_Dn)");
    hSX3FvsB_g->GetYaxis()->SetTitle("Measured Back E");
    hSX3FvsB_g->Draw("colz");

    c1->cd(2);
    hCorrectedFvB->SetTitle("After Correction");
    hCorrectedFvB->Draw("colz");
    TF1 *diag = new TF1("diag", "x", 0, 40000);
    diag->SetLineColor(kRed);
    diag->SetLineWidth(2);
    diag->Draw("same");

    std::cout << "Terminate() completed successfully." << std::endl;
}