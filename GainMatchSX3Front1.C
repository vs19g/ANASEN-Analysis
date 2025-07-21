#define GainMatchSX3_cxx

#include "GainMatchSX3.h"
#include "Armory/ClassSX3.h"
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TApplication.h>
#include <map>
#include <vector>
#include <tuple>
#include <fstream>
#include <iostream>
#include <algorithm>

// Constants
const int MAX_DET = 24;
const int MAX_BK = 4;
const int MAX_UP = 4;
const int MAX_DOWN = 4;

// Gain arrays
double backGain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool backGainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

double frontGain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool frontGainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

// Data container
std::map<std::tuple<int, int, int, int>, std::vector<std::tuple<double, double, double>>> dataPoints;

// Load back gains
void LoadBackGains(const std::string &filename)
{
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

// Front gain matching function
Bool_t GainMatchSX3::Process(Long64_t entry)
{
    // Link SX3 branches
    
    b_sx3Multi->GetEntry(entry);
    b_sx3ID->GetEntry(entry);
    b_sx3Ch->GetEntry(entry);
    b_sx3E->GetEntry(entry);
    b_sx3T->GetEntry(entry);

    sx3.CalIndex();

    Long64_t nentries = tree->GetEntries(Long64_t entry);
    std::cout << "Total entries: " << nentries << std::endl;

    TH2F *hBefore = new TH2F("hBefore", "Before Correction;E_Up+E_Dn;Back Energy", 400, 0, 40000, 400, 0, 40000);
    TH2F *hAfter = new TH2F("hAfter", "After Correction;E_Up+E_Dn;Corrected Back Energy", 400, 0, 40000, 400, 0, 40000);

    for (Long64_t entry = 0; entry < nentries; ++entry)
    {
        tree->GetEntry(entry);
        sx3.CalIndex();

        std::vector<std::pair<int, int>> ID;

        for (int i = 0; i < sx3.multi; i++)
        {
            if (sx3.e[i] > 100)
            {
                ID.push_back({sx3.id[i], i});
            }
        }

        if (ID.empty())
            continue;

        // Sort by id
        std::sort(ID.begin(), ID.end(), [](auto &a, auto &b) { return a.first < b.first; });

        std::vector<std::pair<int, int>> sx3ID;
        sx3ID.push_back(ID[0]);
        bool found = false;

        for (size_t i = 1; i < ID.size(); i++)
        {
            if (ID[i].first == sx3ID.back().first)
            {
                sx3ID.push_back(ID[i]);
                if (sx3ID.size() >= 3)
                    found = true;
            }
            else if (!found)
            {
                sx3ID.clear();
                sx3ID.push_back(ID[i]);
            }
        }

        if (!found)
            continue;

        int sx3ChUp = -1, sx3ChDn = -1, sx3ChBk = -1;
        float sx3EUp = 0.0, sx3EDn = 0.0, sx3EBk = 0.0;
        int sx3id = sx3ID[0].first;

        for (auto &[id, idx] : sx3ID)
        {
            if (sx3.ch[idx] < 8)
            {
                if (sx3.ch[idx] % 2 == 0)
                {
                    sx3ChDn = sx3.ch[idx] / 2;
                    sx3EDn = sx3.e[idx];
                }
                else
                {
                    sx3ChUp = sx3.ch[idx] / 2;
                    sx3EUp = sx3.e[idx];
                }
            }
            else
            {
                sx3ChBk = sx3.ch[idx] - 8;
                sx3EBk = sx3.e[idx];
            }
        }

        if (sx3ChUp < 0 || sx3ChDn < 0 || sx3ChBk < 0)
            continue;

        if (!backGainValid[sx3id][sx3ChBk][sx3ChUp][sx3ChDn])
            continue;

        double corrBk = sx3EBk * backGain[sx3id][sx3ChBk][sx3ChUp][sx3ChDn];

        hBefore->Fill(sx3EUp + sx3EDn, sx3EBk);
        hAfter->Fill(sx3EUp + sx3EDn, corrBk);

        dataPoints[{sx3id, sx3ChBk, sx3ChUp, sx3ChDn}].emplace_back(corrBk, sx3EUp, sx3EDn);
    }

    // === Fit front gains ===
    std::ofstream outFile("sx3_GainMatchfront.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening sx3_GainMatchfront.txt!" << std::endl;
        return;
    }

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
    std::cout << "Front gain matching complete." << std::endl;

    // === Draw diagnostic plots ===
    gStyle->SetOptStat(1110);
    TCanvas *c = new TCanvas("c", "Gain Matching Diagnostics", 1200, 600);
    c->Divide(2, 1);

    c->cd(1);
    hBefore->Draw("colz");
    TF1 *diag1 = new TF1("diag1", "x", 0, 40000);
    diag1->SetLineColor(kRed);
    diag1->Draw("same");

    c->cd(2);
    hAfter->Draw("colz");
    TF1 *diag2 = new TF1("diag2", "x", 0, 40000);
    diag2->SetLineColor(kRed);
    diag2->Draw("same");
}

int main(int argc, char **argv)
{
    TApplication app("app", &argc, argv);

    // Load back gains
    LoadBackGains("sx3_GainMatchback.txt");

    // Open tree
    TFile *f = TFile::Open("input_tree.root"); // <<< Change file name
    if (!f || f->IsZombie())
    {
        std::cerr << "Cannot open input_tree.root!" << std::endl;
        return 1;
    }
    TTree *tree = (TTree *)f->Get("tree");
    if (!tree)
    {
        std::cerr << "Tree not found!" << std::endl;
        return 1;
    }

    // Run front gain matching
    GainMatchSX3(tree);

    app.Run();
    return 0;
}
