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

#include "Armory/ClassSX3.h"

#include "TVector3.h"

TH2F *hSX3FvsB;
TH2F *hQQQFVB;

int padID = 0;

SX3 sx3_contr;
TCutG *cut;
std::map<std::tuple<int, int, int>, std::vector<std::pair<double, double>>> dataPoints;

void GainMatchSX3::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);

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
        ID.push_back(std::pair<int, int>(sx3.id[i], i));
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
                        sx3ChDn = sx3.ch[index]/2;
                        sx3EDn = sx3.e[index];
                    }
                    else
                    {
                        sx3ChUp = sx3.ch[index]/2;
                        sx3EUp = sx3.e[index];
                    }
                }
                else
                {
                    sx3ChBk = sx3.ch[index];
                    sx3EBk = sx3.e[index];
                }
            }
            hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);

            for (int i = 0; i < sx3.multi; i++)
            {
                TString histName = Form("hSX3FVB_id%d_F%d_L+R%d", sx3.id[i], sx3ChUp, sx3ChBk);
                TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                if (!hist2d)
                {
                    hist2d = new TH2F(histName, Form("hSX3FVB_id%d_F%d_L+R%d", sx3.id[i], sx3ChUp, sx3ChBk), 400, 0, 16000, 400, 0, 16000);
                }

                hist2d->Fill(sx3EUp + sx3EDn, sx3EBk);
                if (cut && cut->IsInside(sx3EUp + sx3EDn, sx3EBk))
                {
                    // Accumulate data for gain matching
                    dataPoints[{sx3.id[i], sx3ChUp, sx3ChBk}].emplace_back(sx3EBk, sx3EUp + sx3EDn);
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

    std::ofstream outFile("sx3_GainMatch.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    for (const auto &kv : dataPoints)
    {
        auto [id, ud,bk] = kv.first;
        const auto &pts = kv.second;
        if (pts.size() < 5)
            continue;

        std::vector<double> bkE, udE;
        for (const auto &pr : pts)
        {
            bkE.push_back(pr.first);
            udE.push_back(pr.second);
        }

        TGraph g(bkE.size(), bkE.data(), udE.data());
        TF1 f("f", "[0]*x", 0, 16000);
        g.Fit(&f, "QNR");
        gainArray[id][ud][bk] = f.GetParameter(0);
        gainValid[id][ud][bk] = true;
    }

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

    outFile.close();
    std::cout << "Gain matching complete." << std::endl;

    // === Plot all gain-matched QQQ points together with a 2D histogram ===
    TH2F *hAll = new TH2F("hAll", "All SX3 Gain-Matched;Corrected Back E;Up+dn E",
                          400, 0, 16000, 400, 0, 16000);

    // Fill the combined TH2F with corrected data
    for (auto &kv : dataPoints)
    {
        int id, ud, bk;
        std::tie(id, ud, bk) = kv.first;
        if (!gainValid[id][ud][bk])
            continue;
        auto &pts = kv.second;
        for (auto &pr : pts)
        {
            double corrBack = pr.first * gainArray[id][ud][bk];
            double udE = pr.second;
            hAll->Fill(corrBack, udE);
        }
    }
}
