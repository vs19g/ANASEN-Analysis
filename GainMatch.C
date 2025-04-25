#define GainMatch_cxx

#include "GainMatch.h"
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
#include "Armory/ClassPW.h"

#include "TVector3.h"

TH2F *hSX3FvsB;
TH2F *hQQQFVB;

int padID = 0;

SX3 sx3_contr;
PW pw_contr;
TVector3 hitPos;
bool HitNonZero;
TCutG *cut;
std::map<std::tuple<int, int, int>, std::vector<std::pair<double, double>>> dataPoints;

void GainMatch::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hQQQFVB = new TH2F("hQQQFVB", "number of good QQQ vs QQQ id", 400, 0, 16000, 400, 0, 16000);

    sx3_contr.ConstructGeo();
    pw_contr.ConstructGeo();

    // Load the TCutG object
    TFile *cutFile = TFile::Open("qqqcorr.root");
    if (!cutFile || cutFile->IsZombie())
    {
        std::cerr << "Error: Could not open qqqcorr.root" << std::endl;
        return;
    }
    cut = dynamic_cast<TCutG *>(cutFile->Get("qqqcorr"));
    if (!cut)
    {
        std::cerr << "Error: Could not find TCutG named 'qqqcorr' in qqqcorr.root" << std::endl;
        return;
    }
    cut->SetName("qqqcorr"); // Ensure the cut has the correct name
}

Bool_t GainMatch::Process(Long64_t entry)
{
    hitPos.Clear();
    HitNonZero = false;

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
                    sx3ChBk = sx3.ch[index];
                    sx3EBk = sx3.e[index];
                }
            }
            hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);
        }
    }

    for (int i = 0; i < qqq.multi; i++)
    {
        for (int j = i + 1; j < qqq.multi; j++)
        {
            if (qqq.id[i] == qqq.id[j])
            {
                int chWedge = -1;
                int chRing = -1;
                float eWedge = 0.0;
                float eRing = 0.0;
                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16)
                {
                    chWedge = qqq.ch[i];
                    eWedge = qqq.e[i];
                    chRing = qqq.ch[j] - 16;
                    eRing = qqq.e[j];
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16)
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j];
                    chRing = qqq.ch[i] - 16;
                    eRing = qqq.e[i];
                }
                else
                    continue;

                hQQQFVB->Fill(eWedge, eRing);

                TString histName = Form("hQQQFVB_id%d_r%d_w%d", qqq.id[i], chRing, chWedge);
                TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                if (!hist2d)
                {
                    hist2d = new TH2F(histName, Form("QQQ Det%d R%d W%d;Wedge E;Ring E", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000);
                }

                hist2d->Fill(eWedge, eRing);
                if (cut && cut->IsInside(eWedge, eRing))
                {
                    // Accumulate data for gain matching
                    dataPoints[{qqq.id[i], chRing, chWedge}].emplace_back(eWedge, eRing);
                }
            }
        }
    }

    return kTRUE;
}

void GainMatch::Terminate()
{
    const int MAX_DET = 4;
    const int MAX_RING = 16;
    const int MAX_WEDGE = 16;

    double gainArray[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    bool gainValid[MAX_DET][MAX_RING][MAX_WEDGE] = {{{false}}};

    std::ofstream outFile("qqq_gainmatch.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    for (const auto &kv : dataPoints)
    {
        auto [id, ring, wedge] = kv.first;
        const auto &pts = kv.second;
        if (pts.size() < 5)
            continue;

        std::vector<double> wE, rE;
        for (const auto &pr : pts)
        {
            wE.push_back(pr.first);
            rE.push_back(pr.second);
        }

        TGraph g(wE.size(), wE.data(), rE.data());
        TF1 f("f", "[0]*x", 0, 16000);
        g.Fit(&f, "QNR");
        gainArray[id][ring][wedge] = f.GetParameter(0);
        gainValid[id][ring][wedge] = true;
    }

    for (int id = 0; id < MAX_DET; ++id)
    {
        for (int ring = 0; ring < MAX_RING; ++ring)
        {
            for (int wedge = 0; wedge < MAX_WEDGE; ++wedge)
            {
                if (gainValid[id][ring][wedge])
                {
                    outFile << id << " " << wedge << " " << ring << " " << gainArray[id][ring][wedge] << std::endl;
                    printf("Gain match Det%d Ring%d Wedge%d → %.4f \n", id, ring, wedge, gainArray[id][ring][wedge]);
                }
            }
        }
    }

    outFile.close();
    std::cout << "Gain matching complete." << std::endl;

    // === Plot all gain-matched QQQ points together with a 2D histogram ===
    TH2F *hAll = new TH2F("hAll", "All QQQ Gain-Matched;Corrected Wedge E;Ring E",
                          400, 0, 16000, 400, 0, 16000);

    // Fill the combined TH2F with corrected data
    for (auto &kv : dataPoints)
    {
        int id, ring, wedge;
        std::tie(id, ring, wedge) = kv.first;
        if (!gainValid[id][ring][wedge])
            continue;
        auto &pts = kv.second;
        for (auto &pr : pts)
        {
            double corrWedge = pr.first * gainArray[id][ring][wedge];
            double ringE = pr.second;
            hAll->Fill(corrWedge, ringE);
        }
    }
}
