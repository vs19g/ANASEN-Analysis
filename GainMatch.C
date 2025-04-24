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

const int MAX_DET = 4;
const int MAX_WEDGE = 16;
const int MAX_RING = 16;

bool gainValid[MAX_DET][MAX_RING][MAX_WEDGE] = {{{false}}};
double gainArray[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};

void GainMatch::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hQQQFVB = new TH2F("hQQQFVB", "number of good QQQ vs QQQ id", 400, 0, 16000, 400, 0, 16000);

    sx3_contr.ConstructGeo();
    pw_contr.ConstructGeo();

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
                    hist2d = new TH2F(histName, Form("QQQ GainMatch Det%d R%d W%d;Wedge E;Ring E", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000);
                }

                hist2d->Fill(eWedge, eRing);
            }
        }
    }

    return kTRUE;
}

void GainMatch::Terminate()
{

    TFile *cutFile = TFile::Open("qqqcorr.root");
    TCutG *cut = (TCutG *)cutFile->Get("qqqcorr");
    const int MAX_DET = 4;
    const int MAX_RING = 16;
    const int MAX_WEDGE = 16;

    // Store gains and validity
    static double gainArray[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    static bool gainValid[MAX_DET][MAX_RING][MAX_WEDGE] = {{{false}}};

    std::ofstream outFile("qqq_gainmatch.txt");
    if (!outFile.is_open())
    {
        printf("Error opening output file!");
        return;
    }

    // Collect all (wedgeE, ringE) pairs per detector/ring/wedge
    std::map<std::tuple<int, int, int>, std::vector<std::pair<double, double>>> dataPoints;
    TIter next(gDirectory->GetList());
    TObject *obj;
    while ((obj = next()))
    {
        if (!obj->InheritsFrom("TH2F"))
            continue;
        TH2F *hist2d = static_cast<TH2F *>(obj);
        TString name = hist2d->GetName();
        if (!name.BeginsWith("hQQQFVB_id"))
            continue;
        if (hist2d->GetEntries() < 100)
            continue;

        int id, ring, wedge;
        sscanf(name.Data(), "hQQQFVB_id%d_r%d_w%d", &id, &ring, &wedge);

        for (int binX = 1; binX <= hist2d->GetNbinsX(); ++binX)
        {
            TH1D *projY = hist2d->ProjectionY("_py", binX, binX + 1);
            if (projY->GetEntries() < 30)
            {
                delete projY;
                continue;
            }

            double mean = projY->GetMean();
            TF1 *fit = new TF1("fit", "gaus", mean - 100, mean + 100);
            projY->Fit(fit, "QNR");

            double wedgeE = hist2d->GetXaxis()->GetBinCenter(binX);
            double ringE = fit->GetParameter(1);
            dataPoints[{id, ring, wedge}].emplace_back(wedgeE, ringE);

            delete fit;
            delete projY;
        }
    }

    // Fit slopes with sigma-clipping
    for (auto &kv : dataPoints)
    {
        auto [id, ring, wedge] = kv.first;
        auto &pts = kv.second;
        if (pts.size() < 5)
            continue;

        // Build vectors
        std::vector<double> wE, rE;
        for (auto &pr : pts)
        {
            wE.push_back(pr.first);
            rE.push_back(pr.second);
        }

        // Initial fit
        TGraph *g0 = new TGraph(wE.size(), wE.data(), rE.data());
        TF1 *f0 = new TF1("f0", "[0]*x + [1]", 0, 16000);
        g0->Fit(f0, "QNR");
        double m0 = f0->GetParameter(0), b0 = f0->GetParameter(1);

        // Clip to cut
        std::vector<double> wC, rC;
        for (size_t i = 0; i < wE.size(); ++i)
        {
            if (cut->IsInside(wE[i], rE[i]))
            {
                wC.push_back(wE[i]);
                rC.push_back(rE[i]);
            }
        }

        // Final fit
        if (wC.size() >= 5)
        {
            TGraph *g1 = new TGraph(wC.size(), wC.data(), rC.data());
            TF1 *f1 = new TF1("f1", "[0]*x + [1]", 0, 16000);
            g1->Fit(f1, "QNR");
            gainArray[id][ring][wedge] = f1->GetParameter(0);
            gainValid[id][ring][wedge] = true;
            delete g1;
            delete f1;
        }
        delete g0;
        delete f0;
    }

    // Write out valid entries
    for (int id = 0; id < MAX_DET; ++id)
    {
        for (int ring = 0; ring < MAX_RING; ++ring)
        {
            for (int wedge = 0; wedge < MAX_WEDGE; ++wedge)
            {
                if (!gainValid[id][ring][wedge])
                    continue;
                outFile << id << " " << wedge << " " << ring << " " << gainArray[id][ring][wedge] << std::endl;
                // printf("Gain match Det%d Ring%d Wedge%d → %.4f \n", id, ring, wedge, gainArray[id][ring][wedge]);
            }
        }
    }
    outFile.close();
    printf("Gain matching complete (with sigma-clipping)");

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
