#define Calibration_cxx

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
#include <TVector3.h>
#include "Armory/ClassSX3.h"
#include "Armory/ClassPW.h"
#include "TGraphErrors.h"
#include "Calibration.h"

int padID = 0;

SX3 sx3_contr;
PW pw_contr;
PW pwinstance;
TVector3 hitPos;
// TVector3 anodeIntersection;
std::map<int, std::pair<double, double>> slopeInterceptMap;

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

TH2F *hSX3FvsB;
TH2F *hSX3FvsB_g;
TH2F *hSX3;
TH1F *hZProj;
TH2F *hsx3IndexVE;
TH2F *hsx3IndexVE_gm;
TH2F *hqqqIndexVE;
TH2F *hqqqIndexVE_gm;
TH2F *hsx3Coin;
TH2F *hqqqCoin;
TH2F *hqqqPolar;
TH1F *hsx3E_raw;
TH1F *hsx3E_calib;

TCutG *cut;
TCutG *cut1;

// Gain arrays

const int MAX_SX3 = 24;
const int MAX_UP = 4;
const int MAX_DOWN = 4;
const int MAX_BK = 4;
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
double backGain[MAX_SX3][MAX_BK] = {{0}};
bool backGainValid[MAX_SX3][MAX_BK] = {{false}};
double frontGain[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool frontGainValid[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};
double uvdslope[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
double qqqGain[MAX_QQQ][MAX_BK][MAX_UP] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_BK][MAX_UP] = {{{false}}};
TH1F *hSX3Spectra[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN];
TH1F *hQQQSpectra[MAX_QQQ][MAX_RING][MAX_WEDGE];

void Calibration::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hSX3FvsB_g = new TH2F("hSX3FvsB_g", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hSX3 = new TH2F("hSX3", "SX3 Front v Back; Fronts; Backs", 8, 0, 8, 4, 0, 4);
    hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);
    hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hsx3IndexVE_gm = new TH2F("hsx3IndexVE_cal", "SX3 index vs Energy (calibrated); SX3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hsx3E_raw = new TH1F("hsx3E_raw", "SX3 Back Energy (raw); Energy (arb); Counts", 4000, 0, 16000);
    hsx3E_calib = new TH1F("hsx3E_calib", "SX3 Back Energy (gm); Energy (kev); Counts", 4000, 0, 16000);
    hqqqIndexVE = new TH2F("hqqqIndexVE", "QQQ index vs Energy; QQQ index ; Energy", 4 * 2 * 16, 0, 4 * 2 * 16, 400, 0, 5000);
    hqqqIndexVE_gm = new TH2F("hqqqIndexVE_cal", "QQQ index vs Energy (calibrated); QQQ index ; Energy", 4 * 2 * 16, 0, 4 * 2 * 16, 400, 0, 5000);
    hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);
    hqqqCoin = new TH2F("hqqqCoin", "QQQ Coincident", 4 * 2 * 16, 0, 4 * 2 * 16, 4 * 2 * 16, 0, 4 * 2 * 16);

    hqqqPolar = new TH2F("hqqqPolar", "QQQ Polar ID", 16 * 4, -TMath::Pi(), TMath::Pi(), 16, 10, 50);

    sx3_contr.ConstructGeo();
    pw_contr.ConstructGeo();
    // ----------------------- Load Back Gains
    {
        std::string filename = "sx3_GainMatchback.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int id, bk;
            double gain;
            while (infile >> id >> bk >> gain)
            {
                backGain[id][bk] = gain;
                backGainValid[id][bk] = (gain > 0);
            }
            infile.close();
            std::cout << "Loaded back gains from " << filename << std::endl;
        }
    }

    // ----------------------- Load Front Gains
    {
        std::string filename = "sx3_GainMatchfront.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int id, bk, u, d;
            double gain;
            while (infile >> id >> bk >> u >> d >> gain)
            {
                frontGain[id][bk][u][d] = gain;
                frontGainValid[id][bk][u][d] = (gain > 0);
            }
            infile.close();
            std::cout << "Loaded front gains from " << filename << std::endl;
        }
    }

    // ----------------------- Load QQQ Gains
    {
        std::string filename = "qqq_GainMatch.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int det, ring, wedge;
            double gain;
            while (infile >> det >> ring >> wedge >> gain)
            {
                qqqGain[det][ring][wedge] = gain;
                qqqGainValid[det][ring][wedge] = (gain > 0);
            }
            infile.close();
            std::cout << "Loaded QQQ gains from " << filename << std::endl;
        }
    }

    for (int id = 0; id < MAX_SX3; id++)
    {
        for (int bk = 0; bk < MAX_BK; bk++)
        {
            for (int up = 0; up < MAX_UP; up++)
            {
                for (int dn = 0; dn < MAX_DOWN; dn++)
                {
                    TString hname = Form("hCal_id%d_bk%d_up%d_dn%d", id, bk, up, dn);
                    TString htitle = Form("SX3 id%d bk%d up%d dn%d; Energy (arb); Counts", id, bk, up, dn);
                    hSX3Spectra[id][bk][up][dn] = new TH1F(hname, htitle, 4000, 0, 16000);
                }
            }
        }
    }
    for (int det = 0; det < MAX_QQQ; det++)
    {
        for (int ring = 0; ring < MAX_RING; ring++)
        {
            for (int wedge = 0; wedge < MAX_WEDGE; wedge++)
            {
                TString hname = Form("hCal_qqq%d_ring%d_wedge%d", det, ring, wedge);
                TString htitle = Form("QQQ det%d ring%d wedge%d; Energy (arb); Counts", det, ring, wedge);
                hQQQSpectra[det][ring][wedge] = new TH1F(hname, htitle, 4000, 0, 16000);
            }
        }
    }

    SX3 sx3_contr;
}
Bool_t Calibration::Process(Long64_t entry)
{
    hitPos.Clear();
    HitNonZero = false;

    // Load branches
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

    sx3.CalIndex();
    qqq.CalIndex();
    pc.CalIndex();

    // ########################################################### Raw data
    sx3ecut = false;
    std::vector<std::pair<int, int>> ID; // first = id, 2nd = index
    for (int i = 0; i < sx3.multi; i++)
    {
        ID.emplace_back(sx3.id[i], i);
        hsx3IndexVE->Fill(sx3.index[i], sx3.e[i]);

        if (sx3.e[i] > 100)
            sx3ecut = true;

        for (int j = i + 1; j < sx3.multi; j++)
            hsx3Coin->Fill(sx3.index[i], sx3.index[j]);
    }

    // --- SX3 safe handling ---
    if (!ID.empty())
    {
        std::sort(ID.begin(), ID.end(), [](auto &a, auto &b)
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
                    found = true;
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
            float sx3EUp = 0.0f, sx3EDn = 0.0f, sx3EBk = 0.0f;

            for (size_t i = 0; i < sx3ID.size(); i++)
            {
                int index = sx3ID[i].second;

                if (sx3.ch[index] < 8)
                {
                    if ((sx3.ch[index] % 2) == 0) // even -> down
                    {
                        sx3ChDn = sx3.ch[index];
                        sx3EDn = sx3.e[index];
                    }
                    else // odd -> up
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

                bool haveFrontPair = (sx3ChUp >= 0 || sx3ChDn >= 0);
                bool haveBack = (sx3ChBk >= 0);
                double GM_EUp = 0.0, GM_EDn = 0.0, calibEBack = 0.0;

                if (haveBack)
                {
                    // --- ALWAYS fill raw ADC for diagnostics
                    // (temporarily use the existing spectrum to confirm fills)
                    // If you don't want raw values mixed with calibrated later, create a separate _raw array.
                    hSX3Spectra[sx3ID[i].first][sx3ChBk][sx3ChUp][sx3ChDn]->Fill(sx3EUp);

                    // --- If gain is available, also fill calibrated energy
                    if (frontGainValid[sx3ID[i].first][sx3ChBk][sx3ChUp][sx3ChDn])
                    {
                        GM_EUp = frontGain[sx3ID[i].first][sx3ChBk][sx3ChUp][sx3ChDn] * sx3EUp;
                        if (GM_EUp > 50.0)
                            hSX3Spectra[sx3ID[i].first][sx3ChBk][sx3ChUp][sx3ChDn]->Fill(GM_EUp); // optional: mixes raw+calib
                    }
                    // --- If back gain is available, also fill calibrated energy
                    hsx3E_raw->Fill(sx3EBk);

                    if (backGainValid[sx3ID[i].first][sx3ChBk])
                    {
                        calibEBack = backGain[sx3ID[i].first][sx3ChBk] * sx3EBk;
                        if (calibEBack > 50.0)
                            hsx3E_calib->Fill(calibEBack); // optional: mixes raw+calib
                    }

                    // Keep the other diagnostic plots
                    hsx3IndexVE_gm->Fill(sx3.index[sx3ID[i].second], GM_EUp);
                    hSX3->Fill(sx3ChDn + 4, sx3ChBk);
                    hSX3->Fill(sx3ChUp, sx3ChBk);
                    hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);

                    if (GM_EUp > 50.0 && sx3EBk > 50.0)
                    {
                        sx3_contr.CalSX3Pos(sx3ID[i].first, sx3ChUp, sx3ChDn, sx3ChBk, GM_EUp, sx3EDn);
                        hitPos = sx3_contr.GetHitPos();
                        HitNonZero = true;
                    }
                }
            }
        }
    }

    // ======================= QQQ =======================
    for (int i = 0; i < qqq.multi; i++)
    {
        int det = qqq.id[i];
        if (qqq.e[i] > 100)
            qqqEcut = true;

        for (int j = 0; j < qqq.multi; j++)
        {
            if (j == i)
                continue;
            hqqqCoin->Fill(qqq.index[i], qqq.index[j]);
        }

        for (int j = i + 1; j < qqq.multi; j++)
        {
            if (qqq.id[i] == qqq.id[j])
            {
                int chWedge = -1, chRing = -1;
                if (qqq.ch[i] < qqq.ch[j])
                {
                    chRing = qqq.ch[j] - 16;
                    chWedge = qqq.ch[i];
                }
                else
                {
                    chRing = qqq.ch[i];
                    chWedge = qqq.ch[j] - 16;
                }

                double Ecal = qqq.e[i];
                if (det >= 0 && det < MAX_QQQ &&
                    chRing >= 0 && chRing < MAX_RING &&
                    chWedge >= 0 && chWedge < MAX_WEDGE)
                {
                    // ALWAYS fill raw energy for diagnostics
                    hQQQSpectra[det][chRing][chWedge]->Fill(qqq.e[i]);

                    // If calibrated gain is present, also fill calibrated energy
                    if (qqqGainValid[det][chRing][chWedge])
                    {
                        double Ecal = qqq.e[i] * qqqGain[det][chRing][chWedge];
                        hQQQSpectra[det][chRing][chWedge]->Fill(Ecal); // optional: mixes raw+calib
                    }
                }

                hqqqIndexVE_gm->Fill(qqq.index[i], Ecal);
                hqqqIndexVE->Fill(qqq.index[i], qqq.e[i]);

                double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
                double rho = 50. + 40. / 16. * (chRing + 0.5);
                hqqqPolar->Fill(theta, rho);

                if (!HitNonZero)
                {
                    double x = rho * TMath::Cos(theta);
                    double y = rho * TMath::Sin(theta);
                    hitPos.SetXYZ(x, y, 23 + 75 + 30);
                    HitNonZero = true;
                }
            }
        }
    }

    return kTRUE;
}
void Calibration::Terminate()
{
    const double AM241_ALPHA = 5486.0; // keV

    // ----------------------- Summary Plots
    TH2F *hSX3Summary = new TH2F("hSX3Summary", "SX3 Channel Means;Channel Index;Mean (ADC)",
                                 MAX_SX3 * MAX_BK * MAX_UP * MAX_DOWN, 0, MAX_SX3 * MAX_BK * MAX_UP * MAX_DOWN,
                                 200, 0, 10000);

    TH2F *hQQQSummary = new TH2F("hQQQSummary", "QQQ Channel Means;Channel Index;Mean (ADC)",
                                 MAX_QQQ * MAX_RING * MAX_WEDGE, 0, MAX_QQQ * MAX_RING * MAX_WEDGE,
                                 200, 0, 10000);

    // ----------------------- SX3 Calibration (quick check with mean)
    for (int id = 0; id < MAX_SX3; id++)
    {
        for (int bk = 0; bk < MAX_BK; bk++)
        {
            for (int up = 0; up < MAX_UP; up++)
            {
                for (int dn = 0; dn < MAX_DOWN; dn++)
                {
                    TH1F *hSpec = hSX3Spectra[id][bk][up][dn];
                    if (!hSpec || hSpec->GetEntries() < 200)
                        continue;

                    double mean = hSpec->GetMean();

                    int sx3Index = (((id * MAX_BK + bk) * MAX_UP + up) * MAX_DOWN + dn);
                    hSX3Summary->Fill(sx3Index, mean);

                    std::cout << Form("SX3 id%d bk%d up%d dn%d → mean %.1f",
                                      id, bk, up, dn, mean)
                              << std::endl;
                }
            }
        }
    }

    // ----------------------- QQQ Calibration (quick check with mean)
    for (int det = 0; det < MAX_QQQ; det++)
    {
        for (int ring = 0; ring < MAX_RING; ring++)
        {
            for (int wedge = 0; wedge < MAX_WEDGE; wedge++)
            {
                TH1F *hSpec = hQQQSpectra[det][ring][wedge];
                if (!hSpec || hSpec->GetEntries() < 200)
                    continue;

                double mean = hSpec->GetMean();

                int qqqIndex = ((det * MAX_RING + ring) * MAX_WEDGE + wedge);
                hQQQSummary->Fill(qqqIndex, mean);

                std::cout << Form("QQQ det%d ring%d wedge%d → mean %.1f",
                                  det, ring, wedge, mean)
                          << std::endl;
            }
        }
    }

    // ----------------------- Draw Summary
    TCanvas *cSum = new TCanvas("cSum", "Calibration Summary (Means)", 1200, 600);
    cSum->Divide(2, 1);

    cSum->cd(1);
    hSX3Summary->Draw("COLZ");

    cSum->cd(2);
    hQQQSummary->Draw("COLZ");

    cSum->Update();
}
