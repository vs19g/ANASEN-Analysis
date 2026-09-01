#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCutG.h"
#include "TKey.h"
#include "TString.h"
#include "TSystem.h"
#include <map>
#include <vector>
#include <utility>
#include <cstdlib>

namespace
{
    constexpr double kAdcMax = 16384; // 14-bit-ish ADC range for axes

    struct HitData
    {
        double e;
        unsigned long long t;
    };
}

void MiscEdE(const char *f0,
             const char *f1 = "", const char *f2 = "", const char *f3 = "",
             const char *f4 = "", const char *f5 = "", const char *f6 = "",
             const char *f7 = "", const char *f8 = "", const char *f9 = "")
{
    std::vector<TString> inputs;
    for (const char *f : {f0, f1, f2, f3, f4, f5, f6, f7, f8, f9})
        if (f && f[0])
            inputs.emplace_back(f);
    if (inputs.empty())
    {
        printf("MiscEdE: no input files given.\n");
        return;
    }

    TChain *tree = new TChain("tree");
    for (auto &in : inputs)
    {
        if (gSystem->AccessPathName(in))
        {
            printf("MiscEdE: WARNING input not found, skipping: %s\n", in.Data());
            continue;
        }
        tree->Add(in);
    }
    if (tree->GetNtrees() == 0)
    {
        printf("MiscEdE: no readable inputs.\n");
        return;
    }

    // ---------------------------------------------------------
    // Bash Environment Variable Extraction
    // ---------------------------------------------------------
    int runNumber = -1;
    if (gSystem->Getenv("CURRENT_RUN")) {
        runNumber = std::atoi(gSystem->Getenv("CURRENT_RUN"));
    }
    printf("MiscEdE: Run Number loaded from Bash: %d\n", runNumber);

    TString cutDir = "~/ANASEN_analysis"; // Fallback if not run via bash
    if (gSystem->Getenv("CUT_DIR")) {
        cutDir = gSystem->Getenv("CUT_DIR");
    }

    // ---------------------------------------------------------
    // Run-Dependent Channel Routing
    // ---------------------------------------------------------
    // Defaults (Run <= 180)
    int snLolli = 405, chLolli = 9;
    int snSi    = 405, chSi    = 11;
    int snRF    = 405, chRF    = 15;
    int snMCP   = 405, chMCP   = 14;

    if (runNumber > 180 && runNumber <= 282) {
        chRF = 10;
    } 
    else if (runNumber > 282 && runNumber <= 322) {
        snRF = 89; chRF = 0;
        snMCP = 89; chMCP = 1;
    } 
    else if (runNumber > 322) {
        snRF = 89; chRF = 0;
        snMCP = 89; chMCP = 2;
    }

    // Composite Unique IDs (UID = Board * 100 + Channel) to prevent collisions
    int uidLollipop = snLolli * 100 + chLolli;
    int uidSiMon    = snSi * 100 + chSi;
    int uidRF       = snRF * 100 + chRF;
    int uidMCP      = snMCP * 100 + chMCP;

    printf("  -> Routed Lollipop : Board %d, Ch %d (UID: %d)\n", snLolli, chLolli, uidLollipop);
    printf("  -> Routed Si Monitor: Board %d, Ch %d (UID: %d)\n", snSi, chSi, uidSiMon);
    printf("  -> Routed RF       : Board %d, Ch %d (UID: %d)\n", snRF, chRF, uidRF);
    printf("  -> Routed MCP      : Board %d, Ch %d (UID: %d)\n", snMCP, chMCP, uidMCP);

    // ---------------------------------------------------------
    // Exact schema matching EventBuilder.cpp
    // ---------------------------------------------------------
    const int MAX_MULTI = 2000;
    unsigned int multi = 0;
    unsigned short sn[MAX_MULTI];
    unsigned short ch[MAX_MULTI];
    unsigned short e[MAX_MULTI];
    unsigned long long e_t[MAX_MULTI]; // Timestamp branch is "e_t"

    tree->SetBranchAddress("multi", &multi);
    tree->SetBranchAddress("sn", sn);
    tree->SetBranchAddress("ch", ch);
    tree->SetBranchAddress("e", e);
    tree->SetBranchAddress("e_t", e_t);

    TString stem = gSystem->BaseName(inputs[0].Data());
    if (stem.EndsWith(".root"))
        stem.Remove(stem.Length() - 5);
    TString outName = "MiscEdE_" + stem + ".root";
    TFile *out = new TFile(outName, "recreate");

    // ---------------------------------------------------------
    // Load the 17F and 16O TCutGs from file
    // ---------------------------------------------------------
    TCutG *cut17F = nullptr;
    TString path17F = gSystem->ExpandPathName(Form("%s/17FCut.root", cutDir.Data()));
    if (!gSystem->AccessPathName(path17F))
    {
        TFile *fCut = TFile::Open(path17F);
        if (fCut && !fCut->IsZombie())
        {
            for (auto keyObj : *fCut->GetListOfKeys())
            {
                TKey *key = (TKey *)keyObj;
                if (TString(key->GetClassName()) == "TCutG")
                {
                    cut17F = (TCutG *)key->ReadObj()->Clone("cut17F");
                    break;
                }
            }
            fCut->Close();
        }
    }

    TCutG *cut16O = nullptr;
    TString path16O = gSystem->ExpandPathName(Form("%s/16OCut.root", cutDir.Data()));
    if (!gSystem->AccessPathName(path16O))
    {
        TFile *fCut = TFile::Open(path16O);
        if (fCut && !fCut->IsZombie())
        {
            for (auto keyObj : *fCut->GetListOfKeys())
            {
                TKey *key = (TKey *)keyObj;
                if (TString(key->GetClassName()) == "TCutG")
                {
                    cut16O = (TCutG *)key->ReadObj()->Clone("cut16O");
                    break;
                }
            }
            fCut->Close();
        }
    }

    out->cd();

    // ---------------------------------------------------------
    // Prepare Histograms
    // ---------------------------------------------------------
    std::map<int, TH2D*> hEvsChMap;
    auto getEvsCh = [&](int boardSN) -> TH2D* {
        if (hEvsChMap.count(boardSN)) return hEvsChMap[boardSN];
        TH2D *h = new TH2D(Form("h2_E_vs_ch_bd%d", boardSN), Form("Board %d: E vs channel;board channel;E [ADC]", boardSN), 16, -0.5, 15.5, 800, 0, kAdcMax);
        hEvsChMap[boardSN] = h;
        return h;
    };

    TH2D *h2_EvsdT_LolliSi = new TH2D("EvsdT_Lollipop_Si", "Energy Si vs dT (Lollipop - Si);dT [ticks];Si E [ADC]", 500, 0, 2000, 400, 0, kAdcMax);
    TH2D *h2_TOF_SiRF = new TH2D("TOF_Si_RF", "Time of Flight: Energy Si vs (T_Si - T_RF);T_Si - T_RF [ticks];Si E [ADC]", 500, 0, 2000, 400, 0, kAdcMax);

    // Dynamic Pairwise E-dE mapped by UID
    std::map<std::pair<int, int>, TH2D *> hEdE;
    auto getEdE = [&](int uidA, int uidB) -> TH2D * {
        if (uidA > uidB) std::swap(uidA, uidB);
        auto key = std::make_pair(uidA, uidB);
        auto it = hEdE.find(key);
        if (it != hEdE.end()) return it->second;
        
        TString name = Form("h2_ede_bd%d_ch%d_vs_bd%d_ch%d", uidA/100, uidA%100, uidB/100, uidB%100);
        TString title = Form("E-dE: Bd %d Ch %d (x) vs Bd %d Ch %d (y);X E [ADC];Y E [ADC]", uidA/100, uidA%100, uidB/100, uidB%100);
        
        TH2D *h = new TH2D(name, title, 400, 0, kAdcMax, 400, 0, kAdcMax);
        hEdE[key] = h;
        return h;
    };

    // 17F Gated Plots
    TH2D *hEdE_17F = nullptr;
    TH2D *hEdE_inv = nullptr;
    TH2D *h2_EvsdT_LolliSi_17F = nullptr;
    TH2D *h2_TOF_SiRF_17F = nullptr;
    TH1D *h1_dT_17F = nullptr;

    if (cut17F)
    {
        hEdE_17F = new TH2D("EdE_lollipopIC_vs_SiMonitor_17FCut", "lollipop IC vs Si monitor (17F Gated);Si monitor E [ADC];lollipop IC E [ADC]", 400, 0, kAdcMax, 400, 0, kAdcMax);
        hEdE_inv = new TH2D("EdE_lollipopIC_vs_SiMonitor_17FCut_Inverse", "lollipop IC vs Si monitor (Inverse 17F Gated);Si monitor E [ADC];lollipop IC E [ADC]", 400, 0, kAdcMax, 400, 0, kAdcMax);
        h2_EvsdT_LolliSi_17F = new TH2D("EvsdT_Lollipop_Si_17FCut", "Energy Si vs dT (17F Gated);dT (Lollipop - Si) [ticks];Si E [ADC]", 500, 0, 2000, 400, 0, kAdcMax);
        h2_TOF_SiRF_17F = new TH2D("TOF_Si_RF_17FCut", "Si TOF vs Energy (17F Gated);T_Si - T_RF [ticks];Si E [ADC]", 500, 0, 2000, 400, 0, kAdcMax);
        h1_dT_17F = new TH1D("dT_Lollipop_Si_17FCut", "dT Lollipop-Si (17F Gated);T_Lollipop - T_Si [ticks]", 500, 0, 2000);
    }

    // 16O Gated Plots
    TH2D *hEdE_16O = nullptr;
    TH2D *h2_EvsdT_LolliSi_16O = nullptr;
    TH2D *h2_TOF_SiRF_16O = nullptr;
    TH1D *h1_dT_16O = nullptr;

    if (cut16O)
    {
        hEdE_16O = new TH2D("EdE_lollipopIC_vs_SiMonitor_16OCut", "lollipop IC vs Si monitor (16O Gated);Si monitor E [ADC];lollipop IC E [ADC]", 400, 0, kAdcMax, 400, 0, kAdcMax);
        h2_EvsdT_LolliSi_16O = new TH2D("EvsdT_Lollipop_Si_16OCut", "Energy Si vs dT (16O Gated);dT (Lollipop - Si) [ticks];Si E [ADC]", 500, 0, 2000, 400, 0, kAdcMax);
        h2_TOF_SiRF_16O = new TH2D("TOF_Si_RF_16OCut", "Si TOF vs Energy (16O Gated);T_Si - T_RF [ticks];Si E [ADC]", 500, 0, 2000, 800, 0, kAdcMax);
        h1_dT_16O = new TH1D("dT_Lollipop_Si_16OCut", "dT Lollipop-Si (16O Gated);T_Lollipop - T_Si [ticks]", 500, 0, 2000);
    }

    Long64_t nEnt = tree->GetEntries();
    printf("MiscEdE: %lld entries across %d file(s) -> %s\n", nEnt, tree->GetNtrees(), outName.Data());

    Long64_t used = 0;
    for (Long64_t i = 0; i < nEnt; ++i)
    {
        tree->GetEntry(i);

        // We map hits by UID instead of standard channel to prevent overlaps between Bd405 and Bd89
        std::map<int, HitData> hit;
        for (unsigned int j = 0; j < multi && j < (unsigned)MAX_MULTI; ++j)
        {
            // Only care about MISC boards 405 and 89
            if (sn[j] != 405 && sn[j] != 89)
                continue;

            int c = ch[j];
            int uid = (sn[j] * 100) + ch[j]; 
            double en = e[j];
            unsigned long long ts = e_t[j];

            auto it = hit.find(uid);
            if (it == hit.end() || en > it->second.e)
                hit[uid] = {en, ts};
                
            getEvsCh(sn[j])->Fill(c, en);
        }

        if (hit.empty())
            continue;
        ++used;

        // ---------------------------------------------------------
        // Generic Dynamic Pairwise E-dE Logic across all UIDs
        // ---------------------------------------------------------
        for (auto a = hit.begin(); a != hit.end(); ++a)
        {
            for (auto b = std::next(a); b != hit.end(); ++b)
            {
                if (a->second.e > 500 && b->second.e > 500)
                {
                    getEdE(a->first, b->first)->Fill(a->second.e, b->second.e);
                }
            }
        }

        // ---------------------------------------------------------
        // Core Timing and PID Gating Logic (Lollipop-Si specific)
        // ---------------------------------------------------------
        if (hit.count(uidLollipop) && hit.count(uidSiMon))
        {
            double e_lol = hit[uidLollipop].e;
            double e_si = hit[uidSiMon].e;

            // Calculate delta T safely across potential CAEN rollovers
            double dt_lol_si = static_cast<double>(hit[uidLollipop].t) - static_cast<double>(hit[uidSiMon].t);

            if (e_lol > 500 && e_si > 500)
            {
                h2_EvsdT_LolliSi->Fill(dt_lol_si, e_si);

                bool is17F = cut17F && cut17F->IsInside(e_si, e_lol);
                bool is16O = cut16O && cut16O->IsInside(e_si, e_lol);

                if (is17F)
                {
                    hEdE_17F->Fill(e_si, e_lol);
                    h2_EvsdT_LolliSi_17F->Fill(dt_lol_si, e_si);
                    h1_dT_17F->Fill(dt_lol_si);
                }
                else if (cut17F)
                {
                    hEdE_inv->Fill(e_si, e_lol);
                }

                if (is16O)
                {
                    hEdE_16O->Fill(e_si, e_lol);
                    h2_EvsdT_LolliSi_16O->Fill(dt_lol_si, e_si);
                    h1_dT_16O->Fill(dt_lol_si);
                }

                // Generate absolute TOF if the RF channel fired
                if (hit.count(uidRF))
                {
                    double tof = static_cast<double>(hit[uidSiMon].t) - static_cast<double>(hit[uidRF].t);
                    h2_TOF_SiRF->Fill(tof, e_si);

                    if (is17F)
                        h2_TOF_SiRF_17F->Fill(tof, e_si);
                    if (is16O)
                        h2_TOF_SiRF_16O->Fill(tof, e_si);
                }
            }
        }
    }

    // ---------------------------------------------------------
    // Purge Sparse E-dE Histograms (< 500 counts)
    // ---------------------------------------------------------
    for (auto it = hEdE.begin(); it != hEdE.end();)
    {
        if (it->second->GetEntries() < 500)
        {
            delete it->second; // Free memory and remove from ROOT file
            it = hEdE.erase(it);
        }
        else
        {
            ++it;
        }
    }

    // Convenience alias: Transpose generic plot to correctly place stopping E on X
    if (hEdE.count({uidLollipop, uidSiMon}))
    {
        TH2D *src = hEdE[{uidLollipop, uidSiMon}];
        TH2D *alias = new TH2D("EdE_lollipopIC_vs_SiMonitor_Ungated",
                               Form("lollipop IC vs Si monitor (Ungated);Bd %d Ch %d (Si) E [ADC];Bd %d Ch %d (Lolli) E [ADC]", snSi, chSi, snLolli, chLolli),
                               src->GetNbinsY(), src->GetYaxis()->GetXmin(), src->GetYaxis()->GetXmax(),
                               src->GetNbinsX(), src->GetXaxis()->GetXmin(), src->GetXaxis()->GetXmax());
        for (int ix = 1; ix <= src->GetNbinsX(); ++ix)
            for (int iy = 1; iy <= src->GetNbinsY(); ++iy)
                alias->SetBinContent(iy, ix, src->GetBinContent(ix, iy)); // (x,y) -> (y,x)
    }

    out->Write();
    out->Close();
    printf("MiscEdE: Processed %lld valid MISC events. Wrote %s\n", used, outName.Data());
}