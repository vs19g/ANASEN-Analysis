#define GainMatchQQQ_cxx

#include "GainMatchQQQ.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <fstream>
#include <utility>
#include <algorithm>
#include "Armory/HistPlotter.h"
#include "TVector3.h"

TH2F *hQQQFVB;
HistPlotter *plotter;

int padID = 0;

TCutG *cut;
std::map<std::tuple<int, int, int>, std::vector<std::pair<double, double>>> dataPoints;

void GainMatchQQQ::Begin(TTree * /*tree*/)
{
    plotter = new HistPlotter("GainQQQ.root", "TFILE");
    TString option = GetOption();

    hQQQFVB = new TH2F("hQQQFVB", "QQQ Front vs Back; Front E; Back E", 800, 0, 16000, 800, 0, 16000);


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

Bool_t GainMatchQQQ::Process(Long64_t entry)
{
    b_qqqMulti->GetEntry(entry);
    b_qqqID->GetEntry(entry);
    b_qqqCh->GetEntry(entry);
    b_qqqE->GetEntry(entry);
    b_qqqT->GetEntry(entry);

    qqq.CalIndex();

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

                // TString histName = Form("hQQQFVB_id%d_r%d_w%d", qqq.id[i], chRing, chWedge);
                // TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                // if (!hist2d)
                // {
                //     hist2d = new TH2F(histName, Form("QQQ Det%d R%d W%d;Wedge E;Ring E", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000);
                // }

                // hist2d->Fill(eWedge, eRing);
                // if (cut && cut->IsInside(eWedge, eRing))
                double ratio = eRing / eWedge;
                double maxslope=1.5;
                //gate gets rid of noisy off diagonal events forming  a 'V' about the center
                //TODO: These are very likely nearest-neighbor charge-sharing events, that will go away if appropriately summed

                // if(ratio < maxslope && ratio > 1./maxslope || eWedge<200 || eRing<200) //method adopted from Sudarshan's approach
                bool validPoint = false;
                if(ratio < maxslope && ratio > 1./maxslope)// || eWedge<200 || eRing<200) //method adopted from Sudarshan's approach
                {
                    // Accumulate data for gain matching
                    dataPoints[{qqq.id[i], chRing, chWedge}].emplace_back(eWedge, eRing);
                    plotter->Fill2D("hAll_in", 4000, 0, 16000, 4000, 0, 16000, eWedge, eRing);
                    validPoint = true;
                }
                
                if(!validPoint){
                    plotter->Fill2D("hAll_out", 4000, 0, 16000, 4000, 0, 16000, eWedge, eRing);
                }
            }
        }
    }

    return kTRUE;
}

void GainMatchQQQ::Terminate()
{
    const int MAX_DET = 4;
    const int MAX_RING = 16;
    const int MAX_WEDGE = 16;

    double gainArray[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    bool gainValid[MAX_DET][MAX_RING][MAX_WEDGE] = {{{false}}};

    std::ofstream outFile("qqq_GainMatch.txt");
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
                          4000, 0, 16000, 4000, 0, 16000);

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
    
    plotter->FlushToDisk();
}
