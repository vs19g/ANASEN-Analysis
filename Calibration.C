
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
#include "Armory/HistPlotter.h"
#include "TVector3.h"
#include "Calibration.h"

TH2F *hQQQFVB;
HistPlotter *plotter;
int padID = 0;

TCutG *cut;
std::map<std::tuple<int, int, int>, std::vector<std::pair<double, double>>> dataPoints;

bool qqqEcut = false;

// Gain Arrays
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
double qqqGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

void Calibration::Begin(TTree * /*tree*/)
{
    plotter = new HistPlotter("Calib.root", "TFILE");
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
    for (int det = 0; det < MAX_QQQ; det++)
    {
        for (int ring = 0; ring < MAX_RING; ring++)
        {
            for (int wedge = 0; wedge < MAX_WEDGE; wedge++)
            {
                TString hname = Form("hCal_qqq%d_ring%d_wedge%d", det, ring, wedge);
                TString htitle = Form("QQQ det%d ring%d wedge%d; Energy (arb); Counts", det, ring, wedge);
                // hQQQSpectra[det][ring][wedge] = new TH1F(hname, htitle, 4000, 0, 16000);
            }
        }
    }
}

Bool_t Calibration::Process(Long64_t entry)
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
            if (qqq.e[i] > 100)
                qqqEcut = true;
            if (qqq.id[i] == qqq.id[j])
            {
                int chWedge = -1;
                int chRing = -1;
                float eWedgeRaw = 0.0;
                float eWedge = 0.0;
                float eRingRaw = 0.0;
                float eRing = 0.0;
                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
                {
                    chWedge = qqq.ch[i];
                    eWedgeRaw = qqq.e[i];
                    eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    printf("Wedge E: %.2f  Gain: %.4f \n", eWedge, qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16]);
                    chRing = qqq.ch[j] - 16;
                    eRingRaw = qqq.e[j];
                    eRing = qqq.e[j]; //*qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i]-16];
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
                    eWedgeRaw = qqq.e[j];

                    chRing = qqq.ch[i] - 16;
                    eRing = qqq.e[i];// * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    eRingRaw = qqq.e[i];
                }
                else
                    continue;

                // hQQQFVB->Fill(eWedge, eRing);
                plotter->Fill2D(Form("hRaw_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge),400,0,16000,400,0,16000, eWedgeRaw, eRingRaw,"ERaw");
                plotter->Fill2D(Form("hGM_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge),400,0,16000,400,0,16000, eWedge, eRing,"EGM");
                plotter->Fill2D("hRawQQQ", 4000, 0, 16000, 4000, 0, 16000, eWedgeRaw, eRingRaw);
                plotter->Fill2D("hGMQQQ", 4000, 0, 16000, 4000, 0, 16000, eWedge, eRing);

                TString histName = Form("hQQQFVB_id%d_r%d_w%d", qqq.id[i], chRing, chWedge);
                // TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                // if (!hist2d)
                // {
                //     hist2d = new TH2F(histName, Form("QQQ Det%d R%d W%d;Wedge E;Ring E", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000);
                // }

                // hist2d->Fill(eWedge, eRing);
                // if (cut && cut->IsInside(eWedge, eRing))
                {
                    // Accumulate data for gain matching
                    dataPoints[{qqq.id[i], chRing, chWedge}].emplace_back(eWedge, eRing);
                }
            }
        }
    }

    return kTRUE;
}

void Calibration::Terminate()
{
    const double AM241_PEAK = 5485.56; // keV
    double calibArray[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
    bool calibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

    std::ofstream outFile("qqq_Calib.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

   // ======================= Loop over all channels =======================
for (auto &kv : dataPoints) {

    int det, ring, wedge;
    std::tie(det, ring, wedge) = kv.first;

    const std::vector<std::pair<double,double>> &pts = kv.second;

    if (pts.size() < 5)
        continue;

    // Build TGraph from stored (wedgeGM, ringE)
    std::vector<double> wedgeGM, ringE;
    wedgeGM.reserve(pts.size());
    ringE.reserve(pts.size());

    for (auto &p : pts) {
        wedgeGM.push_back(p.first);   // gain-matched wedge energy (ADC)
        ringE.push_back(p.second);    // ring energy (ADC)
    }

    TGraph g(pts.size(), wedgeGM.data(), ringE.data());
    g.SetTitle(Form("QQQ Det%d Ring%d Wedge%d", det, ring, wedge));

    // Fit a line through origin: E_ring = a * E_wedge
    TF1 f("f","[0]*x",0,16000);
    g.Fit(&f,"QNR");     // Quiet, No draw, use Range

    double slope_raw = f.GetParameter(0);

    if (slope_raw <= 0)
        continue;

    // Convert raw slope into keV calibration:
    // Use the Am241 peak expected position:
    //   E_keV = ADC * slope_keV
    double slope_keV = AM241_PEAK / (AM241_PEAK / slope_raw); 
    // Simplifies to:
    // slope_keV = slope_raw; // slope now directly converts ADC → keV

    calibArray[det][ring][wedge] = slope_keV;
    calibValid[det][ring][wedge] = true;

    outFile << det << " " << ring << " " << wedge << " "
            << slope_keV << "\n";

    printf("Calib Det=%d Ring=%d Wedge=%d  slope=%.5f\n",
           det, ring, wedge, slope_keV);
}

    outFile.close();
    std::cout << "Gain matching complete." << std::endl;

    // === Plot all gain-matched QQQ points together with a 2D histogram ===
    TH2F *hAll = new TH2F("hAll", "All QQQ Gain-Matched;Corrected Wedge E;Ring E",
                          800, 0, 16000, 800, 0, 16000);

    // Fill the combined TH2F with corrected data
    for (auto &kv : dataPoints)
    {
        int id, ring, wedge;
        std::tie(id, ring, wedge) = kv.first;
        if (!calibValid[id][ring][wedge])
            continue;
        auto &pts = kv.second;
        for (auto &pr : pts)
        {
            double corrWedge = pr.first * calibArray[id][ring][wedge];
            double corrRing = pr.second * calibArray[id][ring][wedge];
            hAll->Fill(corrWedge, corrRing);
            plotter->Fill2D("hAll", 4000, 0, 16000, 4000, 0, 16000, corrWedge, corrRing); // Create the histogram in the plotter
        }
    }
    plotter->FlushToDisk();
}
