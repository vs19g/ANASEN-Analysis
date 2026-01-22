
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
double qqqwGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
// double qqqrGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqwGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
// bool qqqrGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

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
            double gainw, gainr;
            while (infile >> det >> ring >> wedge >> gainw >> gainr)
            {
                qqqwGain[det][ring][wedge] = gainw;
                // qqqrGain[det][ring][wedge] = gainr;
                qqqwGainValid[det][ring][wedge] = (gainw > 0);
                // qqqrGainValid[det][ring][wedge] = (gainr > 0);
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
                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && /*qqqrGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16] &&*/ qqqwGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
                {
                    chWedge = qqq.ch[i];
                    eWedgeRaw = qqq.e[i];
                    eWedge = qqq.e[i] * qqqwGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    // printf("Wedge E: %.2f  Gain: %.4f \n", eWedge, qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16]);
                    chRing = qqq.ch[j] - 16;
                    eRingRaw = qqq.e[j];
                    eRing = qqq.e[j];// * qqqrGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && /*qqqrGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16] &&*/ qqqwGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j] * qqqwGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
                    eWedgeRaw = qqq.e[j];

                    chRing = qqq.ch[i] - 16;
                    eRing = qqq.e[i];// * qqqrGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    eRingRaw = qqq.e[i];
                }
                else
                    continue;

                // hQQQFVB->Fill(eWedge, eRing);
                plotter->Fill2D(Form("hRaw_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000, eWedgeRaw, eRingRaw, "ERaw");
                plotter->Fill2D(Form("hGM_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000, eWedge, eRing, "EGM");
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
                const double MIN_ADC = 1500.0;
                const double MAX_ADC = 3000.0;

                // if (eWedge >= MIN_ADC && eWedge <= MAX_ADC &&
                //     eRing >= MIN_ADC && eRing <= MAX_ADC)
                double ratio = (eWedge > 0.0) ? (eRing / eWedge) : 0.0;

                double maxslope = 1.5;

                bool validPoint = false;
                if (ratio < maxslope && ratio > 1. / maxslope)
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
    const double AM241_PEAK = 5485.56;
    const double P_PEAK = 7000; // keV

    double calibArray[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
    bool calibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

    std::ofstream outFile("qqq_Calib.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening qqq_Calib.txt!" << std::endl;
        return;
    }

    //----------------------------------------------------------------------
    // 1. Create per–channel 1D spectra in ADC from stored gain-matched data
    //----------------------------------------------------------------------

    std::map<std::tuple<int, int, int>, TH1F *> spectra;

    for (auto &kv : dataPoints)
    {
        int det, ring, wedge;
        std::tie(det, ring, wedge) = kv.first;

        TString hname = Form("hSpec_d%d_r%d_w%d", det, ring, wedge);
        TH1F *h = new TH1F(hname, hname, 4000, 0, 16000);

        for (auto &p : kv.second)
        {
            double eWedge = p.first; // already gain-matched ADC
            double eRing = p.second;

            // Use ring ADC for calibration (cleaner alpha peak)
            h->Fill(eRing);
        }

        spectra[kv.first] = h;
    }

    //----------------------------------------------------------------------
    // 2. Fit Am-241 peak and extract keV/ADC calibration slope
    //----------------------------------------------------------------------

    for (auto &kv : spectra)
    {
        int det, ring, wedge;
        std::tie(det, ring, wedge) = kv.first;
        TH1F *h = kv.second;

        if (!h || h->GetEntries() < 50)
            continue;

        int binMax = h->GetMaximumBin();
        double adcPeak = h->GetXaxis()->GetBinCenter(binMax);

        if (adcPeak <= 0)
            continue;

        // double slope_keV = AM241_PEAK / adcPeak; // keV per ADC
        double slope_keV = P_PEAK / adcPeak; // keV per ADC

        calibArray[det][ring][wedge] = slope_keV;
        calibValid[det][ring][wedge] = true;

        outFile << det << " " << wedge << " " << ring << " "
                << slope_keV << "\n";

        // printf("QQQ DET=%d R=%d W=%d  ADCpeak=%.1f  slope_keV=%.6f\n",det, ring, wedge, adcPeak, slope_keV);
    }

    outFile.close();
    std::cout << "Wrote QQQ calibration file qqq_Calib.txt\n";

    //----------------------------------------------------------------------
    // 3. Build fully calibrated 2D combined histogram
    //----------------------------------------------------------------------

    TH2F *hCal = new TH2F("hCal",
                          "All QQQ Calibrated;Wedge Energy (keV);Ring Energy (keV)",
                          800, 0, 7000,
                          800, 0, 7000);

    for (auto &kv : dataPoints)
    {
        int det, ring, wedge;
        std::tie(det, ring, wedge) = kv.first;

        if (!calibValid[det][ring][wedge])
            continue;

        double slope = calibArray[det][ring][wedge];

        for (auto &p : kv.second)
        {
            double eWGM = p.first; // gain matched ADC
            double eRGM = p.second;

            double eWkeV = eWGM * slope / 1000;
            double eRkeV = eRGM * slope / 1000;

            hCal->Fill(eWkeV, eRkeV);
            plotter->Fill2D("hCalQQQ", 4000, 0, 10, 4000, 0, 10, eWkeV, eRkeV);
            plotter->Fill2D(Form("hRCal_qqq%d", det), 16, 0, 15, 400, 0, 24, ring, eRkeV, "RingCal");
            plotter->Fill2D(Form("hWCal_qqq%d", det), 16, 0, 15, 400, 0, 24, wedge, eWkeV, "WedgeCal");
        }
    }

    plotter->FlushToDisk();
    std::cout << "Calibrated 2D QQQ histogram saved.\n";
}
