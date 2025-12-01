
#define QQQ_Calcheck_cxx

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
#include "QQQ_Calcheck.h"

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
double qqqrGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqwGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
bool qqqrGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
double qqqCalib[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqCalibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

void QQQ_Calcheck::Begin(TTree * /*tree*/)
{
    plotter = new HistPlotter("Cal_checkQQQ.root", "TFILE");
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
            double gainw,gainr;
            while (infile >> det >> ring >> wedge >> gainw>> gainr)
            {
                qqqwGain[det][ring][wedge] = gainw;
                qqqrGain[det][ring][wedge] = gainr;
                qqqwGainValid[det][ring][wedge] = (gainw > 0);
                qqqrGainValid[det][ring][wedge] = (gainr > 0);
            }
            infile.close();
            std::cout << "Loaded QQQ gains from " << filename << std::endl;
        }
    }
    // ----------------------- Load QQQ Calibrations
    {
        std::string filename = "qqq_Calib.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int det, ring, wedge;
            double slope;
            while (infile >> det >> ring >> wedge >> slope)
            {
                qqqCalib[det][ring][wedge] = slope;
                qqqCalibValid[det][ring][wedge] = (slope > 0);
            }
            infile.close();
            std::cout << "Loaded QQQ calibrations from " << filename << std::endl;
        }
    }
}

Bool_t QQQ_Calcheck::Process(Long64_t entry)
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
                float eWedgeMeV = 0.0;
                float eRingRaw = 0.0;
                float eRing = 0.0;
                float eRingMeV = 0.0;
                // plug in gains
                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqrGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16] && qqqwGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
                {
                    chWedge = qqq.ch[i];
                    eWedgeRaw = qqq.e[i];
                    eWedge = qqq.e[i] * qqqwGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    // printf("Wedge E: %.2f  Gain: %.4f \n", eWedge, qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16]);
                    chRing = qqq.ch[j] - 16;
                    eRingRaw = qqq.e[j];
                    eRing = qqq.e[j] * qqqrGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i]-16];
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqrGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16] && qqqwGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j] * qqqwGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
                    eWedgeRaw = qqq.e[j];

                    chRing = qqq.ch[i] - 16;
                    eRing = qqq.e[i] * qqqrGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    eRingRaw = qqq.e[i];
                }
                else
                    continue;
                // plug in calibrations
                if (qqqCalibValid[qqq.id[i]][chRing][chWedge])
                {
                    eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chRing][chWedge] / 1000;
                    eRingMeV = eRing * qqqCalib[qqq.id[i]][chRing][chWedge] / 1000;
                }
                else
                    continue;

                // hQQQFVB->Fill(eWedge, eRing);
                plotter->Fill2D(Form("hRaw_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 400, 0, 8000, 400, 0, 8000, eWedgeRaw, eRingRaw, "ERaw");
                plotter->Fill2D(Form("hGM_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 400, 0, 16000, 400, 0, 16000, eWedge, eRing, "EGM");
                plotter->Fill2D(Form("hCal_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 400, 0, 10, 400, 0, 10, eWedgeMeV, eRingMeV, "ECal");
                if(eWedgeRaw >1500 && eRingRaw>1500 )
                plotter->Fill2D(Form("hCal_cut_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 400, 0, 10, 400, 0, 10, eWedgeMeV, eRingMeV, "ECal_cut");
                
                plotter->Fill2D(Form("hRCal_qqq%d", qqq.id[i]), 16, 0, 15, 1000, 0, 30, chRing, eRingMeV, "RingCal");
                plotter->Fill2D(Form("hWCal_qqq%d", qqq.id[i]), 16, 0, 15, 1000, 0, 30, chWedge, eWedgeMeV, "WedgeCal");
                plotter->Fill2D("hRawQQQ", 4000, 0, 8000, 4000, 0, 8000, eWedgeRaw, eRingRaw);
                plotter->Fill2D("hGMQQQ", 4000, 0, 16000, 4000, 0, 16000, eWedge, eRing);
                plotter->Fill2D("hCalQQQ", 4000, 0, 10, 4000, 0, 10, eWedgeMeV, eRingMeV);
            }
        }
    }

    return kTRUE;
}

void QQQ_Calcheck::Terminate()
{
    plotter->FlushToDisk();
    std::cout << "Calibration check file for 2D QQQ histogram saved.\n";
}
