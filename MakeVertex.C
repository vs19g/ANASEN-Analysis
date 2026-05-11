#define MakeVertex_cxx

Int_t colors[40] = {
    kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange,
    kSpring, kTeal, kAzure, kViolet, kPink, kGray, kWhite,
    kRed + 2, kGreen + 2, kBlue + 2, kYellow + 2, kMagenta + 2, kCyan + 2, kOrange + 2,
    kSpring + 2, kTeal + 2, kAzure + 2, kViolet + 2, kPink + 2,
    kRed - 7, kGreen - 7, kBlue - 7, kYellow - 7, kMagenta - 7, kCyan - 7, kOrange - 7,
    kSpring - 7, kTeal - 7, kAzure - 7, kViolet - 7, kPink - 7, kGray + 2};

#include "MakeVertex.h"
#include "Armory/ClassPW.h"
#include "Armory/HistPlotter.h"
#include "Armory/SX3Geom.h"
#include "Armory/PC_StepLadder_Correction.h"
#include "Armory/Kinematics.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TBranch.h>
#include <TGraph2D.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TH3D.h>
#include <TVector3.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <algorithm>

bool realtime = true;
bool process_alpha_proton_scattering = true;
double source_vertex = 53; // 53
const double qqq_z = 100.0;
const double anode_gain = 1.5146e-5; // channels --> MeV
std::string dataset;

TF1 pcfix_func("func", model_invert, -200, 200);
TGraph *MeV_to_cm = NULL, *cm_to_MeV = NULL;
TGraph *MeV_to_cm_p = NULL, *cm_to_MeVp = NULL;

TApplication *app = NULL;
TH1F *hha = NULL, *hhc = NULL;
TH3D *frame = NULL;
TCanvas *can1 = NULL, *can2 = NULL;

TPolyLine3D *pla[24] = {NULL};
TPolyLine3D *plc[24] = {NULL};
TPolyLine3D *qqqw[16][4] = {NULL};
TPolyLine3D *trajectory = NULL;
TGraph2D *qqqg = NULL, *crossoverg = NULL, *guessg = NULL;

double z_to_crossover_rho(double z)
{
    return 9.20645e-5 * z * z + 34.1973;
}

double z_to_crossover_rho_cathode(double z)
{
    return 9.20645e-5 * z * z + 34.1973;
}

// Global instances
PW pw_contr;
PW pwinstance;
TVector3 hitPos;
double qqqenergy, qqqtimestamp;
class Event
{
public:
    Event(TVector3 p, double e1, double e2, double t1, double t2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2) {}
    Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2) {}
    // Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2, int m1, int m2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2), multi1(m1), multi2(m2) {}
    Event(TVector3 p, double e1, double e2, double t1, double t2, int a, int c, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), Anodech(a), Cathodech(c), ch1(c1), ch2(c2) {}

    TVector3 pos;
    int ch1 = -1;        // int(ch1/16) gives qqq id, ch1%16 gives ring#
    int ch2 = -1;        // int(ch2/16) gives qqq id, ch2%16 gives wedge#
    double Energy1 = -1; // Front for QQQ, Anode for PC
    double Energy2 = -1; // Back for QQQ, Cathode for PC
    double Time1 = -1;
    double Time2 = -1;
    int Anodech = -1;
    int Cathodech = -1;

    // misc elements;
    int multi1 = -1, multi2 = -1;
};

// Calibration globals
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
double qqqGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
double qqqCalib[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqCalibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

double sx3BackGain[24][4][4] = {{{1.}}};
double sx3FrontGain[24][4] = {{1.}};
double sx3FrontOffset[24][4] = {{0.}};
double sx3RightGain[24][4] = {{1.}};

// PC Arrays
double pcSlope[48];
double pcIntercept[48];

HistPlotter *plotter;

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

void protonAlphaHistograms(HistPlotter *plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);

void MakeVertex::Begin(TTree * /*tree*/)
{
    pcfix_func.SetNpx(100000);
    TString option = GetOption();
    if (option != "")
        plotter = new HistPlotter(option.Data(), "TFILE");
    else
        plotter = new HistPlotter("Analyzer_SX3.root", "TFILE");

    pw_contr.ConstructGeo();
    pwinstance.ConstructGeo();
    if (gROOT->IsBatch())
        realtime = false;

    // ---------------------------------------------------------
    // 1. CRITICAL FIX: Initialize PC Arrays to Default (Raw)
    // ---------------------------------------------------------
    for (int i = 0; i < 48; i++)
    {
        pcSlope[i] = 1.0;     // Default slope = 1 (preserves Raw energy)
        pcIntercept[i] = 0.0; // Default intercept = 0
    }

    if (getenv("DATASET"))
        dataset = std::string(getenv("DATASET"));
    if (getenv("source_vertex"))
        source_vertex = (double)std::atof(std::string(getenv("source_vertex")).c_str());
    std::cout << "Dataset set to  " << dataset << std::endl;
    std::cout << "source_vertex set to  " << source_vertex << std::endl;

    if (getenv("flipa"))
    {
        int flip_offset = std::atoi(getenv("anode_offset"));
        int yes_to_flip = std::atoi(getenv("flipa"));
        if (yes_to_flip && flip_offset)
        {
            std::cout << "Flipping anodes and offseting by " << flip_offset << " wires." << std::endl;
        }
        else if (flip_offset)
        {
            std::cout << "Offseting anodes without flip by " << flip_offset << " wires." << std::endl;
        }
    }

    fflush(stdout);
    // usleep(4e5);
    //  Load PC Calibrations
    std::ifstream inputFile("slope_intercept_results_" + dataset + ".dat");
    if (inputFile.is_open())
    {
        std::string line;
        int index;
        double slope, intercept;
        while (std::getline(inputFile, line))
        {
            std::stringstream ss(line);
            ss >> index >> slope >> intercept;
            if (index >= 0 && index <= 47)
            {
                pcSlope[index] = slope;
                pcIntercept[index] = intercept;
            }
        }
        inputFile.close();
    }
    else
    {
        std::cerr << "Error opening slope_intercept.dat" << std::endl;
    }

    // ... (Load QQQ Gains and Calibs - same as before) ...
    {
        std::string filename = "qqq_GainMatch.dat";
        std::ifstream infile(filename);
        if (infile.is_open())
        {
            int det, ring, wedge;
            double gainw, gainr;
            while (infile >> det >> wedge >> ring >> gainw >> gainr)
            {
                qqqGain[det][wedge][ring] = gainw;
                qqqGainValid[det][wedge][ring] = (gainw > 0);
                // std::cout << "QQQ Gain Loaded: Det " << det << " Ring " << ring << " Wedge " << wedge << " GainW " << gainw << " GainR " << gainr << std::endl;
            }
            infile.close();
        }
    }

    {
        std::string filename = "qqq_Calib.dat";
        std::ifstream infile(filename);
        if (infile.is_open())
        {
            int det, ring, wedge;
            double slope;
            while (infile >> det >> wedge >> ring >> slope)
            {
                qqqCalib[det][wedge][ring] = slope;
                qqqCalibValid[det][wedge][ring] = (slope > 0);
                // std::cout << "QQQ Calib Loaded: Det " << det << " Ring " << ring << " Wedge " << wedge << " Slope " << slope << std::endl;
            }
            infile.close();
        }
    }

    {
        std::ifstream infile("sx3cal/" + dataset + "/backgains.dat");
        std::string temp;
        int backpos, frontpos, clkpos;
        if (infile.is_open())
            while (infile >> clkpos >> temp >> frontpos >> temp >> backpos >> sx3BackGain[clkpos][frontpos][backpos])
                ; // std::cout << sx3BackGain[clkpos][frontpos][backpos] << std::endl;
        infile.close();

        infile.open("sx3cal/" + dataset + "/frontgains.dat");
        if (infile.is_open())
            while (infile >> clkpos >> temp >> temp >> frontpos >> sx3FrontOffset[clkpos][frontpos] >> sx3FrontGain[clkpos][frontpos])
                ; // std::cout << sx3FrontOffset[clkpos][frontpos] << " " << sx3FrontGain[clkpos][frontpos] << std::endl;
        infile.close();

        infile.open("sx3cal/" + dataset + "/rightgains.dat");
        if (infile.is_open())
            while (infile >> clkpos >> frontpos >> temp >> sx3RightGain[clkpos][frontpos])
            {
                sx3RightGain[clkpos][frontpos] = TMath::Abs(sx3RightGain[clkpos][frontpos]);
            }
        infile.close();
    }
    //    MeV_to_cm = new TGraph("eloss_calculations/alphas_in_250torr_mix_filtered_6MeV.txt","%lf %*lf %lf");
    MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_20MeV.dat", "%lf %*lf %lf");
    cm_to_MeV = new TGraph(MeV_to_cm->GetN(), MeV_to_cm->GetY(), MeV_to_cm->GetX());

    MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV.dat", "%lf %*lf %lf");
    cm_to_MeVp = new TGraph(MeV_to_cm_p->GetN(), MeV_to_cm_p->GetY(), MeV_to_cm_p->GetX());

    // cm_to_MeV.Eval(MeV_to_cm.Eval(detectedE)-PathLength) gives energy of particle before it traversed 'path length'

    if (realtime)
    {
        can1 = new TCanvas("wireindex", "c1", 0, 0, 640, 480);
        can2 = new TCanvas("3d", "c2", 650, 0, 640, 480);
        can1->cd();
        // can2->SetFillColor(30);
        frame = new TH3D("frame", "frame", 1000, -100, 100, 1000, -100, 100, 1000, -200, 200);
        hha = new TH1F("hha", "Anode Ecal vs wire#", 48, -12, 36);
        hhc = new TH1F("hhc", "Cathode Ecal vs wire#", 48, -12, 36);
        hha->SetLineColor(kRed);
        hha->GetYaxis()->SetRangeUser(0, 16384);
        hha->GetXaxis()->SetTitle("press any key, interrupt/refresh or double click to continue..");
        hha->Draw();
        hhc->Draw("SAME");
        can1->Modified();
        can1->Update();
        can1->BuildLegend();
        can2->cd();
        frame->Draw();
        for (int i = 0; i < 24; i++)
        {
            plc[i] = new TPolyLine3D(2);
            pla[i] = new TPolyLine3D(2);
            pla[i]->SetPoint(0, pwinstance.An[i].first.X(), pwinstance.An[i].first.Y(), pwinstance.An[i].first.Z());
            pla[i]->SetPoint(1, pwinstance.An[i].second.X(), pwinstance.An[i].second.Y(), pwinstance.An[i].second.Z());
            plc[i]->SetPoint(0, pwinstance.Ca[i].first.X(), pwinstance.Ca[i].first.Y(), pwinstance.Ca[i].first.Z());
            plc[i]->SetPoint(1, pwinstance.Ca[i].second.X(), pwinstance.Ca[i].second.Y(), pwinstance.Ca[i].second.Z());
            plc[i]->SetLineStyle(kDotted);
            pla[i]->SetLineStyle(kDotted);
            pla[i]->SetLineWidth(1.);
            plc[i]->SetLineWidth(1.);
            plc[i]->Draw("same");
            pla[i]->Draw("same");
            plc[i]->SetLineColor(colors[i]);
            pla[i]->SetLineColor(colors[i]);
        }
        crossoverg = new TGraph2D(1);
        crossoverg->SetName("crossoverg");
        crossoverg->SetMarkerStyle(20);
        crossoverg->SetMarkerColor(kBlue + 3);
        qqqg = new TGraph2D(1);
        qqqg->SetName("qqqg");
        qqqg->SetMarkerColor(kRed);
        qqqg->SetMarkerStyle(42);

        crossoverg->SetPoint(0, 0, 0, 0);
        qqqg->SetPoint(0, 0, 0, qqq_z);
        crossoverg->Draw("P same");
        qqqg->Draw("P same");

        trajectory = new TPolyLine3D(2);
        trajectory->SetPoint(0, 0, 0, 0);
        trajectory->SetPoint(1, 0, 0, 0);
        trajectory->Draw("same");

        can2->Modified();
        can2->Update();
    }
}

Bool_t MakeVertex::Process(Long64_t entry)
{
    hitPos.Clear();
    qqqenergy = -1;
    qqqtimestamp = -1;
    HitNonZero = false;
    bool qqq1000cut = false;
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

    double timecut_low = getenv("timecut_low") ? std::atof(getenv("timecut_low")) : 0;
    double timecut_high = getenv("timecut_high") ? std::atof(getenv("timecut_high")) : 1e15;

    if (pc.multi > 0)
    {
        for (int i = 0; i < pc.multi; i++)
        {
            if (pc.t[i] * 1e-9 < timecut_high && pc.t[i] * 1e-9 >= timecut_low)
            {
                // good, keep it moving
            }
            else
            {
                return kTRUE;
            }
        }
    }

    sx3.CalIndex();
    qqq.CalIndex();
    pc.CalIndex();

    std::vector<Event> SX3_Events;
    if (sx3.multi > 1)
    {
        std::array<sx3det, 24> Fsx3;
        // std::cout << "-----" << std::endl;
        bool found_upstream_sx3 = 0;
        for (int i = 0; i < sx3.multi; i++)
        {
            int id = sx3.id[i];
            if (id >= 12)
                continue;
            if (sx3.ch[i] >= 8)
            {
                int sx3ch = sx3.ch[i] - 8;
                sx3ch = (sx3ch + 3) % 4;
                if (id >= 12)
                {
                    found_upstream_sx3 = 1;
                    // std::cout << Form("f%d(",id) << sx3ch << "," << sx3.e[i] << ") " << std::flush;
                }
                // if(sx3ch==0 || sx3ch==3) continue;
                double value = sx3.e[i];
                int gch = sx3.id[i] * 4 + (sx3.ch[i] - 8);
                if (id < 12)
                    Fsx3.at(id).fillevent("BACK", sx3ch, value);
                Fsx3.at(id).ts = static_cast<double>(sx3.t[i]);
                plotter->Fill2D("sx3backs_all_raw", 100, 0, 100, 800, 0, 4096, gch, sx3.e[i]);
            }
            else
            {
                int sx3ch = sx3.ch[i] / 2;
                double value = sx3.e[i];
                if (id >= 12)
                {
                    found_upstream_sx3 = 1;
                    // std::cout << Form("b%d(",id) << sx3ch << "," << value << ") " << std::flush;
                }
                if (sx3.ch[i] % 2 == 0)
                {
                    Fsx3.at(id).fillevent("FRONT_L", sx3ch, value * sx3RightGain[id][sx3ch]);
                }
                else
                {
                    Fsx3.at(id).fillevent("FRONT_R", sx3ch, value);
                }
            }
        } // end for (i in sx3.multi)
        // if(found_upstream_sx3) std::cout << std::endl;

        for (int id = 0; id < 24; id++)
        {
            // std::cout << id << " " << Fsx3.at(id).valid_front_chans.size() << " " << Fsx3.at(id).valid_back_chans.size() << std::endl;;
            try
            {
                Fsx3.at(id).validate();
            }
            catch (std::exception exc)
            {
                std::cout << "oops! anyway " << std::endl;
                continue;
            }
            auto det = Fsx3.at(id);
            bool no_charge_sharing_strict = det.valid_front_chans.size() == 1 && det.valid_back_chans.size() == 1;
            if (det.valid)
            {
                // std::cout << det.frontEL << " " << det.frontEL*sx3RightGain[id][det.stripF] << std::endl;
                // plotter->Fill2D("be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_b"+std::to_string(det.stripB),200,-1,1,800,0,8192,det.frontX,det.backE,"evsx");

                plotter->Fill2D("matched_be_vs_x_sx3_id_" + std::to_string(id) + "_f" + std::to_string(det.stripF), 200, -60, 60, 800, 0, 8192,
                                det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx_matched");
                // plotter->Fill2D("fe_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_"+std::to_string(det.stripB),200,-1,1,800,0,4096,det.frontX,det.backE,"evsx");
                // plotter->Fill2D("l_vs_r_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF),800,0,4096,800,0,4096,det.frontEL,det.frontER,"l_vs_r");
            }
            if (det.valid && (id == 9 || id == 7 || id == 1 || id == 3) && det.stripF != DEFAULT_NULL && det.stripB != DEFAULT_NULL)
            {
                double z = det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF];
                z = z + (75.0 / 2.0) - 3.0; // convert local sx3z to detector global coordinate system as indicated by measurements.
                // Note that this will be different for the upstream barrel, when it gets implemented
                double backE = det.backE * sx3BackGain[id][det.stripF][det.stripB];
                // if(backE<2000) continue;
                det.stripF = 3 - det.stripF;
                double beta_n = 15.0 + TMath::ATan2((2 * det.stripF - 3) * 40.30, 8.0 * 88.0 * TMath::Cos(15.0 * M_PI / 180.0)) * 180. / M_PI; // how much to add per strip to the starting position
                double phi_n = ((-id + 0.5) * 30 + beta_n);
                phi_n += 45;

                // if(getenv("flip180")) {
                //	if(std::string(getenv("flip180"))=="1") {
                // if(dataset=="17F")
                //	phi_n+=180;//run 37 in 17F-->
                // }
                //}
                phi_n *= M_PI / 180.; // starting-position phi + strip contribution
                Event sx3ev(TVector3(88.0 * TMath::Cos(phi_n), 88.0 * TMath::Sin(phi_n), z), backE * 0.001, -1, det.ts, -1, det.stripB + 4 * id, det.stripF + 4 * id);
                SX3_Events.push_back(sx3ev);
                plotter->Fill2D("sx3backs_gm", 100, 0, 100, 800, 0, 8192, det.stripB + 4 * id, backE);

                // plotter->Fill2D("SX3CartesianPlot", 200, -100, 100, 200, -100, 100, 88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n), "hCalSX3");
                plotter->Fill2D("SX3CartesianPlot" + std::to_string(id), 200, -100, 100, 200, -100, 100, 88.0 * TMath::Cos(phi_n), 88.0 * TMath::Sin(phi_n), "hCalSX3");
            }
            if (det.valid && det.stripF != DEFAULT_NULL && det.stripB != DEFAULT_NULL)
            {
                plotter->Fill2D("sx3backs_raw", 100, 0, 100, 800, 0, 8192, det.stripB + 4 * id, det.backE);
            }
        }
    }
    // return kTRUE;
    //  QQQ Processing

    int qqqCount = 0;
    int qqqAdjCh = 0;
    // REMOVE WHEN RERUNNING USING THE NEW CALIBRATION FILE
    std::vector<Event> QQQ_Events, PC_Events;
    std::vector<Event> QQQ_Events_Raw, PC_Events_Raw;
    std::vector<Event> QQQ_Events2; // clustering done

    std::unordered_map<int, std::tuple<int, int, double, double>> qvecr[4], qvecw[4];
    if (qqq.multi > 1)
    {
        // if(qqq.multi>=3) std::cout << "-----" << std::endl;
        for (int i = 0; i < qqq.multi; i++)
        {
            if (qqq.ch[i] / 16)
            {
                if (qvecr[qqq.id[i]].find(qqq.ch[i]) != qvecr[qqq.id[i]].end())
                    std::cout << "mayday!" << std::endl;
                qvecr[qqq.id[i]][qqq.ch[i]] = std::tuple(qqq.id[i], qqq.ch[i], qqq.e[i], qqq.t[i]);
            }
            else
            {
                if (qvecw[qqq.id[i]].find(qqq.ch[i]) != qvecw[qqq.id[i]].end())
                    std::cout << "mayday!" << std::endl;
                qvecw[qqq.id[i]][qqq.ch[i]] = std::tuple(qqq.id[i], qqq.ch[i], qqq.e[i], qqq.t[i]);
            }
        }
    }

    bool PCSX3TimeCut = false;
    bool PCASX3TimeCut = false;
    bool PCCSX3TimeCut = false;

    bool PCQQQTimeCut = false;
    bool PCAQQQTimeCut = false;
    bool PCCQQQTimeCut = false;
    for (int i = 0; i < qqq.multi; i++)
    {
        plotter->Fill2D("QQQ_Index_Vs_Energy", 16 * 8, 0, 16 * 8, 2000, 0, 16000, qqq.index[i], qqq.e[i], "hRawQQQ");

        for (int j = 0; j < qqq.multi; j++)
        {
            if (j == i)
                continue;
            plotter->Fill2D("QQQ_Coincidence_Matrix", 16 * 8, 0, 16 * 8, 16 * 8, 0, 16 * 8, qqq.index[i], qqq.index[j], "hRawQQQ");
        }

        for (int k = 0; k < pc.multi; k++)
        {
            if (pc.index[k] < 24 && pc.e[k] > 10)
            {
                plotter->Fill2D("QQQ_Vs_Anode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
                plotter->Fill2D("QQQ_Vs_PC_Index", 16 * 8, 0, 16 * 8, 24, 0, 24, qqq.index[i], pc.index[k], "hRawQQQ");
            }
            else if (pc.index[k] >= 24 && pc.e[k] > 10)
            {
                plotter->Fill2D("QQQ_Vs_Cathode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
            }
        }

        for (int j = i + 1; j < qqq.multi; j++)
        {
            if (qqq.id[i] == qqq.id[j])
            {
                qqqCount++;

                int chWedge = -1;
                int chRing = -1;
                double eWedge = 0.0;
                double eWedgeMeV = 0.0;
                double eRing = 0.0;
                double eRingMeV = 0.0;
                double tRing = 0.0;
                double tWedge = 0.0;

                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
                {
                    chWedge = qqq.ch[i];
                    eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    chRing = qqq.ch[j] - 16;
                    eRing = qqq.e[j];
                    tRing = static_cast<double>(qqq.t[j]);
                    tWedge = static_cast<double>(qqq.t[i]);
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
                    chRing = qqq.ch[i] - 16;
                    eRing = qqq.e[i];
                    tRing = static_cast<double>(qqq.t[i]);
                    tWedge = static_cast<double>(qqq.t[j]);
                }
                else
                    continue;

                plotter->Fill1D("Wedgetime_Vs_Ringtime", 100, -1000, 1000, tWedge - tRing, "hTiming");
                plotter->Fill2D("RingE_vs_Index", 16 * 4, 0, 16 * 4, 1000, 0, 16000, chRing + qqq.id[i] * 16, eRing, "hRawQQQ");
                plotter->Fill2D("WedgeE_vs_Index", 16 * 4, 0, 16 * 4, 1000, 0, 16000, chWedge + qqq.id[i] * 16, eWedge, "hRawQQQ");
                plotter->Fill2D("WedgeE_Vs_RingECal", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");

                if (qqqCalibValid[qqq.id[i]][chWedge][chRing])
                {
                    eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;
                    eRingMeV = eRing * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;

                    if (eRingMeV / eWedgeMeV > 3.0 || eRingMeV / eWedgeMeV < 1.0 / 3.0)
                        continue;
                    // if(eRingMeV<1.2 || eWedgeMeV<1.2) continue;

                    double theta = 2 * TMath::Pi() * (-qqq.id[i] * 16 + (15 - chWedge) + 0.5) / (16 * 4);
                    double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
                    // z used to be 75+30+23=128
                    // we found a 12mm shift towards the vertex later --> 116
                    Event qqqevent(TVector3(rho * TMath::Cos(theta), rho * TMath::Sin(theta), qqq_z), eRingMeV, eWedgeMeV, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);
                    Event qqqeventr(TVector3(rho * TMath::Cos(theta), rho * TMath::Sin(theta), qqq_z), eRing, eWedge, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);
                    if (qqq.id[i] >= 1)
                    {
                        QQQ_Events.push_back(qqqevent);
                        QQQ_Events_Raw.push_back(qqqeventr);
                        plotter->Fill2D("WedgeE_Vs_RingECal_selected", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");
                    }
                    plotter->Fill2D("QQQCartesianPlot", 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
                    plotter->Fill2D("QQQCartesianPlot" + std::to_string(qqq.id[i]), 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
                    plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hPCQQQ");
                }
                else
                    continue;

                for (int k = 0; k < pc.multi; k++)
                {
                    plotter->Fill2D("RingCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
                    plotter->Fill2D("WedgeCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
                    plotter->Fill2D("WedgeCh_vs_Anode_Index" + std::to_string(qqq.id[i]), 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k]);
                    plotter->Fill2D("RingCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
                    plotter->Fill2D("WedgeCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");

                    if (pc.index[k] < 24 && pc.e[k] > 10)
                    {
                        plotter->Fill2D("Timing_Difference_QQQ_PC", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
                        plotter->Fill2D("DelT_Vs_QQQRingECal", 500, -2000, 2000, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hTiming");
                        // if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
                        if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
                        {
                            PCAQQQTimeCut = true;
                            plotter->Fill2D("CalibratedQQQEvsPCE_R", 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
                            plotter->Fill2D("CalibratedQQQEvsPCE_W", 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
                        }
                    }

                    if (pc.index[k] >= 24 && pc.e[k] > 10)
                    {
                        if (tRing - static_cast<double>(pc.t[k]) < -200)
                            PCCQQQTimeCut = true;
                        // if (tRing - static_cast<double>(pc.t[k]) > 200) PCCQQQTimeCut = true;
                        plotter->Fill2D("Timing_Difference_QQQ_PC_Cathode", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
                    }
                } // end of pc k loop

                if (!HitNonZero)
                {
                    // double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
                    // double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
                    double theta = 2 * TMath::Pi() * (-qqq.id[i] * 16 + (15 - chWedge) + 0.5) / (16 * 4);
                    double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
                    double x = rho * TMath::Cos(theta);
                    double y = rho * TMath::Sin(theta);
                    hitPos.SetXYZ(x, y, qqq_z);
                    // if(realtime) qqqg->SetPoint(0,hitPos.X(),hitPos.Y(),hitPos.Z());
                    if (realtime)
                        qqqg->AddPoint(hitPos.X(), hitPos.Y(), hitPos.Z());
                    qqqenergy = eRingMeV;
                    qqqtimestamp = tRing;
                    HitNonZero = true;
                }
            } // if j==i
        } // j loop end
    } // i loop end

    PCQQQTimeCut = PCAQQQTimeCut && PCCQQQTimeCut;
    plotter->Fill1D("QQQ_Multiplicity", 10, 0, 10, qqqCount, "hRawQQQ");

    typedef std::unordered_map<int, std::tuple<int, double, double>> WireEvent; // this stores nearest neighbour wire events, or a 'cluster'
    WireEvent aWireEvents, cWireEvents;                                         // naming for book keeping

    aWireEvents.clear();
    aWireEvents.reserve(24);
    if (realtime)
    {
        hha->Reset();
        hhc->Reset();
    }
    // PC Gain Matching and Filling
    double anodeT = -99999;
    double cathodeT = 99999;
    int anodeIndex = -1;
    int cathodeIndex = -1;
    for (int i = 0; i < pc.multi; i++)
    {
        // std::cout << pc.index[i] << " " << pc.e[i] << " " << std::endl;
        if (pc.e[i] > 10)
        {
            plotter->Fill2D("PC_Index_Vs_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], static_cast<double>(pc.e[i]), "hRawPC");
        }
        else
        {
            continue;
        }

        if (pc.index[i] < 48)
        {
            pc.e[i] = pcSlope[pc.index[i]] * pc.e[i] + pcIntercept[pc.index[i]];
            plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], pc.e[i], "hGMPC");
        }

        if (pc.index[i] < 24)
        {
            anodeT = static_cast<double>(pc.t[i]);
            anodeIndex = pc.index[i];

            if (getenv("flipa"))
            {
                int flip_offset = std::atoi(getenv("anode_offset"));
                int yes_to_flip = std::atoi(getenv("flipa"));
                if (yes_to_flip && flip_offset)
                {
                    int flipped_index = (23 - anodeIndex + flip_offset) % 24;
                    aWireEvents[flipped_index] = std::tuple(flipped_index, pc.e[i], static_cast<double>(pc.t[i]));
                    // std::cout << "Flipping anodes and offseting by " << flip_offset << " wires." << std::endl;
                }
                else if (flip_offset)
                {
                    int offset_index = (anodeIndex + flip_offset) % 24;
                    aWireEvents[pc.index[i]] = std::tuple(offset_index, pc.e[i], static_cast<double>(pc.t[i]));
                    // std::cout << "Offseting anodes without flip by " << offset_index << " wires." << std::endl;
                }
                else
                    aWireEvents[pc.index[i]] = std::tuple(pc.index[i], pc.e[i], static_cast<double>(pc.t[i]));
            }
            else
                aWireEvents[pc.index[i]] = std::tuple(pc.index[i], pc.e[i], static_cast<double>(pc.t[i]));
            if (realtime)
                hha->SetBinContent(hha->FindFixBin(anodeIndex), pc.e[i]);
        }
        else
        {
            cathodeT = static_cast<double>(pc.t[i]);
            cathodeIndex = pc.index[i] - 24;
            if (getenv("flipc"))
            {
                int flip_offset = std::atoi(getenv("flipc"));
                int flipped_index = (cathodeIndex + flip_offset) % 24;
                cWireEvents[flipped_index] = std::tuple(flipped_index, pc.e[i], static_cast<double>(pc.t[i]));
            }
            else
            {
                cWireEvents[pc.index[i] - 24] = std::tuple(pc.index[i] - 24, pc.e[i], static_cast<double>(pc.t[i]));
            }

            if (realtime)
                hhc->SetBinContent(hhc->FindFixBin(cathodeIndex), pc.e[i]);
        }

        if (anodeT != -99999 && cathodeT != 99999)
        {
            for (int j = 0; j < qqq.multi; j++)
            {
                plotter->Fill1D("PC_Time_qqq", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
                plotter->Fill2D("PC_Time_Vs_QQQ_ch", 200, -2000, 2000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, qqq.ch[j], "hTiming");
                plotter->Fill2D("PC_Time_vs_AIndex", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
                plotter->Fill2D("PC_Time_vs_CIndex", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
                // plotter->Fill1D("PC_Time_A" + std::to_string(anodeIndex) + "_C" + std::to_string(cathodeIndex), 200, -1000, 1000, anodeT - cathodeT, "TimingPC");
            }

            for (int j = 0; j < sx3.multi; j++)
            {
                plotter->Fill1D("PC_Time_sx3", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
            }
            for (auto sx3event : SX3_Events)
            {
                bool TCC = sx3event.Time1 - cathodeT < 0;
                bool TCA = sx3event.Time1 - anodeT < 0;
                // plotter->Fill2D("sx3_z_phi_awire"+std::to_string(anodeIndex)+"_TC"+std::to_string(TCA), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
                // plotter->Fill2D("sx3_z_phi_cwire"+std::to_string(cathodeIndex)+"_TC"+std::to_string(TCC), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
            }

            plotter->Fill1D("PC_Time", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
        }

        for (int j = i + 1; j < pc.multi; j++)
        {
            plotter->Fill2D("PC_Coincidence_Matrix", 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
            plotter->Fill2D("PC_Coincidence_Matrix_anodeMinusCathode_lt_-200_" + std::to_string(anodeT - cathodeT < -200), 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
            plotter->Fill2D("Anode_V_Anode", 24, 0, 24, 24, 0, 24, pc.index[i], pc.index[j], "hGMPC");
        }
    }
    anodeHits.clear();
    cathodeHits.clear();
    corrcatMax.clear();

    int aID = 0;
    int cID = 0;
    double aE = 0;
    double cE = 0;
    double aESum = 0;
    double cESum = 0;
    double aEMax = 0;
    double cEMax = 0;
    int aIDMax = 0;
    int cIDMax = 0;

    for (int i = 0; i < pc.multi; i++)
    {
        // if (pc.e[i] > 100)
        {
            if (pc.index[i] < 24)
            {
                anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
            }
            else if (pc.index[i] >= 24)
            {
                cathodeHits.push_back(std::pair<int, double>(pc.index[i] - 24, pc.e[i]));
            }
        }
    }

    std::sort(anodeHits.begin(), anodeHits.end(), [](std::pair<int, double> a, std::pair<int, double> b)
              { return a.first < b.first; });

    std::sort(cathodeHits.begin(), cathodeHits.end(), [](std::pair<int, double> a, std::pair<int, double> b)
              { return a.first < b.first; });

    // clusters = collection of (collection of wires) where each wire is (index, energy, timestamp)
    std::vector<std::vector<std::tuple<int, double, double>>> aClusters = pwinstance.Make_Clusters(aWireEvents);
    std::vector<std::vector<std::tuple<int, double, double>>> cClusters = pwinstance.Make_Clusters(cWireEvents);

    std::vector<std::pair<double, double>> sumE_AC;
    for (auto aCluster : aClusters)
    {
        for (auto cCluster : cClusters)
        {
            if (aCluster.size() == 0)
                continue;
            if (cCluster.size() == 0)
                continue;
            // both have at least 1, here. Keep the a1, c1 events
            auto [crossover, alpha, apSumE, cpSumE, apMaxE, cpMaxE, apTSMaxE, cpTSMaxE] = pwinstance.FindCrossoverProperties(aCluster, cCluster);
            if (alpha != 9999999 && apSumE != -1)
            {
                // Event PCEvent(crossover,apMaxE,cpMaxE,apTSMaxE,cpTSMaxE);
                // Event PCEvent(crossover,apSumE,cpSumE,apTSMaxE,cpTSMaxE);
                Event PCEvent(crossover, apSumE, cpMaxE, apTSMaxE, cpTSMaxE); // run12 shows cathode-max and anode-sum provide best dE signals.
                // std::cout << apSumE << " " << crossover.Perp() << " " << apMaxE << " " << apTSMaxE << std::endl;
                PCEvent.multi1 = aCluster.size();
                PCEvent.multi2 = cCluster.size();
                PCEvent.Anodech = std::get<0>(aCluster[0]);
                PCEvent.Cathodech = std::get<0>(cCluster[0]);
                PC_Events.push_back(PCEvent);
                sumE_AC.push_back(std::pair(apSumE, cpSumE));
            }
            else
            {
                ; // std::cout << "AAAA " << std::endl;
            }
        }
    }

    if (process_alpha_proton_scattering)
    {
        protonAlphaHistograms(plotter, QQQ_Events, SX3_Events, PC_Events);
        // return kTRUE;
    } // end if(process_alpha_proton_scattering)

    if (QQQ_Events.size() && PC_Events.size())
        plotter->Fill2D("PCEv_vs_QQQEv", 20, 0, 20, 20, 0, 20, QQQ_Events.size(), PC_Events.size());

    plotter->Fill2D("ac_vs_cc", 20, 0, 20, 20, 0, 20, aClusters.size(), cClusters.size(), "wiremult");
    for (auto cluster : aClusters)
    {
        plotter->Fill1D("aClusters" + std::to_string(aClusters.size()), 20, -5, 15, cluster.size(), "wiremult");
    }
    for (auto cluster : cClusters)
    {
        plotter->Fill1D("cClusters" + std::to_string(cClusters.size()), 20, -5, 15, cluster.size(), "wiremult");
    }

    if (cClusters.size() && aClusters.size())
    {
        plotter->Fill2D("ac_vs_cc_ign0", 20, 0, 20, 20, 0, 20, aClusters.size(), cClusters.size(), "wiremult");
    }

    for (auto sx3event : SX3_Events)
    {
        for (int i = 0; i < 24; i++)
        {
            if (aWireEvents.find(i) != aWireEvents.end())
            {
                auto awire = aWireEvents[i];
                if (sx3event.Time1 - (double)std::get<2>(awire) < -150)
                {
                    // plotter->Fill2D("sx3_z_phi2_awire"+std::to_string(std::get<0>(awire)), 400,-100,100, 100, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
                    // plotter->Fill2D("sx3_z_strip#_awire"+std::to_string(std::get<0>(awire)), 400,-100,100, 100, -50,50,sx3event.pos.Z(), sx3event.ch2);
                    plotter->Fill2D("onewire_dEa_Esx3_TC1_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, sx3event.Energy1, std::get<1>(awire), "1wire");
                    plotter->Fill2D("onewire_aNum_sx3Phi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -360, 360, i, sx3event.pos.Phi() * 180. / M_PI, "1wire");
                }
            }

            if (cWireEvents.find(i) != cWireEvents.end())
            {
                auto cwire = cWireEvents[i];
                if (sx3event.Time1 - (double)std::get<2>(cwire) < -150)
                {
                    // plotter->Fill2D("sx3_z_phi2_cwire"+std::to_string(std::get<0>(cwire)),400,-100,100, 100, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
                    // plotter->Fill2D("sx3_z_strip#_cwire"+std::to_string(std::get<0>(cwire)),400,-100,100, 100, -50,50,sx3event.pos.Z(), sx3event.ch2 );
                    plotter->Fill2D("onewire_dEc_Esx3_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, sx3event.Energy1, std::get<1>(cwire), "1wire");
                    plotter->Fill2D("onewire_cNum_sx3Phi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -360, 360, i, sx3event.pos.Phi() * 180. / M_PI, "1wire");
                }
            }
        } // for 'i' loop

        for (const auto acluster : aClusters)
        {
            auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster, "ANODE");
            int a_number = acluster.size();
            TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire, sx3event.pos.Phi());
            if (sx3event.Time1 - apTSMaxE < -150)
            {
                plotter->Fill2D("dEa_interp_Esx3_TC1_ignC" + std::to_string(acluster.size()), 400, 0, 10, 800, 0, 40000, sx3event.Energy1, apSumE, "ainterp_noc");
                plotter->Fill2D("aPhi_interp_sx3Phi_TC1_ignC" + std::to_string(acluster.size()), 120, -360, 360, 120, -360, 360, pc_closest.Phi() * 180. / M_PI, sx3event.pos.Phi() * 180. / M_PI, "ainterp_noc");
                plotter->Fill2D("aZ_interp_sx3Z_TC1_ignC" + std::to_string(acluster.size()), 400, -200, 200, 300, -100, 200, pc_closest.Z(), sx3event.pos.Z(), "ainterp_noc");
                TVector3 x2(pc_closest), x1(sx3event.pos);
                TVector3 v = x2 - x1;
                double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
                TVector3 vector_closest_to_axis = x1 + t_minimum * v;
                plotter->Fill2D("vertexZ_interp_sx3Z_TC1_ignC" + std::to_string(acluster.size()), 400, -200, 200, 300, -100, 200, vector_closest_to_axis.Z(), sx3event.pos.Z(), "ainterp_noc");
                plotter->Fill2D("vertexXY_interp_TC1_ignC" + std::to_string(acluster.size()), 200, -100, 100, 200, -100, 100, vector_closest_to_axis.X(), vector_closest_to_axis.Y(), "ainterp_noc");
            }
        }
    }

    for (auto qqqevent : QQQ_Events)
    {
        for (int i = 0; i < 24; i++)
        {
            if (aWireEvents.find(i) != aWireEvents.end())
            {
                auto awire = aWireEvents[i];
                if (qqqevent.Time1 - (double)std::get<2>(awire) < -150)
                {
                    plotter->Fill2D("onewire_dEa_Eqqq_TC1_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, std::get<1>(awire), "1wire");
                    plotter->Fill2D("onewire_aNum_QQQPhi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -360, 360, i, qqqevent.pos.Phi() * 180. / M_PI, "1wire");
                }
            }

            if (cWireEvents.find(i) != cWireEvents.end())
            {
                auto cwire = cWireEvents[i];
                if (qqqevent.Time1 - (double)std::get<2>(cwire) < -150)
                {
                    plotter->Fill2D("onewire_dEc_Eqqq_TC1_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, std::get<1>(cwire), "1wire");
                    plotter->Fill2D("onewire_cNum_QQQPhi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -360, 360, i, qqqevent.pos.Phi() * 180. / M_PI, "1wire");
                }
            }
        } // for 'i' loop
    }

    bool PCSX3PhiCut = false;
    for (auto pcevent : PC_Events)
    {
        if (aClusters.size() == 1 && cClusters.size() == 1)
        {
            // plotter->Fill1D("pcz_a"+std::to_string(aClusters.at(0).size())+"_c"+std::to_string(cClusters.at(0).size()),800,-200,200,pcevent.pos.Z(),"wiremult");
            std::string detid = "_+_";
            if (SX3_Events.size())
                detid = "+sx3";
            if (QQQ_Events.size())
                detid = "+qqq";
            // plotter->Fill1D("pcz_a"+std::to_string(aClusters.at(0).size())+"_c"+std::to_string(cClusters.at(0).size())+detid,800,-200,200,pcevent.pos.Z(),"wiremult");
        }

        PCSX3TimeCut = false;
        for (auto sx3event : SX3_Events)
        {
            plotter->Fill1D("dt_pcA_sx3B" + std::to_string(sx3event.ch2), 640, -2000, 2000, sx3event.Time1 - pcevent.Time1, "hTiming");
            plotter->Fill1D("dt_pcC_sx3B" + std::to_string(sx3event.ch2), 640, -2000, 2000, sx3event.Time1 - pcevent.Time2, "hTiming");
            if (sx3event.Time1 - pcevent.Time1 < 0) //-150 for alphas
                PCASX3TimeCut = 1;
            if (sx3event.Time1 - pcevent.Time2 < 0) //-200 for alphas
                PCCSX3TimeCut = 1;
            PCSX3TimeCut = PCASX3TimeCut && PCCSX3TimeCut;

            bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi() + TMath::Pi() / 4. && sx3event.pos.Phi() >= pcevent.pos.Phi() - TMath::Pi() / 4.;
            PCSX3PhiCut = phicut;

            plotter->Fill1D("dt_pcA_sx3B", 640, -2000, 2000, sx3event.Time1 - pcevent.Time1);
            plotter->Fill1D("dt_pcC_sx3B", 640, -2000, 2000, sx3event.Time1 - pcevent.Time2);
            plotter->Fill2D("dt_pcA_vs_sx3RE", 640, -2000, 2000, 400, 0, 30, sx3event.Time1 - pcevent.Time1, sx3event.Energy1);
            plotter->Fill2D("dE_E_Anodesx3B", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1);
            plotter->Fill2D("dE_E_Cathodesx3B", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2);
            if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
                plotter->Fill2D("dE_E_Anodesx3B_a1c2", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1);
            if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
                plotter->Fill2D("dE_E_Cathodesx3B_a1c2", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2);
            if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
                plotter->Fill2D("dE_E_Anodesx3B_a2c1", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1);
            if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
                plotter->Fill2D("dE_E_Cathodesx3B_a2c1", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2);
            plotter->Fill2D("sx3phi_vs_pcphi" + std::to_string(sx3event.Time1 - pcevent.Time1 < -150), 100, -360, 360, 100, -360, 360, sx3event.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI);
            if (PCSX3TimeCut)
            {
                plotter->Fill1D("dt_pcA_sx3B_timecut", 640, -2000, 2000, sx3event.Time1 - pcevent.Time1);
                plotter->Fill1D("dt_pcC_sx3B_timecut", 640, -2000, 2000, sx3event.Time1 - pcevent.Time2);
                plotter->Fill2D("xyplot_sx3" + std::to_string(sx3event.ch2 / 4), 100, -100, 100, 100, -100, 100, sx3event.pos.X(), sx3event.pos.Y());
                plotter->Fill2D("xyplot_sx3" + std::to_string(sx3event.ch2 / 4), 100, -100, 100, 100, -100, 100, pcevent.pos.X(), pcevent.pos.Y());
                plotter->Fill2D("pcz_vs_pcphi_TimeCut", 600, -200, 200, 120, -360, 360, pcevent.pos.Z(), pcevent.pos.Phi() * 180 / M_PI); // x-axis is all Si det, y-axis is PC anode+cathode only
            }

            double sx3rho = 88.0;           // approximate barrel radius
            double sx3z = sx3event.pos.Z(); // w.r.t target origin at 90 for run12
            double pcz = pcevent.pos.Z();
            double calcsx3theta = TMath::ATan2(sx3rho - z_to_crossover_rho(pcz), sx3z - pcz);
            plotter->Fill2D("dE2_E_Anodesx3B", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * TMath::Sin(calcsx3theta));
            plotter->Fill2D("dE2_E_Cathodesx3B", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * TMath::Sin(calcsx3theta));

            double sx3theta = TMath::ATan2(sx3rho, sx3z - source_vertex);
            double pczguess = 37.0 / TMath::Tan(sx3theta) + source_vertex;
            double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z()) / TMath::Tan(sx3theta) + source_vertex;
            double sinTheta = TMath::Sin(sx3theta);

            TVector3 x2(pcevent.pos), x1(sx3event.pos);
            TVector3 v = x2 - x1;
            double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
            TVector3 vector_closest_to_z_sx3 = x1 + t_minimum * v;
            plotter->Fill1D("VertexReconZ_SX3" + std::to_string(PCSX3TimeCut), 600, -1300, 1300, vector_closest_to_z_sx3.Z(), "hPCZSX3");
            plotter->Fill2D("VertexReconXY_SX3" + std::to_string(PCSX3TimeCut), 100, -100, 100, 100, -100, 100, vector_closest_to_z_sx3.X(), vector_closest_to_z_sx3.Y(), "hPCZSX3");

            plotter->Fill2D("pcz_vs_time", 2000, 0, 2000, 600, -200, 200, pcevent.Time1 * 1e-9, pcevent.pos.Z());                   // x-axis is all Si det, y-axis is PC anode+cathode only
            plotter->Fill2D("pcphi_vs_time", 2000, 0, 2000, 180, -360, 360, pcevent.Time1 * 1e-9, pcevent.pos.Phi() * 180. / M_PI); // x-axis is all Si det, y-axis is PC anode+cathode only
            // plotter->Fill2D("pcz_vs_time_strip"+std::to_string(sx3event.ch2),2000,0,2000,600,-200,200,pcevent.Time1*1e-9,pcevent.pos.Z()); //x-axis is all Si det, y-axis is PC anode+cathode only
            plotter->Fill2D("sx3phi_vs_time", 2000, 0, 2000, 180, -360, 360, pcevent.Time1 * 1e-9, sx3event.pos.Phi() * 180. / M_PI); // x-axis is all Si det, y-axis is PC anode+cathode only

            plotter->Fill2D("pcz_vs_sx3pczguess", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z()); // x-axis is all Si det, y-axis is PC anode+cathode only
            if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
            {
                // if(pcevent.multi1==1) {
                plotter->Fill2D("pcz_vs_sx3pczguess_A1C2", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z());
                double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());

                TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
                TVector3 v = x2f - x1;
                double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
                TVector3 r_rhoMin_fix = x1 + t_minimum * v;
                plotter->Fill1D("VertexRecon_pczfix_sx3", 800, -300, 300, r_rhoMin_fix.Z());
                plotter->Fill1D("pczfix_A1C2_1d_sx3", 600, -200, 200, pcz_fix);
                plotter->Fill2D("pczfix_vs_sx3pczguess_A1C2", 600, -200, 200, 600, -200, 200, pczguess, pcz_fix);
                plotter->Fill2D("pcz_vs_sx3pczguess_A1C2_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z());

                double sinTheta_customV = TMath::Sin((sx3event.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta());
                plotter->Fill2D("dE3_E_CathodeSX3_A1C2_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV);
                plotter->Fill2D("dE3_E_AnodeSX3_A1C2_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV);

                if (TMath::Abs(r_rhoMin_fix.Z()) < 200.0)
                {
                    plotter->Fill2D("dE3_E_AnodeSX3B_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV);
                    plotter->Fill2D("dE3_E_CathodeSX3B_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV);
                }

                // ==============================================================================
                // BENCHMARKING: Twisted Wire (1A0C logic) vs Cathode Charge Division (A1C2) for SX3
                // ==============================================================================
                if (aClusters.size() == 1)
                { // Ensure we unambiguously grab the correct anode wire
                    int aWireID = std::get<0>(aClusters.front().front());

                    // 1. Get wire geometry
                    TVector3 a1 = pwinstance.An[aWireID].first;
                    TVector3 wireVec = pwinstance.An[aWireID].first - pwinstance.An[aWireID].second;

                    // 2. Define track plane (Z-axis to SX3 hit)
                    TVector3 planeNormal(-TMath::Sin(sx3event.pos.Phi()), TMath::Cos(sx3event.pos.Phi()), 0.0);
                    double dot_wireVec = wireVec.Dot(planeNormal);

                    if (TMath::Abs(dot_wireVec) > 1e-6)
                    {

                        // 4. Reconstruct Vertex Z using ONLY the SX3 hit and the twisted wire

                        // 3. Find intersection of wire and track plane
                        double t_intersect = -(a1.Dot(planeNormal)) / dot_wireVec;
                        TVector3 pcz_intersect = a1 + t_intersect * wireVec;
                        TVector3 x2(pcz_intersect), x1(sx3event.pos);
                        TVector3 v = x2 - x1;
                        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
                        TVector3 vector_minimisedto_z = x1 + t_minimum * v;
                        double deltaRho = sx3event.pos.Perp() - pcz_intersect.Perp();
                        double deltaZ = sx3event.pos.Z() - pcz_intersect.Z();
                        double vertex_recon = sx3event.pos.Z() - sx3event.pos.Perp() * (deltaZ / deltaRho);

                        // ==============================================================================
                        // 5. FILL BENCHMARK PLOTS (Saved in the "1wire" folder)
                        // ==============================================================================
                        // A. Compare the PC Z-coordinate (Twisted Wire vs Cathodes)
                        plotter->Fill1D("Benchmark_SX3_PCZ_Difference", 400, -50, 50, pcz_intersect.Z() - pcevent.pos.Z(), "1wire");
                        plotter->Fill2D("Benchmark_SX3_PCZ_Twisted_vs_Cathode", 400, -200, 200, 400, -200, 200, pcevent.pos.Z(), pcz_intersect.Z(), "1wire");
                        plotter->Fill2D("Benchmark_SX3_PCZ_Twisted_vs_Cathode_sx3" + std::to_string(sx3event.ch2), 400, -200, 200, 400, -200, 200, pcevent.pos.Z(), pcz_intersect.Z(), "1wire");

                        // B. Compare the Vertex Z-coordinate
                        // r_rhoMin_fix.Z() is your cathode-based closest approach vertex already calculated just above this block
                        plotter->Fill1D("Benchmark_SX3_VertexZ_Difference", 400, -100, 100, vertex_recon - r_rhoMin_fix.Z(), "1wire");
                        plotter->Fill2D("Benchmark_SX3_VertexZ_Twisted_vs_0Cathode", 400, -200, 200, 400, -200, 200, r_rhoMin_fix.Z(), vertex_recon, "1wire");
                        plotter->Fill2D("Benchmark_SX3_VertexZ_Twisted_vs_0Cathodenew", 400, -200, 200, 400, -200, 200, vector_minimisedto_z.Z(), vertex_recon, "1wire");
                        plotter->Fill2D("Benchmark_SX3_VertexZ_Twisted_vs_0Cathode_sx3" + std::to_string(sx3event.ch2), 400, -200, 200, 400, -200, 200, r_rhoMin_fix.Z(), vertex_recon, "1wire");
                        plotter->Fill2D("Benchmark_SX3_VertexZ_Twisted_vs_0Cathode_anode" + std::to_string(aWireID), 400, -200, 200, 400, -200, 200, r_rhoMin_fix.Z(), vertex_recon, "1wire");
                        plotter->Fill2D("Benchmark_SX3XY" + std::to_string(sx3event.ch2), 400, -100, 100, 400, -100, 100, vector_minimisedto_z.X(), vector_minimisedto_z.Y(), "1wire");
                    }
                }
                // ==============================================================================

                // ==============================================================================
                // BENCHMARKING: Numerical Delta-Phi Minimization Scan
                // ==============================================================================
                if (aClusters.size() == 1)
                {
                    int aWireID = std::get<0>(aClusters.front().front());

                    // 1. Get wire geometry
                    TVector3 a1 = pwinstance.An[aWireID].first;  // Top of the wire
                    TVector3 a2 = pwinstance.An[aWireID].second; // Bottom of the wire
                    TVector3 wireVec = a2 - a1;                  // Vector pointing down the wire

                    // Variables to track our minimums during the scan
                    double min_delta_phi = 9999.0;
                    double best_t = -1.0;
                    TVector3 best_pcz_intersect;

                    // 2. THE SCAN: Walk down the wire in 1000 tiny steps
                    // (For a 380mm wire, this is checking every 0.38 mm)
                    int num_steps = 1000;
                    for (int i = 0; i <= num_steps; ++i)
                    {
                        double t_test = (double)i / num_steps;    // Ranges from 0.0 to 1.0
                        TVector3 test_pt = a1 + t_test * wireVec; // The 3D point at this step

                        // Calculate absolute Delta Phi between Si hit and this specific point on the wire
                        double dPhi = TMath::Abs(TVector2::Phi_mpi_pi(sx3event.pos.Phi() - test_pt.Phi()));

                        // If this is the smallest Delta Phi we've seen so far, save it!
                        if (dPhi < min_delta_phi)
                        {
                            min_delta_phi = dPhi;
                            best_t = t_test;
                            best_pcz_intersect = test_pt;
                        }
                    }

                    // 3. Extract the Z coordinate that yielded the minimum Delta Phi
                    double pcz_minimized = best_pcz_intersect.Z();

                    // 4. Reconstruct Vertex Z using our minimized PC Z
                    double anode_rho = best_pcz_intersect.Perp();
                    double deltaRho = sx3event.pos.Perp() - anode_rho;
                    double deltaZ = sx3event.pos.Z() - pcz_minimized;

                    double vertex_recon_minimized = sx3event.pos.Z() - sx3event.pos.Perp() * (deltaZ / deltaRho);

                    // ==============================================================================
                    // 5. FILL PLOTS
                    // ==============================================================================
                    // Look at how close we actually got to the Si Phi.
                    // If min_delta_phi > 0.1 radians, it means the track never truly matched the wire!
                    plotter->Fill1D("Benchmark_SX3_Min_DeltaPhi", 5000, -10, 10, min_delta_phi, "1wire");

                    // Standard benchmarking comparisons against the A1C2 Cathode baseline
                    plotter->Fill1D("Benchmark_SX3_PCZ_Diff_Scan", 800, -180, 180, pcz_minimized - pcevent.pos.Z(), "1wire");
                    plotter->Fill2D("Benchmark_SX3_PCZA12C_vs_minimized", 800, -200, 200, 800, -200, 200, pcevent.pos.Z(), pcz_minimized, "1wire");

                    plotter->Fill1D("Benchmark_SX3_VertexZ_Diff_minimized", 400, -100, 100, vertex_recon_minimized - r_rhoMin_fix.Z(), "1wire");
                    plotter->Fill2D("Benchmark_SX3_VertexZA12C_vs_minimized", 800, -200, 200, 800, -200, 200, r_rhoMin_fix.Z(), vertex_recon_minimized, "1wire");
                }
            }
            if (pcevent.multi1 == 1 && pcevent.multi2 == 3)
            {
                plotter->Fill2D("pcz_vs_sx3pczguess_A1C3", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z());
                plotter->Fill2D("pcz_vs_sx3pczguess_A1C3_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z());
            }

            plotter->Fill2D("pcz_vs_sx3pczguess_int", 600, -200, 200, 600, -200, 200, pcz_guess_int, pcevent.pos.Z()); // x-axis is all Si det, y-axis is PC anode+cathode only
            plotter->Fill2D("pcz_vs_sx3pczguess_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z());
            // plotter->Fill2D("pcz_vs_sx3pczguess_phi"+std::to_string(sx3event.pos.Phi()*180/M_PI),300,0,200,600,-200,200,pczguess,pcevent.pos.Z());

            /*plotter->Fill2D("pcz_vs_sx3z_strip="+std::to_string(sx3event.ch2),300,0,100,600,-200,200,sx3z,pcevent.pos.Z(),"sx3_vs_pc_zcorr");
            plotter->Fill2D("pcz_vs_sx3z_strip="+std::to_string(sx3event.ch2)+"_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),300,0,100,600,-200,200,sx3z,pcevent.pos.Z(),"sx3_vs_pc_zcorr");

            plotter->Fill2D("pcdEC_vs_sx3z_strip="+std::to_string(sx3event.ch2)+"_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),800,0,20000,600,-200,200,pcevent.Energy2,sx3z,"sx3_vs_pc_zcorr");
            plotter->Fill2D("pcdEA_vs_sx3z_strip="+std::to_string(sx3event.ch2)+"_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),800,0,20000,600,-200,200,pcevent.Energy1,sx3z,"sx3_vs_pc_zcorr");*/

            /*for(auto cc: cClusters)
             for(auto ac: aClusters) {
                plotter->Fill2D("pcz_sx3_phicut_a"+std::to_string(ac.size())+"_c"+std::to_string(cc.size())+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
                if(ac.size()==2 && cc.size()==1) {
                    plotter->Fill2D("pcz_sx3_phicut_a("+std::to_string(std::get<0>(ac.at(0)))+","+std::to_string(std::get<0>(ac.at(1)))+")_c"+std::to_string(std::get<0>(cc.at(0)))+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
                    plotter->Fill2D("a2c1_vs_sx3_strip",24,0,24,64,0,64,0.5*(std::get<0>(ac.at(0))+std::get<0>(ac.at(1))),sx3event.ch2,"hPCZSX3");
                    //plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<-150)+"_a("+std::to_string(std::get<0>(ac.at(0)))+","+std::to_string(std::get<0>(ac.at(1)))+")_c"+std::to_string(std::get<0>(cc.at(0))),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
                }
                if(cc.size()==2 && ac.size()==1) {
                    plotter->Fill2D("pcz_sx3_phicut_c("+std::to_string(std::get<0>(cc.at(0)))+","+std::to_string(std::get<0>(cc.at(1)))+")_a"+std::to_string(std::get<0>(ac.at(0)))+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
                    plotter->Fill2D("c2a1_vs_sx3_strip",24,0,24,64,0,64,0.5*(std::get<0>(cc.at(0))+std::get<0>(cc.at(1))),sx3event.ch2,"hPCZSX3");
                    plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<-150)+"_c("+std::to_string(std::get<0>(cc.at(0)))+","+std::to_string(std::get<0>(cc.at(1)))+")_a"+std::to_string(std::get<0>(ac.at(0))),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
                    //plotter->Fill2D("pcz_vs_sx3z_2C1A_phiCut_TC"+std::to_string(PCSX3TimeCut),300,0,200,600,-400,400,sx3z,pcevent.pos.Z());
                }

                if(ac.size()==1 && cc.size()==1) {
                    plotter->Fill2D("pcz_sx3_phicut_a("+std::to_string(std::get<0>(ac.at(0)))+")_c"+std::to_string(std::get<0>(cc.at(0)))+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
                    //plotter->Fill2D("a2c1_vs_sx3_strip",24,0,24,64,0,64,0.5*(std::get<0>(ac.at(0))+std::get<0>(ac.at(1))),sx3event.ch2,"hPCZSX3");
                    //plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<-150)+"_a("+std::to_string(std::get<0>(ac.at(0)))+")_c"+std::to_string(std::get<0>(cc.at(0))),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
                }
            }*/
            // end for

            bool sx3PhiCut = (TMath::Abs(sx3event.pos.Phi() - pcevent.pos.Phi()) < 45.0 * M_PI / 180.);
            PCSX3PhiCut = sx3PhiCut;
            plotter->Fill1D("pcz_sx3Coinc_phiCut" + std::to_string(sx3PhiCut) + "_TC" + std::to_string(PCSX3TimeCut), 300, 0, 200, sx3z);
            plotter->Fill2D("pcz_vs_sx3z_phiCut" + std::to_string(sx3PhiCut) + "_TC" + std::to_string(PCSX3TimeCut), 300, 0, 200, 600, -400, 400, sx3z, pcevent.pos.Z());

            // plotter->Fill2D("sx3E_vs_sx3z"+std::to_string(sx3event.ch2),400,0,30,300,0,200,sx3event.Energy1,sx3z);
            plotter->Fill2D("sx3E_vs_sx3z", 400, 0, 30, 300, 0, 200, sx3event.Energy1, sx3z);

            plotter->Fill2D("pcdEA_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy1, sx3z);
            plotter->Fill2D("pcdEC_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy2, sx3z);

            plotter->Fill2D("pcdEA_vs_sx3z" + std::to_string(sx3event.ch2), 800, 0, 20000, 300, 0, 200, pcevent.Energy1, sx3z, "pcE_vs_sx3pos");
            plotter->Fill2D("pcdEC_vs_sx3z" + std::to_string(sx3event.ch2), 800, 0, 20000, 300, 0, 200, pcevent.Energy2, sx3z, "pcE_vs_sx3pos");

            plotter->Fill2D("pcdE2A_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy1 * sinTheta, sx3z);
            plotter->Fill2D("pcdE2C_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy2 * sinTheta, sx3z);
            plotter->Fill2D("phi_vs_stripnum", 180, -180, 180, 48, 0, 48, pcevent.pos.Phi() * 180. / M_PI, sx3event.ch2);
            plotter->Fill2D("E_theta_AnodeSX3", 400, -20, 180, 300, 0, 15, sx3theta * 180 / M_PI, sx3event.Energy1);
        }
        if (PCSX3TimeCut)
        {
            plotter->Fill1D("PCZ_sx3", 800, -200, 200, pcevent.pos.Z(), "hPCZSX3");
            plotter->Fill1D("PCZ", 800, -200, 200, pcevent.pos.Z(), "phicut");
            /*for(auto cc: cClusters)
            for(auto ac: aClusters) {
                plotter->Fill1D("PCZsx3_phicut_a"+std::to_string(ac.size())+"_c"+std::to_string(cc.size()),800,-200,200,pcevent.pos.Z(),"hPCZSX3");
            }*/
        }
    } // end PC-SX3 coincidence

    /*for(size_t ii=0; ii<QQQ_Events.size(); ii++) {
        for(size_t jj=ii+1; jj<QQQ_Events.size(); jj++) {
            //if(TMath::Abs(QQQ_Events.at(ii).pos.Phi()*180/M_PI-QQQ_Events.at(jj).pos.Phi()*180/M_PI)>20) continue;
            if(QQQ_Events.at(ii).ch1 == QQQ_Events.at(jj).ch1) continue;
            if(QQQ_Events.at(ii).ch2 == QQQ_Events.at(jj).ch2) continue;
            if(QQQ_Events.at(ii).ch1 == QQQ_Events.at(jj).ch1-1) continue;
            if(QQQ_Events.at(ii).ch2 == QQQ_Events.at(jj).ch2-1) continue;
            if(QQQ_Events.at(ii).ch1 == QQQ_Events.at(jj).ch1+1) continue;
            if(QQQ_Events.at(ii).ch2 == QQQ_Events.at(jj).ch2+1) continue;

            double dt = QQQ_Events.at(ii).Time1-QQQ_Events.at(jj).Time1;
            plotter->Fill1D("dt_qqqi_qqqj",800,-2000,2000,dt);
            if(TMath::Abs(dt) > 150) continue;
            plotter->Fill1D("dt_qqqi_qqqj_coinc",800,-2000,2000,dt);
            double sum_e = QQQ_Events.at(ii).Energy1+QQQ_Events.at(jj).Energy1;
            plotter->Fill2D("sum_qqqE",400,0,30,400,0,30,QQQ_Events.at(ii).Energy1,sum_e);
            plotter->Fill2D("qqq_matrix",400,0,30,400,0,30,QQQ_Events.at(ii).Energy1,QQQ_Events.at(jj).Energy1);
            plotter->Fill2D("qqq_matrix",400,0,30,400,0,30,QQQ_Events.at(jj).Energy1,QQQ_Events.at(ii).Energy1);
            plotter->Fill2D("qqq_ch2_ch2",400,0,400,400,0,400,QQQ_Events.at(jj).ch2,QQQ_Events.at(ii).ch2);
            plotter->Fill2D("qqq_ch1_ch1",400,0,400,400,0,400,QQQ_Events.at(jj).ch1,QQQ_Events.at(ii).ch1);

            if(sum_e > 6.50 && sum_e < 7.50) {
                plotter->Fill2D("qqq_ang1_ang2",180,-360,360,180,-360,360,QQQ_Events.at(jj).pos.Phi()*180/M_PI,QQQ_Events.at(ii).pos.Phi()*180/M_PI);
                //if(PC_Events.size()<2) continue;
                for(auto pcevent: PC_Events) {
                    plotter->Fill2D("pcphi_vs_qqqphi_i_esumcut",180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,QQQ_Events.at(ii).pos.Phi()*180/M_PI);
                    plotter->Fill2D("pcphi_vs_qqqphi_j_esumcut",180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,QQQ_Events.at(jj).pos.Phi()*180/M_PI);
                }
            }
        }
    }*/

    ///////////////////Single wire analysis for the anodes///////////////////

    if (aClusters.size() == 1 && cClusters.size() == 0)
    {
        // Extract the primary anode hit properties
        auto anodeHit = aClusters.front().front();
        int aWireID = std::get<0>(anodeHit);
        double aEnergy = std::get<1>(anodeHit);
        double aTime = std::get<2>(anodeHit);

        // Get the 3D endpoints of the fired twisted anode wire from your geometry class
        TVector3 a1 = pwinstance.An[aWireID].first;
        TVector3 wireVec = pwinstance.An[aWireID].first - pwinstance.An[aWireID].second;

        // Loop over SX3_Events directly
        for (auto sx3event : SX3_Events)
        {
            if (sx3event.Time1 - aTime < -150) // Time cut for protons
            {
                // 1. Define the plane of the track (Z-axis to SX3 hit)
                TVector3 planeNormal(-TMath::Sin(sx3event.pos.Phi()), TMath::Cos(sx3event.pos.Phi()), 0.0);

                // 2. Find intersection of the twisted wire with this track plane
                double dot_wireVec = wireVec.Dot(planeNormal);

                // Prevent divide-by-zero if wire is perfectly parallel to the track plane
                if (TMath::Abs(dot_wireVec) < 1e-6)
                    continue;

                double t_intersect = -(a1.Dot(planeNormal)) / dot_wireVec;

                // Calculate the exact 3D coordinate on the wire that matches the SX3 phi
                TVector3 pcz_intersect = a1 + t_intersect * wireVec;

                // 3. Reconstruct Vertex Z
                double deltaRho = sx3event.pos.Perp() - pcz_intersect.Perp();
                double deltaZ = sx3event.pos.Z() - pcz_intersect.Z();

                double vertex_recon = sx3event.pos.Z() - sx3event.pos.Perp() * (deltaZ / deltaRho);

                // 4. Energy Loss Correction in Silicon
                double path_length = (sx3event.pos - TVector3(0, 0, vertex_recon)).Mag() * 0.1;
                double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1) - path_length);

                double theta_recon = (sx3event.pos - TVector3(0, 0, vertex_recon)).Theta();
                double sinTheta = TMath::Sin(theta_recon);

                // 5. Fill Diagnostics
                plotter->Fill1D("1A0C_twisted_pcz_recon_Phi_SX3" + std::to_string(PCSX3PhiCut), 600, -300, 300, pcz_intersect.Z(), "1A0C");
                plotter->Fill1D("1A0C_twisted_vertex_recon_Phi_SX3" + std::to_string(PCSX3PhiCut), 600, -300, 300, vertex_recon, "1A0C");

                plotter->Fill2D("1A0C_sx3_E_vs_theta_raw_Phi_SX3" + std::to_string(PCSX3PhiCut), 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, sx3event.Energy1, "1A0C");
                plotter->Fill2D("1A0C_sx3_E_vs_theta_corr_Phi_SX3" + std::to_string(PCSX3PhiCut), 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, sx3Efix, "1A0C");

                plotter->Fill2D("1A0C_dE_Ecorr_Anode_SX3_Phi" + std::to_string(PCSX3PhiCut), 400, 0, 30, 800, 0, 40000, sx3Efix, aEnergy * sinTheta, "1A0C");

                // Track where on the wire the hit occurred (0 to 1 is inside the physical PC)
                plotter->Fill1D("1A0C_wire_t_parameter_Phi" + std::to_string(PCSX3PhiCut), 200, -0.5, 1.5, t_intersect, "1A0C");
            }
        }

        // Loop over QQQ_Events directly

        for (auto qqqevent : QQQ_Events)
        {
            if (qqqevent.Time1 - aTime < -150) // Time cut for protons
            {
                // 1. Define the plane of the track (Z-axis to SX3 hit)
                TVector3 planeNormal(-TMath::Sin(qqqevent.pos.Phi()), TMath::Cos(qqqevent.pos.Phi()), 0.0);

                // 2. Find intersection of the twisted wire with this track plane
                double dot_wireVec = wireVec.Dot(planeNormal);

                // Prevent divide-by-zero if wire is perfectly parallel to the track plane
                if (TMath::Abs(dot_wireVec) < 1e-6)
                    continue;

                double t_intersect_QQQ = -(a1.Dot(planeNormal)) / dot_wireVec;

                // Calculate the exact 3D coordinate on the wire that matches the SX3 phi
                TVector3 pcz_intersect = a1 + t_intersect_QQQ * wireVec;

                // 3. Reconstruct Vertex Z
                double deltaRho = qqqevent.pos.Perp() - pcz_intersect.Perp();
                double deltaZ = qqqevent.pos.Z() - pcz_intersect.Z();

                double vertex_recon = qqqevent.pos.Z() - qqqevent.pos.Perp() * (deltaZ / deltaRho);

                // 4. Energy Loss Correction in Silicon
                double path_length = (qqqevent.pos - TVector3(0, 0, vertex_recon)).Mag() * 0.1;

                // FIXED: Changed sx3Efix to qqqEfix
                double qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1) - path_length);

                double theta_recon = (qqqevent.pos - TVector3(0, 0, vertex_recon)).Theta();
                double sinTheta = TMath::Sin(theta_recon);

                // 5. Fill Diagnostics
                plotter->Fill1D("1A0C_twisted_pcz_recon_Phi_QQQ" + std::to_string(PCSX3PhiCut), 600, -300, 300, pcz_intersect.Z(), "1A0C");
                plotter->Fill1D("1A0C_twisted_vertex_recon_Phi_QQQ" + std::to_string(PCSX3PhiCut), 600, -300, 300, vertex_recon, "1A0C");

                // FIXED: Changed "sx3" to "qqq" in the histogram names to avoid overwriting your barrel data
                plotter->Fill2D("1A0C_qqq_E_vs_theta_raw_Phi_" + std::to_string(PCSX3PhiCut), 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, qqqevent.Energy1, "1A0C");
                plotter->Fill2D("1A0C_qqq_E_vs_theta_corr_Phi_" + std::to_string(PCSX3PhiCut), 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, qqqEfix, "1A0C");

                plotter->Fill2D("1A0C_dE_Ecorr_Anode_QQQ_Phi" + std::to_string(PCSX3PhiCut), 400, 0, 30, 800, 0, 40000, qqqEfix, aEnergy * sinTheta, "1A0C");

                // Track where on the wire the hit occurred (0 to 1 is inside the physical PC)
                plotter->Fill1D("1A0C_wire_t_parameter_QQQ_Phi" + std::to_string(PCSX3PhiCut), 200, -0.5, 1.5, t_intersect_QQQ, "1A0C");
            }
        }
    }

    for (auto pcevent : PC_Events)
    {
        for (auto qqqevent : QQQ_Events)
        {
            plotter->Fill1D("dt_pcA_qqqR", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1);
            plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE", 640, -2000, 2000, 400, 0, 30, qqqevent.Time1 - pcevent.Time1, qqqevent.Energy1);
            plotter->Fill1D("dt_pcC_qqqW", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2);
            plotter->Fill2D("phiPC_vs_phiQQQ", 180, -360, 360, 180, -360, 360, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI);
            double sinTheta = TMath::Sin((qqqevent.pos - TVector3(0, 0, source_vertex)).Theta()); /// TMath::Sin((TVector3(51.5,0,128.) - TVector3(0,0,85)).Theta());

            TVector3 x2(pcevent.pos);
            TVector3 x1(qqqevent.pos);
            TVector3 v = x2 - x1;
            double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
            TVector3 r_rhoMin = x1 + t_minimum * v;

            // bool timecut = (qqqevent.Time1 - pcevent.Time1 < -150);
            bool timecut = (qqqevent.Time1 - pcevent.Time1 < -150);
            bool lowercut_cath = pcevent.Energy2 * sinTheta < 250 && (qqqevent.Energy2 < 5.0 || qqqevent.Energy1 < 5.0);
            bool phicut = qqqevent.pos.Phi() <= pcevent.pos.Phi() + TMath::Pi() / 4. && qqqevent.pos.Phi() >= pcevent.pos.Phi() - TMath::Pi() / 4.;

            if (lowercut_cath && phicut)
            {
                plotter->Fill1D("dt_pcA_qqqR_pidlow_PC1", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1);
                plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE_pidlow_PC1", 640, -2000, 2000, 400, 0, 30, qqqevent.Time1 - pcevent.Time1, qqqevent.Energy1);
                plotter->Fill1D("dt_pcC_qqqW_pidlow_PC1", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2);
            }
            if (timecut)
            { // && qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4. ) {

                plotter->Fill2D("dE_E_AnodeQQQR", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1);
                plotter->Fill2D("dE_E_CathodeQQQR", 400, 0, 30, 800, 0, 10000, qqqevent.Energy2, pcevent.Energy2);
                if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
                    plotter->Fill2D("dE_E_AnodeQQQR_a1c2", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1);
                if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
                    plotter->Fill2D("dE_E_CathodeQQQR_a1c2", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2);
                if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
                    plotter->Fill2D("dE_E_AnodeQQQR_a2c1", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1);
                if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
                    plotter->Fill2D("dE_E_CathodeQQQR_a2c1", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2);

                if (phicut)
                {
                    plotter->Fill2D("dE2_E_AnodeQQQR_TC1PC1_pidlow" + std::to_string(lowercut_cath), 400, 0, 30, 800, 0, 4000, qqqevent.Energy1, pcevent.Energy1 * sinTheta);
                    plotter->Fill2D("dE2_E_CathodeQQQW_TC1PC1_pidlow" + std::to_string(lowercut_cath), 400, 0, 30, 800, 0, 1000, qqqevent.Energy2, pcevent.Energy2 * sinTheta);
                    // plotter->Fill2D("E_theta_AnodeQQQR_TC1PC1_pidlow"+std::to_string(lowercut_cath),75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
                    plotter->Fill2D("E_theta_zoomin_AnodeQQQR_TC1PC1_pidlow" + std::to_string(lowercut_cath), 60, 0, 30, 300, 0, 15, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, qqqevent.Energy1);
                }

                plotter->Fill2D("dE2_E_AnodeQQQR_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 4000, qqqevent.Energy1, pcevent.Energy1 * sinTheta);
                plotter->Fill2D("dE2_E_CathodeQQQR_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 1000, qqqevent.Energy2, pcevent.Energy2 * sinTheta);
                plotter->Fill2D("dEC_vs_dEA_TC1_PC" + std::to_string(phicut), 800, 0, 40000, 800, 0, 10000, pcevent.Energy1, pcevent.Energy2);
                plotter->Fill2D("qqqphi_vs_time", 2000, 0, 2000, 180, -360, 360, pcevent.Time1 * 1e-9, qqqevent.pos.Phi() * 180. / M_PI); // x-axis is all Si det, y-axis is PC anode+cathode only

                plotter->Fill1D("dt_pcA_qqqR_timecut", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1);
                plotter->Fill1D("dt_pcC_qqqW_timecut", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2);
                plotter->Fill2D("dE_theta_AnodeQQQR", 90, 0, 90, 400, 0, 20000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, pcevent.Energy1);
                plotter->Fill2D("dE2_theta_AnodeQQQR_zoomin", 60, 0, 30, 400, 0, 5000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, pcevent.Energy1 * sinTheta);
                plotter->Fill2D("dE2_theta_AnodeQQQR", 90, 0, 90, 400, 0, 20000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, pcevent.Energy1 * sinTheta);
                plotter->Fill2D("phiPC_vs_phiQQQ_TimeCut", 180, -360, 360, 180, -360, 360, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI);

                // plotter->Fill2D("E_theta_AnodeQQQR_TC1_PC"+std::to_string(phicut),75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
                // plotter->Fill2D("E_theta_zoomin_AnodeQQQR_TC1_PC"+std::to_string(phicut),60,0,30,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
                // plotter->Fill2D("E2_theta_AnodeQQQR",75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
                plotter->Fill2D("Etot2_theta_AnodeQQQR", 75, 0, 90, 300, 0, 15, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, qqqevent.Energy1 + pcevent.Energy1 * anode_gain * sinTheta);

                plotter->Fill2D("dE_theta_CathodeQQQR", 75, 0, 90, 800, 0, 10000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, pcevent.Energy2);
                plotter->Fill2D("dE2_theta_CathodeQQQR", 75, 0, 90, 800, 0, 10000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, pcevent.Energy2 * sinTheta);
                plotter->Fill2D("dE2_theta_CathodeQQQR_zoomin", 60, 0, 30, 800, 0, 3000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, pcevent.Energy2 * sinTheta);

                plotter->Fill2D("dE_phi_AnodeQQQR", 100, -180, 180, 800, 0, 40000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Phi() * 180 / M_PI, pcevent.Energy1);
                plotter->Fill2D("dE_phi_CathodeQQQR", 100, -180, 180, 800, 0, 10000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Phi() * 180 / M_PI, pcevent.Energy2);
                plotter->Fill1D("PCZ", 800, -200, 200, pcevent.pos.Z(), "phicut");
                // plotter->Fill1D("PCZ_phicut_a"+std::to_string(aClusters.at(0).size())+"_c"+std::to_string(cClusters.at(0).size()),800,-200,200,pcevent.pos.Z(),"wiremult");

                double pcz_guess_37 = 37. / TMath::Tan((qqqevent.pos - TVector3(0, 0, source_vertex)).Theta()) + source_vertex;
                plotter->Fill2D("pczguess_vs_pc_37", 180, 0, 200, 150, 0, 200, pcz_guess_37, pcevent.pos.Z(), "phicut");

                double pcz_guess_42 = 42. / TMath::Tan((qqqevent.pos - TVector3(0, 0, source_vertex)).Theta()) + source_vertex;
                plotter->Fill2D("pczguess_vs_pc_42", 180, 0, 200, 150, 0, 200, pcz_guess_42, pcevent.pos.Z(), "phicut");

                double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z()) / TMath::Tan((qqqevent.pos - TVector3(0, 0, source_vertex)).Theta()) + source_vertex;
                // plotter->Fill2D("pczguess_vs_pc_int",180,0,200,150,0,200,pcz_guess_int,pcevent.pos.Z(),"phicut");
                plotter->Fill2D("pczguess_vs_pc_int", 400, -200, 200, 600, -400, 400, pcz_guess_int, pcevent.pos.Z(), "phicut");
                if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
                {
                    double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
                    TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
                    TVector3 v = x2f - x1;
                    double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
                    TVector3 r_rhoMin_fix = x1 + t_minimum * v;

                    double sinTheta_customV = TMath::Sin((qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta());
                    plotter->Fill2D("dE3_E_CathodeQQQW_A1C2_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, qqqevent.Energy2, pcevent.Energy2 * sinTheta_customV);
                    plotter->Fill2D("dE3_E_AnodeQQQR_A1C2_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV);

                    plotter->Fill1D("VertexRecon_pczfix_qqq", 800, -400, 400, r_rhoMin_fix.Z());
                    plotter->Fill1D("VertexRecon_pczfix_qqq_PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 800, -400, 400, r_rhoMin_fix.Z());

                    if (TMath::Abs(r_rhoMin_fix.Z()) < 200.0)
                    {
                        plotter->Fill2D("dE3_E_AnodeQQQR_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV);
                        plotter->Fill2D("dE3_E_CathodeQQQR_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV);
                    }

                    plotter->Fill1D("pczfix_A1C2_1d_qqq", 600, -200, 200, pcz_fix);
                    plotter->Fill2D("pczfix_vs_qqqpczguess_A1C2", 600, -200, 200, 600, -200, 200, pcz_guess_int, pcz_fix);
                    plotter->Fill2D("pczguess_vs_pc_int_A1C2", 400, -200, 200, 600, -400, 400, pcz_guess_int, pcevent.pos.Z(), "phicut");

                    double path_length = (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Mag() * 0.1;
                    // std::cout << path_length << std::endl;
                    double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length);
                    double qqqEfix_p = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1) - path_length);

                    plotter->Fill2D("E_thetaf_AnodeQQQR_TC1_PC" + std::to_string(phicut), 180, 0, 180, 600, 0, 15, (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta() * 180 / M_PI, qqqevent.Energy1);
                    if (lowercut_cath)
                        plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 180, 0, 180, 600, 0, 15, (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta() * 180 / M_PI, qqqEfix_p);
                    else
                    {
                        std::string zcut = "_" + std::to_string((TMath::Abs(r_rhoMin_fix.Z()) < 180));
                        plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath) + zcut, 180, 0, 180, 600, 0, 15, (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta() * 180 / M_PI, qqqEfix);
                    }

                    std::string morecuts = "_pidlow" + std::to_string(lowercut_cath) + "_vertexfix=" + std::to_string(floor(r_rhoMin_fix.Z() / 20) * 20 + 10);
                    // plotter->Fill2D("E_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut)+morecuts,180,0,180,800,0,8,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqevent.Energy1,"morecuts");

                    // plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut)+morecuts,180,0,180,800,0,8,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqEfix,"morecuts");

                    plotter->Fill2D("dE3_Ef_AnodeQQQR_TC1" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 600, 0, 15, 800, 0, 40000, qqqEfix, pcevent.Energy1 * sinTheta_customV);
                    plotter->Fill2D("dE3_Ef_CathodeQQQR_TC1PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 600, 0, 15, 800, 0, 10000, qqqEfix, pcevent.Energy2 * sinTheta_customV);

                    // ==============================================================================
                    // BENCHMARKING: Twisted Wire (1A0C logic) vs Cathode Charge Division (A1C2)
                    // ==============================================================================
                    if (aClusters.size() == 1)
                    { // Ensure we unambiguously grab the correct anode wire
                        int aWireID = std::get<0>(aClusters.front().front());

                        // 1. Get wire geometry
                        TVector3 a1 = pwinstance.An[aWireID].first;
                        TVector3 wireVec = pwinstance.An[aWireID].first - pwinstance.An[aWireID].second;

                        // 2. Define track plane (Z-axis to QQQ hit)
                        TVector3 planeNormal(-TMath::Sin(qqqevent.pos.Phi()), TMath::Cos(qqqevent.pos.Phi()), 0.0);
                        double dot_wireVec = wireVec.Dot(planeNormal);

                        if (TMath::Abs(dot_wireVec) > 1e-6)
                        {
                            // 3. Find intersection of wire and track plane
                            double t_intersect = -(a1.Dot(planeNormal)) / dot_wireVec;
                            TVector3 pcz_intersect = a1 + t_intersect * wireVec;

                            // 4. Reconstruct Vertex Z using ONLY the QQQ hit and the twisted wire (Ignoring Cathodes)
                            double deltaRho = qqqevent.pos.Perp() - pcz_intersect.Perp();
                            double deltaZ = qqqevent.pos.Z() - pcz_intersect.Z();
                            double vertex_recon_twisted = qqqevent.pos.Z() - qqqevent.pos.Perp() * (deltaZ / deltaRho);

                            // ==============================================================================
                            // 5. FILL BENCHMARK PLOTS (Saved in the "1wire" folder)
                            // ==============================================================================
                            // A. Compare the PC Z-coordinate (Twisted Wire vs Cathodes)
                            plotter->Fill1D("Benchmark_PCZ_Difference", 400, -50, 50, pcz_intersect.Z() - pcevent.pos.Z(), "1wire");
                            plotter->Fill2D("Benchmark_PCZ_Twisted_vs_Cathode", 400, -200, 200, 400, -200, 200, pcevent.pos.Z(), pcz_intersect.Z(), "1wire");

                            // B. Compare the Vertex Z-coordinate
                            plotter->Fill1D("Benchmark_VertexZ_Difference", 400, -100, 100, vertex_recon_twisted - r_rhoMin_fix.Z(), "1wire");
                            plotter->Fill1D("Benchmark_VertexZ_Difference", 400, -100, 100, vertex_recon_twisted - r_rhoMin_fix.Z(), "1wire");
                            plotter->Fill2D("Benchmark_VertexZ_Twisted_vs_Cathode", 400, -200, 200, 400, -200, 200, r_rhoMin_fix.Z(), vertex_recon_twisted, "1wire");

                            // C. Diagnostic: Where on the wire did it hit?
                            plotter->Fill1D("Benchmark_TwistedWire_t", 200, -0.5, 1.5, t_intersect, "1wire");
                        }
                    }
                    // ==============================================================================
                }
                double qqqrho = qqqevent.pos.Perp();
                double qqqz = (qqqevent.pos - TVector3(0, 0, source_vertex)).Z();
                double tan_theta = qqqrho / qqqz;
                double pcz_guess_int2 = z_to_crossover_rho(pcevent.pos.Z()) / tan_theta + source_vertex;
                plotter->Fill2D("pczguess_vs_pc_int2", 180, 0, 200, 150, 0, 200, pcz_guess_int2, pcevent.pos.Z(), "phicut");

                double qqqz2 = (qqqevent.pos - r_rhoMin).Z();
                double tan_theta2 = qqqrho / qqqz2;
                double pcz_guess_int3 = z_to_crossover_rho(pcevent.pos.Z()) / tan_theta2 + r_rhoMin.Z();
                plotter->Fill2D("pczguess_vs_pc_int3", 180, 0, 200, 150, 0, 200, pcz_guess_int3, pcevent.pos.Z(), "phicut");
                // plotter->Fill2D("pczguess_vs_pc_int2_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),180,0,200,150,0,200,pcz_guess_int2,pcevent.pos.Z(),"phicut");

                double pcz_guess = pcz_guess_int;
                plotter->Fill2D("pctheta_vs_qqqtheta_sv", 180, -360, 360, 180, -360, 360, (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, (pcevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, "phicut");
                plotter->Fill2D("pctheta_vs_qqqtheta_rmz", 180, -360, 360, 180, -360, 360, (qqqevent.pos - TVector3(0, 0, r_rhoMin.Z())).Theta() * 180 / M_PI, (pcevent.pos - TVector3(0, 0, r_rhoMin.Z())).Theta() * 180 / M_PI, "phicut");
                plotter->Fill2D("pctheta_vs_qqqtheta_rm", 180, -360, 360, 180, -360, 360, (qqqevent.pos - r_rhoMin).Theta() * 180 / M_PI, (pcevent.pos - r_rhoMin).Theta() * 180 / M_PI, "phicut");
                plotter->Fill2D("pczguess_vs_pc_phi=" + std::to_string(qqqevent.pos.Phi() * 180. / M_PI), 300, 0, 200, 150, 0, 200, pcz_guess, pcevent.pos.Z(), "phicut");
            }
        }
    } // end PC QQQ coincidence
    // HALFTIME! Can stop here in future versions
    // return kTRUE;
    if (anodeHits.size() >= 1 && cathodeHits.size() >= 1)
    {
        // 2. CRITICAL FIX: Define reference vector 'a'
        // In Analyzer.cxx, 'a' was left over from the loop. We use the first anode wire as reference here.
        // (Assuming pwinstance.An is populated and wires are generally parallel).
        TVector3 refAnode = pwinstance.An[0].first - pwinstance.An[0].second;

        {
            for (const auto &anode : anodeHits)
            {
                aID = anode.first;
                aE = anode.second;
                aESum += aE;
                if (aE > aEMax)
                {
                    aEMax = aE;
                    aIDMax = aID;
                }
            }

            for (const auto &cathode : cathodeHits)
            {
                cID = cathode.first;
                cE = cathode.second;
                plotter->Fill2D("AnodeMax_Vs_Cathode_Coincidence_Matrix", 24, 0, 24, 24, 0, 24, aIDMax, cID, "hRawPC");
                plotter->Fill2D("Anode_Vs_Cathode_Coincidence_Matrix", 24, 0, 24, 24, 0, 24, aID, cID, "hRawPC");
                plotter->Fill2D("Anode_Vs_Cathode_Coincidence_Matrix_qqq" + std::to_string(HitNonZero), 24, 0, 24, 24, 0, 24, aID, cID, "hRawPC");
                plotter->Fill2D("Anode_vs_CathodeE", 2000, 0, 30000, 2000, 0, 30000, aE, cE, "hGMPC");
                plotter->Fill2D("CathodeMult_V_CathodeE", 6, 0, 6, 2000, 0, 30000, cathodeHits.size(), cE, "hGMPC");
                /*for (int j = -4; j < 3; j++)
                {
                  if ((aIDMax + 24 + j) % 24 == 23 - cID)
                  {
                    corrcatMax.push_back(std::pair<int, double>(cID, cE));
                    cESum += cE;
                  }
                }*/
                if ((aIDMax + cID) % 24 == 22 || (aIDMax + cID) % 24 == 23 || (aIDMax + cID) % 24 >= 0 || (aIDMax + cID) % 24 <= 3)
                {
                    corrcatMax.push_back(std::pair<int, double>(cID, cE));
                    cESum += cE;
                    if (cE > cEMax)
                    {
                        cEMax = cE;
                        cIDMax = cID;
                    }
                }
            }
        }
    }

    TVector3 anodeIntersection, vector_closest_to_z;
    anodeIntersection.Clear();
    vector_closest_to_z.Clear();
    if (corrcatMax.size() > 0)
    {
        double x = 0, y = 0, z = 0;
        for (const auto &corr : corrcatMax)
        {
            if (pwinstance.Crossover[aIDMax][corr.first][0].z > 9000000)
                continue;
            if (cESum > 0)
            {
                x += (corr.second) / cESum * pwinstance.Crossover[aIDMax][corr.first][0].x;
                y += (corr.second) / cESum * pwinstance.Crossover[aIDMax][corr.first][0].y;
                z += (corr.second) / cESum * pwinstance.Crossover[aIDMax][corr.first][0].z;
            }
        }
        if (x == 0 && y == 0 && z == 0)
            ;
        // to ignore events with no valid crossover points
        else
        {
            anodeIntersection = TVector3(x, y, z);
            if (realtime)
            {
                // crossoverg->SetPoint(0,x,y,z);
                crossoverg->AddPoint(x, y, z);
            }
            // std::cout << "Anode Intersection: " << anodeIntersection.X() << ", " << anodeIntersection.Y() << ", " << anodeIntersection.Z() << " " << aIDMax <<  std::endl;
        }
    }
    bool PCQQQPhiCut = false;
    // flip the algorithm for cathode 1 multi anode events
    if ((hitPos.Phi() > (anodeIntersection.Phi() - TMath::PiOver4())) && (hitPos.Phi() < (anodeIntersection.Phi() + TMath::PiOver4())))
    {
        PCQQQPhiCut = true;
    }

    if (anodeIndex != -1 && cathodeIndex != -1 && hitPos.Perp() != 0 && anodeIntersection.Perp() != 0 && realtime && PCQQQPhiCut && PCQQQTimeCut)
    {
        // can1->Modified();
        // can1->Update();
        TVector3 x2(anodeIntersection);
        TVector3 x1(hitPos);
        TVector3 v = x2 - x1;
        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
        TVector3 r_rhoMin = x1 + t_minimum * v;

        trajectory->SetPoint(0, x1.X(), x1.Y(), x1.Z());
        trajectory->SetPoint(1, r_rhoMin.X(), r_rhoMin.Y(), r_rhoMin.Z());

        for (auto cath : corrcatMax)
        {
            plc[cath.first]->SetLineWidth(3);
            // plc[cath.first]->SetLineStyle(kLine);
        }
        for (auto anodeW : anodeHits)
        {
            pla[anodeW.first]->SetLineWidth(3);
            // pla[anodeW.first]->SetLineStyle(kLine);
        }
        // can2->Modified();
        // can2->Update();
        // while(can1->WaitPrimitive());

        // pla[anodeIndex]->SetLineWidth(1);
        // pla[anodeIndex]->SetLineStyle(kDotted);
        for (auto anodeW : anodeHits)
        {
            pla[anodeW.first]->SetLineWidth(1);
            pla[anodeW.first]->SetLineStyle(kDotted);
        }
        for (auto cath : corrcatMax)
        {
            plc[cathodeIndex]->SetLineStyle(kDotted);
            plc[cath.first]->SetLineWidth(1);
        }
    }

    if (anodeIntersection.Z() != 0 && anodeIntersection.Perp() > 0 && HitNonZero)
    {
        plotter->Fill1D("PC_Z_Projection", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
        plotter->Fill2D("Z_Proj_VsDelTime", 600, -300, 300, 200, -2000, 2000, anodeIntersection.Z(), anodeT - cathodeT, "hPCzQQQ");
        plotter->Fill2D("IntPhi_vs_QQQphi", 100, -200, 200, 80, -200, 200, anodeIntersection.Phi() * 180. / TMath::Pi(), hitPos.Phi() * 180. / TMath::Pi(), "hPCQQQ");
        // plotter->Fill2D("Inttheta_vs_QQQtheta", 90, 0, 180, 20, 0, 45, anodeIntersection.Theta() * 180. / TMath::Pi(), hitPos.Theta() * 180. / TMath::Pi(), "hPCQQQ");
        // plotter->Fill2D("Inttheta_vs_QQQtheta_TC" + std::to_string(PCQQQTimeCut)+ "_PC"+std::to_string(PCQQQPhiCut), 90, 0, 180, 20, 0, 45, anodeIntersection.Theta() * 180. / TMath::Pi(), hitPos.Theta() * 180. / TMath::Pi(), "hPCQQQ");
        plotter->Fill2D("IntPhi_vs_QQQphi_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 100, -200, 200, 80, -200, 200, anodeIntersection.Phi() * 180. / TMath::Pi(), hitPos.Phi() * 180. / TMath::Pi(), "hPCQQQ");
    }

    if (anodeIntersection.Z() != 0 && anodeIntersection.Perp() > 0 && PCSX3TimeCut)
    {
        plotter->Fill1D("PC_Z_Projection_sx3", 600, -200, 200, anodeIntersection.Z(), "hPCZSX3");
    }
    if (anodeIntersection.Z() != 0 && cathodeHits.size() >= 2)
        plotter->Fill1D("PC_Z_Projection_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");

    if (anodeIntersection.Z() != 0 && cathodeHits.size() == 1)
    {
        plotter->Fill1D("PC_Z_proj_1C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
        plotter->Fill2D("IntersectionPhi_vs_AnodeZ_1C", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hPCzQQQ");
    }

    if (anodeIntersection.Z() != 0 && cathodeHits.size() == 2)
    {
        plotter->Fill1D("PC_Z_proj_2C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
        plotter->Fill2D("IntersectionPhi_vs_AnodeZ_2C", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hGMPC");
    }
    if (anodeIntersection.Z() != 0 && cathodeHits.size() > 2)
    {
        plotter->Fill1D("PC_Z_proj_nC", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
        plotter->Fill2D("IntersectionPhi_vs_AnodeZ_nC", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hGMPC");
    }
    if (anodeHits.size() > 0 && cathodeHits.size() > 0)
        plotter->Fill2D("AHits_vs_CHits", 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");

    // make another plot with nearest neighbour constraint
    bool hasNeighbourAnodes = false;
    bool hasNeighbourCathodes = false;

    // 1. Check Anodes for neighbours (including wrap-around 0-23)
    for (size_t i = 0; i < anodeHits.size(); i++)
    {
        for (size_t j = i + 1; j < anodeHits.size(); j++)
        {
            int diff = std::abs(anodeHits[i].first - anodeHits[j].first);
            if (diff == 1 || diff == 23)
            { // 23 handles the cylindrical wrap
                hasNeighbourAnodes = true;
                break;
            }
        }
        if (hasNeighbourAnodes)
            break;
    }

    // 2. Check Cathodes for neighbours (including wrap-around 0-23)
    for (size_t i = 0; i < cathodeHits.size(); i++)
    {
        for (size_t j = i + 1; j < cathodeHits.size(); j++)
        {
            int diff = std::abs(cathodeHits[i].first - cathodeHits[j].first);
            if (diff == 1 || diff == 23)
            {
                hasNeighbourCathodes = true;
                break;
            }
        }
        if (hasNeighbourCathodes)
            break;
    }

    // ---------------------------------------------------------
    // FILL PLOTS
    // ---------------------------------------------------------
    if (anodeHits.size() > 0 && cathodeHits.size() > 0)
    {
        plotter->Fill2D("AHits_vs_CHits_NA" + std::to_string(hasNeighbourAnodes), 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");
        plotter->Fill2D("AHits_vs_CHits_NC" + std::to_string(hasNeighbourCathodes), 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");

        // Constraint Plot: Only fill if BOTH planes have adjacent hits
        // This effectively removes events with only isolated single-wire hits (noise)
        if (hasNeighbourAnodes && hasNeighbourCathodes)
        {
            plotter->Fill2D("AHits_vs_CHits_NN", 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");
        }
    }

    if (HitNonZero && anodeIntersection.Z() != 0)
    {
        pw_contr.CalTrack2(hitPos, anodeIntersection);
        plotter->Fill1D("VertexRecon", 600, -1300, 1300, pw_contr.GetZ0());
        plotter->Fill1D("VertexRecon_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, pw_contr.GetZ0());

        if (cathodeHits.size() == 2)
            plotter->Fill1D("VertexRecon_2c_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, pw_contr.GetZ0());

        TVector3 x2(anodeIntersection), x1(hitPos);

        TVector3 v = x2 - x1;
        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
        vector_closest_to_z = x1 + t_minimum * v;

        plotter->Fill1D("VertexRecon_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");

        if (qqqenergy < 4.0)
            plotter->Fill1D("VertexRecon_Z(qqqE<4.0MeV)_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");

        if (vector_closest_to_z.Perp() < 20)
        {
            plotter->Fill1D("VertexRecon_RadialCut_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");
        }

        plotter->Fill2D("VertexRecon_XY_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 100, -100, 100, 100, -100, 100, vector_closest_to_z.X(), vector_closest_to_z.Y(), "customVertex");
        if (cathodeHits.size() == 2)
        {
            plotter->Fill1D("VertexRecon2C_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");
            if (vector_closest_to_z.Perp() < 20)
            {
                plotter->Fill1D("VertexRecon2C_RadialCut_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");
            }
            plotter->Fill2D("VertexRecon2C_XY_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 100, -100, 100, 100, -100, 100, vector_closest_to_z.X(), vector_closest_to_z.Y(), "customVertex");
            plotter->Fill2D("VertexRecon2C_RhoZ_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 100, -100, 100, 600, -1300, 1300, vector_closest_to_z.Perp(), vector_closest_to_z.Z(), "customVertex");
            plotter->Fill2D("VertexRecon2C_Z_vs_QQQE_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, 800, 0, 20, vector_closest_to_z.Z(), qqqenergy, "customVertex");
        }
    }

    for (int i = 0; i < qqq.multi; i++)
    {
        if (anodeIntersection.Perp() > 0)
        { // suppress x,y=0,0 events
            if (PCQQQTimeCut)
            {
                plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");
                plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, hitPos.X(), hitPos.Y(), "hPCQQQ");
            }
            plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");
        }
        for (int j = i + 1; j < qqq.multi; j++)
        {
            if (qqq.id[i] == qqq.id[j])
            {
                int chWedge = -1;
                int chRing = -1;
                double eWedge = 0.0;
                double eWedgeMeV = 0.0;
                double eRing = 0.0;
                double eRingMeV = 0.0;
                double tRing = 0.0;
                int qqqID = -1;
                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
                {
                    chWedge = qqq.ch[i];
                    eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
                    chRing = qqq.ch[j] - 16;
                    eRing = qqq.e[j];
                    tRing = static_cast<double>(qqq.t[j]);
                    qqqID = qqq.id[i];
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
                    chRing = qqq.ch[i] - 16;
                    tRing = static_cast<double>(qqq.t[i]);
                    eRing = qqq.e[i];
                    qqqID = qqq.id[i];
                }
                else
                    continue;

                if (qqqCalibValid[qqq.id[i]][chRing][chWedge])
                {
                    eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chRing][chWedge] / 1000;
                    eRingMeV = eRing * qqqCalib[qqq.id[i]][chRing][chWedge] / 1000;
                }
                else
                    continue;

                // if (anodeIntersection.Z() != 0)
                {
                    plotter->Fill2D("PC_Z_vs_QQQRing", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCzQQQ");
                    plotter->Fill2D("PC_Z_vs_QQQRho", 600, -300, 300, 40, 40, 110, anodeIntersection.Z(), hitPos.Perp(), "hPCzQQQ");
                }

                if (anodeIntersection.Z() != 0 && cathodeHits.size() == 2)
                {
                    plotter->Fill2D("PC_Z_vs_QQQRing_2C", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCzQQQ");
                    plotter->Fill2D("PC_Z_vs_QQQRing_2C" + std::to_string(qqq.id[i]), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCzQQQ");
                    plotter->Fill2D("PC_Z_vs_QQQWedge_2C", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chWedge, "hPCzQQQ");
                }
                plotter->Fill2D("VertexRecon_QQQRingTC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, 16, 0, 16, vector_closest_to_z.Z(), chRing, "hPCQQQ");
                double phi = TMath::ATan2(anodeIntersection.Y(), anodeIntersection.X()) * 180. / TMath::Pi();
                plotter->Fill2D("PolarAngle_Vs_QQQWedge" + std::to_string(qqqID), 360, -360, 360, 16, 0, 16, phi, chWedge, "hPCQQQ");
                // plotter->Fill2D("EdE_PC_vs_QQQ_timegate_ls1000"+std::to_string())

                plotter->Fill2D("PC_Z_vs_QQQRing_Det" + std::to_string(qqqID), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCQQQ");
                // double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
                // double rho = 50. + 40. / 16. * (chRing + 0.5);

                for (int k = 0; k < pc.multi; k++)
                {
                    if (pc.index[k] >= 24)
                        continue;

                    //          double sinTheta = TMath::Sin((hitPos-vector_closest_to_z).Theta());
                    double sinTheta = TMath::Sin((anodeIntersection - TVector3(0, 0, 90.0)).Theta());
                    //           double sinTheta = TMath::Sin((anodeIntersection-vector_closest_to_z).Theta());
                    //          double sinTheta = TMath::Sin((hitPos-TVector3(0,0,30.0)).Theta());
                    //          double sinTheta = TMath::Sin(hitPos.Theta());

                    if (cathodeHits.size() == 2 && PCQQQPhiCut)
                    {
                        plotter->Fill2D("CalibratedQQQE_RvsCPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eRingMeV, pc.e[k] * sinTheta, "hPCQQQ");
                        plotter->Fill2D("CalibratedQQQE_WvsCPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eWedgeMeV, pc.e[k] * sinTheta, "hPCQQQ");
                        plotter->Fill2D("CalibratedQQQE_RvsPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
                        plotter->Fill2D("CalibratedQQQE_WvsPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
                        plotter->Fill2D("PCQQQ_dTimevsdPhi", 200, -2000, 2000, 80, -200, 200, tRing - static_cast<double>(pc.t[k]), (hitPos.Phi() - anodeIntersection.Phi()) * 180. / TMath::Pi(), "hTiming");
                    }
                }
            } /// qqq i==j case end
        } // j loop end
    } // qqq i loop end

    TVector3 guessVertex(0, 0, source_vertex); // for run12, subtract anodeIntersection.Z() by ~74.0 seems to work
    // rho=40.0 mm is halfway between the cathodes(rho=42) and anodes(rho=37)
    // double pcz_guess = 37.0/TMath::Tan((hitPos-guessVertex).Theta())  + guessVertex.Z(); //this is ideally kept to be all QQQ+userinput for calibration of pcz
    double pcz_guess = z_to_crossover_rho(anodeIntersection.Z()) / TMath::Tan((hitPos - guessVertex).Theta()) + guessVertex.Z(); // this is ideally kept to be all QQQ+userinput for calibration of pcz
    if (PCQQQTimeCut && PCQQQPhiCut && hitPos.Perp() > 0 && anodeIntersection.Perp() > 0 && cathodeHits.size() >= 2)
    {
        plotter->Fill2D("pczguess_vs_qqqE", 100, 0, 200, 800, 0, 20, pcz_guess, qqqenergy, "pczguess");
        double pczoffset = 0.0;
        // plotter->Fill2D("pczguess_vs_pcz_rad="+std::to_string(hitPos.Perp()),100,0,200,150,0,200,pcz_guess,anodeIntersection.Z(),"pczguess"); //entirely qqq-derived position vs entirely PC derived position
        plotter->Fill2D("pczguess_vs_pcz_phi=" + std::to_string(hitPos.Phi() * 180. / M_PI), 200, 0, 200, 200, 0, 200, pcz_guess, anodeIntersection.Z() + pczoffset, "pczguess"); // entirely qqq-derived position vs entirely PC derived position
        plotter->Fill2D("pczguess_vs_pcz", 200, 0, 200, 200, 0, 200, pcz_guess, anodeIntersection.Z() + pczoffset);
        plotter->Fill2D("pcz_vs_pcPhi_rad=" + std::to_string(hitPos.Perp()), 360, 0, 360, 150, 0, 200, anodeIntersection.Phi() * 180. / M_PI, anodeIntersection.Z() + pczoffset, "pczguess");
    }
    for (int i = 0; i < sx3.multi; i++)
    {
        // plotting sx3 strip hits vs anode phi
        if (sx3.ch[i] < 8 && anodeIntersection.Perp() > 0)
            plotter->Fill2D("PCPhi_vs_SX3Strip", 100, -200, 200, 8 * 24, 0, 8 * 24, anodeIntersection.Phi() * 180. / TMath::Pi(), sx3.id[i] * 8 + sx3.ch[i]);
    }

    if (anodeIntersection.Z() != 0 && cathodeHits.size() == 3)
    {
        plotter->Fill1D("PC_Z_proj_3C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    }

    if (anodeIntersection.Perp() != 0)
    {
        plotter->Fill2D("AnodeMaxE_Vs_Cathode_Sum_Energy", 2000, 0, 20000, 2000, 0, 10000, aEMax, cESum, "hGMPC");
        plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy", 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
        plotter->Fill2D("AnodeMaxE_Vs_Cathode_Max_Energy", 800, 0, 20000, 800, 0, 10000, aEMax, cEMax, "hGMPC");
        // double sinTheta = TMath::Sin((anodeIntersection - TVector3(0,0,source_vertex)).Theta());///TMath::Sin((TVector3(51.5,0,128.) - TVector3(0,0,85)).Theta());
        // plotter->Fill2D("AnodeMaxE_Vs_Cathode_Max_Energy_path_corrected", 800, 0, 20000, 800, 0, 10000, aEMax*sinTheta, cEMax*sinTheta, "hGMPC");
        plotter->Fill2D("AnodeSumE_Vs_Cathode_Sum_Energy", 800, 0, 20000, 800, 0, 10000, aESum, cESum, "hGMPC");
        plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_TC" + std::to_string(PCQQQTimeCut) + "_PC" + std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
        // plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_path_corrected"+std::to_string(PCQQQTimeCut)+"_PC"+std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cEMax*sinTheta, "hGMPC");
        // plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_path_corrected", 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cEMax*sinTheta, "hGMPC");

        if (PCQQQTimeCut && PCQQQPhiCut)
        {
            plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_TC" + std::to_string(PCQQQTimeCut) + "_PC" + std::to_string(PCQQQPhiCut) + "_cMax" + std::to_string(cIDMax), 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
        }
        // plotter->Fill2D("AnodeSumE_Vs_CathodeSum_Energy_path_corrected", 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cESum*sinTheta, "hGMPC");
        // plotter->Fill2D("AnodeSumE_Vs_CathodeSum_Energy_path_corrected_TC"+std::to_string(PCQQQTimeCut)+"_PC"+std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cESum*sinTheta, "hGMPC");        */
    }
    plotter->Fill1D("Correlated_Cathode_MaxAnode", 6, 0, 5, corrcatMax.size(), "hGMPC");
    plotter->Fill2D("Correlated_Cathode_VS_MaxAnodeEnergy", 6, 0, 5, 2000, 0, 30000, corrcatMax.size(), aEMax, "hGMPC");
    plotter->Fill1D("AnodeHits", 12, 0, 11, anodeHits.size(), "hGMPC");
    plotter->Fill2D("AnodeMaxE_vs_AnodeHits", 12, 0, 11, 2000, 0, 30000, anodeHits.size(), aEMax, "hGMPC");

    if (anodeHits.size() < 1)
    {
        plotter->Fill1D("NoAnodeHits_CathodeHits", 6, 0, 5, cathodeHits.size(), "hGMPC");
    }

    for (auto cwevent : cWireEvents)
    {
        // plotter->Fill1D("cwdtqqq_vs_cw"+std::to_string(PCQQQTimeCut),800,-2000,2000,24,0,24,std::get<2>(cwevent)-qqqtimestamp,std::get<0>(cwevent));
        for (auto awevent : aWireEvents)
        {
            plotter->Fill2D("aw_vs_cw", 24, 0, 24, 24, 0, 24, std::get<0>(awevent), std::get<0>(cwevent));
            plotter->Fill2D("aw_vs_cw_dtq" + std::to_string(PCQQQTimeCut), 24, 0, 24, 24, 0, 24, std::get<0>(awevent), std::get<0>(cwevent));
        }
    }
    for (auto awevent : aWireEvents)
    {
        // plotter->Fill1D("awdtqqq_vs_aw"+std::to_string(PCQQQTimeCut),800,-2000,2000,24,0,24,std::get<2>(awevent)-qqqtimestamp,std::get<0>(awevent));
    }

    return kTRUE;
}

void MakeVertex::Terminate()
{
    plotter->FlushToDisk(10);
    /*    can1->Modified();
        can1->Update();
        can2->Modified();
        can2->Update();
        while(can1->WaitPrimitive());
        while(can2->WaitPrimitive());*/
}

void protonAlphaHistograms(HistPlotter *plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events)
{

    // Sidetrack for a(p,p)
    std::string aplabel = "a(p,p)";
    Kinematics apkin_p(1.008664916, 4.002603254, 1.008664916, 4.002603254, 7.0); // m3 is proton
    Kinematics apkin_a(1.008664916, 4.002603254, 4.002603254, 1.008664916, 7.0); // m3 is alpha

    for (auto qqqevent : QQQ_Events)
    {
        for (auto sx3event : SX3_Events)
        {
            plotter->Fill1D("ap_qqq_sx3_dt", 800, -2000, 2000, qqqevent.Time1 - sx3event.Time1, aplabel);
            if (TMath::Abs(qqqevent.Time1 - sx3event.Time1) > 300)
                continue;
            // sx3event.pos.SetZ(sx3event.pos.Z()+5.0);
            plotter->Fill1D("ap_qqq_sx3_dt_timecut", 800, -2000, 2000, qqqevent.Time1 - sx3event.Time1, aplabel);
            plotter->Fill1D("ap_qqq_sx3_dphi", 180, -360, 360, qqqevent.pos.Phi() * 180 / M_PI - sx3event.pos.Phi() * 180 / M_PI, aplabel);
            plotter->Fill2D("ap_qqq_sx3_dphi_vs_qqqphi", 180, -360, 360, 180, -360, 360, qqqevent.pos.Phi() * 180 / M_PI - sx3event.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, aplabel);
            plotter->Fill2D("ap_qqq_sx3_matrix", 400, 0, 10, 400, 0, 10, qqqevent.Energy1, sx3event.Energy1, aplabel);

            for (auto pcevent : PC_Events)
            {

                double pcz_fix = pcfix_func.Eval(pcevent.pos.Z()) - 5.0;
                TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
                TVector3 x1(qqqevent.pos);
                TVector3 v = x2f - x1;
                double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
                TVector3 r_rhoMin_fix = x1 + t_minimum * v;
                double vertex_z = r_rhoMin_fix.Z();
                double theta_q = (qqqevent.pos - TVector3(0, 0, vertex_z)).Theta();
                // double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
                double sinTheta_customV = TMath::Sin(theta_q);
                double theta_s = (sx3event.pos - TVector3(0, 0, vertex_z)).Theta();
                // double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
                double sinTheta_s = TMath::Sin(theta_s);
                // if(vertex_z<0 || vertex_z>100) continue;

                // double sinTheta = TMath::Sin((qqqevent.pos - pcevent.pos).Theta());
                // plotter->Fill2D("sinTheta2_vs_sinTheta",80,-2,2,80,-2,2,sinTheta,sinTheta_customV,aplabel);

                plotter->Fill2D("ap_dE_E_Anodesx3B", 400, 0, 10, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, aplabel);
                plotter->Fill2D("ap_dE_E_Cathodesx3B", 400, 0, 10, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, aplabel);
                plotter->Fill2D("ap_dE_E_AnodeQQQ", 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, aplabel);
                plotter->Fill2D("ap_dE_E_CathodeQQQ", 400, 0, 10, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, aplabel);
                plotter->Fill2D("ap_dE3_E_AnodeQQQ", 400, 0, 10, 400, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV, aplabel);
                plotter->Fill2D("ap_dE3_E_CathodeQQQ", 400, 0, 10, 400, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV, aplabel);

                plotter->Fill2D("ap_dPhi_QQQ_PC", 180, -360, 360, 180, -360, 360, pcevent.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, aplabel);
                plotter->Fill2D("ap_dPhi_SX3_PC", 180, -360, 360, 180, -360, 360, pcevent.pos.Phi() * 180 / M_PI, sx3event.pos.Phi() * 180 / M_PI, aplabel);
                plotter->Fill1D("ap_dt_Anode_QQQ", 600, -2000, 2000, pcevent.Time1 - qqqevent.Time1, aplabel);
                plotter->Fill1D("ap_dt_Cathode_QQQ", 600, -2000, 2000, pcevent.Time2 - qqqevent.Time1, aplabel);
                plotter->Fill1D("ap_dt_Anode_SX3", 600, -2000, 2000, pcevent.Time1 - sx3event.Time1, aplabel);
                plotter->Fill1D("ap_dt_Cathode_SX3", 600, -2000, 2000, pcevent.Time2 - sx3event.Time1, aplabel);
                plotter->Fill1D("ap_pczfix", 600, -300, 300, pcz_fix, aplabel);
                plotter->Fill1D("ap_pcz", 600, -300, 300, pcevent.pos.Z(), aplabel);

                double path_length_q = (qqqevent.pos - TVector3(0, 0, vertex_z)).Mag() * 0.1;
                double path_length_s = (sx3event.pos - TVector3(0, 0, vertex_z)).Mag() * 0.1;
                // double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
                // double path_length_s = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;

                // We know that alphas predominantly are detected in QQQs, and protons in SX3s, and that protons don't leave much of a trace in dE layer.
                // Using the estimated path lengths, we correct alpha eloss in qqq, and protons in sx3. The result should (hopefully be) vertex independent.

                double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length_q);
                double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1) - path_length_s);
                // plotter->Fill2D("qqqEf_sx3E_matrix_all",400,0,10,400,0,10,qqqEfix,sx3event.Energy1,aplabel);
                plotter->Fill2D("ap_qqqEf_sx3Ef_matrix", 400, 0, 10, 400, 0, 10, qqqEfix, sx3Efix, aplabel);

                plotter->Fill2D("ap_Ef_vs_theta_qqq", 100, 0, 180, 400, 0, 10, theta_q * 180 / M_PI, qqqEfix, aplabel);
                plotter->Fill2D("ap_Ef_vs_theta_sx3", 100, 0, 180, 400, 0, 10, theta_s * 180 / M_PI, sx3Efix, aplabel);
                plotter->Fill2D("ap_theta_vs_theta_qqq_sx3", 100, 0, 180, 100, 0, 180, theta_q * 180 / M_PI, theta_s * 180 / M_PI, aplabel);
                plotter->Fill1D("ap_VertexReconZ", 400, -200, 200, vertex_z, aplabel);
                plotter->Fill2D("ap_VertexReconXY", 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), aplabel);
                plotter->Fill1D("ap_Ex_from_protons", 200, -10, 10, apkin_p.getExc(sx3Efix, theta_s * 180 / M_PI), aplabel);
                plotter->Fill1D("ap_Ex_from_alpha", 200, -10, 10, apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI), aplabel);

                if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
                { // one-anode, two-cathode events, as originally intended
                    // std::cout << "Test" << std::endl;
                    plotter->Fill1D("ap_VertexReconZ_a1c2", 400, -200, 200, vertex_z, aplabel);
                    plotter->Fill2D("ap_VertexReconXY_a1c2", 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), aplabel);
                    plotter->Fill2D("ap_theta_vs_theta_qqq_sx3_a1c2", 100, 0, 180, 100, 0, 180, theta_q * 180 / M_PI, theta_s * 180 / M_PI, aplabel);
                    plotter->Fill2D("ap_Ef_vs_theta_qqq_a1c2", 100, 0, 180, 400, 0, 10, theta_q * 180 / M_PI, qqqEfix, aplabel);
                    plotter->Fill1D("ap_Ex_from_protons_a1c2", 200, -10, 10, apkin_p.getExc(sx3Efix, theta_s * 180 / M_PI), aplabel);
                    plotter->Fill1D("ap_Ex_from_alpha_a1c2", 200, -10, 10, apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI), aplabel);

                    // std::cout << apkin_p.getExc(sx3Efix,theta_s*180/M_PI) << " " << apkin_a.getExc(qqqEfix,theta_q*180/M_PI)<<  std::endl;
                    plotter->Fill2D("ap_Ef_vs_theta_sx3_a1c2", 100, 0, 180, 400, 0, 10, theta_s * 180 / M_PI, sx3Efix, aplabel);

                    // plotter->Fill2D("qqqEf_sx3E_matrix",400,0,10,400,0,10,qqqEfix,sx3event.Energy1,aplabel);
                    plotter->Fill2D("ap_qqq_sx3_matrix_a1c2", 400, 0, 10, 400, 0, 10, qqqevent.Energy1, sx3event.Energy1, aplabel);
                    plotter->Fill2D("ap_qqqEf_sx3Ef_matrix_a1c2", 400, 0, 10, 400, 0, 10, qqqEfix, sx3Efix, aplabel);
                    // std::cout << sx3event.Energy1 	<< " " <<  path_length_s << " " << sx3Efix << std::endl;

                    // plotter->Fill2D("dE3_Ef_AnodeQQQ_a1c2",400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV,aplabel);
                    // plotter->Fill2D("dE3_Ef_CathodeQQQ_a1c2",400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,aplabel);

                } // end if(a1c2) loop
            } // end PC_Events for loop

        } // end SX3_Events for loop
    } // end QQQ_Events for loop, end sidetrack a(p,p)

    return;
}
