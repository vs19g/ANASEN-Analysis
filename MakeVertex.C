#define MakeVertex_cxx

#include "MakeVertex.h"
#include "Armory/ClassPW.h"
#include "Armory/HistPlotter.h"
#include "Armory/SX3Geom.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TBranch.h>
#include "TVector3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <algorithm>

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

  TVector3 pos;
  int ch1 = -1;        // int(ch1/16) gives qqq id, ch1%16 gives ring#
  int ch2 = -1;        // int(ch2/16) gives qqq id, ch2%16 gives wedge#
  double Energy1 = -1; // Front for QQQ, Anode for PC
  double Energy2 = -1; // Back for QQQ, Cathode for PC
  double Time1 = -1;
  double Time2 = -1;
};

// Calibration globals
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
const double qqqpos = 100.0;
const double vertexpos = 14.2;
const double pcrad = 37.0;
double qqqGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
double qqqCalib[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqCalibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
// TCutg *cutQQQ;

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

void MakeVertex::Begin(TTree * /*tree*/)
{
  TString option = GetOption();
  plotter = new HistPlotter("Analyzer_SX3.root", "TFILE");
  pw_contr.ConstructGeo();
  pwinstance.ConstructGeo();

  // ---------------------------------------------------------
  // 1. CRITICAL FIX: Initialize PC Arrays to Default (Raw)
  // ---------------------------------------------------------
  for (int i = 0; i < 48; i++)
  {
    pcSlope[i] = 1.0;     // Default slope = 1 (preserves Raw energy)
    pcIntercept[i] = 0.0; // Default intercept = 0
  }

  // Calculate Crossover Geometry ONCE
  TVector3 a, c, diff;
  double a2, ac, c2, adiff, cdiff, denom, alpha;

  for (size_t i = 0; i < pwinstance.An.size(); i++)
  {
    a = pwinstance.An[i].first - pwinstance.An[i].second;

    for (size_t j = 0; j < pwinstance.Ca.size(); j++)
    {
      c = pwinstance.Ca[j].first - pwinstance.Ca[j].second;
      diff = pwinstance.An[i].first - pwinstance.Ca[j].first;
      a2 = a.Dot(a);
      c2 = c.Dot(c);
      ac = a.Dot(c);
      adiff = a.Dot(diff);
      cdiff = c.Dot(diff);
      denom = a2 * c2 - ac * ac;
      alpha = (ac * cdiff - c2 * adiff) / denom;

      Crossover[i][j][0].x = pwinstance.An[i].first.X() + alpha * a.X();
      Crossover[i][j][0].y = pwinstance.An[i].first.Y() + alpha * a.Y();
      Crossover[i][j][0].z = pwinstance.An[i].first.Z() + alpha * a.Z();

      if (Crossover[i][j][0].z < -190 || Crossover[i][j][0].z > 190 || (i + j) % 24 == 12)
      {
        Crossover[i][j][0].z = 9999999;
      }

      Crossover[i][j][1].x = alpha;
      Crossover[i][j][1].y = 0;
    }
  }

  // Load PC Calibrations
  std::ifstream inputFile("slope_intercept_results.dat");
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

  // Load QQQ Cuts from file
  //   {
  //     std::string filename = "QQQ_PCCut.root";
  //     TFile *cutFile = TFile::Open(filename.c_str(), "READ");
  //     if (cutFile && !cutFile->IsZombie())
  //     {
  //       cutQQQ = (TCutg *)cutFile->Get("cutQQQPC");
  //       if (cutQQQ)
  //       {
  //         std::cout << "Loaded QQQ PC cut from " << filename << std::endl;
  //       }
  //       else
  //       {
  //         std::cerr << "Error: cutQQQPC not found in " << filename << std::endl;
  //       }
  //       cutFile->Close();
  //   }
  // }

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
    std::ifstream infile("sx3cal/backgains.dat");
    std::string temp;
    int backpos, frontpos, clkpos;
    std::cout << "foo" << std::endl;
    if (infile.is_open())
      while (infile >> clkpos >> temp >> frontpos >> temp >> backpos >> sx3BackGain[clkpos][frontpos][backpos])
        std::cout << sx3BackGain[clkpos][frontpos][backpos] << std::endl;
    infile.close();

    infile.open("sx3cal/frontgains.dat");
    if (infile.is_open())
      while (infile >> clkpos >> temp >> temp >> frontpos >> sx3FrontOffset[clkpos][frontpos] >> sx3FrontGain[clkpos][frontpos])
        std::cout << sx3FrontOffset[clkpos][frontpos] << " " << sx3FrontGain[clkpos][frontpos] << std::endl;
    infile.close();

    infile.open("sx3cal/rightgains.dat");
    if (infile.is_open())
      while (infile >> clkpos >> frontpos >> temp >> sx3RightGain[clkpos][frontpos])
      {
        sx3RightGain[clkpos][frontpos] = TMath::Abs(sx3RightGain[clkpos][frontpos]);
      }
    infile.close();
  }
  std::cout << "aaa" << std::endl;
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

  sx3.CalIndex();
  qqq.CalIndex();
  pc.CalIndex();

  std::vector<Event> sx3Events;
  // if(sx3.multi>1) {
  // 	std::array<sx3det,12> Fsx3;
  // 	//std::cout << "-----" << std::endl;
  // 	for(int i=0; i<sx3.multi; i++) {
  // 		if(sx3.id[i]>=12) continue;
  // 		int id = sx3.id[i];
  // 		if(sx3.ch[i]>=8) {
  // 			int sx3ch=sx3.ch[i]-8;
  // 			sx3ch=(sx3ch+3)%4;
  // 			if(sx3ch==0 || sx3ch==3) continue;
  // 			float value=sx3.e[i];
  // 			int gch = sx3.id[i]*4+(sx3.ch[i]-8);
  // 			Fsx3.at(id).fillevent("BACK",sx3ch,value);
  // 			Fsx3.at(id).ts = static_cast<double>(sx3.t[i]);
  // 			plotter->Fill2D("sx3backs_raw",100,0,100,800,0,4096,gch,sx3.e[i]);
  // 		} else {
  // 			int sx3ch=sx3.ch[i]/2;
  // 			double value=sx3.e[i];
  // 			if(sx3.ch[i]%2==0) {
  // 	  		    Fsx3.at(id).fillevent("FRONT_L",sx3ch,value*sx3RightGain[id][sx3ch]);
  // 			} else {
  // 				Fsx3.at(id).fillevent("FRONT_R",sx3ch,value);
  // 			}
  // 		}
  // 	}
  // 	for(int id=0; id<12; id++) {
  // 		Fsx3.at(id).validate();
  // 		auto det = Fsx3.at(id);
  // 		bool no_charge_sharing_strict = det.valid_front_chans.size()==1 && det.valid_back_chans.size()==1;
  // 		if(det.valid) {
  // 		    //std::cout << det.frontEL << " " << det.frontEL*sx3RightGain[id][det.stripF] << std::endl;
  // 			plotter->Fill2D("be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_b"+std::to_string(det.stripB),200,-1,1,800,0,8192,
  // 							det.frontX,det.backE,"evsx");
  // 			//std::cout << sx3BackGain[id][det.stripF][det.stripB] << " " << sx3FrontGain[id][det.stripF] << std::endl;
  // 			plotter->Fill2D("matched_be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_"+std::to_string(id*4+det.stripF),200,-30,30,800,0,8192,
  // 							det.frontX*sx3FrontGain[id][det.stripF]+sx3FrontOffset[id][det.stripF],det.backE*sx3BackGain[id][det.stripF][det.stripB],"evsx_matched");
  // 			//plotter->Fill2D("fe_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_"+std::to_string(det.stripB),200,-1,1,800,0,4096,det.frontX,det.backE,"evsx");
  // 			plotter->Fill2D("l_vs_r_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF),800,0,4096,800,0,4096,det.frontEL,det.frontER,"l_vs_r");
  // 		}
  // 		if(det.valid && (id ==9 || id==7 || id == 1 || id==3) && det.stripF!=DEFAULT_NULL && det.stripB!=DEFAULT_NULL) {
  // 			double z = det.frontX*sx3FrontGain[id][det.stripF]+sx3FrontOffset[id][det.stripF];
  // 			double backE = det.backE*sx3BackGain[id][det.stripF][det.stripB];
  // 			Event sx3ev(TVector3(0,0,z),backE,-1,det.ts,-1,det.stripB+4*id,det.stripF+4*id);
  // 			sx3Events.push_back(sx3ev);
  // 		}
  // 	}
  // }
  // return kTRUE;
  // QQQ Processing

  int qqqCount = 0;
  int qqqAdjCh = 0;
  // REMOVE WHEN RERUNNING USING THE NEW CALIBRATION FILE
  // for (int i = 0; i < qqq.multi; i++)
  // {
  //   //if ((qqq.id[i] == 3 || qqq.id[i] == 1) && qqq.ch[i] < 16)
  //   if (qqq.id[i] == 1 && qqq.ch[i] < 16) //for run 12, 26Al
  //   {
  //     qqq.ch[i] = 16 - qqq.ch[i];
  //   }
  // }
  // for (int i = 0; i < qqq.multi; i++)
  // {
  //   if (qqq.id[i] == 0 && qqq.ch[i] >= 16)
  //   {
  //     qqq.ch[i] = 31 - qqq.ch[i] + 16;
  //   }
  // }

  std::vector<Event> QQQ_Events, PC_Events;
  std::vector<Event> QQQ_Events_Raw, PC_Events_Raw;
  std::vector<Event> QQQ_Events2; // clustering done

  std::unordered_map<int, std::tuple<int, int, double, double>> qvecr[4], qvecw[4];
  if (qqq.multi > 1)
  {
    // if(qqq.multi>=3) std::cout << "-----" << std::endl;
    for (int i = 0; i < qqq.multi; i++)
    {
      // if(qqq.multi>=3) std::cout << std::setprecision(16) << "qqq"<< qqq.id[i] << " " << std::string(qqq.ch[i]/16?"ring":"wedge") << qqq.ch[i]%16 << " " << qqq.e[i] << " " << qqq.t[i] - qqq.t[0] << std::endl;
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

  bool PCQQQTimeCut = false;
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
      if (pc.index[k] < 24 && pc.e[k] > 50)
      {
        plotter->Fill2D("QQQ_Vs_Anode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
        plotter->Fill2D("QQQ_Vs_PC_Index", 16 * 8, 0, 16 * 8, 24, 0, 24, qqq.index[i], pc.index[k], "hRawQQQ");
      }
      else if (pc.index[k] >= 24 && pc.e[k] > 50)
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

        if (qqqCalibValid[qqq.id[i]][chWedge][chRing])
        {
          eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;
          eRingMeV = eRing * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;

          if (eRingMeV / eWedgeMeV > 3.0 || eRingMeV / eWedgeMeV < 1.0 / 3.0)
            continue;
          // if(eRingMeV<4.0 || eWedgeMeV<4.0) continue;

          double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
          double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
                                                           // z used to be 75+30+23=128
          // we found a 12mm shift towards the vertex later --> 116
          Event qqqevent(TVector3(rho * TMath::Cos(theta), rho * TMath::Sin(theta), qqqpos), eRingMeV, eWedgeMeV, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);
          Event qqqeventr(TVector3(rho * TMath::Cos(theta), rho * TMath::Sin(theta), qqqpos), eRing, eWedge, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);
          QQQ_Events.push_back(qqqevent);
          QQQ_Events_Raw.push_back(qqqeventr);
          plotter->Fill2D("QQQCartesianPlot", 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
          plotter->Fill2D("QQQCartesianPlot" + std::to_string(qqq.id[i]), 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
          if (PCQQQTimeCut)
          {
            plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hPCQQQ");
          }
          plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hPCQQQ");
        }
        else
          continue;

        plotter->Fill2D("WedgeE_Vs_RingECal", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");
        plotter->Fill2D("WedgeE_Vs_RingECal_selected", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");

        for (int k = 0; k < pc.multi; k++)
        {
          plotter->Fill2D("RingCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Anode_Index" + std::to_string(qqq.id[i]), 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k]);
          plotter->Fill2D("RingCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");

          if (pc.index[k] < 24 && pc.e[k] > 50)
          {
            plotter->Fill2D("Timing_Difference_QQQ_PC", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
            plotter->Fill2D("DelT_Vs_QQQRingECal", 500, -2000, 2000, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hTiming");
            plotter->Fill2D("CalibratedQQQEvsPCE_R", 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
            plotter->Fill2D("CalibratedQQQEvsPCE_W", 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
            if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
            // if (tRing - static_cast<double>(pc.t[k]) < -150 && tRing - static_cast<double>(pc.t[k]) > -450) // 27Al
            // if (tRing - static_cast<double>(pc.t[k]) < -70 && tRing - static_cast<double>(pc.t[k]) > -150) // 17F
            {
              PCQQQTimeCut = true;
            }
          }

          if (pc.index[k] >= 24 && pc.e[k] > 50)
          {
            plotter->Fill2D("Timing_Difference_QQQ_PC_Cathode", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
          }
        } // end of pc k loop

        if (!HitNonZero)
        {
          double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
          double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
          double x = rho * TMath::Cos(theta);
          double y = rho * TMath::Sin(theta);
          hitPos.SetXYZ(x, y, (qqqpos));
          qqqenergy = eRingMeV;
          qqqtimestamp = tRing;
          HitNonZero = true;
        }
      } // if j==i
    } // j loop end
  } // i loop end

  plotter->Fill1D("QQQ_Multiplicity", 10, 0, 10, qqqCount, "hRawQQQ");

  /*if(QQQ_Events.size()>=1) {
    std::cout<< " ---->" << std::endl;
    for(auto qe: QQQ_Events) {
      std::cout << qe.ch1/16 << " " <<qe.ch2/16 << " " << qe.ch1%16 << " "<< qe.ch2%16 << " " << qe.Energy1 << " " << qe.Energy2 << " " << std::endl;
    }
  }*/

  typedef std::unordered_map<int, std::tuple<int, double, double>> WireEvent; // this stores nearest neighbour wire events, or a 'cluster'
  WireEvent aWireEvents, cWireEvents;                                         // naming for book keeping
  aWireEvents.clear();
  aWireEvents.reserve(24);

  // PC Gain Matching and Filling
  double anodeT = -99999;
  double cathodeT = 99999;
  int anodeIndex = -1;
  int cathodeIndex = -1;
  for (int i = 0; i < pc.multi; i++)
  {
    if (pc.e[i] > 50)
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
      aWireEvents[pc.index[i]] = std::tuple(pc.index[i], pc.e[i], static_cast<double>(pc.t[i]));
    }
    else
    {
      cathodeT = static_cast<double>(pc.t[i]);
      cathodeIndex = pc.index[i] - 24;
      cWireEvents[pc.index[i] - 24] = std::tuple(pc.index[i] - 24, pc.e[i], static_cast<double>(pc.t[i]));
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
  int aIDMax = 0;

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
      // if (aCluster.size() <= 1 && cCluster.size() <= 1)
      //   continue;
      auto [crossover, alpha, apSumE, cpSumE, apMaxE, cpMaxE, apTSMaxE, cpTSMaxE] = pwinstance.FindCrossoverProperties(aCluster, cCluster);
      if (alpha != 9999999 && apSumE != -1)
      {
        // Event PCEvent(crossover,apMaxE,cpMaxE,apTSMaxE,cpTSMaxE);
        // Event PCEvent(crossover,apSumE,cpSumE,apTSMaxE,cpTSMaxE);
        // Event PCEvent(crossover, apSumE, cpMaxE, apTSMaxE, cpTSMaxE); // run12 shows cathode-max and anode-sum provide best dE signals.
        Event PCEvent(crossover, apSumE, cpMaxE, apTSMaxE, cpTSMaxE, aCluster.size(), cCluster.size()); // changed to include cluster size info --VS
        // std::cout << apSumE << " " << crossover.Perp() << " " << apMaxE << " " << apTSMaxE << std::endl;
        PC_Events.push_back(PCEvent);
        sumE_AC.push_back(std::pair(apSumE, cpSumE));
      }
    }
  }
  if (QQQ_Events.size() && PC_Events.size())
    plotter->Fill2D("PCEv_vs_QQQEv", 20, 0, 20, 20, 0, 20, QQQ_Events.size(), PC_Events.size());

  for (auto pcevent : PC_Events)
  {
    for (auto sx3event : sx3Events)
    {
      plotter->Fill1D("dt_pcA_sx3B" + std::to_string(sx3event.ch2), 640, -2000, 2000, sx3event.Time1 - pcevent.Time1);
      plotter->Fill1D("dt_pcC_sx3B" + std::to_string(sx3event.ch2), 640, -2000, 2000, sx3event.Time1 - pcevent.Time2);
      plotter->Fill2D("dE_E_Anodesx3B", 400, 0, 10, 800, 0, 40000, sx3event.Energy1 * 0.001, pcevent.Energy1);

      plotter->Fill2D("dE_E_Cathodesx3B", 400, 0, 10, 800, 0, 10000, sx3event.Energy1 * 0.001, pcevent.Energy2);
      double sx3z = sx3event.pos.Z() + (75.0 / 2.0) - 3.0; // w.r.t target origin at 90 for run12
      double sx3rho = 88.0;                                // approximate barrel radius
      double sx3theta = TMath::ATan2(sx3rho, sx3z - vertexpos);
      double pczguess = pcrad / TMath::Tan(sx3theta) + vertexpos;
      plotter->Fill2D("pcz_vs_sx3pczguess", 300, 0, 200, 150, 0, 200, pczguess, pcevent.pos.Z());
      plotter->Fill2D("pcz_vs_sx3pczguess" + std::to_string(sx3event.ch2), 300, 0, 200, 150, 0, 200, pczguess, pcevent.pos.Z());
      plotter->Fill2D("pcz_vs_sx3z", 300, 0, 200, 150, 0, 200, sx3z, pcevent.pos.Z());
    }
  }

  for (auto pcevent : PC_Events)
  {
    for (auto qqqevent : QQQ_Events)
    {
      plotter->Fill1D("dt_pcA_qqqR", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1);
      plotter->Fill1D("dt_pcC_qqqW", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2);
      plotter->Fill2D("dE_E_AnodeQQQR", 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1);
      plotter->Fill2D("dE_E_CathodeQQQR", 400, 0, 10, 800, 0, 10000, qqqevent.Energy2, pcevent.Energy2);
      double sinTheta = TMath::Sin((qqqevent.pos - TVector3(0, 0, vertexpos)).Theta()) / TMath::Sin((TVector3(51.5, 0, qqqpos) - TVector3(0, 0, vertexpos)).Theta());
      plotter->Fill2D("dE2_E_AnodeQQQR", 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta);
      plotter->Fill2D("dE2_E_CathodeQQQR", 400, 0, 10, 800, 0, 10000, qqqevent.Energy2, pcevent.Energy2 * sinTheta);

      if (qqqevent.pos.Phi() <= pcevent.pos.Phi() + TMath::Pi() / 4. && qqqevent.pos.Phi() >= pcevent.pos.Phi() - TMath::Pi() / 4.)
      {
        plotter->Fill1D("PCZ", 800, -200, 200, pcevent.pos.Z(), "phicut");
        double pcz_guess = pcrad / TMath::Tan((qqqevent.pos - TVector3(0, 0, vertexpos)).Theta()) + vertexpos; // this is ideally kept to be all QQQ+userinput for calibration of pcz
        plotter->Fill2D("pczguess_vs_pc", 300, 0, 200, 150, 0, 200, pcz_guess, pcevent.pos.Z(), "phicut");
        plotter->Fill2D("pczguess_vs_pc_phi=" + std::to_string(qqqevent.pos.Phi() * 180. / M_PI), 300, 0, 200, 150, 0, 200, pcz_guess, pcevent.pos.Z(), "phicut");
        // plotter->Fill1D("PCZ",800,-200,200,pcevent.pos.Z(),"phicut");
      }

      int aSize = pcevent.ch1;
      int cSize = pcevent.ch2;

      if (cSize == 1)
      {
        if (aSize == 1)
          plotter->Fill1D("pcz_a1c1Cluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
        else if (aSize == 2)
          plotter->Fill1D("pcz_a2c1Cluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
        else if (aSize >= 3)
          plotter->Fill1D("pcz_aNc1Cluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
      }
      else if (cSize == 2)
      {
        if (aSize == 1)
          plotter->Fill1D("pcz_a1c2Cluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
        else if (aSize == 2)
          plotter->Fill1D("pcz_a2c2Cluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
        else if (aSize >= 3)
          plotter->Fill1D("pcz_aNc2Cluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
      }
      else if (cSize >= 3)
      {
        if (aSize == 1)
          plotter->Fill1D("pcz_a1cNCluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
        else if (aSize == 2)
          plotter->Fill1D("pcz_a2cNCluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
        else if (aSize >= 3)
          plotter->Fill1D("pcz_aNcNCluster", 600, -300, 300, pcevent.pos.Z(), "hPCzQQQ");
      }
    }
  }
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
        plotter->Fill2D("Anode_vs_CathodeE", 2000, 0, 30000, 2000, 0, 30000, aE, cE, "hGMPC");
        plotter->Fill2D("CathodeMult_V_CathodeE", 6, 0, 6, 2000, 0, 30000, cathodeHits.size(), cE, "hGMPC");
        for (int j = -4; j < 3; j++)
        {
          if ((aIDMax + 24 + j) % 24 == 23 - cID)
          {
            corrcatMax.push_back(std::pair<int, double>(cID, cE));
            cESum += cE;
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
      if (Crossover[aIDMax][corr.first][0].z > 9000000)
        continue;
      if (cESum > 0)
      {
        x += (corr.second) / cESum * Crossover[aIDMax][corr.first][0].x;
        y += (corr.second) / cESum * Crossover[aIDMax][corr.first][0].y;
        z += (corr.second) / cESum * Crossover[aIDMax][corr.first][0].z;
      }
    }
    if (x == 0 && y == 0 && z == 0)
      ;
    // to ignore events with no valid crossover points
    else
      anodeIntersection = TVector3(x, y, z);
    // << "Anode Intersection: " << anodeIntersection.X() << ", " << anodeIntersection.Y() << ", " << anodeIntersection.Z() << std::endl;
  }
  bool PCQQQPhiCut = false;
  // flip the algorithm for cathode 1 multi anode events
  if ((hitPos.Phi() > (anodeIntersection.Phi() - TMath::PiOver4())) && (hitPos.Phi() < (anodeIntersection.Phi() + TMath::PiOver4())))
  {
    PCQQQPhiCut = true;
  }

  for (double Tz = 60; Tz <= 100; Tz += 1.0)
  {
    TVector3 TargetPos(0, 0, Tz);
    if (PCQQQPhiCut && anodeIntersection.Perp() > 0 && anodeIntersection.Z() != 0 && cathodeHits.size() >= 2)
    {
      plotter->Fill2D("Inttheta_vs_QQQtheta_TC" + std::to_string(PCQQQTimeCut) + "_TZ" + std::to_string(Tz), 400, 0, 180, 90, 0, 90, (anodeIntersection - TargetPos).Theta() * 180. / TMath::Pi(), (hitPos - TargetPos).Theta() * 180. / TMath::Pi(), "TPosVariation");
      // plotter->Fill2D("R_ratio_to_Z_ratio" + std::to_string(PCQQQTimeCut) + "_TZ" + std::to_string(Tz), 100, -2, 2, 100, -2, 2, (anodeIntersection - TargetPos).Z()/(hitPos-TargetPos).Z(), ((anodeIntersection - TargetPos).Perp()+2.5)/(hitPos-TargetPos).Perp(), "TPosVariation");
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
      plotter->Fill2D("VertexRecon2C_Z_vs_QQQE_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, 800, 0, 20000, vector_closest_to_z.Z(), qqqenergy, "customVertex");
    }
  }

  for (int i = 0; i < qqq.multi; i++)
  {
    if (anodeIntersection.Perp() > 0)
    { // suppress x,y=0,0 events
      if (PCQQQTimeCut)
      {
        plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");
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

  TVector3 guessVertex(0, 0, 90.); // for run12, subtract anodeIntersection.Z() by ~74.0 seems to work
  // rho=40.0 mm is halfway between the cathodes(rho=42) and anodes(rho=37)
  double pcz_guess = 42.0 / TMath::Tan((hitPos - guessVertex).Theta()) + guessVertex.Z(); // this is ideally kept to be all QQQ+userinput for calibration of pcz
  if (PCQQQTimeCut && PCQQQPhiCut && hitPos.Perp() > 0 && anodeIntersection.Perp() > 0 && cathodeHits.size() >= 2)
  {
    plotter->Fill2D("pczguess_vs_qqqE", 100, 0, 200, 800, 0, 20, pcz_guess, qqqenergy, "pczguess");
    double pczoffset = 30.0;
    // plotter->Fill2D("pczguess_vs_pcz_rad="+std::to_string(hitPos.Perp()),100,0,200,150,0,200,pcz_guess,anodeIntersection.Z(),"pczguess"); //entirely qqq-derived position vs entirely PC derived position
    plotter->Fill2D("pczguess_vs_pcz_phi=" + std::to_string(hitPos.Phi() * 180. / M_PI), 100, 0, 200, 150, 0, 200, pcz_guess, anodeIntersection.Z() + pczoffset, "pczguess"); // entirely qqq-derived position vs entirely PC derived position
    plotter->Fill2D("pczguess_vs_pcz", 100, 0, 200, 150, 0, 200, pcz_guess, anodeIntersection.Z() + pczoffset);
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

  plotter->Fill2D("AnodeMaxE_Vs_Cathode_Sum_Energy", 2000, 0, 30000, 2000, 0, 30000, aEMax, cESum, "hGMPC");
  plotter->Fill1D("Correlated_Cathode_MaxAnode", 6, 0, 5, corrcatMax.size(), "hGMPC");
  plotter->Fill2D("Correlated_Cathode_VS_MaxAnodeEnergy", 6, 0, 5, 2000, 0, 30000, corrcatMax.size(), aEMax, "hGMPC");
  plotter->Fill1D("AnodeHits", 12, 0, 11, anodeHits.size(), "hGMPC");
  plotter->Fill2D("AnodeMaxE_vs_AnodeHits", 12, 0, 11, 2000, 0, 30000, anodeHits.size(), aEMax, "hGMPC");

  if (anodeHits.size() < 1)
  {
    plotter->Fill1D("NoAnodeHits_CathodeHits", 6, 0, 5, cathodeHits.size(), "hGMPC");
  }

  return kTRUE;
}

void MakeVertex::Terminate()
{
  plotter->FlushToDisk();
}
