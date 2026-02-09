#define TrackRecon_cxx

#include "TrackRecon.h"
#include "Armory/ClassPW.h"
#include "Armory/HistPlotter.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "TVector3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>

// Global instances
PW pw_contr;
PW pwinstance;
TVector3 hitPos;

// Calibration globals
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
double qqqGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
double qqqCalib[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqCalibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
// TCutg *cutQQQ;

// PC Arrays
double pcSlope[48];
double pcIntercept[48];

HistPlotter *plotter;

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

void TrackRecon::Begin(TTree *tree)
{ // get file name
  std::cout << tree->GetCurrentFile()->GetName() << std::endl;
  // get substring from file name to identify run number

  TString option = GetOption();
  std::string treefilename(tree->GetCurrentFile()->GetName());

  plotter = new HistPlotter(treefilename.substr(0, treefilename.length() - std::string(".root").length()) + "_histograms.root", "TFILE");
  // plotter = new HistPlotter("Analyzer.root", "TFILE");

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
    std::cerr << "Error opening slope_intercept.txt" << std::endl;
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
}

Bool_t TrackRecon::Process(Long64_t entry)
{
  hitPos.Clear();
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

  // QQQ Processing

  int qqqCount = 0;
  int qqqAdjCh = 0;
  // REMOVE WHEN RERUNNING USING THE NEW CALIBRATION FILE
  // for (int i = 0; i < qqq.multi; i++)
  // {
  //   if ((qqq.id[i] == 3 || qqq.id[i] == 1) && qqq.ch[i] < 16)
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
        }
        else
          continue;

        plotter->Fill2D("WedgeE_Vs_RingECal", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");

        for (int k = 0; k < pc.multi; k++)
        {
          plotter->Fill2D("RingCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Anode_Index" + std::to_string(qqq.id[i]), 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k]);
          plotter->Fill2D("RingCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");

          if (pc.index[k] < 24 && pc.e[k] > 50)
          {
            // plotter->Fill2D("QQQ_CalibW_Vs_PC_Energy", 1000, 0, 16, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hCalQQQ");
            // plotter->Fill2D("QQQ_CalibR_Vs_PC_Energy", 1000, 0, 16, 2000, 0, 30000, eRingMeV, pc.e[k], "hCalQQQ");

            // if (tRing - static_cast<double>(pc.t[k]) < 0 && tRing - static_cast<double>(pc.t[k]) > -600)
            // // {
            // //   plotter->Fill2D("QQQ_CalibW_Vs_PC_Energy_Tight", 1000, 0, 16, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hCalQQQ");
            // //   plotter->Fill2D("QQQ_CalibR_Vs_PC_Energy_Tight", 1000, 0, 16, 2000, 0, 30000, eRingMeV, pc.e[k], "hCalQQQ");
            // // }
            // // else
            // // {
            // //   plotter->Fill2D("QQQ_CalibW_Vs_PC_Energy_OffTime", 1000, 0, 16, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hCalQQQ");
            // //   plotter->Fill2D("QQQ_CalibR_Vs_PC_Energy_OffTime", 1000, 0, 16, 2000, 0, 30000, eRingMeV, pc.e[k], "hCalQQQ");
            // // }
            plotter->Fill2D("Timing_Difference_QQQR_Anode_vRing", 1250, -2500, 2500, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
            plotter->Fill2D("DelT_Vs_QQQRingECal", 500, -2500, 2500, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hTiming");
            plotter->Fill2D("CalibratedQQQEvsAnodeE_R", 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
            plotter->Fill2D("CalibratedQQQEvsAnodeE_W", 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
            if (tRing - static_cast<double>(pc.t[k]) < -150) // 27Al
            // if (tRing - static_cast<double>(pc.t[k]) < -75 && tRing - static_cast<double>(pc.t[k]) > -145) // 17F
            {
              PCQQQTimeCut = true;
            }
          }

          if (pc.index[k] >= 24 && pc.e[k] > 50)
          {
            plotter->Fill2D("Timing_Difference_QQQR_Cathode_vRing", 1250, -2500, 2500, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
          }
        }

        double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
        double rho = 50. + 50. / 16. * (chRing + 0.5);

        plotter->Fill2D("QQQPolarPlot", 16 * 4, -TMath::Pi(), TMath::Pi(), 32, 40, 100, theta, rho, "hCalQQQ");
        plotter->Fill2D("QQQCartesianPlot", 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
        plotter->Fill2D("QQQCartesianPlot" + std::to_string(qqq.id[i]), 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
        if (PCQQQTimeCut)
        {
          plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hPCQQQ");
        }
        plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hPCQQQ");

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

  plotter->Fill1D("QQQ_Multiplicity", 10, 0, 10, qqqCount, "hRawQQQ");

  // PC Gain Matching and Filling
  double anodeT = -99999;
  double cathodeT = 99999;
  int anodeIndex = -1;
  int cathodeIndex = -1;
  for (int i = 0; i < pc.multi; i++)
  {
    if (pc.e[i] > 10)
    {
      plotter->Fill2D("PC_Index_Vs_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], static_cast<double>(pc.e[i]), "hRawPC");
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
    }
    else
    {
      cathodeT = static_cast<double>(pc.t[i]);
      cathodeIndex = pc.index[i] - 24;
    }

    if (anodeT != -99999 && cathodeT != 99999)
    {
      for (int j = 0; j < qqq.multi; j++)
      {
        plotter->Fill1D("AC_Time_qqq_coin", 500, -2500, 2500, anodeT - cathodeT, "hTiming");
        plotter->Fill2D("AC_Time_Vs_QQQ_ch", 500, -2500, 2500, 16 * 8, 0, 16 * 8, anodeT - cathodeT, qqq.ch[j], "hTiming");
        plotter->Fill2D("AC_Time_vs_AIndex", 500, -2500, 2500, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
        plotter->Fill2D("AC_Time_vs_CIndex", 500, -2500, 2500, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
        // plotter->Fill1D("AC_Time_A" + std::to_string(anodeIndex) + "_C" + std::to_string(cathodeIndex), 200, -1000, 1000, anodeT - cathodeT, "TimingPC");
      }

      for (int j = 0; j < sx3.multi; j++)
      {
        plotter->Fill1D("AC_Time_sx3_coinc", 500, -2500, 2500, anodeT - cathodeT, "hTiming");
      }

      plotter->Fill1D("AC_Time", 500, -2500, 2500, anodeT - cathodeT, "hTiming");
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
        anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
      else if (pc.index[i] >= 24)
        cathodeHits.push_back(std::pair<int, double>(pc.index[i] - 24, pc.e[i]));
    }
  }

  // std::sort(anodeHits.begin(), anodeHits.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
  //           { return a.second > b.second; });
  // std::sort(cathodeHits.begin(), cathodeHits.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
  //           { return a.second > b.second; });

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

  TVector3 anodeIntersection;
  anodeIntersection.Clear();
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
    // std::cout << "Anode Intersection: " << anodeIntersection.X() << ", " << anodeIntersection.Y() << ", " << anodeIntersection.Z() << std::endl;
  }
  bool PCQQQPhiCut = false;
  // flip the algorithm for cathode 1 multi anode events
  if ((hitPos.Phi() > (anodeIntersection.Phi() - TMath::PiOver4())) && (hitPos.Phi() < (anodeIntersection.Phi() + TMath::PiOver4())))
  {
    PCQQQPhiCut = true;
  }

  for (double AIz = 20; AIz <= 100; AIz += 5.0)
  {
    TVector3 TargetPos(0, 0, AIz);
    if (PCQQQPhiCut && anodeIntersection.Perp() != 0 && cathodeHits.size() >= 2)
      // TVector3 anodePosAtZ(anodeIntersection.X() * (AIz / anodeIntersection.Z()), anodeIntersection.Y() * (AIz / anodeIntersection.Z()), AIz);
      // TVector3 anodePosAtZ(anodeIntersection.X(), anodeIntersection.Y(),anodeIntersection.Z() + AIz);
      plotter->Fill2D("Inttheta_vs_QQQtheta_TC" + std::to_string(PCQQQTimeCut) + "_TZ" + std::to_string(AIz), 180, 0, 180, 90, 0, 90, (anodeIntersection - TargetPos).Theta() * 180. / TMath::Pi(),
                      (hitPos - TargetPos).Theta() * 180. / TMath::Pi(), "TPosVariation");
  }

  if (anodeIntersection.Perp() != 0)
  {
    plotter->Fill1D("PC_Z_Projection", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("Z_Proj_VsDelTime", 600, -300, 300, 200, -2000, 2000, anodeIntersection.Z(), anodeT - cathodeT, "hPCzQQQ");
    plotter->Fill2D("IntPhi_vs_QQQphi", 100, -200, 200, 80, -200, 200, anodeIntersection.Phi() * 180. / TMath::Pi(), hitPos.Phi() * 180. / TMath::Pi(), "hPCQQQ");
    plotter->Fill1D("IntRho", 200, 0, 100, anodeIntersection.Perp(), "hRawPC");
    plotter->Fill2D("Inttheta_vs_QQQtheta", 90, 0, 180, 20, 0, 45, anodeIntersection.Theta() * 180. / TMath::Pi(), hitPos.Theta() * 180. / TMath::Pi(), "hPCQQQ");
    plotter->Fill2D("Inttheta_vs_QQQtheta_TC" + std::to_string(PCQQQTimeCut), 90, 0, 180, 20, 0, 45, anodeIntersection.Theta() * 180. / TMath::Pi(), hitPos.Theta() * 180. / TMath::Pi(), "hPCQQQ");
    plotter->Fill2D("IntPhi_vs_QQQphi_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 100, -200, 200, 80, -200, 200, anodeIntersection.Phi() * 180. / TMath::Pi(), hitPos.Phi() * 180. / TMath::Pi(), "hPCQQQ");
  }
  if (anodeIntersection.Perp() != 0 && cathodeHits.size() >= 2)
    plotter->Fill1D("PC_Z_Projection_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");

  if (anodeIntersection.Perp() != 0 && cathodeHits.size() == 1)
  {
    plotter->Fill1D("PC_Z_proj_1C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("IntersectionPhi_vs_AnodeZ_1C", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hPCzQQQ");
  }

  if (anodeIntersection.Perp() != 0 && cathodeHits.size() == 2)
  {
    plotter->Fill1D("PC_Z_proj_2C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("IntersectionPhi_vs_AnodeZ_2C", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hGMPC");
  }
  if (anodeIntersection.Perp() != 0 && cathodeHits.size() > 2)
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
    plotter->Fill1D("VertexRecon", 600, -300, 300, pw_contr.GetZ0());

    if (PCQQQPhiCut && PCQQQTimeCut)
    {
      if (cathodeHits.size() == 2)
        plotter->Fill1D("VertexRecon_TC_PhiC_2C", 600, -300, 300, pw_contr.GetZ0());
    }
    plotter->Fill1D("VertexRecon_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -300, 300, pw_contr.GetZ0());
  }

  for (int i = 0; i < qqq.multi; i++)
  {
    if (PCQQQTimeCut)
    {
      plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");
    }
    plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");

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
        plotter->Fill2D("Vertex_V_QQQRing", 600, -300, 300, 16, 0, 16, pw_contr.GetZ0(), chRing, "hPCQQQ");
        double phi = TMath::ATan2(anodeIntersection.Y(), anodeIntersection.X()) * 180. / TMath::Pi();
        // while (phi > 180)
        //   phi -= 180;
        // while (phi < -180)
        //   phi += 180;
        plotter->Fill2D("PolarAngle_Vs_QQQWedge" + std::to_string(qqqID), 360, -360, 360, 16, 0, 16, phi, chWedge, "hPCQQQ");
        // plotter->Fill2D("EdE_PC_vs_QQQ_timegate_ls1000"+std::to_string())

        plotter->Fill2D("PC_Z_vs_QQQRing_Det" + std::to_string(qqqID), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCQQQ");
        for (int k = 0; k < pc.multi; k++)
        {
          if (pc.index[k] >= 24)
            continue;
          plotter->Fill2D("CalibratedQQQE_RvsAnodeE_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
          plotter->Fill2D("CalibratedQQQE_WvsAnodeE_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
          plotter->Fill2D("AnodeQQQ_dTimevsdPhi", 200, -2000, 2000, 80, -200, 200, tRing - static_cast<double>(pc.t[k]), (hitPos.Phi() - anodeIntersection.Phi()) * 180. / TMath::Pi(), "hTiming");
          plotter->Fill1D("AnodeQQQ_Time", 200, -2000, 2000, tRing - static_cast<double>(pc.t[k]));
        }
      }
    }
  }
  for (int i = 0; i < sx3.multi; i++)
  {
    // plotting sx3 strip hits vs anode phi
    if (sx3.ch[i] < 8)
      plotter->Fill2D("AnodePhi_vs_SX3Strip", 100, -200, 200, 8 * 24, 0, 8 * 24, anodeIntersection.Phi() * 180. / TMath::Pi(), sx3.id[i] * 8 + sx3.ch[i]);
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

void TrackRecon::Terminate()
{
  plotter->FlushToDisk();
}