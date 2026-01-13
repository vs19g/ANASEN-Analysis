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

void TrackRecon::Begin(TTree * /*tree*/)
{
  TString option = GetOption();
  plotter = new HistPlotter("Analyzer_QQQ.root", "TFILE");

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

      if (Crossover[i][j][0].z < -190 || Crossover[i][j][0].z > 190)
      {
        Crossover[i][j][0].z = 9999999;
      }

      Crossover[i][j][1].x = alpha;
      Crossover[i][j][1].y = 0;
    }
  }

  // Load PC Calibrations
  std::ifstream inputFile("slope_intercept_results.txt");
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
    std::string filename = "qqq_GainMatch.txt";
    std::ifstream infile(filename);
    if (infile.is_open())
    {
      int det, ring, wedge;
      double gainw, gainr;
      while (infile >> det >> ring >> wedge >> gainw >> gainr)
      {
        qqqGain[det][ring][wedge] = gainw;
        qqqGainValid[det][ring][wedge] = (gainw > 0);
      }
      infile.close();
    }
  }
  {
    std::string filename = "qqq_Calib.txt";
    std::ifstream infile(filename);
    if (infile.is_open())
    {
      int det, ring, wedge;
      double slope;
      while (infile >> det >> ring >> wedge >> slope)
      {
        qqqCalib[det][ring][wedge] = slope;
        qqqCalibValid[det][ring][wedge] = (slope > 0);
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
  qqq1000cut = false;
  int qqqCount = 0;
  int qqqAdjCh = 0;

  for (int i = 0; i < qqq.multi; i++)
  {
    plotter->Fill2D("QQQ_Index_Vs_Energy", 16 * 8, 0, 16 * 8, 2000, 0, 16000, qqq.index[i], qqq.e[i], "hRawQQQ");

    if (qqq.e[i] > 100)
    {
      qqqEcut = true;
    }
    if (qqq.e[i] > 1000)
      qqq1000cut = true;

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
        plotter->Fill2D("QQQ_Vs_PC_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k]);
        plotter->Fill2D("QQQ_Index_Vs_PC_Index", 16 * 8, 0, 16 * 8, 24, 0, 24, qqq.index[i], pc.index[k]);
      }
      else if(pc.index[k] >= 24 && pc.e[k] > 50)
      {
        plotter->Fill2D("QQQ_Vs_PC_Energy_Cathode", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k]);
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

        plotter->Fill1D("Wedgetime_Vs_Ringtime", 2000, -1000, 1000, tWedge - tRing, "hCalQQQ");

        if (qqqCalibValid[qqq.id[i]][chRing][chWedge])
        {
          eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chRing][chWedge] / 1000;
          eRingMeV = eRing * qqqCalib[qqq.id[i]][chRing][chWedge] / 1000;
        }
        else
          continue;

        plotter->Fill2D("WedgeE_Vs_RingECal", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");

        for (int k = 0; k < pc.multi; k++)
        {
          if (pc.index[k] < 24 && pc.e[k] > 50)
          {
            plotter->Fill2D("QQQ_CalibW_Vs_PC_Energy", 1000, 0, 16, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hCalQQQ");
            plotter->Fill2D("QQQ_CalibR_Vs_PC_Energy", 1000, 0, 16, 2000, 0, 30000, eRingMeV, pc.e[k], "hCalQQQ");
            if (tRing - static_cast<double>(pc.t[k]) < 0 && tRing - static_cast<double>(pc.t[k]) > -600)
            {
              plotter->Fill2D("QQQ_CalibW_Vs_PC_Energy_Tight", 1000, 0, 16, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hCalQQQ");
              plotter->Fill2D("QQQ_CalibR_Vs_PC_Energy_Tight", 1000, 0, 16, 2000, 0, 30000, eRingMeV, pc.e[k], "hCalQQQ");
            }
            else
            {
              plotter->Fill2D("QQQ_CalibW_Vs_PC_Energy_OffTime", 1000, 0, 16, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hCalQQQ");
              plotter->Fill2D("QQQ_CalibR_Vs_PC_Energy_OffTime", 1000, 0, 16, 2000, 0, 30000, eRingMeV, pc.e[k], "hCalQQQ");
            }
            plotter->Fill2D("Timing_Difference_QQQ_PC", 20000, -1000, 1000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hCalQQQ");
            plotter->Fill2D("DelT_Vs_QQQRingECal", 20000, -1000, 1000, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hCalQQQ");
          }
        }

        double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
        double rho = 50. + 40. / 16. * (chRing + 0.5);

        plotter->Fill2D("QQQPolarPlot", 16 * 4, -TMath::Pi(), TMath::Pi(), 32, 40, 100, theta, rho, "hCalQQQ");

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

  plotter->Fill1D("QQQ_Multiplicity", 10, 0, 10, qqqCount, "hCalQQQ");

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
      plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 24, 0, 24, 2000, 0, 30000, pc.index[i], pc.e[i], "hGMPC");
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
        plotter->Fill1D("PC_Time_qqq", 200, -1000, 1000, anodeT - cathodeT, "hGMPC");
        plotter->Fill2D("PC_Time_Vs_QQQ_ch", 200, -1000, 1000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, qqq.ch[j], "hGMPC");
        plotter->Fill2D("PC_Time_vs_AIndex", 200, -1000, 1000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hGMPC");
        plotter->Fill2D("PC_Time_vs_CIndex", 200, -1000, 1000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hGMPC");
        plotter->Fill1D("PC_Time_A" + std::to_string(anodeIndex) + "_C" + std::to_string(cathodeIndex), 200, -1000, 1000, anodeT - cathodeT, "TimingPC");
      }

      for (int j = 0; j < sx3.multi; j++)
      {
        plotter->Fill1D("PC_Time_sx3", 200, -1000, 1000, anodeT - cathodeT, "hGMPC");
      }

      plotter->Fill1D("PC_Time", 200, -1000, 1000, anodeT - cathodeT, "hGMPC");
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

  if (anodeHits.size() == 2 && cathodeHits.size() >= 1)
  {
    // 2. CRITICAL FIX: Define reference vector 'a'
    // In Analyzer.cxx, 'a' was left over from the loop. We use the first anode wire as reference here.
    // (Assuming pwinstance.An is populated and wires are generally parallel).
    TVector3 refAnode = pwinstance.An[0].first - pwinstance.An[0].second;

    if (((TMath::TanH(hitPos.Y() / hitPos.X())) > (TMath::TanH(refAnode.Y() / refAnode.X()) - TMath::PiOver4())) ||
        ((TMath::TanH(hitPos.Y() / hitPos.X())) < (TMath::TanH(refAnode.Y() / refAnode.X()) + TMath::PiOver4())))
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
  if (qqq1000cut)
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
    anodeIntersection = TVector3(x, y, z);
  }

  if (anodeIntersection.Z() != 0)
  {
    plotter->Fill1D("PC_Z_Projection", 600, -300, 300, anodeIntersection.Z(), "hGMPC");
    plotter->Fill2D("Z_Proj_VsDelTime", 600, -300, 300, 200, -1000, 1000, anodeIntersection.Z(), anodeT - cathodeT, "hGMPC");
  }

  if (anodeIntersection.Z() != 0 && cathodeHits.size() == 1)
  {
    plotter->Fill1D("PC_Z_proj_1C", 600, -300, 300, anodeIntersection.Z(), "hGMPC");
  }
  if (anodeIntersection.Z() != 0 && cathodeHits.size() == 2)
  {
    plotter->Fill1D("PC_Z_proj_2C", 600, -300, 300, anodeIntersection.Z(), "hGMPC");
  }
if(anodeHits.size()>0 && cathodeHits.size()>0)
  plotter->Fill2D("AHits_vs_CHits", 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");


  for (int i = 0; i < qqq.multi; i++)
  {
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
        int qqqID = -1;
        if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
        {
          chWedge = qqq.ch[i];
          eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
          chRing = qqq.ch[j] - 16;
          eRing = qqq.e[j];
          qqqID = qqq.id[i];
        }
        else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
        {
          chWedge = qqq.ch[j];
          eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
          chRing = qqq.ch[i] - 16;
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
          plotter->Fill2D("PC_Z_vs_QQQRing", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hGMPC");
        }

        // plotter->Fill2D("EdE_PC_vs_QQQ_timegate_ls1000"+std::to_string())

        plotter->Fill2D("PC_Z_vs_QQQRing_Det" + std::to_string(qqqID), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hGMPC");
      }
    }
  }

  if (anodeIntersection.Z() != 0 && cathodeHits.size() == 3)
  {
    plotter->Fill1D("PC_Z_proj_3C", 600, -300, 300, anodeIntersection.Z(), "hGMPC");
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

  if (HitNonZero && anodeIntersection.Z() != 0)
  {
    pw_contr.CalTrack2(hitPos, anodeIntersection);
    plotter->Fill1D("VertexRecon", 600, -300, 300, pw_contr.GetZ0(), "hGMPC");
  }

  return kTRUE;
}

void TrackRecon::Terminate()
{
  plotter->FlushToDisk();
}