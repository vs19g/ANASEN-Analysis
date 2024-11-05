#define Analyzer_cxx

#include "Analyzer.h"
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <algorithm>
#include <utility>

#include "Armory/ClassPW.h"
#include "Armory/ClassSX3.h"

#include "TVector3.h"

TH2F *hsx3IndexVE;
TH2F *hqqqIndexVE;
TH2F *hpcIndexVE;

TH2F *hsx3Coin;
TH2F *hqqqCoin;
TH2F *hpcCoin;

TH2F *hqqqPolar;
TH2F *hsx3VpcIndex;
TH2F *hqqqVpcIndex;
TH2F *hqqqVpcE;
TH2F *hsx3VpcE;
TH2F *hanVScatsum;
TH2F *hICvsSi;
TH2F *hAnodeHits;
TH2F *hSiEvsMCPt;
TH2F *hRfvsMCPt;
TH1F *hAnodeHits1d;
TH1F *hPCMultiplicity;
TH1F *hRFtime;
TH1F *hSi;
TH1F *hSi_gated;
TH1F *hSiMCPt;

int padID = 0;

SX3 sx3_contr;
PW pw_contr;
TVector3 hitPos;
bool HitNonZero;

TH1F *hZProj;
TCutG *PCCoinc;

TCutG *alpha_cut_up;
TCutG *alpha_cut_down;
TCutG *cutg;
bool inCut;
bool inCutUp;
bool inCutDown;
bool inCutG;

void Analyzer::Begin(TTree * /*tree*/) {
  TString option = GetOption();

  hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
  hsx3IndexVE->SetNdivisions(-612, "x");
  hqqqIndexVE = new TH2F("hqqqIndexVE", "QQQ index vs Energy; QQQ index ; Energy", 4 * 2 * 16, 0, 4 * 2 * 16, 400, 0, 5000);
  hqqqIndexVE->SetNdivisions(-1204, "x");
  hpcIndexVE = new TH2F("hpcIndexVE", "PC index vs Energy; PC index ; Energy", 2 * 24, 0, 2 * 24, 6400, 0, 30000);
  hpcIndexVE->SetNdivisions(-1204, "x");

  hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);
  hqqqCoin = new TH2F("hqqqCoin", "QQQ Coincident", 4 * 2 * 16, 0, 4 * 2 * 16, 4 * 2 * 16, 0, 4 * 2 * 16);
  hpcCoin = new TH2F("hpcCoin", "PC Coincident", 2 * 24, 0, 2 * 24, 2 * 24, 0, 2 * 24);

  hqqqPolar = new TH2F("hqqqPolar", "QQQ Polar ID", 16 * 4, -TMath::Pi(), TMath::Pi(), 16, 10, 50);

  hsx3VpcIndex = new TH2F("hsx3Vpcindex", "sx3 vs pc; sx3 index; pc index", 24 * 12, 0, 24 * 12, 48, 0, 48);
  hsx3VpcIndex->SetNdivisions(-612, "x");
  hsx3VpcIndex->SetNdivisions(-12, "y");
  hqqqVpcIndex = new TH2F("hqqqVpcindex", "qqq vs pc; qqq index; pc index", 4 * 2 * 16, 0, 4 * 2 * 16, 48, 0, 48);
  hqqqVpcIndex->SetNdivisions(-612, "x");
  hqqqVpcIndex->SetNdivisions(-12, "y");

  hqqqVpcE = new TH2F("hqqqVpcEnergy", "qqq vs pc; qqq energy; pc energy", 400, 0, 5000, 6400, 0, 30000);
  hqqqVpcE->SetNdivisions(-612, "x");
  hqqqVpcE->SetNdivisions(-12, "y");

  hsx3VpcE = new TH2F("hsx3VpcEnergy", "sx3 vs pc; sx3 energy; pc energy", 400, 0, 5000, 6400, 0, 30000);
  hsx3VpcE->SetNdivisions(-612, "x");
  hsx3VpcE->SetNdivisions(-12, "y");

  hZProj = new TH1F("hZProj", "ZProjection", 600, -600, 600);
  hAnodeHits1d = new TH1F("hAnodeHits1d", "Anode Hits", 24, 0, 24);
  hAnodeHits = new TH2F("hAnodeHits", "Anode vs Anode Energy, Anode ID; Anode E", 24, 0, 23, 400, 0, 30000);
  hPCMultiplicity = new TH1F("hPCMultiplicity", "Number of PC/Event", 40, 0, 40);
  hanVScatsum = new TH2F("hanVScatsum", "Anode vs Cathode Sum; Anode E; Cathode E", 6400, 0, 30000, 6400, 0, 30000);
  hICvsSi = new TH2F("hICvsSi", "IC vs Si; Si E; IC E", 800, 0, 20000, 400, 0, 8000);
  hSi = new TH1F("hSi", "Si E", 800, 0, 20000);
  hSi_gated = new TH1F("hSi_gated", "Si E", 800, 0, 20000);
  hRFtime = new TH1F("hRFtime", "Rf-MCP time(ns)", 3000, -3000, 3000);
  hSiEvsMCPt = new TH2F("hSiEsMCPt", "Si E vs MCP time; SI E; MCP time", 800, 0, 20000, 3000, -3000, 3000);
  hSiMCPt = new TH1F("hSiMCPt", "Si vs MCP time", 1500, -3000, 3000);
  hRfvsMCPt = new TH2F("hRfvsMCPt", "RF vs MCP time; RF(ns) ; MCP time(ns)", 1000, -2000, 2000, 1000, -2000, 2000);

  sx3_contr.ConstructGeo();
  pw_contr.ConstructGeo();
  TFile *f1 = new TFile("PCCoinc.root");
  PCCoinc = (TCutG *)f1->Get("PCCoinc");
  TFile *f2 = new TFile("alpha_cut_up.root");
  alpha_cut_up = (TCutG *)f2->Get("alpha_cut_up");
  TFile *f3 = new TFile("alpha_cut_down.root");
  alpha_cut_down = (TCutG *)f3->Get("alpha_cut_down");
  TFile *f4 = new TFile("CUTG.root");
  cutg = (TCutG *)f4->Get("CUTG");

  // TFile *f1 = new TFile("AnCatSum.root");
  // AnCatSum= (TCutG*)f1->Get("AnCatSum");
}

Bool_t Analyzer::Process(Long64_t entry) {

  // if ( entry > 100 ) return kTRUE;

  hitPos.Clear();
  HitNonZero = false;

  // if( entry > 1) return kTRUE;
  // printf("################### ev : %llu \n", entry);

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
  b_miscCh->GetEntry(entry);
  b_miscE->GetEntry(entry);
  b_miscID->GetEntry(entry);
  b_miscMulti->GetEntry(entry);
  b_miscT->GetEntry(entry);
  b_miscTf->GetEntry(entry);

  sx3.CalIndex();
  qqq.CalIndex();
  pc.CalIndex();

  // sx3.Print();

  // ########################################################### Raw data
  //  //======================= SX3

  std::vector<std::pair<int, int>> ID; // first = id, 2nd = index
  for (int i = 0; i < sx3.multi; i++) {
    if (sx3.e[i] > 50) {
      ID.push_back(std::pair<int, int>(sx3.id[i], i));

      hsx3IndexVE->Fill(sx3.index[i], sx3.e[i]);

      for (int j = i + 1; j < sx3.multi; j++) {
        hsx3Coin->Fill(sx3.index[i], sx3.index[j]);
      }

      // for( int j = 0; j < pc.multi; j++){
      // hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
      //  if( sx3.ch[index] > 8 ){
      //    hsx3VpcE->Fill( sx3.e[i], pc.e[j] );
      //   }
      // }
    }
  }

  if (ID.size() > 0) {
    std::sort(ID.begin(), ID.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b) { return a.first < b.first; });
    // printf("##############################\n");
    // for( size_t i = 0; i < ID.size(); i++) printf("%zu | %d %d \n", i, ID[i].first, ID[i].second );

    std::vector<std::pair<int, int>> sx3ID;
    sx3ID.push_back(ID[0]);
    bool found = false;
    for (size_t i = 1; i < ID.size(); i++) {
      if (ID[i].first == sx3ID.back().first) {
        sx3ID.push_back(ID[i]);
        if (sx3ID.size() >= 3) {
          found = true;
        }
      } else {
        if (!found) {
          sx3ID.clear();
          sx3ID.push_back(ID[i]);
        }
      }
    }

    // printf("---------- sx3ID Multi : %zu \n", sx3ID.size());

    if (found) {
      int sx3ChUp, sx3ChDn, sx3ChBk;
      float sx3EUp, sx3EDn;
      // printf("------ sx3 ID : %d, multi: %zu\n", sx3ID[0].first, sx3ID.size());
      for (size_t i = 0; i < sx3ID.size(); i++) {
        int index = sx3ID[i].second;
        // printf(" %zu | index %d | ch : %d, energy : %d \n", i, index, sx3.ch[index], sx3.e[index]);

        if (sx3.ch[index] < 8) {
          if (sx3.ch[index] % 2 == 0) {
            sx3ChDn = sx3.ch[index];
            sx3EDn = sx3.e[index];
          } else {
            sx3ChUp = sx3.ch[index];
            sx3EUp = sx3.e[index];
          }
        } else {
          sx3ChBk = sx3.ch[index];
        }
        for (int j = 0; j < pc.multi; j++) {
          // hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
          if (sx3.ch[i] > 8 && pc.index[j] < 24 && pc.e[j] > 50) {
            hsx3VpcE->Fill(sx3.e[i], pc.e[j]);
            // printf(" sx3 Ch: %d, pc Ch: %d , : %d\n", sx3.index[i], pc.index[j], sx3.t[i] - pc.t[j]);
            //  hpcIndexVE->Fill( pc.index[i], pc.e[i] );
          }
        }
      }

      sx3_contr.CalSX3Pos(sx3ID[0].first, sx3ChUp, sx3ChDn, sx3ChBk, sx3EUp, sx3EDn);
      hitPos = sx3_contr.GetHitPos();
      HitNonZero = true;
      // hitPos.Print();
    }
  }

  // //======================= QQQ
  for (int i = 0; i < qqq.multi; i++) {
    if (qqq.e[i] > 50) {

      // for( int j = 0; j < pc.multi; j++){
      // if(pc.index[j]==4){
      hqqqIndexVE->Fill(qqq.index[i], qqq.e[i]);
      // }
      // }
      for (int j = 0; j < qqq.multi; j++) {
        if (j == i)
          continue;
        hqqqCoin->Fill(qqq.index[i], qqq.index[j]);
      }

      for (int j = 0; j < pc.multi; j++) {
        if (pc.index[j] < 24 && pc.e[j] > 50) {
          hqqqVpcE->Fill(qqq.e[i], pc.e[j]);
          //  hpcIndexVE->Fill( pc.index[i], pc.e[i] );
          hqqqVpcIndex->Fill(qqq.index[i], pc.index[j]);
        }
      }
      // }
      // if( qqq.used[i] == true ) continue;
      for (int j = i + 1; j < qqq.multi; j++) {

        // if( qqq.id[i] == qqq.id[j] && (16 - qqq.ch[i]) * (16 - qqq.ch[j]) < 0  ){ // must be same detector and wedge and ring
        if (qqq.id[i] == qqq.id[j]) { // must be same detector

          int chWedge = -1;
          int chRing = -1;
          if (qqq.ch[i] < qqq.ch[j]) {
            chRing = qqq.ch[j] - 16;
            chWedge = qqq.ch[i];
          } else {
            chRing = qqq.ch[i];
            chWedge = qqq.ch[j] - 16;
          }

          // printf(" ID : %d , chWedge : %d, chRing : %d \n", qqq.id[i], chWedge, chRing);

          double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
          double rho = 10. + 40. / 16. * (chRing + 0.5);
          // if(qqq.e[i]>50){
          hqqqPolar->Fill(theta, rho);
          // }
          // qqq.used[i] = true;
          // qqq.used[j] = true;

          if (!HitNonZero) {
            double x = rho * TMath::Cos(theta);
            double y = rho * TMath::Sin(theta);
            hitPos.SetXYZ(x, y, 23 + 75 + 30);
            HitNonZero = true;
          }
        }
      }
    }
  }
  // //======================= PC

  std::vector<std::pair<int, double>> anodeHits;
  std::vector<std::pair<int, double>> cathodeHits;
  int aID = 0;
  int cID = 0;
  float cEMax = 0;
  int cIDMax = 0;
  float cEnextMax = 0;
  int cIDnextMax = 0;
  float aE = 0;
  float cE = 0;

  // Define the excluded SX3 and QQQ channels
  std::unordered_set<int> excludeSX3 = {34, 35, 36, 37, 61, 62, 67, 73, 74, 75, 76, 77, 78, 79, 80, 93, 97, 100, 103, 108, 109, 110, 111, 112};
  std::unordered_set<int> excludeQQQ = {0, 17, 109, 110, 111, 112, 113, 119, 127, 128};

  for (int i = 0; i < pc.multi; i++) {
    // for(int j=0; j<pc.multi;j++){
    //   if(pc.id[j]==0){
    //     anodeCount++;
    //     }
    //   }

    if (pc.e[i] > 100 & pc.multi < 7) {
      // hpcIndexVE->Fill( pc.index[i], pc.e[i] );
      // for( int j = i+1; j < pc.multi; j++){
      //   hpcCoin->Fill( pc.index[i], pc.index[j]);
      // }

      // for (int j=0;j<sx3.multi;j++){
      //   if(excludeSX3.find(sx3.index[j]) == excludeSX3.end()){

      // hpcIndexVE->Fill( pc.index[i], pc.e[i] );

      for (int j = i + 1; j < pc.multi; j++) {
        inCut = false;
        if (PCCoinc->IsInside(pc.index[i], pc.index[j])) {
          inCut = true;
        }
        // hpcCoin->Fill( pc.index[i], pc.index[j]);
      }
      // if(pc.e[i]>100){
      if (pc.index[i] < 24) {
        anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
        // anodeCount++;
      } else if (pc.index[i] >= 24) {
        cathodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
      }
      // }
      //   }
      // }
      // hpcIndexVE->Fill( pc.index[i], pc.e[i] );
    }
  }
  hPCMultiplicity->Fill(pc.multi);

  float aESum = 0;
  float cESum = 0;
  if (anodeHits.size() == 1 && cathodeHits.size() >= 1) {

    inCutDown = false;
    inCutUp = false;
    for (const auto &anode : anodeHits) {

      // for(int l=0; l<sx3.multi; l++){
      //   if (sx3.index[l]==80){

      aID = anode.first;
      aE = anode.second;
      aESum += aE;
      // printf("aID : %d, aE : %f, cE : %f\n", aID, aE, cE);
    }

    for (const auto &cathode : cathodeHits) {
      cID = cathode.first;
      cE = cathode.second;
      cESum += cE;

      if (cE > cEMax) {
        cEMax = cE;
        cIDMax = cID;
      }
      if (cE > cEnextMax && cE < cEMax) {
        cEnextMax = cE;
        cIDnextMax = cID;
      }
    }

    if (alpha_cut_down->IsInside(aE, cESum)) {
      inCutDown = true;
    }
    if (alpha_cut_up->IsInside(aE, cESum)) {
      inCutUp = true;
    }

    // if (inCutUp)
    // {
    for (int i = 0; i < pc.multi; i++) {
      for (int j = i + 1; j < pc.multi; j++) {
        hpcCoin->Fill(pc.index[i], pc.index[j]);
        hpcIndexVE->Fill(pc.index[i], pc.e[i]);
      }
    }
    // }
    // if (inCut) {
    hanVScatsum->Fill(aE, cESum);
    hAnodeHits->Fill(aID, aE);
    hAnodeHits1d->Fill(anodeHits.size());
    // }
    // }
  }

  // Miscellaneous channels including the Lollipop IC and Si detectors and hot needle IC
  // Misc ch 0,1, 2, 3, 4 in order are the LIC, LSi, HNIC-difference, MCP, and Rf
  bool timing = false;
  inCutG = false;
  double SiE = 0;
  double SiT = 0;
  double MCPt = 0;
  double MCPE = 0;
  double Rft = 0;
  double ICt = 0;
  double ICe = 0;
  double SiCFDt = 0;
  for (int i = 0; i < misc.multi; i++) {
    // if (misc.ch[i] == 1 && misc.e[i] > 10000 && misc.e[i] < 15000) {
    // if(misc.e[i] > 7500 && misc.e[i]<15000) hSi->Fill(misc.e[i]);

    if (misc.ch[i] == 1) {
      // hSi->Fill(misc.e[i]);
      SiE = misc.e[i];
      SiT = misc.t[i] + misc.tf[i] * 4. / 1000;
      // hSi->Fill(misc.e[i]);
    }
    if (misc.ch[i] == 2) {
      ICt = misc.t[i] + misc.tf[i] * 4. / 1000;
      ICe = misc.e[i];
      hSi->Fill(misc.e[i]);
    }
    if (misc.ch[i] == 3) {
      // only analyze the first MCP in any event
      if (MCPt == 0) {
        MCPt = misc.t[i] + misc.tf[i] * 4. / 1000;
        MCPE = misc.e[i];
      }
    }
    if (misc.ch[i] == 4) {
      // only analyze the first RF in any event
      if (Rft == 0) {
        Rft = misc.t[i] + misc.tf[i] * 4. / 1000;
      }
    }
    if (misc.ch[i] == 5) {
      if (SiCFDt == 0) {
        SiCFDt = misc.t[i] + misc.tf[i] * 4. / 1000;
      }
    }

    // hSiEvsMCPt1->Fill(SiE, Rft-MCPt);
    // hSiEvsMCPt->Fill(ICe, MCPt - Rft);
    if (MCPt != 0 && Rft != 0) {
      // if (SiE > 10200 && SiE < 12200) {
      // hRfvsMCPt->Fill(Rft - ICt, MCPt - ICt);

      hSiMCPt->Fill(MCPt - ICt);
      // if(misc.ch[i] == 2 && misc.e[i] > 1000 && misc.e[i]<2000)
      hRFtime->Fill(Rft - ICt);
      // }
      // printf("RF time : %lld %lld %lld\n", Rft, MCPt, (MCPt - Rft));
      // }
    }
    // inCutG = true;
    // if (misc.ch[i] == 1) hSi->Fill(misc.e[i]);

    // for (int j = 0; j < qqq.multi; j++) {
      // if (pc.id[j] == 0) {
        hRfvsMCPt->Fill(Rft-ICt, MCPt -ICt);
        hSiEvsMCPt->Fill(ICe, MCPt - ICt);
      // }
    // }
    for (int j = i + 1; j < misc.multi; j++) {
      //   if (cutg->IsInside(misc.e[i], misc.e[j])) {
      //     inCutG = true;
      //   })

      if (misc.ch[j] == 4 && misc.ch[i] == 3) {
        // hRFtime->Fill(misc.t[j]*1. + misc.tf[j] * 4. / 1000 - (misc.t[i]*1. + misc.tf[i] * 4. / 1000));

        if (misc.t[j] + misc.tf[j] * 4. / 1000 - (misc.t[i] + misc.tf[i] * 4. / 1000) > 20 && misc.t[j] + misc.tf[j] * 4. / 1000 - (misc.t[i] + misc.tf[i] * 4. / 1000) < 100) {
          timing = true;
        }
        // printf("RF time : %lld %lld %lld %lld %lld\n", misc.t[i], misc.t[j], misc.tf[i], misc.tf[j], (misc.t[j]*1000 + misc.tf[j]*4 - (misc.t[i]*1000 + misc.tf[i]*4)));
      }
    }

    // for (int j = i + 1; j < misc.multi; j++) {
    if (timing == true) {
      // hICvsSi->Fill(misc.e[i], misc.e[j]);
      if (misc.ch[i] == 1) {
        hSi_gated->Fill(misc.e[i]);
        // }
      }
      // }
    }
  }

  if (HitNonZero) {
    // pw_contr.CalTrack1( hitPos, aID, cIDMax, cIDnextMax, cEMax, cEnextMax,1);

    // pw_contr.CalTrack(hitPos, aID, cID);
    hZProj->Fill(pw_contr.GetZ0());
  }

  // ########################################################### Track constrcution

  // ############################## DO THE KINEMATICS

  return kTRUE;
}

void Analyzer::Terminate() {

  gStyle->SetOptStat("neiou");
  TCanvas *canvas = new TCanvas("cANASEN", "ANASEN", 2000, 2000);
  // TCanvas *a = new TCanvas("aANASEN", "ANASEN", 800, 600);
  canvas->Divide(3, 3);
  // hRFtime->Draw();
  // TCanvas *b = new TCanvas("bANASEN", "ANASEN", 800, 600);
  // // hICvsSi->Draw("colz");
  // hSi->Draw();

  // =============================================== pad-1
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  hsx3IndexVE->Draw("colz");

  //=============================================== pad-2
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  hqqqIndexVE->Draw("colz");

  //=============================================== pad-3
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  hpcIndexVE->Draw("colz");

  //=============================================== pad-4
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  hsx3Coin->Draw("colz");

  //=============================================== pad-5
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  canvas->cd(padID)->SetLogz(true);

  hqqqCoin->Draw("colz");

  //=============================================== pad-6
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  hpcCoin->Draw("colz");

  //=============================================== pad-7
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  // hsx3VpcIndex ->Draw("colz");
  hsx3VpcE->Draw("colz");

  //=============================================== pad-8
  padID++;
  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);

  // hqqqVpcIndex ->Draw("colz");

  hqqqVpcE->Draw("colz");
  //=============================================== pad-9
  padID++;

  // canvas->cd(padID)->DrawFrame(-50, -50, 50, 50);
  // hqqqPolar->Draw("same colz pol");

  canvas->cd(padID);
  canvas->cd(padID)->SetGrid(1);
  //  hZProj->Draw();
  hanVScatsum->Draw("colz");
  // hAnodeHits->Draw("colz");
  // // hAnodeMultiplicity->Draw();
}
