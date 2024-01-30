#define Analyzer_cxx

#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

TH2F * hsx3IndexVE;
TH2F * hqqqIndexVE;
TH2F * hpcIndexVE;

TH2F * hsx3Coin;
TH2F * hqqqCoin;
TH2F * hpcCoin;

TH2F * hqqqPolar;
TH2F * hsx3VpcIndex;

int padID = 0;


void Analyzer::Begin(TTree * /*tree*/){
  TString option = GetOption();

  hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy",  24*12, 0,  24*12, 400, 0, 500); hsx3IndexVE->SetNdivisions( -612, "x");
  hqqqIndexVE = new TH2F("hqqqIndexVE", "QQQ index vs Energy; QQQ index ; Energy", 4*2*16, 0, 4*2*16, 400, 0, 500); hqqqIndexVE->SetNdivisions( -1204, "x");
  hpcIndexVE = new TH2F("hpcIndexVE",   "PC index vs Energy; PC index ; Energy",     2*24, 0,   2*24, 400, 0, 4000); hpcIndexVE->SetNdivisions( -1204, "x");


  hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident",  24*12, 0,  24*12,  24*12, 0, 24*12);
  hqqqCoin = new TH2F("hqqqCoin", "QQQ Coincident", 4*2*16, 0, 4*2*16, 4*2*16, 0, 4*2*16);
  hpcCoin  = new TH2F("hpcCoin",  "PC Coincident",    2*24, 0,   2*24,   2*24, 0, 2*24);

  hqqqPolar = new TH2F("hqqqPolar", "QQQ Polar ID", 16*4, -TMath::Pi(), TMath::Pi(),16, 10, 50);

  hsx3VpcIndex = new TH2F("hsx3Vpcindex", "sx3 vs pc; sx3 index; pc index",  24*12, 0, 24*12, 48, 0, 48);
  hsx3VpcIndex->SetNdivisions( -612, "x");
  hsx3VpcIndex->SetNdivisions( -12, "y");

}

Bool_t Analyzer::Process(Long64_t entry){

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

  sx3.CalIndex();
  qqq.CalIndex();
  pc.CalIndex();

  // sx3.Print();

  //########################################################### Raw data
  // //======================= SX3
  for( int i = 0; i < sx3.multi; i ++){

    hsx3IndexVE->Fill( sx3.index[i], sx3.e[i] );

    for( int j = i+1; j < sx3.multi; j++){
      hsx3Coin->Fill( sx3.index[i], sx3.index[j]);
    }

    for( int j = 0; j < pc.multi; j++){
      hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
    }
  }

  // //======================= QQQ
  for( int i = 0; i < qqq.multi; i ++){
    hqqqIndexVE->Fill( qqq.index[i], qqq.e[i] );

    for( int j = i + 1; j < qqq.multi; j++){
      hqqqCoin->Fill( qqq.index[i], qqq.index[j]);

      if( qqq.used[i] == true ) continue;
      
      if( qqq.id[i] == qqq.id[j] && (16 - qqq.ch[i]) * (16 - qqq.ch[j]) < 0 ){ // must be same detector and wedge and ring
        
        int chWedge = qqq.ch[i];
        int chRing  = qqq.ch[j] - 16;
        if( qqq.ch[i] >= 16 ) {
          chWedge = qqq.ch[j];
          chRing  = qqq.ch[i] - 16;
        }

        //printf(" ID : %d , chWedge : %d, chRing : %d \n", qqq.id[i], chWedge, chRing);

        double theta = -TMath::Pi() + 2*TMath::Pi()/16/4.*(qqq.id[i]*16 + chWedge +0.5);
        double rho   = 10.+40./16.*(chRing+0.5);

        hqqqPolar->Fill( theta, rho);

        qqq.used[i] = true;
        qqq.used[j] = true;

      }
    }
  }

  // //======================= PC
  for( int i = 0; i < pc.multi; i ++){
    hpcIndexVE->Fill( pc.index[i], pc.e[i] );

    for( int j = i+1; j < pc.multi; j++){
      hpcCoin->Fill( pc.index[i], pc.index[j]);
    }
  }

  //########################################################### Track constrcution



  return kTRUE;
}

void Analyzer::Terminate(){

  gStyle->SetOptStat("neiou");
  TCanvas * canvas = new TCanvas("cANASEN", "ANASEN", 1000, 1000); 
  canvas->Divide(3,3);

  //hsx3VpcIndex->Draw("colz");

  //=============================================== pad-1
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hsx3IndexVE->Draw("colz");

  //=============================================== pad-2
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hqqqIndexVE->Draw("colz");

  //=============================================== pad-3
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hpcIndexVE->Draw("colz");

  //=============================================== pad-4
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hsx3Coin->Draw("colz");

  //=============================================== pad-5
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hqqqCoin->Draw("colz");

  //=============================================== pad-6
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hpcCoin->Draw("colz");

  //=============================================== pad-7
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  canvas->cd(padID)->DrawFrame(-50, -50, 50, 50);
  hqqqPolar->Draw("same colz pol");  

}
