#define Analyzer_cxx

#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <utility>
#include <algorithm>

#include "Armory/ClassSX3.h"
#include "Armory/ClassPW.h"

#include "TVector3.h"

TH2F * hsx3IndexVE;
TH2F * hqqqIndexVE;
TH2F * hpcIndexVE;

TH2F * hsx3Coin;
TH2F * hqqqCoin;
TH2F * hpcCoin;

TH2F * hqqqPolar;
TH2F * hsx3VpcIndex;
TH2F * hqqqVpcIndex;
TH2F * hqqqVpcE;
TH2F * hsx3VpcE;
TH2F * hanVScatsum;
TH2F * hAnodeHits;
TH1F * hAnodeMultiplicity;


int padID = 0;

SX3 sx3_contr;
PW pw_contr;
TVector3 hitPos;
bool HitNonZero;

TH1F * hZProj;
TCutG *AnCatSum;
bool inCut;

void Analyzer::Begin(TTree * /*tree*/){
  TString option = GetOption();

  hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy",  24*12, 0,  24*12, 400, 0, 5000); hsx3IndexVE->SetNdivisions( -612, "x");
  hqqqIndexVE = new TH2F("hqqqIndexVE", "QQQ index vs Energy; QQQ index ; Energy", 4*2*16, 0, 4*2*16, 400, 0, 5000); hqqqIndexVE->SetNdivisions( -1204, "x");
  hpcIndexVE = new TH2F("hpcIndexVE",   "PC index vs Energy; PC index ; Energy",     2*24, 0,   2*24, 400, 0, 4000); hpcIndexVE->SetNdivisions( -1204, "x");


  hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident",  24*12, 0,  24*12,  24*12, 0, 24*12);
  hqqqCoin = new TH2F("hqqqCoin", "QQQ Coincident", 4*2*16, 0, 4*2*16, 4*2*16, 0, 4*2*16);
  hpcCoin  = new TH2F("hpcCoin",  "PC Coincident",    2*24, 0,   2*24,   2*24, 0, 2*24);

  hqqqPolar = new TH2F("hqqqPolar", "QQQ Polar ID", 16*4, -TMath::Pi(), TMath::Pi(),16, 10, 50);

  hsx3VpcIndex = new TH2F("hsx3Vpcindex", "sx3 vs pc; sx3 index; pc index",  24*12, 0, 24*12, 48, 0, 48);
  hsx3VpcIndex->SetNdivisions( -612, "x");
  hsx3VpcIndex->SetNdivisions( -12, "y");
  hqqqVpcIndex = new TH2F("hqqqVpcindex", "qqq vs pc; qqq index; pc index",  4*2*16, 0, 4*2*16, 48, 0, 48);
  hqqqVpcIndex->SetNdivisions( -612, "x");
  hqqqVpcIndex->SetNdivisions( -12, "y");

  hqqqVpcE = new TH2F("hqqqVpcEnergy", "qqq vs pc; qqq energy; pc energy", 400, 0, 5000, 400, 0, 5000);
  hqqqVpcE->SetNdivisions( -612, "x");
  hqqqVpcE->SetNdivisions( -12, "y");
  
  hsx3VpcE = new TH2F("hsx3VpcEnergy", "sx3 vs pc; sx3 energy; pc energy", 400, 0, 5000, 400, 0, 5000);
  hsx3VpcE->SetNdivisions( -612, "x");
  hsx3VpcE->SetNdivisions( -12, "y");

  hZProj = new TH1F("hZProj", "ZProjection", 600,-600, 600);
  hAnodeHits = new TH2F("hAnodeHits", "Anode vs Anode Energy, Anode ID; Anode E", 24,0 , 23, 400, 0 , 20000);
  hAnodeMultiplicity = new TH1F("hAnodeMultiplicity", "Number of Anodes/Event", 40, 0, 40); 
  hanVScatsum = new TH2F("hanVScatsum", "Anode vs Cathode Sum; Anode E; Cathode E", 400,0 , 10000, 400, 0 , 16000);

  sx3_contr.ConstructGeo();
  pw_contr.ConstructGeo();
  TFile *f1 = new TFile("AnCatSum.root");
  AnCatSum= (TCutG*)f1->Get("AnCatSum");
}




Bool_t Analyzer::Process(Long64_t entry){

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

  sx3.CalIndex();
  qqq.CalIndex();
  pc.CalIndex();

  // sx3.Print();

  //########################################################### Raw data
  // //======================= SX3

  std::vector<std::pair<int, int>> ID; // first = id, 2nd = index
  for( int i = 0; i < sx3.multi; i ++){
    ID.push_back(std::pair<int, int>(sx3.id[i], i));

    hsx3IndexVE->Fill( sx3.index[i], sx3.e[i] );

    for( int j = i+1; j < sx3.multi; j++){
      hsx3Coin->Fill( sx3.index[i], sx3.index[j]);
    }

    for( int j = 0; j < pc.multi; j++){
      hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
      //  if( sx3.ch[index] > 8 ){
      //    hsx3VpcE->Fill( sx3.e[i], pc.e[j] );
      //   }
    }
  }


  if( ID.size() > 0 ){
    std::sort(ID.begin(), ID.end(),  [](const std::pair<int, int> & a, const std::pair<int, int> & b) {
      return a.first < b.first;
    } );
    // printf("##############################\n");
    // for( size_t i = 0; i < ID.size(); i++) printf("%zu | %d %d \n", i, ID[i].first, ID[i].second );

    std::vector<std::pair<int, int>> sx3ID; 
    sx3ID.push_back(ID[0]);
    bool found = false;
    for( size_t i = 1; i < ID.size(); i++){
      if( ID[i].first == sx3ID.back().first) {
        sx3ID.push_back(ID[i]);
        if( sx3ID.size() >= 3) {
          found = true;
        } 
      }else{
        if( !found ){
          sx3ID.clear();
          sx3ID.push_back(ID[i]);
        }
      }
    }

    // printf("---------- sx3ID Multi : %zu \n", sx3ID.size());

    if( found ){
      int sx3ChUp, sx3ChDn, sx3ChBk;
      float sx3EUp, sx3EDn;
      // printf("------ sx3 ID : %d, multi: %zu\n", sx3ID[0].first, sx3ID.size());
      for( size_t i = 0; i < sx3ID.size(); i++ ){
        int index = sx3ID[i].second;
        // printf(" %zu | index %d | ch : %d, energy : %d \n", i, index, sx3.ch[index], sx3.e[index]);


        if( sx3.ch[index] < 8 ){
          if( sx3.ch[index] % 2 == 0) {
            sx3ChDn = sx3.ch[index];
            sx3EDn = sx3.e[index];
          }else{
            sx3ChUp = sx3.ch[index];
            sx3EUp = sx3.e[index];
          }
        }else{
          sx3ChBk = sx3.ch[index];
        }
         for( int j = 0; j < pc.multi; j++){
      // hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
            if( sx3.ch[index] > 8 ){
              hsx3VpcE->Fill( sx3.e[i], pc.e[j] );
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
  for( int i = 0; i < qqq.multi; i ++){
    // for( int j = 0; j < pc.multi; j++){
        // if(pc.index[j]==4){
          hqqqIndexVE->Fill( qqq.index[i], qqq.e[i] );
        // }
    // }
    for( int j = 0; j < qqq.multi; j++){
      if ( j == i ) continue;
      hqqqCoin->Fill( qqq.index[i], qqq.index[j]);
    }


    for( int j = i + 1; j < qqq.multi; j++){
      for( int k = 0; k < pc.multi; k++){
        if(pc.index[k]<24 && pc.e[k]>50 ){
          hqqqVpcE->Fill( qqq.e[i], pc.e[k] );
          //  hpcIndexVE->Fill( pc.index[i], pc.e[i] );
                hqqqVpcIndex->Fill( qqq.index[i], pc.index[j] );

        }
      // }
      }
      // if( qqq.used[i] == true ) continue;
      
      //if( qqq.id[i] == qqq.id[j] && (16 - qqq.ch[i]) * (16 - qqq.ch[j]) < 0  ){ // must be same detector and wedge and ring
      if( qqq.id[i] == qqq.id[j]  ){ // must be same detector 
        
        int chWedge = -1;
        int chRing  =  -1;
        if( qqq.ch[i] < qqq.ch[j]){
          chRing = qqq.ch[j] - 16;
          chWedge = qqq.ch[i];
        }else{
          chRing = qqq.ch[i];
          chWedge = qqq.ch[j] - 16;
        }

        // printf(" ID : %d , chWedge : %d, chRing : %d \n", qqq.id[i], chWedge, chRing);

        double theta = -TMath::Pi()/2 + 2*TMath::Pi()/16/4.*(qqq.id[i]*16 + chWedge +0.5);
        double rho   = 10.+40./16.*(chRing+0.5);
        // if(qqq.e[i]>50){
        hqqqPolar->Fill( theta, rho);
        // }
        // qqq.used[i] = true;
        // qqq.used[j] = true;

        if( !HitNonZero ){
        double x = rho * TMath::Cos(theta);
        double y = rho * TMath::Sin(theta);
        hitPos.SetXYZ(x, y, 23 + 75 + 30);
        HitNonZero = true;
        }
      }
    }
  
   
  }
  // //======================= PC

  std::vector<std::pair<int, double>> anodeHits;
  std::vector<std::pair<int, double>> cathodeHits;  
  int aID = 0;
  int cID = 0;
  int anodeCount = 0;
  float cEMax = 0;
  int cIDMax = 0;
  float cEnextMax = 0;
  int cIDnextMax = 0;
  float aE = 0;
  float cE = 0;

// Define the excluded SX3 and QQQ channels
  std::unordered_set<int> excludeSX3 = {34, 35, 36, 37, 61, 62, 67, 73, 74, 75, 76, 77, 78, 79, 80, 93, 97, 100, 103, 108, 109, 110, 111, 112};
  std::unordered_set<int> excludeQQQ = {0, 17, 109, 110, 111, 112, 113, 119, 127, 128};
  inCut=false;
  for( int i = 0; i < pc.multi; i ++){
    for(int j=0; j<pc.multi;j++){
      if(pc.id[j]==0){
        anodeCount++;
        }
      }
    
    if(pc.e[i]>50 && anodeCount==1){
      // hpcIndexVE->Fill( pc.index[i], pc.e[i] );
      // for( int j = i+1; j < pc.multi; j++){
      //   hpcCoin->Fill( pc.index[i], pc.index[j]);
      // }

      // for (int j=0;j<sx3.multi;j++){
      //   if(excludeSX3.find(sx3.index[j]) == excludeSX3.end()){

        hpcIndexVE->Fill( pc.index[i], pc.e[i] );
        for( int j = i+1; j < pc.multi; j++){
          hpcCoin->Fill( pc.index[i], pc.index[j]);
        }
          // if(pc.e[i]>100){
            if (pc.index[i] < 24 ){
                anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
                // anodeCount++;
            } else if (pc.index[i] >= 24){
                cathodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
            }
          // }  
      //   }
      // }
    // hpcIndexVE->Fill( pc.index[i], pc.e[i] );
      hAnodeMultiplicity->Fill(pc.multi);

      float aESum = 0;
      if (anodeHits.size()>=1 && cathodeHits.size() >= 1){
        for (const auto& anode : anodeHits) {
          float cESum = 0;
          // for(int l=0; l<sx3.multi; l++){
          //   if (sx3.index[l]==80){

              int aID = anode.first;
              float aE = anode.second;
              aESum += aE;

              for (const auto& cathode : cathodeHits) {
                int cID = cathode.first;
                float cE = cathode.second;
                if(cE>cEMax){
                  cEMax = cE;
                  cIDMax = cID;
                }
                if(cE>cEnextMax && cE<cEMax){
                  cEnextMax = cE;
                  cIDnextMax = cID;
                }

                cESum += cE;
              }

           if( AnCatSum->IsInside(aE, cESum)){
              inCut=true;
            }
              // if(inCut){
              hanVScatsum->Fill(aE,cESum);
              hAnodeHits->Fill(aID, aE);
              //  }
          //   }
          // }
        }
      }
    }    
  }
    // hanVScatsum->Fill(aE,cE);
    

    if( HitNonZero){
      pw_contr.CalTrack1( hitPos, aID, cIDMax, cIDnextMax, cEMax, cEnextMax,1);
      hZProj->Fill(pw_contr.GetZ0());
    }

  

  


  //########################################################### Track constrcution


  //############################## DO THE KINEMATICS


  return kTRUE;
}

void Analyzer::Terminate(){

  gStyle->SetOptStat("neiou");
  TCanvas * canvas = new TCanvas("cANASEN", "ANASEN", 2000, 2000); 
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

  canvas->cd(padID)->SetLogz(true);

  hqqqCoin->Draw("colz");

  //=============================================== pad-6
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hpcCoin->Draw("colz");

  //=============================================== pad-7
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  // hsx3VpcIndex ->Draw("colz"); 
  hsx3VpcE->Draw("colz") ;

  //=============================================== pad-8
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  // hqqqVpcIndex ->Draw("colz");  

  hqqqVpcE ->Draw("colz");
  //=============================================== pad-9
  padID ++; 

  // canvas->cd(padID)->DrawFrame(-50, -50, 50, 50);
  // hqqqPolar->Draw("same colz pol");

 canvas->cd(padID); canvas->cd(padID)->SetGrid(1);
//  hZProj->Draw();
  hanVScatsum->Draw("colz");
  // hAnodeHits->Draw("colz");
// hAnodeMultiplicity->Draw();
}
