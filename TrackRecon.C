#define Analyzer_cxx

#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>

#include <utility>
#include <algorithm>

#include "Armory/ClassSX3.h"
#include "Armory/ClassPC1An.h"

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
int padID = 0;

SX3 sx3_contr;
PC pw_contr;
TVector3 hitPos;
bool HitNonZero;

TH1F * hZProj;

void Analyzer::Begin(TTree * /*tree*/){
  TString option = GetOption();

  hZProj = new TH1F("hZProj", "Z Projection", 200, -600, 600);

  sx3_contr.ConstructGeo();
  pw_contr.ConstructGeo();

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
    for( int j = i + 1; j < qqq.multi; j++){
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
 PCHit_1An hitInfo;

  ID.clear();
  int counter=0;
  std::vector<std::pair<int, double>> E; 
  E.clear();

  if( E.size()==3 ){
    float aE = 0;
    float cE = 0;
    bool multi_an =false;
      for(int l=0;l<E.size();l++){
        if(E[l].first<24 && E[l].first!=20 && E[l].first!=12){
          if(!multi_an){
            aE = E[l].second;
          }
          multi_an=true;
        }
        else {
          cE = E[l].second + cE;
        }
      }
    // printf("anode= %d, cathode = %d\n", aID, cID);
  // }
    if( ID[0].first < 1 ) {
      aID = pc.ch[ID[0].second];
      cID = pc.ch[ID[1].second];
    }else{
      cID = pc.ch[ID[0].second];
      aID = pc.ch[ID[1].second];
    }

    hanVScatsum->Fill(aE,cE);
      
    if( HitNonZero){
      pw_contr.CalTrack3( hitPos, hitinfo, cID);
      hZProj->Fill(pw_contr.GetZ0());
    }

  // }
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

  hqqqCoin->Draw("colz");

  //=============================================== pad-6
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hpcCoin->Draw("colz");

  //=============================================== pad-7
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hsx3VpcIndex ->Draw("colz"); 
  // hsx3VpcE->Draw("colz") ;

  //=============================================== pad-8
  padID ++; canvas->cd(padID); canvas->cd(padID)->SetGrid(1);

  hqqqVpcIndex ->Draw("colz");  

  // hqqqVpcE ->Draw("colz");
  //=============================================== pad-9
  padID ++; 

  // canvas->cd(padID)->DrawFrame(-50, -50, 50, 50);
  // hqqqPolar->Draw("same colz pol");

 canvas->cd(padID); canvas->cd(padID)->SetGrid(1);
 hZProj->Draw();
  // hanVScatsum->Draw("colz");

}
