#define PreAnalysis_cxx

#include "PreAnalysis.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "mapping.h"

TH2F ** hRate;
TH2F ** hEnergy;


int padID = 0;

void PreAnalysis::Begin(TTree * /*tree*/){
  TString option = GetOption();

  //Find the first and last timestamp, to get the run duration

  hRate = new TH2F * [nBd];
  hEnergy = new TH2F * [nBd];
  for( int i = 0; i < nBd; i++){
    if( board.at(i) > 1000) {
      hRate[i] = new TH2F(Form("hRate%d", board.at(i)), Form("Digi-%d; ch; sec", board.at(i)), 64, 0, 64, 100, 0, 100);
      hEnergy[i] = new TH2F(Form("hEnergy%d", board.at(i)), Form("Digi-%d; ch; raw E", board.at(i)), 64, 0, 64, 400, 0, 500);
    }else{
      hRate[i] = new TH2F(Form("hRate%d", board.at(i)), Form("Digi-%d; ch; sec", board.at(i)), 16, 0, 16, 100, 0, 100);
      hEnergy[i] = new TH2F(Form("hEnergy%d", board.at(i)), Form("Digi-%d; ch; raw E", board.at(i)), 16, 0, 16, 400, 0, 5000);
    }    
  }   

  printf("================================ done creating histograms.\n");

}

Bool_t PreAnalysis::Process(Long64_t entry){

  b_multi->GetEntry(entry);
  b_sn->GetEntry(entry);
  b_ch->GetEntry(entry);
  b_e->GetEntry(entry);
  b_e_t->GetEntry(entry);

  // printf("------------- multi: %lu\n", multi);

  for( unsigned int i = 0; i < multi; i++){

    for( int j = 0; j < nBd; j++ ){

      if( sn[i] == board.at(j) ) {
        hRate[j]->Fill(ch[i], e_t[i]/1e9);

        hEnergy[j]->Fill(ch[i], e[i]);
        break;
      }
    }
  } 

  return kTRUE;
}

void PreAnalysis::Terminate(){

  printf("================================ %s\n", __func__);

  gStyle->SetOptStat("neiou");
  TCanvas * canvas = new TCanvas("cANASEN", "Pre-Analysis, ANASEN", 4000, 800); 
  canvas->Divide(10,2);

  for( int i = 0; i < nBd; i++){
    padID++; canvas->cd(padID);
    hRate[i]->Draw("colz");
  }

  for( int i = 0; i < nBd; i++){
    padID++; canvas->cd(padID);
    hEnergy[i]->Draw("colz");
  }

}
