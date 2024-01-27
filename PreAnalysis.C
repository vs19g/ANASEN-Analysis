#define PreAnalysis_cxx

#include "PreAnalysis.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>

#include "mapping.h"

TH1F *** rate;

void PreAnalysis::Begin(TTree * /*tree*/){
   TString option = GetOption();

   //Find the first and last timestamp, to get the run duration

   rate = new TH1F ** [nBd];
   for( int i = 0; i < nBd; i++){
      if( board.at(i) > 1000) {
         rate[i] = new TH1F * [64];
         for(int j = 0; j < 64; j ++){
            rate[i][j] = new TH1F(Form("h%d_%d", board.at(i), j), Form("Digi-%d, ch-%d", board.at(i), j), 100, 0, 100);
         }
      }else{
         rate[i] = new TH1F * [16];
         for(int j = 0; j < 16; j ++){
            rate[i][j] = new TH1F(Form("h%d_%d", board.at(i), j), Form("Digi-%d, ch-%d", board.at(i), j), 100, 0, 100);
         }
      }
   }
   
}

Bool_t PreAnalysis::Process(Long64_t entry){


   return kTRUE;
}

void PreAnalysis::Terminate(){

}
