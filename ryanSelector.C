#define ryanSelector_cxx

#include "ryanSelector.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include <utility>

const std::map<int, unsigned short> board = {
   {0, 17122},  // id, sn
   {1, 17123},
   {2, 22320},
   {3, 22130},
   {4, 22129},
   {5, 15529},
   {6, 15528},
   {7,   379},
   {8,   409},
   {9,   405}
};
const int nBd = board.size();


void ryanSelector::Begin(TTree * /*tree*/){

   TString option = GetOption();


}

void ryanSelector::SlaveBegin(TTree * /*tree*/){

   TString option = GetOption();

}

Bool_t ryanSelector::Process(Long64_t entry){

   b_sn->GetEntry(entry);
   b_e_t->GetEntry(entry);


   return kTRUE;
}

void ryanSelector::SlaveTerminate(){

}

void ryanSelector::Terminate(){

   // TCanvas * canvas = new TCanvas("c1", "c1", 800, 600);

}
