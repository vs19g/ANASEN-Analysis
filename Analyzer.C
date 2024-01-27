#define Analyzer_cxx

#include "Analyzer.h"
#include <TH2.h>
#include <TStyle.h>

TH2F * hsx3IDVch;
// TH2F * hqqqIDVch;
// TH2F * hpcIDVch;

TH2F * hsx3VpcIndex;

// TH2F * hsx3EVIndex;

void Analyzer::Begin(TTree * /*tree*/){
   TString option = GetOption();

   hsx3IDVch = new TH2F("hsx3IDVch", "sx3 ID vs ch; ch ; ID", 24, 0, 24, 12, 0, 12);
   // hqqqIDVch = new TH2F("hqqqIDVch", "qqq ID vs ch; ch ; ID",  4, 0,  4, 32, 0, 32);
   // hpcIDVch  = new TH2F("hpcIDVch",  "pc ID vs ch; ch ; ID",   2, 0,  2, 24, 0, 24);

   hsx3VpcIndex = new TH2F("hsx3Vpcindex", "sx3 vs pc",  24*12, 0, 24*12, 48, 0, 48);

}

Bool_t Analyzer::Process(Long64_t entry){

   if( entry > 1) return kTRUE;
   printf("################### ev : %llu \n", entry);

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

   sx3.Print();

   // //======================= SX3
   // for( int i = 0; i < sx3.multi; i ++){
   //    for( int j = 0; j < pc.multi; j++){
   //       hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
   //    }
   // }


   return kTRUE;
}

void Analyzer::Terminate(){


}
