#define ApplyMapping_cxx

#include "ApplyMapping.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <vector>
#include <utility>


void ApplyMapping::Begin(TTree * /*tree*/){

  TString option = GetOption();

  run = option.Atoi();

}

void ApplyMapping::SlaveBegin(TTree * /*tree*/){

  TString option = GetOption();

}

Bool_t ApplyMapping::Process(Long64_t entry){

  b_event_ID->GetEntry(entry);
  b_multi->GetEntry(entry);
  b_sn->GetEntry(entry);
  b_ch->GetEntry(entry);
  b_e->GetEntry(entry);
  b_e_t->GetEntry(entry);

  eventID = evID;

  // printf("======================== %llu, %u\n", evID, multi);

  sx3.multi = 0;
  qqq.multi = 0;
  pc.multi = 0;

  for( unsigned int i = 0; i < multi; i++){

    //printf("%10u/%10u| %5d, %2u, %6u, %14llu\n", i, multi, sn[i], ch[i], e[i], e_t[i] );
    int globalCh = ch[i];
    for( int j = 0; j < nBd; j++){
      if( board.at(j) == sn[i]){
        globalCh += sn[i] > 1000 ? j * 64 :  7*64 + (j-7) * 16;
      }
    }
    unsigned short ID = mapping[globalCh];

    //=================================== sx3
    if( ID < 10000 ) {

      sx3.id[sx3.multi] = ID / 100;
      sx3.ch[sx3.multi] = ID % 100;
      sx3.e[sx3.multi] = e[i];
      sx3.t[sx3.multi] = e_t[i];
      sx3.multi ++;

    }

    //=================================== qqq
    if( 10000 <= ID && ID < 20000 ) {

      qqq.id[qqq.multi] = (ID - 10000) / 100;
      qqq.ch[qqq.multi] = (ID - 10000) % 100;
      qqq.e[qqq.multi] = e[i];
      qqq.t[qqq.multi] = e_t[i];

      qqq.multi ++;

    }

    //=================================== pc
    if( 20000 <= ID && ID < 30000 ) {

      pc.id[pc.multi] = (ID - 20000) / 100;
      pc.ch[pc.multi] = (ID - 20000) % 100;
      pc.e[pc.multi] = e[i];
      pc.t[pc.multi] = e_t[i];

      pc.multi ++;
    }

  }

  //if(entry == 0 ) sx3.Print();

  /************************************************************************/
  saveFile->cd(); //set focus on this file
  newTree->Fill();
  newTree->Write();
  
  printf("%6lu/%6u [%2d%%]\n\033[A\r", eventID, totnumEntry, TMath::Nint((eventID+1)*100./totnumEntry));

  // clock.Stop("timer");
  // Double_t time = clock.GetRealTime("timer");
  // clock.Start("timer");

  // if ( !shown ) {
  //   if (fmod(time, 10) < 1 ){
  //     printf( "%10lu[%2d%%]|%3.0f min %5.2f sec | expect:%5.2f min\n", 
  //           eventID, 
  //           TMath::Nint((eventID+1)*100./totnumEntry),
  //           TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.,
  //           totnumEntry*time/(eventID+1.)/60.);
  //           shown = true;
  //   }
  // }else{
  //   if (fmod(time, 10) > 9 ){
  //     shown = false;
  //   }
  // }


  return kTRUE;
}

void ApplyMapping::SlaveTerminate(){

}

void ApplyMapping::Terminate(){

  // Double_t time = clock.GetRealTime("timer");

  // printf( "%10lu[%2d%%]|%3.0f min %5.2f sec | expect:%5.2f min\n", 
  //           eventID, 
  //           TMath::Nint((eventID+1)*100./totnumEntry),
  //           TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.,
  //           totnumEntry*time/(eventID+1.)/60.);
  //           shown = true;

  saveFile->cd(); //set focus on this file
  newTree->Write(); 

  UInt_t eventNumber = newTree->GetEntries();

  saveFile->Close();

  printf("-------------- done, saved in %s, %u\n", saveFileName.Data(), eventNumber);
  
}
