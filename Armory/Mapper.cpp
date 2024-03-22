#include <string>
#include <cstdio>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TBenchmark.h>

#include "../mapping_alpha.h"
#include "ClassDet.h"

//===============================
int main(int argc, char **argv){

  printf("=========================================\n");
  printf("===               Mapper              ===\n");
  printf("=========================================\n");  
  if (argc != 2)    {
    printf("Incorrect number of arguments:\n");
    printf("%s [inFile]\n", argv[0]);
    printf("\n\n");
    return 1;
  }

  ///============= read input
  std::string inFileName = argv[1];

  PrintMapping();

  TFile * inFile = new TFile(inFileName.c_str(), "READ");
  TTree * tree = (TTree*) inFile->Get("tree");
  unsigned long long totnumEntry = tree->GetEntries();

  ULong64_t       evID;
  UInt_t          multi;
  UShort_t         sn[MAXMULTI];
  UShort_t         ch[MAXMULTI];
  UShort_t          e[MAXMULTI];
  UShort_t         e2[MAXMULTI];
  ULong64_t       e_t[MAXMULTI];
  UShort_t        e_f[MAXMULTI];

  tree->SetBranchAddress("evID",  &evID);
  tree->SetBranchAddress("multi", &multi);
  tree->SetBranchAddress("sn",    sn);
  tree->SetBranchAddress("ch",    ch);
  tree->SetBranchAddress("e",     e);
  tree->SetBranchAddress("e2",    e2);
  tree->SetBranchAddress("e_t",   e_t);
  tree->SetBranchAddress("e_f",   e_f);

  ///================== new tree
  TString outFileName = inFileName;
  TString runStr = outFileName; 
  int pos = outFileName.Last('/');
  pos = outFileName.Index("_", pos+1); // find next "_"
  runStr.Remove(0, pos+1);
  runStr.Remove(3);
  pos = outFileName.Index("_", pos+1); // find next "_"
  outFileName.Remove(pos); // remove the rest
  outFileName += "_mapped.root";

  ULong_t eventID;
  UInt_t run = runStr.Atoi();
  
  Det sx3;
  Det qqq;
  Det pc ;

  printf(" Raw root file : %s\n", inFileName.c_str());
  printf("           Run : %03d\n", run);
  printf("   total Entry : %lld \n", totnumEntry);
  printf(" Out file name : %s \n", outFileName.Data());

  TFile * saveFile = new TFile( outFileName,"RECREATE");
  TTree * newTree =  new TTree("tree","tree");
  

  newTree->Branch("evID",   &eventID,"eventID/l"); 
  newTree->Branch("run",        &run,"run/i");  

  newTree->Branch("sx3Multi", &sx3.multi, "sx3Multi/s");
  newTree->Branch("sx3ID",    &sx3.id,    "sx3ID[sx3Multi]/s");
  newTree->Branch("sx3Ch",    &sx3.ch,    "sx3Ch[sx3Multi]/s");
  newTree->Branch("sx3E",     &sx3.e,     "sx3Energy[sx3Multi]/s");
  newTree->Branch("sx3T",     &sx3.t,     "sx3Time[sx3Multi]/l");

  newTree->Branch("qqqMulti", &qqq.multi, "qqqMulti/s");
  newTree->Branch("qqqID",    &qqq.id,    "qqqID[qqqMulti]/s");
  newTree->Branch("qqqCh",    &qqq.ch,    "qqqCh[qqqMulti]/s");
  newTree->Branch("qqqE",     &qqq.e,     "qqqEnergy[qqqMulti]/s");
  newTree->Branch("qqqT",     &qqq.t,     "qqqTime[qqqMulti]/l");

  newTree->Branch("pcMulti", &pc.multi, "pcMulti/s");
  newTree->Branch("pcID",    &pc.id,    "pcID[pcMulti]/s");
  newTree->Branch("pcCh",    &pc.ch,    "pcCh[pcMulti]/s");
  newTree->Branch("pcE",     &pc.e,     "pcEnergy[pcMulti]/s");
  newTree->Branch("pcT",     &pc.t,     "pcTime[pcMulti]/l");

  ///================== looping old tree and apply mapping

  //clock   
  // TBenchmark clock;
  // Bool_t shown;

  for( unsigned long long ev = 0; ev < totnumEntry; ev++){
    tree->GetEntry(ev);

    eventID = evID;
    sx3.multi = 0;
    qqq.multi = 0;
    pc.multi = 0;

    qqq.Clear();

    for( unsigned int i = 0; i < multi; i++){

      //printf("%10u/%10u| %5d, %2u, %6u, %14llu\n", i, multi, sn[i], ch[i], e[i], e_t[i] );
      
      //globalCh = digi-ID * nCh(digi-iD) + ch
      int globalCh = -1;

      for( int j = 0; j < nBd; j++){
        if( board.at(j) == sn[i]){
          globalCh = (sn[i] > 1000 ? j * 64 :  7*64 + (j-7) * 16) + ch[i]; //& = number V1740
          break;
        }
      }
      
      if( globalCh == -1) printf("ev %llu\n", ev);

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
    
    saveFile->cd(); //set focus on this file
    newTree->Fill();
    
    if( eventID % 100 == 0 ) printf("%6lu/%6llu [%2d%%]\n\033[A\r", eventID, totnumEntry, TMath::Nint((eventID+1)*100./totnumEntry));

  }

  inFile->Close();

  saveFile->cd(); //set focus on this file
  newTree->Write(); 
  UInt_t eventNumber = newTree->GetEntries();
  saveFile->Close();
  printf("-------------- done,  %u\n", eventNumber);
  
  return 0;

}