//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 22 14:02:44 2024 by ROOT version 6.26/04
// from TTree tree/pulsertest_018_10000_noTrace.root
// found on file: pulsertest_018_10000_noTrace.root
//////////////////////////////////////////////////////////

#ifndef ApplyMapping_h
#define ApplyMapping_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TMath.h>
#include <TBenchmark.h>

#include "mapping.h"

#define MAXMULTI 1000
// Header file for the classes stored in the TTree if any.

class Det{

public:
  Det(): multi(0) {Clear(); }

  UShort_t   multi;
  UShort_t   id[MAXMULTI];
  UShort_t   ch[MAXMULTI];
  UShort_t    e[MAXMULTI];
  ULong64_t   t[MAXMULTI];

  void Clear(){

    multi = 0;
    for( int i = 0; i < MAXMULTI; i++){
      id[i] = 0;
      ch[i] = 0;
       e[i] = 0;
       t[i] = 0;
    }

  }

  void Print(){
    printf("=============================== multi : %u\n", multi);
    for( int i = 0; i < multi; i++) {
      printf(" %3d | %2d-%2d %5u %llu \n", i, id[i], ch[i], e[i], t[i]);
    }
  }

};

class ApplyMapping : public TSelector {
public :
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  ULong64_t       evID;
  UInt_t          multi;
  UShort_t         sn[MAXMULTI];   //[multi]
  UShort_t         ch[MAXMULTI];   //[multi]
  UShort_t          e[MAXMULTI];   //[multi]
  UShort_t         e2[MAXMULTI];   //[multi]
  ULong64_t       e_t[MAXMULTI];   //[multi]
  UShort_t        e_f[MAXMULTI];   //[multi]

  // List of branches
  TBranch        *b_event_ID;   //!
  TBranch        *b_multi;   //!
  TBranch        *b_sn;   //!
  TBranch        *b_ch;   //!
  TBranch        *b_e;   //!
  TBranch        *b_e2;   //!
  TBranch        *b_e_t;   //!
  TBranch        *b_e_f;   //!

  ApplyMapping(TTree * /*tree*/ =0) : fChain(0) { }
  virtual ~ApplyMapping() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  void SetRunID(UInt_t runID) { run = runID; }

  ClassDef(ApplyMapping,0);


  //=============================== 
  TFile * saveFile;
  TTree * newTree;
  TString saveFileName;
  int totnumEntry; // of original root

  //tree  
  ULong_t eventID;
  UInt_t run;
  
  Det sx3;
  Det qqq;
  Det pc ;

  //clock   
  TBenchmark clock;
  Bool_t shown;
   
};

#endif

#ifdef ApplyMapping_cxx
void ApplyMapping::Init(TTree *tree){

  if (!tree) return;

  //======================================   
  totnumEntry = tree->GetEntries();
  printf( "===================================== \n");
  printf( "====== total Entry : %d \n", totnumEntry);
  printf( "===================================== \n");


  fChain = (TChain *) tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evID", &evID, &b_event_ID);
  fChain->SetBranchAddress("multi", &multi, &b_multi);
  fChain->SetBranchAddress("sn", sn, &b_sn);
  fChain->SetBranchAddress("ch", ch, &b_ch);
  fChain->SetBranchAddress("e", e, &b_e);
  fChain->SetBranchAddress("e2", e2, &b_e2);
  fChain->SetBranchAddress("e_t", e_t, &b_e_t);
  fChain->SetBranchAddress("e_f", e_f, &b_e_f);

  PrintMapping();

  //================= Formation of file name
  // TString expName = "";

  // fChain->GetListOfFiles()->Print();
  // int numFile = fChain->GetListOfFiles()->GetLast() + 1;   
  // if( numFile > 0 ) {
  //   int oldRunNum = -100;
  //   bool contFlag = false; // is runNumber continue;
  //   for( int i = 0; i < numFile ; i++){
  //     TString name = fChain->GetListOfFiles()->At(i)->GetTitle();

  //     TString prefix = name;

  //     int found = name.Last('/');        
  //     found = name.Index("_", found+1); // find next "_"
  //     name.Remove(0, found+1);
  //     name.Remove(3); // name should be only runID
  //     int runNum = name.Atoi(); // this should give the 3 digit run number 

  //     if( i == 0 ) {
  //       int found = prefix.Last('/');
  //       prefix.Remove(0, found+1);
  //       found = name.Index("_");
  //       prefix.Remove(0, found+1); // this should give the expName;
  //       expName = prefix;
  //       saveFileName = expName + "_" + prefix +  "_run";
  //     }

  //     if( runNum == oldRunNum + 1 ){
  //       int kk = saveFileName.Sizeof();
  //       if( contFlag == false ){
  //         saveFileName.Remove(kk-2); //remove the "-"
  //         saveFileName += "-";
  //       }else{
  //         saveFileName.Remove(kk-5); //remove the runNum and "-"
  //       }
  //       contFlag = true;
  //     }
  //     if( runNum > oldRunNum + 1) contFlag = false;
      
  //     saveFileName += Form("%03d_", runNum);
  //     oldRunNum = runNum;
  //   }
  //   int kk = saveFileName.Sizeof();
  //   saveFileName.Remove(kk-2); // remove the last "-"
  //   saveFileName += ".root";
  // }else{
  //   saveFileName.Form("%s_default.root", expName.Data());
  // }
  saveFileName = "test.root";

  printf("save file name : %s \n", saveFileName.Data());
  printf("---------------------------------------------\n");
  if( saveFileName == ".root" ) gROOT->ProcessLine(".q");
  
  saveFile = new TFile( saveFileName,"recreate");

  //====================== Create tree   
  newTree =  new TTree("tree","tree");
  
  eventID = 0;
  run = 0;
  sx3.Clear();
  qqq.Clear();
  pc.Clear();

  newTree->Branch("eventID",&eventID,"eventID/l"); 
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

  shown = false;
}

Bool_t ApplyMapping::Notify(){

  return kTRUE;
}

#endif // #ifdef ApplyMapping_cxx
