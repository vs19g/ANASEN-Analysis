#ifndef MakeVertex_h
#define MakeVertex_h

#include <TROOT.h>
#include <TChain.h>
#include <TApplication.h>
#include <TFile.h>
#include <TSelector.h>
#include <iomanip>
#include <vector>         // Required for vectors
#include <utility>        // Required for std::pair

#include "Armory/ClassDet.h"
#include "Armory/ClassPW.h" // YOU ADDED THIS (Correct! Defines Coord)

class MakeVertex : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Det sx3;
   Det qqq;
   Det pc ;
   Det misc;

   ULong64_t       evID;
   UInt_t          run;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_run;   //!
   TBranch        *b_sx3Multi;   //!
   TBranch        *b_sx3ID;   //!
   TBranch        *b_sx3Ch;   //!
   TBranch        *b_sx3E;   //!
   TBranch        *b_sx3T;   //!
   TBranch        *b_qqqMulti;   //!
   TBranch        *b_qqqID;   //!
   TBranch        *b_qqqCh;   //!
   TBranch        *b_qqqE;   //!
   TBranch        *b_qqqT;   //!
   TBranch        *b_pcMulti;   //!
   TBranch        *b_pcID;   //!
   TBranch        *b_pcCh;   //!
   TBranch        *b_pcE;   //!
   TBranch        *b_pcT;   //!
   TBranch        *b_miscMulti;   //!
   TBranch        *b_miscID;   //!
   TBranch        *b_miscCh;   //!
   TBranch        *b_miscE;   //!
   TBranch        *b_miscT;   //!
   TBranch        *b_miscTf;   //!

   // 1. Geometry Cache
   Coord Crossover[24][24][2];

   // 2. Persistent Vectors (REQUIRED for the optimized .cxx to work)
   std::vector<std::pair<int, double>> anodeHits;
   std::vector<std::pair<int, double>> cathodeHits;
   std::vector<std::pair<int, double>> corrcatMax;
   std::vector<std::pair<int, double>> corranoMax;
   std::vector<double> cathodeTimes;
   std::vector<double> anodeTimes;
   
   MakeVertex(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~MakeVertex() { }
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

   ClassDef(MakeVertex,0);
};

#endif

#ifdef MakeVertex_cxx
void MakeVertex::Init(TTree *tree){

   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evID", &evID, &b_eventID);
   fChain->SetBranchAddress("run", &run, &b_run);

   sx3.SetDetDimension(24,12);
   qqq.SetDetDimension(4,32);
   pc.SetDetDimension(2,24);

   fChain->SetBranchAddress("sx3Multi", &sx3.multi, &b_sx3Multi);
   fChain->SetBranchAddress("sx3ID",    &sx3.id, &b_sx3ID);
   fChain->SetBranchAddress("sx3Ch",    &sx3.ch, &b_sx3Ch);
   fChain->SetBranchAddress("sx3E",     &sx3.e, &b_sx3E);
   fChain->SetBranchAddress("sx3T",     &sx3.t, &b_sx3T);
   fChain->SetBranchAddress("qqqMulti", &qqq.multi, &b_qqqMulti);
   fChain->SetBranchAddress("qqqID",    &qqq.id, &b_qqqID);
   fChain->SetBranchAddress("qqqCh",    &qqq.ch, &b_qqqCh);
   fChain->SetBranchAddress("qqqE",     &qqq.e, &b_qqqE);
   fChain->SetBranchAddress("qqqT",     &qqq.t, &b_qqqT);
   fChain->SetBranchAddress("pcMulti",  &pc.multi, &b_pcMulti);
   fChain->SetBranchAddress("pcID",     &pc.id, &b_pcID);
   fChain->SetBranchAddress("pcCh",     &pc.ch, &b_pcCh);
   fChain->SetBranchAddress("pcE",      &pc.e, &b_pcE);
   fChain->SetBranchAddress("pcT",      &pc.t, &b_pcT);
   fChain->SetBranchAddress("miscMulti", &misc.multi, &b_miscMulti);
   fChain->SetBranchAddress("miscID",    &misc.id, &b_miscID);
   fChain->SetBranchAddress("miscCh",    &misc.ch, &b_miscCh);
   fChain->SetBranchAddress("miscE",     &misc.e, &b_miscE);
   fChain->SetBranchAddress("miscT",     &misc.t, &b_miscT);
   fChain->SetBranchAddress("miscf",    &misc.tf, &b_miscTf);
}

Bool_t MakeVertex::Notify(){
   return kTRUE;
}

void MakeVertex::SlaveBegin(TTree * /*tree*/){
   // TString option = GetOption();
}

void MakeVertex::SlaveTerminate(){
}
#endif // #ifdef MakeVertex_cxx
