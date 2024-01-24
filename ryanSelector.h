//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 22 14:02:44 2024 by ROOT version 6.26/04
// from TTree tree/pulsertest_018_10000_noTrace.root
// found on file: pulsertest_018_10000_noTrace.root
//////////////////////////////////////////////////////////

#ifndef ryanSelector_h
#define ryanSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TMath.h>

#define MAXMULTI 500
// Header file for the classes stored in the TTree if any.

class ryanSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

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

   ryanSelector(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~ryanSelector() { }
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

   ClassDef(ryanSelector,0);


   //=============================== 
   TFile * saveFile;
   TTree * newTree;
   TString saveFileName;
   int totnumEntry; // of original root
  
   //tree  
   ULong_t eventID;
   UShort_t run;
   UInt_t multi;
   UShort_t  snC[MAXMULTI];
   UShort_t  chC[MAXMULTI];
   Float_t    eC[MAXMULTI];
   ULong64_t eC_t[MAXMULTI];
   
};

#endif

#ifdef ryanSelector_cxx
void ryanSelector::Init(TTree *tree){

   if (!tree) return;

   //======================================   
   totnumEntry = tree->GetEntries();
   printf( "===================================== \n");
   printf( "====== total Entry : %d \n", totnumEntry);
   printf( "===================================== \n");


   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evID", &evID, &b_event_ID);
   fChain->SetBranchAddress("multi", &multi, &b_multi);
   fChain->SetBranchAddress("sn", sn, &b_sn);
   fChain->SetBranchAddress("ch", ch, &b_ch);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("e2", e2, &b_e2);
   fChain->SetBranchAddress("e_t", e_t, &b_e_t);
   fChain->SetBranchAddress("e_f", e_f, &b_e_f);



   //====================== Create tree   
   newTree =  new TTree("tree","tree");
   
   eventID = -1;
   run = 0;
   
   newTree->Branch("eventID",&eventID,"eventID/l"); 
   newTree->Branch("run",        &run,"run/i"); 
   newTree->Branch("mutli",      &multi,"mutli/s"); 
   newTree->Branch("sn" ,      snC, "sn/s");
   newTree->Branch("ch" ,      chC, "ch/s");
   newTree->Branch("e"  ,       eC, "energy/F");
   newTree->Branch("e_t" ,    eC_t, "timestamp/l");
}

Bool_t ryanSelector::Notify(){

   return kTRUE;
}

#endif // #ifdef ryanSelector_cxx
