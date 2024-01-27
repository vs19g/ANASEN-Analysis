#ifndef PreAnalysis_h
#define PreAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#define MAXMULTI 1000

class PreAnalysis : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   ULong64_t       evID;
   UInt_t          multi;
   UShort_t        sn[MAXMULTI];   //[multi]
   UShort_t        ch[MAXMULTI];   //[multi]
   UShort_t        e[MAXMULTI];   //[multi]
   UShort_t        e2[MAXMULTI];   //[multi]
   ULong64_t       e_t[MAXMULTI];   //[multi]
   UShort_t        e_f[MAXMULTI];   //[multi]
   Bool_t          pileUp[MAXMULTI];   //[multi]

   // List of branches
   TBranch        *b_event_ID;   //!
   TBranch        *b_multi;   //!
   TBranch        *b_sn;   //!
   TBranch        *b_ch;   //!
   TBranch        *b_e;   //!
   TBranch        *b_e2;   //!
   TBranch        *b_e_t;   //!
   TBranch        *b_e_f;   //!
   TBranch        *b_pileUp;   //!

   PreAnalysis(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~PreAnalysis() { }
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

   ClassDef(PreAnalysis,0);
};

#endif

#ifdef PreAnalysis_cxx
void PreAnalysis::Init(TTree *tree){

   // Set branch addresses and branch pointers
   if (!tree) return;
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
   fChain->SetBranchAddress("pileUp", pileUp, &b_pileUp);
}

Bool_t PreAnalysis::Notify(){

   return kTRUE;
}

void PreAnalysis::SlaveBegin(TTree * /*tree*/){

   TString option = GetOption();

}

void PreAnalysis::SlaveTerminate(){

}

#endif // #ifdef PreAnalysis_cxx
