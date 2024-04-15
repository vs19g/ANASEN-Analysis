#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TTree.h"


void processRun(TString runFile){

    TFile * file = new  TFile(runFile);

    TTree * tree = (TTree *) file->Get("tree");

    tree->Process("Analyzer.C+");

}