#include "GainMatchQQQ.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <fstream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "Armory/HistPlotter.h"
#include "TVector3.h"

void make_plots();

TH2F *h1, *h2, *h3, *h4 =nullptr;

int main(){
    TFile* inFile = TFile::Open("Cal_checkQQQ.root");
    
    make_plots();
    return 0;
}
