// #include "/home/splitpole/FSUDAQ_Qt6/Aux/fsuReader.h"
// #include "/home/splitpole/FSUDAQ_Qt6/Aux/fsutsReader.h"
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TString.h>
#include <TMath.h>

#include "mapping.h"

#include "ClassAnasen.h"

class PulserChecker {
public:
  PulserChecker(int sn) : SN(sn){
    t0 = 0;
    t1 = 0;
    dt = 0;
    tStart = -1;
    has1stData = false;
    n = 0;
    mean = 0;
    m2 = 0;
  }

  void addTime(unsigned long long time){
    t0 = t1;
    t1 = time;
    if( has1stData ) {
      dt = t1/1e6 - t0/1e6;
      int mod = 1;
      if( dt >= 2* mean && mean > 0 ){ 
        mod = std::round(dt/mean); 
      }
      if( mod > 1 ) return;
      double new_dt = dt/mod;
      
      n ++;
      double delta = new_dt - mean;
      mean += delta / n;

      double delta2 = new_dt - mean;
      m2 += delta * delta2;
      //if( mod > 1 ) printf("%6d | %llu, %llu (dt %.6f, m : %d) | %.6f\n", SN, t1, t0, dt, mod, new_dt);
    }else{
      tStart = time;
      has1stData = true; 
    }
    // printf("%6d | %13llu, %13llu | %.6f | %13llu\n", SN, t1, t0, dt, tStart);
  }

  unsigned long long getTime0() const {return tStart;}
  int                getSN()    const {return SN;}
  double             getDT()    const {return has1stData ? dt : TMath::QuietNaN() ;}
  double             getMean()  const {return mean;}
  double             getSTD()   const { return n < 2 ? 0.0 :sqrt(m2 / (n - 1)); }

  void Print() { printf("%6d (%2d) | %16llu| n : %6d, mean : %10.6f msec, sd : %10.6f msec\n", 
             SN, SN & 0x1F, tStart, n, getMean(), getSTD());} 

private:
  const int SN;
  unsigned long long tStart, t0, t1;
  double dt;
  bool has1stData;

  int n;
  double mean;
  double m2;

};

void script(TString fileName = "", int maxEvent = -1){

  /*
  //+++++++++++++++++++++++++++++++++++++++++++
  FSUReader * reader = new FSUReader("pulsertest_013_15528_QDC_16_000.fsu", 100);

  reader->ScanNumBlock(1, 1);
  reader->PrintHit(10);

  // for( int i = 0; i < 24 ; i++){
  //   printf("#################### agg : %d \n", i);

  //   reader->ReadNextBlock(false, 0, 1);

  //   reader->GetData()->PrintAllData();

  // }


  // TGraph * graph = new TGraph();
  // int count = 0;
  // unsigned long long tree0 = 0;
  // unsigned long numHit = reader->GetHitCount();
  // for( unsigned long i = 0; i < numHit; i++){
  //   //printf("-------------%6lu \n\033[A\r", i);
  //   unsigned long long t1 = reader->GetHit(i).timestamp;
  //   unsigned short ch = reader->GetHit(i).ch;
  //   unsigned short ee = reader->GetHit(i).energy;
  //   //printf("%6lu | %2u, %6u, %16llu\n", i, ch, ee, t1);
  //   count ++;
  //   graph->SetPoint(i, i, t1/1e9);
  // } 
  // graph->Draw("APL");

  // for( int i = 40; i < 60; i++){
  //   reader->GetHit(i).Print();
  // }
  */

  /*
  //+++++++++++++++++++++++++++++++++++++++++++

  TFile * file0 = new TFile("test_002_1000_noTrace.root"); 
  TFile * f1 = new TFile("test_002_1000.root"); 

  TTree * tree0 = (TTree*) file0->Get("tree");
  TTree * t1 = (TTree*) f1->Get("tree");

  unsigned int m0, m1;
  unsigned short ch0[100], ch1[100];
  unsigned long long ts0[100], ts1[100];


  tree0->SetBranchAddress("multi", &m0);
  tree0->SetBranchAddress("ch", ch0);
  tree0->SetBranchAddress("e_t", ts0);

  t1->SetBranchAddress("multi", &m1);
  t1->SetBranchAddress("ch", ch1);
  t1->SetBranchAddress("e_t", ts1);

  unsigned long long nEntries0 = tree0->GetEntries();
  unsigned long long nEntries1 = t1->GetEntries();

  for( unsigned long long ev = 0; ev < nEntries0; ev++){
    tree0->GetEntry(ev);
    t1->GetEntry(ev);


    if( ch0[0] != ch1[0] || ts0[0] != ts1[0] ) {
      printf("======================= %llu \n", ev);
      printf("ch : %u   %u\n", ch0[0], ch1[0]);
      printf("ts : %llu   %llu\n", ts0[0], ts1[0]);
    }
  }

  file0->Close();
  f1->Close();
  */
  /*
  //+++++++++++++++++++++++++++++++++++++++++++
  printf("######### file : %s \n", fileName.Data());
  TFile * file0 = new TFile(fileName); 

  TTree * tree0 = (TTree*) file0->Get("tree");

  const int MAX_MULTI = 1000;
  unsigned long long                evID = 0;
  unsigned int                     multi = 0;
  unsigned short           sn[MAX_MULTI] = {0};
  unsigned short           ch[MAX_MULTI] = {0};
  unsigned short            e[MAX_MULTI] = {0};
//  unsigned short           e2[MAX_MULTI] = {0};
  unsigned long long      e_t[MAX_MULTI] = {0};
//  unsigned short          e_f[MAX_MULTI] = {0};

  tree0->SetBranchAddress("evID", &evID);
  tree0->SetBranchAddress("multi", &multi);
  tree0->SetBranchAddress("sn", sn);
  tree0->SetBranchAddress("ch", ch);
  tree0->SetBranchAddress("e", e);
  tree0->SetBranchAddress("e_t", e_t);

  unsigned long long nEntries = tree0->GetEntries();

  const std::map<int, int> board = {
    {0, 17122},
    {1, 17123},
    {2, 22320},
    {3, 22130},
    {4, 22129},
    {5, 15529},
    {6, 15528},
    {7,   379},
    {8,   409},
    {9,   405},
  };
  const int nBd = board.size();

  PulserChecker * stat[9];
  for( int i = 0; i < nBd; i++){ stat[i] = new PulserChecker(board.at(i)); }

  printf("Number of Entries : %llu\n", nEntries);

  unsigned long long evEnd = nEntries;
  if( maxEvent >= 0 ) evEnd = maxEvent;

  for( unsigned long long ev = 0; ev < evEnd; ev++){
    tree0->GetEntry(ev);

    // printf("sn : %5d | %14llu \n", sn[0], e_t[0]);

    for( int i = 0; i< nBd; i++){
      for( unsigned int j = 0; j < multi; j++){
        if( sn[j] == stat[i]->getSN() ) {
          stat[i]->addTime(e_t[j]);
          printf(" sn : %5d | %16llu \n", sn[j], e_t[j]);
          break;
        }
      }
    }

    if( ev > 1){
      unsigned long long ts1 = e_t[0];
      tree0->GetEntry(ev-1);
      unsigned long long ts0 = e_t[0];

      if( ts1 <= ts0 )printf("--------------- ev : %llu\n", ev);
    }
  }
  file0->Close();

  printf("=========================\n");

  for( int i = 0; i< nBd; i++){ stat[i]->Print(); }

  unsigned long long time0 = stat[0]->getTime0(); 
  for( int i = 1; i< nBd; i++){ 
    if( stat[i]->getTime0() < time0 ) time0 = stat[i]->getTime0();
  }

  printf("=========================\n");
  printf(" S/N      time Offset ns \n");
  for( int i = 0; i< nBd; i++){ 
    printf(" %5d | %16llu \n", stat[i]->getSN(), stat[i]->getTime0() - time0);
  }
  */
  //+++++++++++++++++++++++++++++++++++++++++++
  // TCanvas * c1 = new TCanvas();
  // ANASEN * kaka = new ANASEN();
  // kaka->DrawAnasen();

  TCanvas * c2 = new TCanvas();
  ANASEN * haha = new ANASEN();

  TVector3 pos (0, 10, 50);
  TVector3 dir (1, 0, 0);

  // dir.SetPhi( 10 * TMath::DegToRad());
  // dir.SetTheta( 90 * TMath::DegToRad());

  std::pair<int, int> wireID = haha->FindWireID(pos, dir, true);
  SX3 sx3 = haha->FindSX3Pos(pos, dir, true);

  haha->CalTrack(sx3.hitPos, wireID.first, wireID.second, true);

  //haha->DrawDeducedTrack(sx3.hitPos, wireID.first, wireID.second);

  // haha->DrawTrack(pos, dir);






}