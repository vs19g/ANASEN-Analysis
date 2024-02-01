#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBenchmark.h"

#include "ClassTransfer.h"
#include "ClassAnasen.h"

//======== Gerneate light particle based on reaction
// find out the CalTrack and the real track
// find out the Q-value uncertaintly

void anasenMS(){

  const int numEvent = 100000;

  //Reaction
  TransferReaction transfer;

  transfer.SetA(12, 6, 0);
  transfer.SetIncidentEnergyAngle(10, 0, 0);

  transfer.Seta( 2, 1);
  transfer.Setb( 1, 1);

  std::vector<float> ExAList = {0};
  std::vector<float> ExList = {0, 1, 2};

  double vertexXRange[2] = { 0, 0}; 
  double vertexYRange[2] = { 0, 0}; 
  double vertexZRange[2] = { 0, 0}; 

  //###################################################
  transfer.CalReactionConstant();

  int nExA = ExAList.size();
  int nEx  = ExList.size();

  ANASEN anasen;

  TString saveFileName = "msAnasen.root";
  printf("\e[32m#################################### building Tree in %s\e[0m\n", saveFileName.Data());
  TFile * saveFile = new TFile(saveFileName, "recreate");
  TTree * tree = new TTree("tree", "tree");

  double KEA;
  tree->Branch("beamKEA",     &KEA, "beamKEA/D");

  double thetaCM, phiCM;
  tree->Branch("thetaCM", &thetaCM, "thetaCM/D");
  tree->Branch("phiCM",     &phiCM, "phiCM/D");

  double thetab, phib, Tb;
  double thetaB, phiB, TB;
  tree->Branch("thetab", &thetab, "thetab/D");
  tree->Branch("phib",     &phib, "phib/D");
  tree->Branch("Tb",         &Tb, "Tb/D");
  tree->Branch("thetaB", &thetaB, "thetaB/D");
  tree->Branch("phiB",     &phiB, "phiB/D");
  tree->Branch("TB",         &TB, "TB/D");

  int ExAID;
  double ExA;
  tree->Branch("ExAID", &ExAID, "ExAID/I");
  tree->Branch("ExA",     &ExA, "ExA/D");

  int ExID;
  double Ex;
  tree->Branch("ExID", &ExID, "ExID/I");
  tree->Branch("Ex",     &Ex, "Ex/D");

  double vertexX, vertexY, vertexZ;
  tree->Branch("vX",     &vertexX, "VertexX/D");
  tree->Branch("vY",     &vertexY, "VertexY/D");
  tree->Branch("vZ",     &vertexZ, "VertexZ/D");

  double sx3X, sx3Y, sx3Z;
  tree->Branch("sx3X",     &sx3X, "sx3X/D");
  tree->Branch("sx3Y",     &sx3Y, "sx3Y/D");
  tree->Branch("sx3Z",     &sx3Z, "sx3Z/D");

  int anodeID, cathodeID;
  tree->Branch("aID", &anodeID, "anodeID/I");
  tree->Branch("cID", &cathodeID, "cathodeID/I");

  int sx3ID, sx3Up, sx3Down, sx3Back;
  double sx3ZFrac;
  tree->Branch("sx3ID",   &sx3ID,   "sx3ID/I");
  tree->Branch("sx3Up",   &sx3Up,   "sx3Up/I");
  tree->Branch("sx3Down", &sx3Down, "sx3Down/I");
  tree->Branch("sx3Back", &sx3Back, "sx3Back/I");
  tree->Branch("sx3ZFrac", &sx3ZFrac, "sx3ZFrac/D");

  double reTheta, rePhi;
  tree->Branch("reTheta", &reTheta, "reconstucted_theta/D");
  tree->Branch("rePhi",     &rePhi, "reconstucted_phi/D");


  //========timer
  TBenchmark clock;
  bool shown ;   
  clock.Reset();
  clock.Start("timer");
  shown = false;

  //================================= Calculate event
  for( int i = 0; i < numEvent ; i++){

    ExAID = gRandom->Integer(nExA);
    ExA = ExAList[ExAID];
    transfer.SetExA(ExA);

    ExID = gRandom->Integer(nEx);
    Ex = ExList[ExID];
    transfer.SetExB(Ex);

    thetaCM = TMath::ACos(2 * gRandom->Rndm() - 1) ; 
    phiCM   = (gRandom->Rndm() - 0.5) * TMath::TwoPi();

    //==== Calculate reaction
    TLorentzVector * output = transfer.Event(thetaCM, phiCM);
    TLorentzVector Pb = output[2];
    TLorentzVector PB = output[3];

    thetab = Pb.Theta() * TMath::RadToDeg();
    thetaB = PB.Theta() * TMath::RadToDeg();

    Tb = Pb.E() - Pb.M();
    TB = PB.E() - PB.M();

    phib = Pb.Phi() * TMath::RadToDeg();
    phiB = PB.Phi() * TMath::RadToDeg();

    vertexX = (vertexXRange[1]- vertexXRange[0])*gRandom->Rndm() +  vertexXRange[0];
    vertexY = (vertexYRange[1]- vertexYRange[0])*gRandom->Rndm() +  vertexYRange[0];
    vertexZ = (vertexZRange[1]- vertexZRange[0])*gRandom->Rndm() +  vertexZRange[0];

    TVector3 vertex(vertexX, vertexY, vertexZ);

    TVector3 dir(1, 0, 0);
    dir.SetTheta(thetab * TMath::DegToRad());
    dir.SetPhi(phib * TMath::DegToRad());

    std::pair<int, int> wireID = anasen.FindWireID(vertex, dir, false);
    SX3 sx3 = anasen.FindSX3Pos(vertex, dir, false);    

    anodeID = wireID.first;
    cathodeID = wireID.second;

    sx3ID = sx3.id;
    if( sx3.id >= 0 ){
      sx3Up = sx3.chUp;
      sx3Down = sx3.chDown;
      sx3Back = sx3.chBack;
      sx3ZFrac = sx3.zFrac;

      sx3X = sx3.hitPos.X();
      sx3Y = sx3.hitPos.Y();
      sx3Z = sx3.hitPos.Z();
      
      // for( int i = 0; i < 12; i++){
      //   sx3Index[i] = -1;
      //   if( i == sx3Up ) sx3Index[i] = sx3ID * 12 + sx3Up; 
      //   if( i == sx3Down ) sx3Index[i] = sx3ID * 12 + sx3Down; 
      //   if( i == sx3Back ) sx3Index[i] = sx3ID * 12 + sx3Back; 
      // }
      anasen.CalTrack(sx3.hitPos, wireID.first, wireID.second, false);
  
      reTheta = anasen.GetTrackTheta() * TMath::RadToDeg();
      rePhi = anasen.GetTrackPhi() * TMath::RadToDeg();

    }else{
      sx3Up = -1;
      sx3Down = -1;
      sx3Back = -1;
      sx3ZFrac = TMath::QuietNaN();
      
      sx3X = TMath::QuietNaN();
      sx3Y = TMath::QuietNaN();
      sx3Z = TMath::QuietNaN();

      // for( int i = 0; i < 12; i++){
      //   sx3Index[i] = -1;
      // }
  
      reTheta = TMath::QuietNaN();
      rePhi = TMath::QuietNaN();
    }



    tree->Fill();

    //#################################################################### Timer  
    clock.Stop("timer");
    Double_t time = clock.GetRealTime("timer");
    clock.Start("timer");

    if ( !shown ) {
      if (fmod(time, 10) < 1 ){
        printf( "%10d[%2d%%]| %8.2f sec | expect: %5.1f min \n", i, TMath::Nint((i+1)*100./numEvent), time , numEvent*time/(i+1)/60);
        shown = 1;
      }
    }else{
      if (fmod(time, 10) > 9 ){
        shown = 0;
      }
    }

  }

  tree->Write();
  int count = tree->GetEntries();
  saveFile->Close();

  printf("=============== done. saved as %s. count(hit==1) : %d\n", saveFileName.Data(), count);

}