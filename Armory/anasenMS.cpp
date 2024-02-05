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

int main(int argc, char **argv){

  printf("=========================================\n");
  printf("===       ANASEN Monte Carlo          ===\n");
  printf("=========================================\n");  
  
  int numEvent = 1000000;
  if( argc >= 2 ) numEvent = atoi(argv[1]);

  //Reaction
  TransferReaction transfer;

  transfer.SetA(12, 6, 0);
  transfer.SetIncidentEnergyAngle(10, 0, 0);

  transfer.Seta( 2, 1);
  transfer.Setb( 1, 1);

  std::vector<float> ExAList = {0};
  std::vector<float> ExList = {0, 1, 2};

  double vertexXRange[2] = { -5,  5}; // mm
  double vertexYRange[2] = { -5,  5}; 
  double vertexZRange[2] = {-70, 70}; 

  double sigmaSX3_W = -1; // mm, < 0 use mid-point
  double sigmaSX3_L = 5; // mm, < 0 use mid-point
  double sigmaPW_A  = 3; // from 0 to 1.
  double sigmaPW_C  = 3; // from 0 to 1.

  //###################################################

  printf("------------ Vertex :\n");
  printf("X : %7.2f - %7.2f mm\n", vertexXRange[0], vertexXRange[1]);
  printf("Y : %7.2f - %7.2f mm\n", vertexYRange[0], vertexYRange[1]);
  printf("Z : %7.2f - %7.2f mm\n", vertexZRange[0], vertexZRange[1]);
  printf("------------ Uncertainty :\n");
  printf(" SX3 horizontal : %.1f\n", sigmaSX3_W);
  printf(" SX3   vertical : %.1f\n", sigmaSX3_L);
  printf("          Anode : %.1f mm\n", sigmaPW_A);
  printf("        Cathode : %.1f mm\n", sigmaPW_C);

  transfer.CalReactionConstant();

  int nExA = ExAList.size();
  int nEx  = ExList.size();

  ANASEN * anasen = new ANASEN();
  SX3 * sx3 = anasen->GetSX3();    
  PW * pw = anasen->GetPW(); 

  TString saveFileName = "SimAnasen.root";
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

  int anodeID[2], cathodeID[2];
  tree->Branch("aID", anodeID, "anodeID/I");
  tree->Branch("cID", cathodeID, "cathodeID/I");

  double anodeDist[2], cathodeDist[2];
  tree->Branch("aDist",   anodeDist, "anodeDist/D");
  tree->Branch("cDist",  cathodeDist, "cathodeDist/D");

  int sx3ID, sx3Up, sx3Dn, sx3Bk;
  double sx3ZFrac;
  tree->Branch("sx3ID",   &sx3ID, "sx3ID/I");
  tree->Branch("sx3Up",   &sx3Up, "sx3Up/I");
  tree->Branch("sx3Dn",   &sx3Dn, "sx3Dn/I");
  tree->Branch("sx3Bk",   &sx3Bk, "sx3Bk/I");
  tree->Branch("sx3ZFrac", &sx3ZFrac, "sx3ZFrac/D");

  double reTheta, rePhi;
  tree->Branch("reTheta", &reTheta, "reconstucted_theta/D");
  tree->Branch("rePhi",     &rePhi, "reconstucted_phi/D");

  double reTheta1, rePhi1;
  tree->Branch("reTheta1", &reTheta1, "reconstucted_theta1/D");
  tree->Branch("rePhi1",     &rePhi1, "reconstucted_phi1/D");


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

    transfer.CalReactionConstant();

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

   
    pw->FindWireID(vertex, dir, false);
    sx3->FindSX3Pos(vertex, dir, false);   

    PWHitInfo hitInfo = pw->GetHitInfo();

    anodeID[0]   = hitInfo.nearestWire.first;
    cathodeID[0] = hitInfo.nearestWire.second;
    anodeID[1]   = hitInfo.nextNearestWire.first;
    cathodeID[1] = hitInfo.nextNearestWire.second;

    anodeDist[0]   = hitInfo.nearestDist.first;
    cathodeDist[0] = hitInfo.nearestDist.second;
    anodeDist[1]   = hitInfo.nextNearestDist.first;
    cathodeDist[1] = hitInfo.nextNearestDist.second;

    sx3ID = sx3->GetID();
    if( sx3ID >= 0 ){
      sx3Up    = sx3->GetChUp();
      sx3Dn    = sx3->GetChDn();
      sx3Bk    = sx3->GetChBk();
      sx3ZFrac = sx3->GetZFrac();

      //Introduce uncertaity
      // TVector3 hitPos = sx3->GetHitPos();
      TVector3 hitPos = sx3->GetHitPosWithSigma(sigmaSX3_W, sigmaSX3_L);

      sx3X = hitPos.X();
      sx3Y = hitPos.Y();
      sx3Z = hitPos.Z();
      
      pw->CalTrack(hitPos, anodeID[0], cathodeID[0], false);
      reTheta = pw->GetTrackTheta() * TMath::RadToDeg();
      rePhi = pw->GetTrackPhi() * TMath::RadToDeg();

      pw->CalTrack2(hitPos, hitInfo, sigmaPW_A, sigmaPW_C, false);
      reTheta1 = pw->GetTrackTheta() * TMath::RadToDeg();
      rePhi1 = pw->GetTrackPhi() * TMath::RadToDeg();

    }else{
      sx3Up = -1;
      sx3Dn = -1;
      sx3Bk = -1;
      sx3ZFrac = TMath::QuietNaN();
      
      sx3X = TMath::QuietNaN();
      sx3Y = TMath::QuietNaN();
      sx3Z = TMath::QuietNaN();

      // for( int i = 0; i < 12; i++){
      //   sx3Index[i] = -1;
      // }
  
      reTheta = TMath::QuietNaN();
      rePhi = TMath::QuietNaN();  

      reTheta1 = TMath::QuietNaN();
      rePhi1 = TMath::QuietNaN();

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

  delete anasen;

  return 0;

}