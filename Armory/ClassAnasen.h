#ifndef ClassAnasen_h
#define ClassAnasen_h

#include <cstdio>
#include <TMath.h>
#include <TVector3.h>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TCanvas.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TRandom.h"

#include "ClassSX3.h"
#include "ClassPW.h"

class ANASEN{
public:
  ANASEN();
  ~ANASEN();

  void DrawTrack(TVector3 pos, TVector3 direction, bool drawEstimatedTrack = false);
  void DrawDeducedTrack(TVector3 sx3Pos, int anodeID, int cathodeID);
  void DrawAnasen(int anodeID1 = -1, 
                  int anodeID2 = -1, 
                  int cathodeID1 = -1, 
                  int cathodeID2 = -1, 
                  int sx3ID = -1, 
                  bool DrawQQQ = false );

private:

  PW pw;
  SX3 sx3;

  const float qqqR1 = 50;
  const float qqqR2 = 100;
  const float qqqZPos = 23 + 75 + 30; 

  void CalGeometry();

  TGeoManager *geom;
  TGeoVolume *worldBox;

  void Construct3DModel(int anodeID1 = -1, 
                        int anodeID2 = -1, 
                        int cathodeID1 = -1, 
                        int cathodeID2 = -1, 
                        int sx3ID = -1, 
                        bool DrawQQQ = true);

};

//!==============================================
inline ANASEN::ANASEN(){

  CalGeometry();

  geom = nullptr;
  worldBox = nullptr;

}

inline ANASEN::~ANASEN(){

  delete geom;

}
//!==============================================
inline void ANASEN::CalGeometry(){

  sx3.ConstructGeo();
  pw.ConstructGeo();

}

inline void ANASEN::Construct3DModel(int anodeID1, int anodeID2, int cathodeID1, int cathodeID2, int sx3ID, bool DrawQQQ ){

  if( geom ) delete geom;

  // Create ROOT manager and master volume
  geom = new TGeoManager("Detector", "ANASEN");

  //--- define some materials
  TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
  //--- define some media
  TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);

  //--- make the top container volume 
  Double_t worldx = 200.; //mm
  Double_t worldy = 200.; //mm
  Double_t worldz = 200.; //mm
  worldBox = geom->MakeBox("ROOT", Vacuum, worldx, worldy, worldz);
  geom->SetTopVolume(worldBox);

  //--- making axis
  TGeoVolume *axisX = geom->MakeTube("axisX", Al, 0, 0.1, 5.);
  axisX->SetLineColor(1);
  worldBox->AddNode(axisX, 1, new TGeoCombiTrans(5, 0, 0., new TGeoRotation("rotA", 90., 90., 0.)));

  TGeoVolume *axisY = geom->MakeTube("axisY", Al, 0, 0.1, 5.);
  axisY->SetLineColor(1);
  worldBox->AddNode(axisY, 1, new TGeoCombiTrans(0, 5, 0., new TGeoRotation("rotB", 0., 90., 0.)));

  TGeoVolume *axisZ = geom->MakeTube("axisZ", Al, 0, 0.1, 5.);
  axisZ->SetLineColor(1);
  worldBox->AddNode(axisZ, 1, new TGeoTranslation(0, 0,  5));

  //.......... convert to wire center dimensions
  TGeoVolume *pcA = geom->MakeTube("tub1", Al, 0, 0.01, pw.GetAnodeLength()/2);
  pcA->SetLineColor(4);  

  int startID = 0;
  int endID = pw.GetNumWire() - 1;

  if( anodeID1 >= 0 && anodeID2 >= 0 ){
    startID = anodeID1;
    endID = anodeID2;
    if( anodeID1 > anodeID2 ) {
      endID = pw.GetNumWire() + anodeID2;
    }
  }

  for( int i = startID; i <= endID; i++){
    TVector3 a = pw.GetAnodneMid(i);
    double wireTheta = pw.GetAnodeTheta(i) * TMath::RadToDeg();
    double wirePhi = pw.GetAnodePhi(i) * TMath::RadToDeg() + 90;

    worldBox->AddNode(pcA, i+1, new TGeoCombiTrans( a.X(), 
                                                    a.Y(), 
                                                    a.Z(), 
                                                    new TGeoRotation("rot1", wirePhi, wireTheta, 0.)));
  }

  TGeoVolume *pcC = geom->MakeTube("tub2", Al, 0, 0.01, pw.GetCathodeLength()/2);
  pcC->SetLineColor(6);

  startID = 0;
  endID = pw.GetNumWire() - 1;

  if( cathodeID1 >= 0 && cathodeID2 >= 0 ){
    startID = cathodeID1;
    endID = cathodeID2;
    if( cathodeID1 > cathodeID2 ) {
      endID = pw.GetNumWire() + cathodeID2;
    }
  }

  for( int i = startID; i <= endID; i++){
    TVector3 a = pw.GetCathodneMid(i);
    double wireTheta = pw.GetCathodeTheta(i) * TMath::RadToDeg();
    double wirePhi = pw.GetCathodePhi(i) * TMath::RadToDeg() + 90;

    worldBox->AddNode(pcC, i+1, new TGeoCombiTrans( a.X(), 
                                                    a.Y(), 
                                                    a.Z(), 
                                                    new TGeoRotation("rot1", wirePhi , wireTheta, 0.)));
  }

  TGeoVolume * sx3Det = geom->MakeBox("box", Al, 0.1, sx3.GetWidth()/2, sx3.GetLength()/2);
  sx3Det->SetLineColor(kGreen+3);

  for( int i = 0; i < sx3.GetNumDet(); i++){
    if( sx3ID != -1 && i != sx3ID ) continue;
    TVector3 aUp = sx3.GetUpMid(i); // center of the SX3 upstream
    TVector3 aDn = sx3.GetDnMid(i); // center of the SX3 Downstream
    double phi = sx3.GetDetPhi(i) * TMath::RadToDeg() + 90;

    worldBox->AddNode(sx3Det, 2*i+1., new TGeoCombiTrans( aUp.X(), 
                                                          aUp.Y(), 
                                                          aUp.Z(), 
                                                          new TGeoRotation("rot1", phi, 0., 0.)));
    worldBox->AddNode(sx3Det, 2*i+1., new TGeoCombiTrans( aDn.X(), 
                                                          aDn.Y(), 
                                                          aDn.Z(), 
                                                          new TGeoRotation("rot1", phi, 0., 0.)));
  }

  if( DrawQQQ ){
    TGeoVolume *qqq = geom->MakeTubs("qqq", Al, qqqR1, qqqR2, 0.5, 5, 85);
    qqq->SetLineColor(7);
    for( int i = 0; i < 4; i++){
      worldBox->AddNode(qqq, i+1, new TGeoCombiTrans( 0, 
                                                      0, 
                                                      qqqZPos, 
                                                      new TGeoRotation("rot1", 360/4 * (i), 0., 0.)));
    }
  }

}

//!============================================== Drawing functions
inline void ANASEN::DrawAnasen(int anodeID1, int anodeID2, int cathodeID1, int cathodeID2, int sx3ID, bool DrawQQQ ){

  Construct3DModel(anodeID1, anodeID2, cathodeID1, cathodeID2, sx3ID, DrawQQQ);

  geom->CloseGeometry();
  geom->SetVisLevel(4);
  worldBox->Draw("ogle");

}

inline  void ANASEN::DrawTrack(TVector3 pos, TVector3 direction, bool drawEstimatedTrack){

  pw.FindWireID(pos, direction);
  sx3.FindSX3Pos(pos, direction);

  std::pair<short, short> wireID = pw.GetNearestID();
  
  Construct3DModel(wireID.first, wireID.first, wireID.second, wireID.second, -1, false);

  double theta = direction.Theta() * TMath::RadToDeg();
  double phi = direction.Phi()  * TMath::RadToDeg();
  // printf("Theta, Phi = %.2f %.2f \n", theta, phi);
  // pos.Print();
  TGeoVolume * Track = geom->MakeTube("track", 0, 0, 0.1, 150.);
  Track->SetLineColor(kRed);
  worldBox->AddNode(Track, 1, new TGeoCombiTrans( pos.X(), pos.Y(), pos.Z(), new TGeoRotation("rotA", phi +  90, theta, 0.)));

  TGeoVolume * startPos = geom->MakeSphere("startPos", 0, 0, 3);
  startPos->SetLineColor(kBlack);
  worldBox->AddNode(startPos, 3, new TGeoCombiTrans( pos.X(), pos.Y(), pos.Z(), new TGeoRotation("rotA", 0, 0, 0.)));

  if( sx3.GetID() >= 0 ){
    TGeoVolume * hit = geom->MakeSphere("hitpos", 0, 0, 3);
    hit->SetLineColor(kRed);

    TVector3 hitPos = sx3.GetHitPos();

    worldBox->AddNode(hit, 2, new TGeoCombiTrans( hitPos.X(), hitPos.Y(), hitPos.Z(), new TGeoRotation("rotA", 0, 0, 0.)));

    if( drawEstimatedTrack ){
      pw.CalTrack(hitPos, wireID.first, wireID.second, true);

      double thetaDeduce = pw.GetTrackTheta() * TMath::RadToDeg();
      double phiDeduce = pw.GetTrackPhi()  * TMath::RadToDeg();

      TGeoVolume * trackDeduce = geom->MakeTube("trackDeduce", 0, 0, 0.1, 100.);
      trackDeduce->SetLineColor(kOrange);
      worldBox->AddNode(trackDeduce, 1, new TGeoCombiTrans( hitPos.X(), hitPos.Y(), hitPos.Z(), new TGeoRotation("rotA", phiDeduce +  90, thetaDeduce, 0.)));

    }

  }
  
  geom->CloseGeometry();
  geom->SetVisLevel(4);
  worldBox->Draw("ogle");

}

inline  void ANASEN::DrawDeducedTrack(TVector3 sx3Pos, int anodeID, int cathodeID){

  pw.CalTrack(sx3Pos, anodeID, cathodeID);

  Construct3DModel(anodeID, anodeID, cathodeID, cathodeID, -1, false);

  double theta = pw.GetTrackTheta() * TMath::RadToDeg();
  double phi = pw.GetTrackPhi()  * TMath::RadToDeg();

  TGeoVolume * Track = geom->MakeTube("axisX", 0, 0, 0.1, 100.);
  Track->SetLineColor(kRed);
  worldBox->AddNode(Track, 1, new TGeoCombiTrans( sx3Pos.X(), sx3Pos.Y(), sx3Pos.Z(), new TGeoRotation("rotA", phi +  90, theta, 0.)));

  TGeoVolume * hit = geom->MakeSphere("hitpos", 0, 0, 3);
  hit->SetLineColor(kRed);
  worldBox->AddNode(hit, 2, new TGeoCombiTrans( sx3Pos.X(), sx3Pos.Y(), sx3Pos.Z(), new TGeoRotation("rotA", 0, 0, 0.)));

  geom->CloseGeometry();
  geom->SetVisLevel(4);
  worldBox->Draw("ogle");

}


#endif 