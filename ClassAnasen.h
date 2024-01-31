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

class ANASEN{
public:
  ANASEN();
  ~ANASEN();

  TVector3 CalSX3Pos(unsigned short ID, unsigned short chUp, unsigned short chDown, unsigned short chBack, float eUp = 0, float eDown = 0 );

  void CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID);

  TVector3 GetTrackPos() const {return trackPos;}
  TVector3 GetTrackVec() const {return trackVec;}

  double GetTrackTheta() const {return trackVec.Theta();}
  double GetTrackPhi()   const {return trackVec.Phi();}

  void DrawAnasen(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1, bool DrawQQQ = false );

private:

  const int nWire = 24;
  const int wireShift = 3;
  const int zLen = 380; //mm
  const int radiusA = 37;
  const int radiusC = 43;

  std::vector<TVector3> P1; // the anode wire position vector in space 
  std::vector<TVector3> P2; // the anode wire position vector in space 
  std::vector<TVector3> Q1; // the cathode wire position vector in space 
  std::vector<TVector3> Q2; // the cathode wire position vector in space 

  std::vector<std::pair<TVector3,TVector3>> S1; // coners of the SX3  0-11, z = mid point
  std::vector<std::pair<TVector3,TVector3>> S2; // coners of the SX3 12-23, z = mid point

  void CalGeometry();

  TVector3 trackPos;
  TVector3 trackVec;

  double trackPosErrorZ;  // mm
  TVector3 tracePosErrorXY; // the mag is the size of the error

  TVector3 trackVecErrorA; // error vector prependicular to the Anode-Pos plan
  TVector3 trackVecErrorC; // error vector prependicular to the Cathode-Pos plan

  const int nSX3 = 12;
  const int sx3Radius = 88;
  const int sx3Width = 40;
  const int sx3Length = 75;
  const int sx3Gap = 46;

  const int qqqR1 = 10;
  const int qqqR2 = 50;
  const int qqqZPos = sx3Gap/2 + sx3Length + 30;

  TGeoManager *geom;
  TGeoVolume *worldBox;

  void Construct3DModel(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1, bool DrawQQQ = true);

};

//==============================================
inline ANASEN::ANASEN(){

  CalGeometry();

  geom = nullptr;
  worldBox = nullptr;

}

inline ANASEN::~ANASEN(){

  delete geom;

}

inline void ANASEN::CalGeometry(){

  TVector3 p1; // anode
  TVector3 p2;
  TVector3 q1; // cathode
  TVector3 q2;
  for(int i = 0; i < nWire; i++ ){

    // Anode rotate right-hand
    p1.SetXYZ( radiusA * TMath::Cos( TMath::TwoPi() * i / nWire ),
               radiusA * TMath::Sin( TMath::TwoPi() * i / nWire ),
               zLen/2);
    p2.SetXYZ( radiusA * TMath::Cos( TMath::TwoPi() * (i + wireShift) / nWire ),
               radiusA * TMath::Sin( TMath::TwoPi() * (i + wireShift) / nWire ),
               -zLen/2);
    P1.push_back(p1);
    P2.push_back(p2);

    // Cathod rotate left-hand
    p1.SetXYZ( radiusC * TMath::Cos( TMath::TwoPi() * i / nWire ),
               radiusC * TMath::Sin( TMath::TwoPi() * i / nWire ),
               zLen/2);
    p2.SetXYZ( radiusC * TMath::Cos( TMath::TwoPi() * (i - wireShift) / nWire ),
               radiusC * TMath::Sin( TMath::TwoPi() * (i - wireShift) / nWire ),
               -zLen/2);
    Q1.push_back(p1);
    Q2.push_back(p2);
  }

  TVector3 sa, sb;
  for(int i = 0; i < nSX3; i++){
    sa.SetXYZ( sx3Radius,  sx3Width/2, sx3Gap/2 + sx3Length/2 );
    sb.SetXYZ( sx3Radius, -sx3Width/2, sx3Gap/2 + sx3Length/2 );

    sa.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    sb.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    S1.push_back(std::pair(sa,sb));

    sa.SetXYZ( sx3Radius,  sx3Width/2, -sx3Gap/2 - sx3Length/2 );
    sb.SetXYZ( sx3Radius, -sx3Width/2, -sx3Gap/2 - sx3Length/2 );

    sa.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    sb.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    S2.push_back(std::pair(sa,sb));
  }

}

inline void ANASEN::Construct3DModel(int anodeID1, int anodeID2, int cathodeID1, int cathodeID2, bool DrawQQQ ){

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
  double dAngle = wireShift * TMath::TwoPi() / nWire;
  double radiusAnew = radiusA * TMath::Cos( dAngle / 2.);
  double wireALength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusA * TMath::Sin(dAngle/2),2) );
  double wireATheta = TMath::ATan2( 2* radiusA * TMath::Sin( dAngle / 2.), zLen);

  // printf("    dAngle : %f\n", dAngle);
  // printf(" newRadius : %f\n", radiusAnew);
  // printf("wireLength : %f\n", wireALength);
  // printf("wire Theta : %f\n", wireATheta);

  TGeoVolume *pcA = geom->MakeTube("tub1", Al, 0, 0.01, wireALength/2);
  pcA->SetLineColor(4);  

  for( int i = 0; i < nWire; i++){
    if( anodeID2 >= 0 &&  (i < anodeID1 || i > anodeID2) ) continue;   
    worldBox->AddNode(pcA, i+1, new TGeoCombiTrans( radiusAnew * TMath::Cos( TMath::TwoPi() / nWire *i + dAngle / 2), 
                                                    radiusAnew * TMath::Sin( TMath::TwoPi() / nWire *i + dAngle / 2), 
                                                    0, 
                                                    new TGeoRotation("rot1", 360/ nWire * (i + wireShift/2.), wireATheta * 180/ TMath::Pi(), 0.)));
  }

  double radiusCnew = radiusC * TMath::Cos( dAngle / 2.);
  double wireCLength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusC * TMath::Sin(dAngle/2),2) );
  double wireCTheta = TMath::ATan2( 2* radiusC * TMath::Sin( dAngle / 2.), zLen);

  TGeoVolume *pcC = geom->MakeTube("tub2", Al, 0, 0.01, wireCLength/2);
  pcC->SetLineColor(6);
  for( int i = 0; i < nWire; i++){
    if( cathodeID2 >= 0 && (i < cathodeID1 || i > cathodeID2) ) continue;   
    worldBox->AddNode(pcC, i+1, new TGeoCombiTrans( radiusCnew * TMath::Cos( TMath::TwoPi() / nWire *i - dAngle/2), 
                                                    radiusCnew * TMath::Sin( TMath::TwoPi() / nWire *i - dAngle/2), 
                                                    0, 
                                                    new TGeoRotation("rot1", 360/ nWire * (i - wireShift/2.), -wireCTheta * 180/ TMath::Pi(), 0.)));
  }

  TGeoVolume * sx3 = geom->MakeBox("box", Al, 0.1, sx3Width/2, sx3Length/2);
  sx3->SetLineColor(kGreen+3);
  for( int i = 0; i < nSX3; i++){
    worldBox->AddNode(sx3, 2*i+1., new TGeoCombiTrans( sx3Radius * TMath::Cos( TMath::TwoPi() / nSX3 * (i + 0.5)), 
                                                       sx3Radius * TMath::Sin( TMath::TwoPi() / nSX3 * (i + 0.5)), 
                                                     sx3Length/2+sx3Gap/2, 
                                                     new TGeoRotation("rot1", 360/nSX3 * (i + 0.5), 0., 0.)));

    worldBox->AddNode(sx3, 2*i+2., new TGeoCombiTrans( sx3Radius * TMath::Cos( TMath::TwoPi() / nSX3 * (i + 0.5)), 
                                                       sx3Radius * TMath::Sin( TMath::TwoPi() / nSX3 * (i + 0.5)), 
                                                       -sx3Length/2-sx3Gap/2, 
                                                       new TGeoRotation("rot1", 360/nSX3 * (i + 0.5), 0., 0.)));
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


inline void ANASEN::DrawAnasen(int anodeID1, int anodeID2, int cathodeID1, int cathodeID2, bool DrawQQQ ){

  
  Construct3DModel(anodeID1, anodeID2, cathodeID1, cathodeID2, DrawQQQ);

  geom->CloseGeometry();
  geom->SetVisLevel(4);
  worldBox->Draw("ogle");

}

inline void ANASEN::CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID){

  trackPos = sx3Pos;

  TVector3 n1 = (P2[anodeID] - P1[anodeID]).Cross((sx3Pos - P1[anodeID]));
  TVector3 n2 = (Q2[anodeID] - Q1[anodeID]).Cross((sx3Pos - Q1[anodeID]));

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n1.Cross(n2)).Unit();

}

inline TVector3 ANASEN::CalSX3Pos(unsigned short ID, unsigned short chUp, unsigned short chDown, unsigned short chBack, float eUp, float eDown){

  TVector3 haha;

  if( (chUp - chDown) != 1 || (chDown % 2) != 0) return haha;

  int reducedID = ID % nSX3;

  TVector3 sa, sb;

  if( ID < nSX3 ){ //down

    sa = S1[reducedID].first;
    sb = S1[reducedID].second;

  }else{

    sa = S2[reducedID].first;
    sb = S2[reducedID].second;

  }

  haha.SetX( (sb.X() - sa.X()) * chUp/8 + sa.X());
  haha.SetY( (sb.Y() - sa.Y()) * chUp/8 + sa.Y());

  if( eUp == 0 || eDown == 0 ){
    haha.SetZ( sa.Z() + (2*(chBack - 7)-1) * sx3Length / 8 );
  }else{
    double frac = (eUp - eDown)/(eUp + eDown); // from +1 (downstream) to -1 (upstream)
    double zPos = sa.Z() +  sx3Length * frac/2;
    haha.SetZ( zPos );
  }

  return haha;

}

#endif 