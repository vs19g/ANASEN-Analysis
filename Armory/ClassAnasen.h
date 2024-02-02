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

struct SX3{
  short id = -1; // -1 when no hit
  short chUp;
  short chDown;
  short chBack;

  double zFrac; // from +1 (downstream) to -1 (upstream)

  double eUp;
  double eDown;
  double eBack;

  TVector3 hitPos;

  void CalZFrac(){
    zFrac = (eUp - eDown)/(eUp + eDown); 
  }

  void Print(){ 
    if( id == -1 ){
      printf("Did not hit any SX3.\n"); 
    }else{
      printf("ID: %d, U,D,B: %d %d %d| zFrac : %.2f\n", id, chUp, chDown, chBack, zFrac); 
      printf("Hit Pos: %.2f, %.2f, %.2f\n", hitPos.X(), hitPos.Y(), hitPos.Z()); 
    }
  }

};

class ANASEN{
public:
  ANASEN();
  ~ANASEN();

  void CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose = false);
  TVector3 CalSX3Pos(unsigned short ID, unsigned short chUp, unsigned short chDown, unsigned short chBack, float eUp = 0, float eDown = 0 );

  TVector3 GetTrackPos() const {return trackPos;}
  TVector3 GetTrackVec() const {return trackVec;}
  double GetTrackTheta() const {return trackVec.Theta();}
  double GetTrackPhi()   const {return trackVec.Phi();}

  void DrawAnasen(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1, int sx3ID = -1, bool DrawQQQ = false );
  void DrawDeducedTrack(TVector3 sx3Pos, int anodeID, int cathodeID);

  //Simulation
  SX3 FindSX3Pos(TVector3 pos, TVector3 direction, bool verbose = false);
  std::pair<int, int> FindWireID(TVector3 pos, TVector3 direction, bool verbose = false);
  void DrawTrack(TVector3 pos, TVector3 direction, bool drawEstimatedTrack = false);

  std::pair<TVector3,TVector3> GetAnode(unsigned short id) const{return An[id];};
  std::pair<TVector3,TVector3> GetCathode(unsigned short id) const{return Ca[id];};

private:

  const int nWire = 24;
  const int wireShift = 3;
  const float zLen = 380; //mm
  const float radiusA = 37;
  const float radiusC = 43;

  std::vector<std::pair<TVector3,TVector3>> An; // the anode wire position vector in space 
  std::vector<std::pair<TVector3,TVector3>> Ca; // the cathode wire position vector in space 

  std::vector<std::pair<TVector3,TVector3>> SDn; // coners of the SX3  0-11, z = mid point
  std::vector<std::pair<TVector3,TVector3>> SUp; // coners of the SX3 12-23, z = mid point
  std::vector<TVector3> SNorml; // normal of the SX3 (outward)

  void CalGeometry();

  TVector3 trackPos;
  TVector3 trackVec;

  double trackPosErrorZ;  // mm
  TVector3 tracePosErrorXY; // the mag is the size of the error

  TVector3 trackVecErrorA; // error vector prependicular to the Anode-Pos plan
  TVector3 trackVecErrorC; // error vector prependicular to the Cathode-Pos plan

  const int nSX3 = 12;
  const float sx3Radius = 88;
  const float sx3Width = 40;
  const float sx3Length = 75;
  const float sx3Gap = 46;

  const float qqqR1 = 50;
  const float qqqR2 = 100;
  const float qqqZPos = sx3Gap/2 + sx3Length + 30;

  // int geomID;
  TGeoManager *geom;
  TGeoVolume *worldBox;

  void Construct3DModel(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1, int sx3ID = -1, bool DrawQQQ = true);

  double Distance(TVector3 a1, TVector3 a2, TVector3 b1, TVector3 b2); 
  std::pair<double, double> Intersect(TVector3 a1, TVector3 a2, TVector3 b1, TVector3 b2, bool verbose = false); 

};

//!==============================================
inline ANASEN::ANASEN(){

  CalGeometry();

  // geomID = 0;
  geom = nullptr;
  worldBox = nullptr;

}

inline ANASEN::~ANASEN(){

  delete geom;

}
//!==============================================
inline void ANASEN::CalGeometry(){

  std::pair<TVector3, TVector3> p1; // anode
  std::pair<TVector3, TVector3> q1; // cathode

  //anode and cathode start at pos-Y axis and count in left-Hand

  //anode wire shift is right-hand.
  //cathode wire shift is left-hand.

  for(int i = 0; i < nWire; i++ ){
    // Anode rotate right-hand
    p1.first.SetXYZ( radiusA * TMath::Cos( TMath::TwoPi() / nWire * (-i)  + TMath::PiOver2()),
                     radiusA * TMath::Sin( TMath::TwoPi() / nWire * (-i)  + TMath::PiOver2()),
                     zLen/2);
    p1.second.SetXYZ( radiusA * TMath::Cos( TMath::TwoPi() / nWire * (-i + wireShift) + TMath::PiOver2()),
                      radiusA * TMath::Sin( TMath::TwoPi() / nWire * (-i + wireShift) + TMath::PiOver2()),
                     -zLen/2);
    An.push_back(p1);

    // Cathod rotate left-hand
    q1.first.SetXYZ( radiusC * TMath::Cos( TMath::TwoPi() / nWire * (-i)  + TMath::PiOver2()),
                     radiusC * TMath::Sin( TMath::TwoPi() / nWire * (-i)  + TMath::PiOver2()),
                     zLen/2);
    q1.second.SetXYZ( radiusC * TMath::Cos( TMath::TwoPi() / nWire * (- i - wireShift) + TMath::PiOver2()),
                      radiusC * TMath::Sin( TMath::TwoPi() / nWire * (- i - wireShift) + TMath::PiOver2()),
                      -zLen/2);
    Ca.push_back(q1);
  }

  // SX3 is couned in left-hand, started at neg-Y axis

  TVector3 sa, sb, sc, sn;
  for(int i = 0; i < nSX3; i++){
    sa.SetXYZ( sx3Radius, -sx3Width/2, sx3Gap/2 + sx3Length/2 );
    sb.SetXYZ( sx3Radius,  sx3Width/2, sx3Gap/2 + sx3Length/2 );

    double rot = TMath::TwoPi() / nSX3 * (-i - 0.5) - TMath::PiOver2();

    sa.RotateZ( rot );
    sb.RotateZ( rot );
    SDn.push_back(std::pair(sa,sb));


    sc.SetXYZ( sx3Radius, -sx3Width/2, sx3Gap/2 );
    sc.RotateZ( rot );

    sn = ((sc-sa).Cross(sb-sa)).Unit();
    SNorml.push_back(sn);

    sa.SetXYZ( sx3Radius, -sx3Width/2, -sx3Gap/2 - sx3Length/2 );
    sb.SetXYZ( sx3Radius,  sx3Width/2, -sx3Gap/2 - sx3Length/2 );

    sa.RotateZ( rot );
    sb.RotateZ( rot );
    SUp.push_back(std::pair(sa,sb));
  }

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
  double dAngle = wireShift * TMath::TwoPi() / nWire;
  double wireALength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusA * TMath::Sin(dAngle/2),2) );
  // double radiusAnew = radiusA * TMath::Cos( dAngle / 2.);
  // double wireATheta = TMath::ATan2( 2* radiusA * TMath::Sin( dAngle / 2.), zLen);

  TGeoVolume *pcA = geom->MakeTube("tub1", Al, 0, 0.01, wireALength/2);
  pcA->SetLineColor(4);  

  int startID = 0;
  int endID = nWire - 1;

  if( anodeID1 >= 0 && anodeID2 >= 0 ){
    startID = anodeID1;
    endID = anodeID2;
    if( anodeID1 > anodeID2 ) {
      endID = nWire + anodeID2;
    }
  }

  for( int i = startID; i <= endID; i++){
    TVector3 a = (An[i].first + An[i].second) * 0.5;
    double wireATheta = (An[i].first - An[i].second).Theta()* TMath::RadToDeg();
    double wireAPhi = (An[i].first - An[i].second).Phi() * TMath::RadToDeg() + 90;

    worldBox->AddNode(pcA, i+1, new TGeoCombiTrans( a.X(), 
                                                    a.Y(), 
                                                    a.Z(), 
                                                    new TGeoRotation("rot1", wireAPhi , wireATheta, 0.)));
  }

  double wireCLength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusC * TMath::Sin(dAngle/2),2) );
  // double radiusCnew = radiusC * TMath::Cos( dAngle / 2.);
  // double wireCTheta = TMath::ATan2( 2* radiusC * TMath::Sin( dAngle / 2.), zLen);

  TGeoVolume *pcC = geom->MakeTube("tub2", Al, 0, 0.01, wireCLength/2);
  pcC->SetLineColor(6);

  startID = 0;
  endID = nWire - 1;

  if( cathodeID1 >= 0 && cathodeID2 >= 0 ){
    startID = cathodeID1;
    endID = cathodeID2;
    if( cathodeID1 > cathodeID2 ) {
      endID = nWire + cathodeID2;
    }
  }

  for( int i = startID; i <= endID; i++){
    TVector3 a = (Ca[i].first + Ca[i].second) * 0.5;
    double wireATheta = (Ca[i].first - Ca[i].second).Theta()* TMath::RadToDeg();
    double wireAPhi = (Ca[i].first - Ca[i].second).Phi() * TMath::RadToDeg() + 90;

    worldBox->AddNode(pcC, i+1, new TGeoCombiTrans( a.X(), 
                                                    a.Y(), 
                                                    a.Z(), 
                                                    new TGeoRotation("rot1", wireAPhi , wireATheta, 0.)));
  }

  TGeoVolume * sx3 = geom->MakeBox("box", Al, 0.1, sx3Width/2, sx3Length/2);
  sx3->SetLineColor(kGreen+3);

  for( int i = 0; i < nSX3; i++){
    if( sx3ID != -1 && i != sx3ID ) continue;
    TVector3 aUp = (SUp[i].first + SUp[i].second)*0.5; // center of the SX3 upstream
    TVector3 aDn = (SDn[i].first + SDn[i].second)*0.5; // center of the SX3 Downstream
    double phi = (SUp[i].second - SUp[i].first).Phi() * TMath::RadToDeg() + 90;

    worldBox->AddNode(sx3, 2*i+1., new TGeoCombiTrans( aUp.X(), 
                                                       aUp.Y(), 
                                                       aUp.Z(), 
                                                       new TGeoRotation("rot1", phi, 0., 0.)));
    worldBox->AddNode(sx3, 2*i+1., new TGeoCombiTrans( aDn.X(), 
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

//!============================================== Aux Functions
inline double ANASEN::Distance(TVector3 a1, TVector3 a2, TVector3 b1, TVector3 b2){
  TVector3 na = a1 - a2;
  TVector3 nb = b1 - b2;
  TVector3 nd = (na.Cross(nb)).Unit();
  return TMath::Abs(nd.Dot(a1-b2));
} 


inline std::pair<double, double> ANASEN::Intersect(TVector3 An, TVector3 p2, TVector3 q1, TVector3 q2, bool verbose){

  //see https://nukephysik101.wordpress.com/2023/12/30/intersect-between-2-line-segments/
  //zero all z-component
  TVector3 a0 = An; a0.SetZ(0);
  TVector3 a1 = p2; a1.SetZ(0);

  TVector3 b0 = q1; b0.SetZ(0);
  TVector3 b1 = q2; b1.SetZ(0);

  double A = ((b0-b1).Cross(a0-a1)).Mag();

  double h = ((b0-a0).Cross(b1-a0)).Z()/ A;
  double k = ((a1-b0).Cross(a0-b0)).Z()/ A;

  if( verbose ) printf(" ----h, k : %f, %f\n", h, k);

  return std::pair(h,k);
} 

//!============================================== Given a position and a direction, find wireID and SX3 position
inline std::pair<int, int> ANASEN::FindWireID(TVector3 pos, TVector3 direction, bool verbose ){

  int anodeID = -1;
  int cathodeID = -1;
  double minAnodeDis = 999999;
  double minCathodeDis = 999999;

  double phi = direction.Phi();
  

  for( int i = 0; i < nWire; i++){

    double disA = 99999999;
    double disC = 99999999;

    double phiS = An[i].first.Phi()  - TMath::PiOver4();
    double phiL = An[i].second.Phi() + TMath::PiOver4();

    // printf("A%2d: %f %f | %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg(), phi * TMath::RadToDeg());

    if( phi > 0 && phiS > phiL ) {
      phiL = phiL + TMath::TwoPi();
      // printf("------ %f %f\n", phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    }
    if( phi < 0 && phiS > phiL ) {
      phiS = phiS - TMath::TwoPi();
      // printf("------ %f %f\n", phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    }

    if( phiS < phi && phi < phiL) {
      disA = Distance( pos, pos + direction, An[i].first, An[i].second);
      if( disA < minAnodeDis ){
        minAnodeDis = disA;
        anodeID = i;
      }
    }

    phiS = Ca[i].second.Phi()- TMath::PiOver4();
    phiL = Ca[i].first.Phi() + TMath::PiOver4();
    // printf("C%2d: %f %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    if( phi > 0 && phiS > phiL ) {
      phiL = phiL + TMath::TwoPi();
      // printf("------ %f %f\n", phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    }
    if( phi < 0 && phiS > phiL ) {
      phiS = phiS - TMath::TwoPi();
      // printf("------ %f %f\n", phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    }

    if(phiS < phi && phi < phiL) {
      disC = Distance( pos, pos + direction, Ca[i].first, Ca[i].second);

      if( disC < minCathodeDis ){
        minCathodeDis = disC;
        cathodeID = i;
      }
    }

    if(verbose) printf(" %2d | %8.2f, %8.2f\n", i, disA, disC);
  }

  if( verbose ) printf("AnodeID %d (%.2f), Cathode %d (%.2f) \n", anodeID, minAnodeDis, cathodeID, minCathodeDis);
  return std::pair(anodeID, cathodeID);
}

inline SX3 ANASEN::FindSX3Pos(TVector3 pos, TVector3 direction, bool verbose){

  SX3 haha;

  haha.id = -1;
  for( int i = 0 ; i < nSX3; i++){

    if(verbose) printf(" %d ", i);
    std::pair<double, double> frac = Intersect( pos, pos + direction, SDn[i].first, SDn[i].second, verbose);


    if( frac.second < 0 || frac.second > 1 ) continue;
    haha.hitPos = pos + frac.first * direction;

    double dis = haha.hitPos.Dot(SNorml[i]);

    if(verbose) {
      printf("reduced distance : %f\n", dis);
      printf(" %d*", (i+1)%nSX3); 
      Intersect( pos, pos + direction, SDn[(i+1)%nSX3].first, SDn[(i+1)%nSX3].second, verbose);
    }

    if( TMath::Abs(dis - sx3Radius) > 0.1 ) continue;

    haha.chDown = 2 * TMath::Floor(frac.second * 4);
    haha.chUp = haha.chDown + 1;

    double zPos = haha.hitPos.Z();
    if( (sx3Gap/2 < zPos && zPos < sx3Gap/2 + sx3Length ) || (-sx3Gap/2 - sx3Length < zPos && zPos < -sx3Gap/2 )  ){

      haha.id = zPos > 0 ? i : i + 12;

      haha.zFrac = zPos > 0 ?  (zPos - sx3Gap/2. - sx3Length/2.)/sx3Length : (zPos - ( - sx3Gap/2. - sx3Length/2.) )/sx3Length ;

      haha.chBack = TMath::Floor( (haha.zFrac + 0.5) * 4 ) + 8;

      if( verbose) haha.Print();

      return haha;
    }else{
      if( verbose ) printf(" zPos out of sensitive region\n");
    }
  }

  if( verbose) haha.Print();
  return haha;
}

//!============================================== Drawing functions
inline void ANASEN::DrawAnasen(int anodeID1, int anodeID2, int cathodeID1, int cathodeID2, int sx3ID, bool DrawQQQ ){

  Construct3DModel(anodeID1, anodeID2, cathodeID1, cathodeID2, sx3ID, DrawQQQ);

  geom->CloseGeometry();
  geom->SetVisLevel(4);
  worldBox->Draw("ogle");

}

inline  void ANASEN::DrawTrack(TVector3 pos, TVector3 direction, bool drawEstimatedTrack){

  std::pair<int, int> id = FindWireID(pos, direction);

  SX3 sx3 = FindSX3Pos(pos, direction);

  int a1 = id.first - 1; if( a1 < 0 ) a1 += nWire;
  int b1 = id.second - 1; if( b1 < 0 ) b1 += nWire;

  //Construct3DModel(a1, id.first+1, b1, id.second+1, false);
  Construct3DModel(id.first, id.first, id.second, id.second, -1, false);

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

  if( sx3.id >= 0 ){
    TGeoVolume * hit = geom->MakeSphere("hitpos", 0, 0, 3);
    hit->SetLineColor(kRed);
    worldBox->AddNode(hit, 2, new TGeoCombiTrans( sx3.hitPos.X(), sx3.hitPos.Y(), sx3.hitPos.Z(), new TGeoRotation("rotA", 0, 0, 0.)));

    if( drawEstimatedTrack ){
      CalTrack(sx3.hitPos, id.first, id.second, true);

      double thetaDeduce = trackVec.Theta() * TMath::RadToDeg();
      double phiDeduce = trackVec.Phi()  * TMath::RadToDeg();

      TGeoVolume * trackDeduce = geom->MakeTube("trackDeduce", 0, 0, 0.1, 100.);
      trackDeduce->SetLineColor(kOrange);
      worldBox->AddNode(trackDeduce, 1, new TGeoCombiTrans( sx3.hitPos.X(), sx3.hitPos.Y(), sx3.hitPos.Z(), new TGeoRotation("rotA", phiDeduce +  90, thetaDeduce, 0.)));

    }

  }
  
  geom->CloseGeometry();
  geom->SetVisLevel(4);
  worldBox->Draw("ogle");

}

inline  void ANASEN::DrawDeducedTrack(TVector3 sx3Pos, int anodeID, int cathodeID){

  CalTrack(sx3Pos, anodeID, cathodeID);

  Construct3DModel(anodeID, anodeID, cathodeID, cathodeID, -1, false);

  double theta = trackVec.Theta() * TMath::RadToDeg();
  double phi = trackVec.Phi()  * TMath::RadToDeg();

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

//!============================================== Duduce trace from experiment
inline void ANASEN::CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose){

  trackPos = sx3Pos;

  TVector3 n1 = (An[anodeID].first - An[anodeID].second).Cross((sx3Pos - An[anodeID].second)).Unit();
  TVector3 n2 = (Ca[cathodeID].first - Ca[cathodeID].second).Cross((sx3Pos - Ca[cathodeID].second)).Unit();

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n2.Cross(n1)).Unit();

  if( verbose ) printf("Theta, Phi = %f, %f \n", trackVec.Theta() *TMath::RadToDeg(), trackVec.Phi()*TMath::RadToDeg()); 

}

inline TVector3 ANASEN::CalSX3Pos(unsigned short ID, unsigned short chUp, unsigned short chDown, unsigned short chBack, float eUp, float eDown){

  TVector3 haha;

  if( (chUp - chDown) != 1 || (chDown % 2) != 0) return haha;

  int reducedID = ID % nSX3;

  TVector3 sa, sb;

  if( ID < nSX3 ){ //down

    sa = SDn[reducedID].second;
    sb = SDn[reducedID].first;

  }else{

    sa = SUp[reducedID].second;
    sb = SUp[reducedID].first;

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