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

  void DrawAnasen(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1, bool DrawQQQ = false );
  void DrawDeducedTrack(TVector3 sx3Pos, int anodeID, int cathodeID);

  //Simulation
  SX3 FindSX3Pos(TVector3 pos, TVector3 direction, bool verbose = false);
  std::pair<int, int> FindWireID(TVector3 pos, TVector3 direction, bool verbose = false);
  void DrawTrack(TVector3 pos, TVector3 direction, bool drawEstimatedTrack = false);

  std::pair<TVector3,TVector3> GetAnode(unsigned short id) const{return P1[id];};
  std::pair<TVector3,TVector3> GetCathode(unsigned short id) const{return Q1[id];};

private:

  const int nWire = 24;
  const int wireShift = 3;
  const float zLen = 380; //mm
  const float radiusA = 37;
  const float radiusC = 43;

  std::vector<std::pair<TVector3,TVector3>> P1; // the anode wire position vector in space 
  std::vector<std::pair<TVector3,TVector3>> Q1; // the cathode wire position vector in space 

  std::vector<std::pair<TVector3,TVector3>> SD; // coners of the SX3  0-11, z = mid point
  std::vector<std::pair<TVector3,TVector3>> SU; // coners of the SX3 12-23, z = mid point
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

  void Construct3DModel(int anodeID1 = -1, int anodeID2 = -1, int cathodeID1 = -1, int cathodeID2 = -1, bool DrawQQQ = true);

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

  for(int i = 0; i < nWire; i++ ){
    // Anode rotate right-hand
    p1.first.SetXYZ( radiusA * TMath::Cos( TMath::TwoPi() / nWire * i ),
                     radiusA * TMath::Sin( TMath::TwoPi() / nWire * i ),
                     zLen/2);
    p1.second.SetXYZ( radiusA * TMath::Cos( TMath::TwoPi() / nWire * (i + wireShift)),
                      radiusA * TMath::Sin( TMath::TwoPi() / nWire * (i + wireShift)),
                     -zLen/2);
    P1.push_back(p1);

    // P1.back().first.Print();
    // P1.back().second.Print();

    // Cathod rotate left-hand
    q1.first.SetXYZ( radiusC * TMath::Cos( TMath::TwoPi() * i / nWire ),
               radiusC * TMath::Sin( TMath::TwoPi() * i / nWire ),
               zLen/2);
    q1.second.SetXYZ( radiusC * TMath::Cos( TMath::TwoPi() * (i - wireShift) / nWire ),
               radiusC * TMath::Sin( TMath::TwoPi() * (i - wireShift) / nWire ),
               -zLen/2);
    Q1.push_back(q1);

    // Q1.back().first.Print();
    // Q1.back().second.Print();

  }

  TVector3 sa, sb, sc, sn;
  for(int i = 0; i < nSX3; i++){
    sa.SetXYZ( sx3Radius, -sx3Width/2, sx3Gap/2 + sx3Length/2 );
    sb.SetXYZ( sx3Radius,  sx3Width/2, sx3Gap/2 + sx3Length/2 );

    sa.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    sb.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    SD.push_back(std::pair(sa,sb));


    sc.SetXYZ( sx3Radius, -sx3Width/2, sx3Gap/2 );
    sc.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );

    sn = ((sc-sa).Cross(sb-sa)).Unit();
    SNorml.push_back(sn);

    sa.SetXYZ( sx3Radius, -sx3Width/2, -sx3Gap/2 - sx3Length/2 );
    sb.SetXYZ( sx3Radius,  sx3Width/2, -sx3Gap/2 - sx3Length/2 );

    sa.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    sb.RotateZ( TMath::TwoPi() / nSX3 * (i + 0.5) );
    SU.push_back(std::pair(sa,sb));
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

//!============================================== Aux Functions
inline double ANASEN::Distance(TVector3 a1, TVector3 a2, TVector3 b1, TVector3 b2){
  TVector3 na = a1 - a2;
  TVector3 nb = b1 - b2;
  TVector3 nd = (na.Cross(nb)).Unit();
  return TMath::Abs(nd.Dot(a1-b2));
} 


inline std::pair<double, double> ANASEN::Intersect(TVector3 p1, TVector3 p2, TVector3 q1, TVector3 q2, bool verbose){

  //see https://nukephysik101.wordpress.com/2023/12/30/intersect-between-2-line-segments/
  //zero all z-component
  TVector3 a0 = p1; a0.SetZ(0);
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

    double phiS = P1[i].first.Phi()  - TMath::PiOver4();
    double phiL = P1[i].second.Phi() + TMath::PiOver4();

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
      disA = Distance( pos, pos + direction, P1[i].first, P1[i].second);
      if( disA < minAnodeDis ){
        minAnodeDis = disA;
        anodeID = i;
      }
    }

    phiS = Q1[i].second.Phi()- TMath::PiOver4();
    phiL = Q1[i].first.Phi() + TMath::PiOver4();
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
      disC = Distance( pos, pos + direction, Q1[i].first, Q1[i].second);

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
    std::pair<double, double> frac = Intersect( pos, pos + direction, SD[i].first, SD[i].second, verbose);


    if( frac.second < 0 || frac.second > 1 ) continue;
    haha.hitPos = pos + frac.first * direction;

    double dis = haha.hitPos.Dot(SNorml[i]);

    if(verbose) {
      printf("reduced distance : %f\n", dis);
      printf(" %d*", (i+1)%nSX3); 
      Intersect( pos, pos + direction, SD[(i+1)%nSX3].first, SD[(i+1)%nSX3].second, verbose);
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
inline void ANASEN::DrawAnasen(int anodeID1, int anodeID2, int cathodeID1, int cathodeID2, bool DrawQQQ ){

  Construct3DModel(anodeID1, anodeID2, cathodeID1, cathodeID2, DrawQQQ);

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
  Construct3DModel(id.first, id.first, id.second, id.second, false);

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

  Construct3DModel(anodeID, anodeID, cathodeID, cathodeID, false);

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

  TVector3 n1 = (P1[anodeID].first - P1[anodeID].second).Cross((sx3Pos - P1[anodeID].second)).Unit();
  TVector3 n2 = (Q1[cathodeID].first - Q1[cathodeID].second).Cross((sx3Pos - Q1[cathodeID].second)).Unit();

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

    sa = SD[reducedID].first;
    sb = SD[reducedID].second;

  }else{

    sa = SU[reducedID].first;
    sb = SU[reducedID].second;

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