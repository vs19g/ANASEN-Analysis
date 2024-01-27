#ifndef ClassAnasen_h
#define ClassAnasen_h

#include <cstdio>
#include <TMath.h>
#include <TVector3.h>

class Anasen{
public:
  Anasen();
  ~Anasen() {}

  void CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID);

  TVector3 GetTrackPos() const {return trackPos;}
  TVector3 GetTrackVec() const {return trackVec;}

  double GetTrackTheta() const {return trackVec.Theta();}
  double GetTrackPhi()   const {return trackVec.Phi();}

  // void DrawTrack();

private:

  const int nWire = 24;
  const int wireShift = 3; // how twisted is the wire.
  const double AnodeRadius = 38 ; // mm
  const double CathodeRadius = 38 ; // mm
  const double PCLength = 100; //mm

  std::vector<TVector3> P1; // the anode wire position vector in space 
  std::vector<TVector3> P2; // the anode wire position vector in space 
  std::vector<TVector3> Q1; // the cathode wire position vector in space 
  std::vector<TVector3> Q2; // the cathode wire position vector in space 

  void CalWireDirection();

  TVector3 trackPos;
  TVector3 trackVec;

  double trackPosErrorZ;  // mm
  TVector3 tracePosErrorXY; // the mag is the size of the error

  TVector3 trackVecErrorA; // error vector prependicular to the Anode-Pos plan
  TVector3 trackVecErrorC; // error vector prependicular to the Cathode-Pos plan

};

inline Anasen::Anasen(){

  CalWireDirection();

}

inline void Anasen::CalWireDirection(){

  TVector3 p1; // anode
  TVector3 p2;
  TVector3 q1; // cathode
  TVector3 q2;
  for(int i = 0; i < nWire; i++ ){

    // Anode rotate right-hand

    p1.SetXYZ( AnodeRadius * TMath::Cos( TMath::TwoPi() * i / nWire ),
               AnodeRadius * TMath::Sin( TMath::TwoPi() * i / nWire ),
               -PCLength/2);
    p2.SetXYZ( AnodeRadius * TMath::Cos( TMath::TwoPi() * (i + wireShift) / nWire ),
               AnodeRadius * TMath::Sin( TMath::TwoPi() * (i + wireShift) / nWire ),
               PCLength/2);
    P1.push_back(p1);
    P2.push_back(p2);

    // Cathod rotate left-hand

    p1.SetXYZ( CathodeRadius * TMath::Cos( TMath::TwoPi() * i / nWire ),
               CathodeRadius * TMath::Sin( TMath::TwoPi() * i / nWire ),
               -PCLength/2);
    p2.SetXYZ( CathodeRadius * TMath::Cos( TMath::TwoPi() * (i - wireShift) / nWire ),
               CathodeRadius * TMath::Sin( TMath::TwoPi() * (i - wireShift) / nWire ),
               PCLength/2);
    Q1.push_back(p1);
    Q2.push_back(p2);
  }

}

inline void Anasen::CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID){

  trackPos = sx3Pos;

  TVector3 n1 = (P2[anodeID] - P1[anodeID]).Cross((sx3Pos - P1[anodeID]));
  TVector3 n2 = (Q2[anodeID] - Q1[anodeID]).Cross((sx3Pos - Q1[anodeID]));

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n1.Cross(n2)).Unit();



}

#endif 