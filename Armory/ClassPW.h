#ifndef ClassPW_h
#define ClassPW_h

#include <cstdio>
#include <TMath.h>
#include <TVector3.h>

struct PWHitInfo{
  std::pair<short, short>   nearestWire; // anode, cathode
  std::pair<double, double> nearestDist; // anode, cathode

  std::pair<short, short>   nextNearestWire; // anode, cathode
  std::pair<double, double> nextNearestDist; // anode, cathode

  void Clear(){
    nearestWire.first = -1;
    nearestWire.second = -1;
    nearestDist.first = 999999999;
    nearestDist.second = 999999999;
    nextNearestWire.first = -1;
    nextNearestWire.second = -1;
    nextNearestDist.first = 999999999;
    nextNearestDist.second = 999999999;
  }
};

//!########################################################
class PW{ // proportional wire
public: 
  PW(){ ClearHitInfo(); };
  ~PW(){};

  PWHitInfo GetHitInfo() const {return hitInfo;}
  std::pair<short, short>   GetNearestID()          const {return hitInfo.nearestWire;}
  std::pair<double, double> GetNearestDistance()    const {return hitInfo.nearestDist;}
  std::pair<short, short>   Get2ndNearestID()       const {return hitInfo.nextNearestWire;}
  std::pair<double, double> Get2ndNearestDistance() const {return hitInfo.nextNearestDist;}

  TVector3 GetTrackPos() const {return trackPos;}
  TVector3 GetTrackVec() const {return trackVec;}
  double GetTrackTheta() const {return trackVec.Theta();}
  double GetTrackPhi()   const {return trackVec.Phi();}

  int GetNumWire() const {return nWire;}
  double GetDeltaAngle() const {return dAngle;}
  double GetAnodeLength() const {return anodeLength;}
  double GetCathodeLength() const {return cathodeLength;}
  TVector3 GetAnodeDn(short id) const {return An[id].first;}
  TVector3 GetAnodeUp(short id) const {return An[id].second;}
  TVector3 GetCathodeDn(short id) const {return Ca[id].first;}
  TVector3 GetCathodeUp(short id) const {return Ca[id].second;}

  TVector3 GetAnodneMid(short id)  const {return (An[id].first + An[id].second) * 0.5; }
  double   GetAnodeTheta(short id) const {return (An[id].first - An[id].second).Theta();}
  double   GetAnodePhi(short id)   const {return (An[id].first - An[id].second).Phi();}

  TVector3 GetCathodneMid(short id)  const {return (Ca[id].first + Ca[id].second) * 0.5; }
  double   GetCathodeTheta(short id) const {return (Ca[id].first - Ca[id].second).Theta();}
  double   GetCathodePhi(short id)   const {return (Ca[id].first - Ca[id].second).Phi();}

  void ClearHitInfo();
  void ConstructGeo();
  void FindWireID(TVector3 pos, TVector3 direction, bool verbose = false);
  void CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose = false);
  void CalTrack2(TVector3 sx3Pos, PWHitInfo hitInfo, bool verbose = false);

  void Print(){
    printf("     The nearest | Anode: %2d(%5.2f) Cathode: %2d(%5.2f)\n", hitInfo.nearestWire.first, 
                                                                         hitInfo.nearestDist.first, 
                                                                         hitInfo.nearestWire.second, 
                                                                         hitInfo.nearestDist.second);

    printf(" The 2nd nearest | Anode: %2d(%5.2f) Cathode: %2d(%5.2f)\n", hitInfo.nextNearestWire.first, 
                                                                         hitInfo.nextNearestDist.first, 
                                                                         hitInfo.nextNearestWire.second, 
                                                                         hitInfo.nextNearestDist.second);
  }

private:

  PWHitInfo hitInfo;

  TVector3 trackPos;
  TVector3 trackVec;

  const int nWire = 24;
  const int wireShift = 3;
  const float zLen = 380; //mm
  const float radiusA = 37;
  const float radiusC = 43;

  double dAngle;
  double anodeLength;
  double cathodeLength;

  std::vector<std::pair<TVector3,TVector3>> An; // the anode wire position vector in space 
  std::vector<std::pair<TVector3,TVector3>> Ca; // the cathode wire position vector in space 

  double Distance(TVector3 a1, TVector3 a2, TVector3 b1, TVector3 b2){
    TVector3 na = a1 - a2;
    TVector3 nb = b1 - b2;
    TVector3 nd = (na.Cross(nb)).Unit();
    return TMath::Abs(nd.Dot(a1-b2));
  } 

};

inline void PW::ClearHitInfo(){
  hitInfo.Clear();
}

inline void PW::ConstructGeo(){

  An.clear();
  Ca.clear();

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

  dAngle = wireShift * TMath::TwoPi() / nWire;
  anodeLength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusA * TMath::Sin(dAngle/2),2) );
  cathodeLength = TMath::Sqrt( zLen*zLen + TMath::Power(2* radiusC * TMath::Sin(dAngle/2),2) );
}

inline void PW::FindWireID(TVector3 pos, TVector3 direction, bool verbose ){

  hitInfo.Clear();
  double phi = direction.Phi();

  for( int i = 0; i < nWire; i++){

    double disA = 99999999;
    double phiS = An[i].first.Phi()  - TMath::PiOver4();
    double phiL = An[i].second.Phi() + TMath::PiOver4();
    // printf("A%2d: %f %f | %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg(), phi * TMath::RadToDeg());
    if( phi > 0 && phiS > phiL )  phiL = phiL + TMath::TwoPi();
    if( phi < 0 && phiS > phiL )  phiS = phiS - TMath::TwoPi();

    if( phiS < phi && phi < phiL) {
      disA = Distance( pos, pos + direction, An[i].first, An[i].second);
      if( disA < hitInfo.nearestDist.first ){
        hitInfo.nearestDist.first = disA;
        hitInfo.nearestWire.first = i;
      }
    }

    double disC = 99999999;
    phiS = Ca[i].second.Phi()- TMath::PiOver4();
    phiL = Ca[i].first.Phi() + TMath::PiOver4();
    // printf("C%2d: %f %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    if( phi > 0 && phiS > phiL ) phiL = phiL + TMath::TwoPi();
    if( phi < 0 && phiS > phiL ) phiS = phiS - TMath::TwoPi();

    if(phiS < phi && phi < phiL) {
      disC = Distance( pos, pos + direction, Ca[i].first, Ca[i].second);
      if( disC < hitInfo.nearestDist.second ){
        hitInfo.nearestDist.second = disC;
        hitInfo.nearestWire.second = i;
      }
    }

    if(verbose) printf(" %2d | %8.2f, %8.2f\n", i, disA, disC);
  }

  //==== find the 2nd nearest wire
  short anode1 = hitInfo.nearestWire.first;

  double haha1 = Distance( pos, pos + direction, An[anode1-1].first, An[anode1-1].second);
  double haha2 = Distance( pos, pos + direction, An[anode1+1].first, An[anode1+1].second);
  if( haha1 < haha2){
    hitInfo.nextNearestWire.first = anode1-1;
    hitInfo.nextNearestDist.first = haha1;
  }else{
    hitInfo.nextNearestWire.first = anode1+1;
    hitInfo.nextNearestDist.first = haha2;
  }

  short cathode1 = hitInfo.nearestWire.second;

  haha1 = Distance( pos, pos + direction, Ca[cathode1-1].first, Ca[cathode1-1].second);
  haha2 = Distance( pos, pos + direction, Ca[cathode1+1].first, Ca[cathode1+1].second);
  if( haha1 < haha2){
    hitInfo.nextNearestWire.second = cathode1-1;
    hitInfo.nextNearestDist.second = haha1;
  }else{
    hitInfo.nextNearestWire.second = cathode1+1;
    hitInfo.nextNearestDist.second = haha2;
  }

  if( verbose ) Print();
}

inline void PW::CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose){

  trackPos = sx3Pos;

  TVector3 n1 = (An[anodeID].first - An[anodeID].second).Cross((sx3Pos - An[anodeID].second)).Unit();
  TVector3 n2 = (Ca[cathodeID].first - Ca[cathodeID].second).Cross((sx3Pos - Ca[cathodeID].second)).Unit();

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n2.Cross(n1)).Unit();

  if( verbose ) printf("Theta, Phi = %f, %f \n", trackVec.Theta() *TMath::RadToDeg(), trackVec.Phi()*TMath::RadToDeg()); 

}

inline void PW::CalTrack2(TVector3 sx3Pos, PWHitInfo hitInfo, bool verbose){

  trackPos = sx3Pos;

  // fraction between the nearest wire and and the 2nd nearest wire by distance
  double totDistA = hitInfo.nearestDist.first + hitInfo.nextNearestDist.first;
  double fracA = hitInfo.nearestDist.first / totDistA;
  short anodeID1 = hitInfo.nearestWire.first;
  short anodeID2 = hitInfo.nextNearestWire.first;
  TVector3 shiftA1 = (An[anodeID2].first - An[anodeID1].first) * fracA;
  TVector3 shiftA2 = (An[anodeID2].second - An[anodeID1].second) * fracA;

  double totDistC = hitInfo.nearestDist.second + hitInfo.nextNearestDist.second;
  double fracC = hitInfo.nearestDist.second / totDistA;
  short cathodeID1 = hitInfo.nearestWire.second;
  short cathodeID2 = hitInfo.nextNearestWire.second;
  TVector3 shiftC1 = (Ca[anodeID2].first - Ca[anodeID1].first) * fracC;
  TVector3 shiftC2 = (Ca[anodeID2].second - Ca[anodeID1].second) * fracC;

  TVector3 a1 = An[anodeID1].first + shiftA1;
  TVector3 a2 = An[anodeID1].second + shiftA2;

  TVector3 c1 = Ca[cathodeID1].first + shiftC1;
  TVector3 c2 = Ca[cathodeID1].second + shiftC2;

  TVector3 n1 = (a1 - a2).Cross((sx3Pos - a2)).Unit();
  TVector3 n2 = (c1 - c2).Cross((sx3Pos - c2)).Unit();

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n2.Cross(n1)).Unit();

  if( verbose ) printf("Theta, Phi = %f, %f \n", trackVec.Theta() *TMath::RadToDeg(), trackVec.Phi()*TMath::RadToDeg()); 

}
#endif 