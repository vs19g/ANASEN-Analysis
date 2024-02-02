#ifndef ClassPW_h
#define ClassPW_h

#include <cstdio>
#include <TMath.h>
#include <TVector3.h>

class PW{ // proportional wire

public: 
  PW(){ Clear(); };
  ~PW(){};

  std::pair<short, short>   GetNearestID()          const {return std::pair(anode1, cathode1);}
  std::pair<double, double> GetNearestDistance()    const {return std::pair(anodeDis1,cathodeDis1);}
  std::pair<short, short>   Get2ndNearestID()       const {return std::pair(anode2,cathode2);}
  std::pair<double, double> Get2ndNearestDistance() const {return std::pair(anodeDis2,cathodeDis2);}

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

  void Clear();
  void ConstructGeo();
  void FindWireID(TVector3 pos, TVector3 direction, bool verbose = false);
  void CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose = false);

  void Print(){
    printf("     The nearest | Anode: %2d(%5.2f) Cathode: %2d(%5.2f)\n", anode1, anodeDis1, cathode1, cathodeDis1);
    printf(" The 2nd nearest | Anode: %2d(%5.2f) Cathode: %2d(%5.2f)\n", anode2, anodeDis2, cathode2, cathodeDis2);
  }

private:

  // the nearest wire
  short anode1;
  short cathode1;
  double anodeDis1;
  double cathodeDis1;

  // the 2nd nearest wire
  short anode2;
  short cathode2;
  double anodeDis2;
  double cathodeDis2;

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

inline void PW::Clear(){
  anode1 = -1;
  cathode1 = -1;
  anodeDis1   = 999999999;
  cathodeDis1 = 999999999;
  anode2 = -1;
  cathode2 = -1;
  anodeDis2   = 999999999;
  cathodeDis2 = 999999999;
}

inline void PW::ConstructGeo(){
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

  Clear();

  double phi = direction.Phi();
  for( int i = 0; i < nWire; i++){

    double disA = 99999999;
    double disC = 99999999;

    double phiS = An[i].first.Phi()  - TMath::PiOver4();
    double phiL = An[i].second.Phi() + TMath::PiOver4();

    // printf("A%2d: %f %f | %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg(), phi * TMath::RadToDeg());
    if( phi > 0 && phiS > phiL )  phiL = phiL + TMath::TwoPi();
    if( phi < 0 && phiS > phiL )  phiS = phiS - TMath::TwoPi();

    if( phiS < phi && phi < phiL) {
      disA = Distance( pos, pos + direction, An[i].first, An[i].second);
      if( disA < anodeDis1 ){
        anodeDis1 = disA;
        anode1 = i;
      }
    }

    phiS = Ca[i].second.Phi()- TMath::PiOver4();
    phiL = Ca[i].first.Phi() + TMath::PiOver4();
    // printf("C%2d: %f %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    if( phi > 0 && phiS > phiL ) phiL = phiL + TMath::TwoPi();
    if( phi < 0 && phiS > phiL ) phiS = phiS - TMath::TwoPi();

    if(phiS < phi && phi < phiL) {
      disC = Distance( pos, pos + direction, Ca[i].first, Ca[i].second);
      if( disC < cathodeDis1 ){
        cathodeDis1 = disC;
        cathode1 = i;
      }
    }

    if(verbose) printf(" %2d | %8.2f, %8.2f\n", i, disA, disC);
  }

  //==== find the 2nd nearest wire
  double haha1 = Distance( pos, pos + direction, An[anode1-1].first, An[anode1-1].second);
  double haha2 = Distance( pos, pos + direction, An[anode1+1].first, An[anode1+1].second);
  if( haha1 < haha2){
    anode2 = anode1-1;
    anodeDis2 = haha1;
  }else{
    anode2 = anode1+1;
    anodeDis2 = haha2;
  }

  haha1 = Distance( pos, pos + direction, Ca[cathode1-1].first, Ca[cathode1-1].second);
  haha2 = Distance( pos, pos + direction, Ca[cathode1+1].first, Ca[cathode1+1].second);
  if( haha1 < haha2){
    cathode2 = cathode1-1;
    cathodeDis2 = haha1;
  }else{
    cathode2 = cathode1+1;
    cathodeDis2 = haha2;
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

#endif 