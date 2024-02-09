#ifndef ClassSX3_h
#define ClassSX3_h

#include <cstdio>
#include <TMath.h>
#include <TVector3.h>

class SX3{
public:
  SX3(){Clear();};
  ~SX3(){}
  
  short GetID() const {return id;}
  short GetChUp() const {return chUp;}
  short GetChDn() const {return chDn;}
  short GetChBk() const {return chBk;}

  TVector3 GetHitPos() const {return hitPos;}
  TVector3 GetHitPosWithSigma(double sigmaY_mm, double sigmaZ_mm);

  double GetZFrac() const {return zFrac;} // range from -0.5 to 0.5

  void Clear();
  void ConstructGeo();
  void FindSX3Pos(TVector3 pos, TVector3 direction, bool verbose = false);
  void CalSX3Pos(unsigned short ID, unsigned short chUp, unsigned short chDown, unsigned short chBack, float eUp, float eDown);

  double GetNumDet() const {return numDet;}
  double GetWidth()  const {return width;}
  double GetLength() const {return length;}
  TVector3 GetDnL(short id) const {return SDn[id].first; } // lower strip ID
  TVector3 GetDnH(short id) const {return SDn[id].second; } // higher strip ID
  TVector3 GetUpL(short id) const {return SUp[id].first; } // lower strip ID
  TVector3 GetUpH(short id) const {return SUp[id].second; } // higher strip ID

  TVector3 GetDnMid(short id) const { return (SDn[id].first + SDn[id].second)*0.5;}
  TVector3 GetUpMid(short id) const { return (SUp[id].first + SUp[id].second)*0.5;}

  double GetDetPhi(short id) const { return (SUp[id].second - SUp[id].first).Phi();}

  void Print(){ 
    if( id == -1 ){
      printf("Did not hit any SX3.\n"); 
    }else{
      printf("ID: %d, U,D,B: %d %d %d| zFrac : %.2f\n", id, chUp, chDn, chBk, zFrac); 
      printf("Hit Pos: %.2f, %.2f, %.2f\n", hitPos.X(), hitPos.Y(), hitPos.Z()); 
    }
  }

  // void CalZFrac(){
  //   zFrac = (eUp - eDn)/(eUp + eDn); 
  // }

private:

  const int numDet = 12;
  const float radius = 88;
  const float width = 40;
  const float length = 75;
  const float gap = 46;

  short id; // -1 when no hit
  short chUp;
  short chDn;
  short chBk;

  double zFrac; // from +1 (downstream) to -1 (upstream)

  double eUp;
  double eDn;
  double eBk;

  TVector3 hitPos;

  std::vector<std::pair<TVector3,TVector3>> SDn; // coners of the SX3  0-11, z = mid point
  std::vector<std::pair<TVector3,TVector3>> SUp; // coners of the SX3 12-23, z = mid point
  std::vector<TVector3> SNorml; // normal of the SX3 (outward)


  std::pair<double, double> Intersect(TVector3 p1, TVector3 p2, TVector3 q1, TVector3 q2, bool verbose){

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

};

inline void SX3::Clear(){
  id = -1;
  chUp = -1;
  chDn = -1;
  chBk = -1;
  zFrac = TMath::QuietNaN();

  eUp = TMath::QuietNaN();
  eDn = TMath::QuietNaN();
  eBk = TMath::QuietNaN();

  SDn.clear();
  SUp.clear();
}

inline void SX3::ConstructGeo(){
  TVector3 sa, sb, sc, sn;

  for(int i = 0; i < numDet; i++){
    sa.SetXYZ( radius, -width/2, gap/2 + length/2 );
    sb.SetXYZ( radius,  width/2, gap/2 + length/2 );

    double rot = TMath::TwoPi() / numDet * (i - 0.5) - TMath::PiOver2();

    sa.RotateZ( rot );
    sb.RotateZ( rot );
    SDn.push_back(std::pair(sa,sb));

    sc.SetXYZ( radius, -width/2, gap/2 );
    sc.RotateZ( rot );

    sn = ((sc-sa).Cross(sb-sa)).Unit();
    SNorml.push_back(sn);

    sa.SetXYZ( radius, -width/2, -gap/2 - length/2 );
    sb.SetXYZ( radius,  width/2, -gap/2 - length/2 );

    sa.RotateZ( rot );
    sb.RotateZ( rot );
    SUp.push_back(std::pair(sa,sb));
  }
}


inline void SX3::FindSX3Pos(TVector3 pos, TVector3 direction, bool verbose){

  id = -1;
  for( int i = 0 ; i < numDet; i++){

    if(verbose) printf(" %d ", i);
    std::pair<double, double> frac = Intersect( pos, pos + direction, SDn[i].first, SDn[i].second, verbose);


    if( frac.second < 0 || frac.second > 1 ) continue;
    hitPos = pos + frac.first * direction;

    double dis = hitPos.Dot(SNorml[i]);

    if(verbose) {
      printf("reduced distance : %f\n", dis);
      printf(" %d*", (i+1)%numDet); 
      Intersect( pos, pos + direction, SDn[(i+1)%numDet].first, SDn[(i+1)%numDet].second, verbose);
    }

    if( TMath::Abs(dis - radius) > 0.1 ) continue;

    chDn = 2 * TMath::Floor(frac.second * 4);
    chUp = chDn + 1;

    double zPos = hitPos.Z();
    if( (gap/2 < zPos && zPos < gap/2 + length ) || (-gap/2 - length < zPos && zPos < -gap/2 )  ){

      id = zPos > 0 ? i : i + 12;

      zFrac = zPos > 0 ?  (zPos - gap/2. - length/2.)/length : (zPos - ( - gap/2. - length/2.) )/length ;

      chBk = TMath::Floor( (zFrac + 0.5) * 4 ) + 8;

      if( verbose) Print();
      return ;

    }else{
      if( verbose ) printf(" zPos out of sensitive region\n");
    }
  }

  if( verbose) Print();
}

inline TVector3 SX3::GetHitPosWithSigma(double sigmaY_mm, double sigmaZ_mm){

  double phi = SNorml[id%numDet].Phi();

  TVector3 haha = hitPos;
  haha.RotateZ(-phi);

  double y = haha.Y() + gRandom->Gaus(0, sigmaY_mm);
  if( sigmaY_mm < 0 ){
    double deltaW = width/4;
    y = TMath::Floor((haha.Y()-deltaW)/deltaW)*deltaW + deltaW*1.5; // when ever land on each strip, set the position to be center of the strip.
    if( y >= 25 ) y = 15;
  }

  double z = haha.Z() + gRandom->Gaus(0, sigmaZ_mm);
  if( sigmaZ_mm < 0 ){
    haha.Z();
    double delta = length/4;
    int sign = z > 0 ? 1 : -1;
    z = TMath::Floor( (abs(z)-gap/2)/delta )*delta + 0.5 * delta + gap/2;
    if( z >= 107.375 ) z = 88.625;
    z = sign * z;
  }

  haha.SetY(y);
  haha.SetZ(z);
  haha.RotateZ(phi);

  return haha;

}


inline void SX3::CalSX3Pos(unsigned short ID, unsigned short chUp, unsigned short chDown, unsigned short chBack, float eUp, float eDown){

  hitPos.Clear();

  if( (chUp - chDown) != 1 || (chDown % 2) != 0) return ;

  int reducedID = ID % numDet;

  TVector3 sa, sb;

  if( ID < numDet ){ //down
    sa = SDn[reducedID].second;
    sb = SDn[reducedID].first;
  }else{
    sa = SUp[reducedID].second;
    sb = SUp[reducedID].first;
  }

  hitPos.SetX( (sb.X() - sa.X()) * chUp/8 + sa.X());
  hitPos.SetY( (sb.Y() - sa.Y()) * chUp/8 + sa.Y());

  if( eUp == 0 || eDown == 0 ){
    hitPos.SetZ( sa.Z() + (2*(chBk - 7)-1) * length / 8 );
  }else{
    double frac = (eUp - eDown)/(eUp + eDown); // from +1 (downstream) to -1 (upstream)
    double zPos = sa.Z() +  length * frac/2;
    hitPos.SetZ( zPos );
  }

}

#endif 