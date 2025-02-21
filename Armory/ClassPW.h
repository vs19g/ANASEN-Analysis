#ifndef ClassPW_h
#define ClassPW_h

#include <cstdio>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom.h>

struct PWHitInfo
{
  std::pair<short, short> nearestWire;   // anode, cathode
  std::pair<double, double> nearestDist; // anode, cathode

  std::pair<short, short> nextNearestWire;   // anode, cathode
  std::pair<double, double> nextNearestDist; // anode, cathode

  void Clear()
  {
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

struct Coord
{
  float x, y, z;
  Coord() : x(0), y(0), z(0) {}
  Coord(const TVector3 &vec)
  {
    x = vec.X(); // TVector3's X() returns the x-coordinate
    y = vec.Y(); // TVector3's Y() returns the y-coordinate
    z = vec.Z(); // TVector3's Z() returns the z-coordinate
  }
};

//! ########################################################
class PW
{ // proportional wire
public:
  PW() { ClearHitInfo(); };
  ~PW() {};

  PWHitInfo GetHitInfo() const { return hitInfo; }
  std::pair<short, short> GetNearestID() const { return hitInfo.nearestWire; }
  std::pair<double, double> GetNearestDistance() const { return hitInfo.nearestDist; }
  std::pair<short, short> Get2ndNearestID() const { return hitInfo.nextNearestWire; }
  std::pair<double, double> Get2ndNearestDistance() const { return hitInfo.nextNearestDist; }

  std::vector<std::pair<TVector3, TVector3>> An; // the anode wire position vector in space
  std::vector<std::pair<TVector3, TVector3>> Ca; // the cathode wire position vector in space

  TVector3 GetTrackPos() const { return trackPos; }
  TVector3 GetTrackVec() const { return trackVec; }
  double GetTrackTheta() const { return trackVec.Theta(); }
  double GetTrackPhi() const { return trackVec.Phi(); }
  double GetZ0();

  int GetNumWire() const { return nWire; }
  double GetDeltaAngle() const { return dAngle; }
  double GetAnodeLength() const { return anodeLength; }
  double GetCathodeLength() const { return cathodeLength; }
  TVector3 GetAnodeDn(short id) const { return An[id].first; }
  TVector3 GetAnodeUp(short id) const { return An[id].second; }
  TVector3 GetCathodeDn(short id) const { return Ca[id].first; }
  TVector3 GetCathodeUp(short id) const { return Ca[id].second; }

  TVector3 GetAnodneMid(short id) const { return (An[id].first + An[id].second) * 0.5; }
  double GetAnodeTheta(short id) const { return (An[id].first - An[id].second).Theta(); }
  double GetAnodePhi(short id) const { return (An[id].first - An[id].second).Phi(); }

  TVector3 GetCathodneMid(short id) const { return (Ca[id].first + Ca[id].second) * 0.5; }
  double GetCathodeTheta(short id) const { return (Ca[id].first - Ca[id].second).Theta(); }
  double GetCathodePhi(short id) const { return (Ca[id].first - Ca[id].second).Phi(); }

  void ClearHitInfo();
  void ConstructGeo();
  void FindWireID(TVector3 pos, TVector3 direction, bool verbose = false);
  void CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose = false);
  void CalTrack2(TVector3 sx3Pos, PWHitInfo hitInfo, double sigmaA = 0, double sigmaC = 0, bool verbose = false);

  void Print()
  {
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
  const float zLen = 380; // mm
  const float radiusA = 37;
  const float radiusC = 43;

  double dAngle;
  double anodeLength;
  double cathodeLength;

  // std::vector<std::pair<TVector3, TVector3>> An; // the anode wire position vector in space
  // std::vector<std::pair<TVector3, TVector3>> Ca; // the cathode wire position vector in space

  double Distance(TVector3 a1, TVector3 a2, TVector3 b1, TVector3 b2)
  {
    TVector3 na = a1 - a2;
    TVector3 nb = b1 - b2;
    TVector3 nd = (na.Cross(nb)).Unit();
    return TMath::Abs(nd.Dot(a1 - b2));
  }
};

inline void PW::ClearHitInfo()
{
  hitInfo.Clear();
}

inline void PW::ConstructGeo()
{

  An.clear();
  Ca.clear();

  std::pair<TVector3, TVector3> p1; // anode
  std::pair<TVector3, TVector3> q1; // cathode

  // anode and cathode start at pos-Y axis and count in right-Hand
  // anode wire shift is right-hand.
  // cathode wire shift is left-hand.

  for (int i = 0; i < nWire; i++)
  {
    // Anode rotate right-hand
    p1.first.SetXYZ(radiusA * TMath::Cos(TMath::TwoPi() / nWire * (i) + TMath::PiOver2()),
                    radiusA * TMath::Sin(TMath::TwoPi() / nWire * (i) + TMath::PiOver2()),
                    zLen / 2);
    p1.second.SetXYZ(radiusA * TMath::Cos(TMath::TwoPi() / nWire * (i + wireShift) + TMath::PiOver2()),
                     radiusA * TMath::Sin(TMath::TwoPi() / nWire * (i + wireShift) + TMath::PiOver2()),
                     -zLen / 2);
    An.push_back(p1);

    // Cathod rotate left-hand with the 3 wire offset accounted for (+1 from the calculated offset from the PC coincidence spectrum)
    q1.first.SetXYZ(radiusC * TMath::Cos(TMath::TwoPi() / nWire * (i + wireShift + 1) + TMath::PiOver2()),
                    radiusC * TMath::Sin(TMath::TwoPi() / nWire * (i + wireShift + 1) + TMath::PiOver2()),
                    zLen / 2);
    q1.second.SetXYZ(radiusC * TMath::Cos(TMath::TwoPi() / nWire * (i + 1) + TMath::PiOver2()),
                     radiusC * TMath::Sin(TMath::TwoPi() / nWire * (i + 1) + TMath::PiOver2()),
                     -zLen / 2);
    Ca.push_back(q1);
  }
  // correcting for the fact that the order of the cathode wires is reversed
  std::reverse(Ca.begin(), Ca.end());
  // adjusting for the 3 wire offset, the rbegin and rend are used as the rotation of the wires is done in the opposite direction i.e. 1,2,3 -> 3,1,2
  // NOT NECESSARY ANY MORE, HAS BEEN IMCORPORATED INTO THE WIREOFFSET IN THE BEGINNING
  // std::rotate(Ca.rbegin(), Ca.rbegin() + 4, Ca.rend());

  dAngle = wireShift * TMath::TwoPi() / nWire;
  anodeLength = TMath::Sqrt(zLen * zLen + TMath::Power(2 * radiusA * TMath::Sin(dAngle / 2), 2));
  cathodeLength = TMath::Sqrt(zLen * zLen + TMath::Power(2 * radiusC * TMath::Sin(dAngle / 2), 2));
}

inline void PW::FindWireID(TVector3 pos, TVector3 direction, bool verbose)
{

  hitInfo.Clear();
  double phi = direction.Phi();

  for (int i = 0; i < nWire; i++)
  {

    double disA = 99999999;
    double phiS = An[i].first.Phi() - TMath::PiOver4();
    double phiL = An[i].second.Phi() + TMath::PiOver4();
    // printf("A%2d: %f %f | %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg(), phi * TMath::RadToDeg());
    if (phi > 0 && phiS > phiL)
      phiL = phiL + TMath::TwoPi();
    if (phi < 0 && phiS > phiL)
      phiS = phiS - TMath::TwoPi();

    if (phiS < phi && phi < phiL)
    {
      disA = Distance(pos, pos + direction, An[i].first, An[i].second);
      if (disA < hitInfo.nearestDist.first)
      {
        hitInfo.nearestDist.first = disA;
        hitInfo.nearestWire.first = i;
      }
    }

    double disC = 99999999;
    phiS = Ca[i].second.Phi() - TMath::PiOver4();
    phiL = Ca[i].first.Phi() + TMath::PiOver4();
    // printf("C%2d: %f %f\n", i, phiS * TMath::RadToDeg(), phiL * TMath::RadToDeg());
    if (phi > 0 && phiS > phiL)
      phiL = phiL + TMath::TwoPi();
    if (phi < 0 && phiS > phiL)
      phiS = phiS - TMath::TwoPi();

    if (phiS < phi && phi < phiL)
    {
      disC = Distance(pos, pos + direction, Ca[i].first, Ca[i].second);
      if (disC < hitInfo.nearestDist.second)
      {
        hitInfo.nearestDist.second = disC;
        hitInfo.nearestWire.second = i;
      }
    }

    if (verbose)
      printf(" %2d | %8.2f, %8.2f\n", i, disA, disC);
  }

  //==== find the 2nd nearest wire
  short anode1 = hitInfo.nearestWire.first;
  short aaa1 = anode1 - 1;
  if (aaa1 < 0)
    aaa1 += nWire;
  short aaa2 = (anode1 + 1) % nWire;

  double haha1 = Distance(pos, pos + direction, An[aaa1].first, An[aaa1].second);
  double haha2 = Distance(pos, pos + direction, An[aaa2].first, An[aaa2].second);
  if (haha1 < haha2)
  {
    hitInfo.nextNearestWire.first = aaa1;
    hitInfo.nextNearestDist.first = haha1;
  }
  else
  {
    hitInfo.nextNearestWire.first = aaa2;
    hitInfo.nextNearestDist.first = haha2;
  }

  short cathode1 = hitInfo.nearestWire.second;
  short ccc1 = cathode1 - 1;
  if (ccc1 < 0)
    ccc1 += nWire;
  short ccc2 = (cathode1 + 1) % nWire;

  haha1 = Distance(pos, pos + direction, Ca[ccc1].first, Ca[ccc1].second);
  haha2 = Distance(pos, pos + direction, Ca[ccc2].first, Ca[ccc2].second);
  if (haha1 < haha2)
  {
    hitInfo.nextNearestWire.second = ccc1;
    hitInfo.nextNearestDist.second = haha1;
  }
  else
  {
    hitInfo.nextNearestWire.second = ccc2;
    hitInfo.nextNearestDist.second = haha2;
  }

  if (verbose)
    Print();
}

inline void PW::CalTrack(TVector3 sx3Pos, int anodeID, int cathodeID, bool verbose)
{

  trackPos = sx3Pos;

  TVector3 n1 = (An[anodeID].first - An[anodeID].second).Cross((sx3Pos - An[anodeID].second)).Unit();
  TVector3 n2 = (Ca[cathodeID].first - Ca[cathodeID].second).Cross((sx3Pos - Ca[cathodeID].second)).Unit();

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n2.Cross(n1)).Unit();

  if (verbose)
    printf("Theta, Phi = %f, %f \n", trackVec.Theta() * TMath::RadToDeg(), trackVec.Phi() * TMath::RadToDeg());
}

inline void PW::CalTrack2(TVector3 sx3Pos, PWHitInfo hitInfo, double sigmaA, double sigmaC, bool verbose)
{

  trackPos = sx3Pos;

  double p1 = TMath::Abs(hitInfo.nearestDist.first + gRandom->Gaus(0, sigmaA));
  double p2 = TMath::Abs(hitInfo.nextNearestDist.first + gRandom->Gaus(0, sigmaA));
  double fracA = p1 / (p1 + p2);
  short anodeID1 = hitInfo.nearestWire.first;
  short anodeID2 = hitInfo.nextNearestWire.first;
  TVector3 shiftA1 = (An[anodeID2].first - An[anodeID1].first) * fracA;
  TVector3 shiftA2 = (An[anodeID2].second - An[anodeID1].second) * fracA;

  double q1 = TMath::Abs(hitInfo.nearestDist.second + gRandom->Gaus(0, sigmaC));
  double q2 = TMath::Abs(hitInfo.nextNearestDist.second + gRandom->Gaus(0, sigmaC));
  double fracC = q1 / (q1 + q2);
  short cathodeID1 = hitInfo.nearestWire.second;
  short cathodeID2 = hitInfo.nextNearestWire.second;
  TVector3 shiftC1 = (Ca[cathodeID2].first - Ca[cathodeID1].first) * fracC;
  TVector3 shiftC2 = (Ca[cathodeID2].second - Ca[cathodeID1].second) * fracC;

  TVector3 a1 = An[anodeID1].first + shiftA1;
  TVector3 a2 = An[anodeID1].second + shiftA2;

  TVector3 c1 = Ca[cathodeID1].first + shiftC1;
  TVector3 c2 = Ca[cathodeID1].second + shiftC2;

  TVector3 n1 = (a1 - a2).Cross((sx3Pos - a2)).Unit();
  TVector3 n2 = (c1 - c2).Cross((sx3Pos - c2)).Unit();

  // if the handiness of anode and cathode revered, it should be n2 cross n1
  trackVec = (n2.Cross(n1)).Unit();

  if (verbose)
    printf("Theta, Phi = %f, %f \n", trackVec.Theta() * TMath::RadToDeg(), trackVec.Phi() * TMath::RadToDeg());
}

inline double PW::GetZ0()
{

  double x = trackPos.X();
  double y = trackPos.Y();
  double rho = TMath::Sqrt(x * x + y * y);
  double theta = trackVec.Theta();

  return trackPos.Z() - rho / TMath::Tan(theta);
}

#endif