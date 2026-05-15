#ifndef ClassPW_h
#define ClassPW_h

#include <cstdio>
#include <iostream>
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

  Coord Crossover[24][24][2];

  inline TVector3 getClosestWirePosAtWirePhi(std::pair<TVector3, TVector3>, double phi);
  inline std::tuple<std::pair<TVector3, TVector3>, double, double, double> GetPseudoWire(const std::vector<std::tuple<int, double, double>> &cluster, std::string type);

  inline std::tuple<TVector3, double, double, double, double, double, double, double>
  FindCrossoverProperties(const std::vector<std::tuple<int, double, double>> &a_cluster, const std::vector<std::tuple<int, double, double>> &c_cluster);

  inline std::vector<std::vector<std::tuple<int, double, double>>>
  Make_Clusters(std::unordered_map<int, std::tuple<int, double, double>> wireEvents);

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
  void CalTrack2(TVector3 sx3Pos, TVector3 anodeInt, bool verbose = false);

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
  const int wireShift = 4;
  // const float zLen = 380; // mm
  //  const float zLen = 348.6; // mm
  const float zLen = 174.3 * 2; // mm
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

  double k = TMath::TwoPi() / 24.; // 48 solder thru holes, wires in every other one
  double offset_a1 = -6 * k - 3 * k;
  double offset_c1 = -4 * k - 2 * k - TMath::TwoPi() / 48; // correct for a half-turn
  // std::cerr << "Here!" << std::endl;
  // #include "../scratch/testing.h"
  double offset_a2 = offset_a1 + wireShift * k;
  double offset_c2 = offset_c1 - wireShift * k;

  for (int i = 0; i < nWire; i++)
  {
    // Anode rotate right-hand coming in towards +z riding with the beam. In this frame, +x is to the right, and +y down
    // updated Feb 2026, Sudarsan B. Photographs indicate that anode wires twist right handed, as one moves from -z to +z with the convention above
    // wire indices increase leftward as one moves to +z (hence -k factor), but wires themselves twist rightward - as indicated by offset_a2 being more +ve w.r.t offset_a1
    //'First' is -z locus, 'second' is +z locus
    p1.first.SetXYZ(radiusA * TMath::Cos(-k * i + offset_a1),
                    radiusA * TMath::Sin(-k * i + offset_a1),
                    -zLen / 2);
    p1.second.SetXYZ(radiusA * TMath::Cos(-k * i + offset_a2),
                     radiusA * TMath::Sin(-k * i + offset_a2),
                     +zLen / 2);

    // Cathodes twist left-hand as indicated by offset_c2 being more negative than offset_c1, under the same system, while wires increase rightward (hence +k factor)
    q1.first.SetXYZ(radiusC * TMath::Cos(k * i + offset_c1),
                    radiusC * TMath::Sin(k * i + offset_c1),
                    -zLen / 2);
    q1.second.SetXYZ(radiusC * TMath::Cos(k * i + offset_c2),
                     radiusC * TMath::Sin(k * i + offset_c2),
                     zLen / 2);
    An.push_back(p1);
    Ca.push_back(q1);
  }

  // Calculate Crossover Geometry ONCE
  TVector3 a, c, diff;
  double a2, ac, c2, adiff, cdiff, denom, alpha;

  for (size_t i = 0; i < An.size(); i++)
  {
    // a = An[i].first - An[i].second;
    a = An[i].second - An[i].first;
    for (size_t j = 0; j < Ca.size(); j++)
    {
      c = Ca[j].second - Ca[j].first;
      diff = An[i].second - Ca[j].second;
      a2 = a.Dot(a);
      c2 = c.Dot(c);
      ac = a.Dot(c);
      adiff = a.Dot(diff);
      cdiff = c.Dot(diff);
      denom = a2 * c2 - ac * ac;
      alpha = (ac * cdiff - c2 * adiff) / denom;

      Crossover[i][j][0].x = An[i].second.X() + alpha * a.X();
      Crossover[i][j][0].y = An[i].second.Y() + alpha * a.Y();
      Crossover[i][j][0].z = An[i].second.Z() + alpha * a.Z();

      if (Crossover[i][j][0].z < -190 || Crossover[i][j][0].z > 190)
      {
        // std::cout << "Weird crossover but ok" << std::endl;
      }
      if ((i + j) % 24 == 12 || Crossover[i][j][0].z < -190 || Crossover[i][j][0].z > 190)
      {
        Crossover[i][j][0].z = 9999999;
        // std::cout << "Weird crossover" << std::endl;
      }

      Crossover[i][j][1].x = alpha;
      Crossover[i][j][1].y = 0;
    }
  }

  dAngle = wireShift * TMath::TwoPi() / nWire;
  anodeLength = TMath::Sqrt(zLen * zLen + TMath::Power(2 * radiusA * TMath::Sin(dAngle / 2), 2));
  cathodeLength = TMath::Sqrt(zLen * zLen + TMath::Power(2 * radiusC * TMath::Sin(dAngle / 2), 2)); // chord length subtending an angle alpha is 2rsin(alpha/2)
}

inline TVector3 PW::getClosestWirePosAtWirePhi(std::pair<TVector3, TVector3> awire, double sx3phi_radian)
{
  // 1. Get wire geometry
  TVector3 a1 = awire.first;  // Top of the wire
  TVector3 a2 = awire.second; // Bottom of the wire
  TVector3 wireVec = a2 - a1; // Vector pointing down the wire

  // Variables to track our minimums during the scan
  double min_delta_phi = 9999.0;
  double best_t = -1.0;
  TVector3 best_pcz_intersect;

  // 2. THE SCAN: Walk down the wire in 1000 tiny steps
  // (For a 380mm wire, this is checking every 0.38 mm)
  int num_steps = 1000;
  for (int i = 0; i <= num_steps; ++i)
  {
    double t_test = (double)i / num_steps;    // Ranges from 0.0 to 1.0
    TVector3 test_pt = a1 + t_test * wireVec; // The 3D point at this step

    // Calculate absolute Delta Phi between Si hit and this specific point on the wire
    if (TMath::IsNaN(sx3phi_radian - test_pt.Phi()))
      continue;
    double dPhi = TMath::Abs(TVector2::Phi_mpi_pi(sx3phi_radian - test_pt.Phi())); // Phi_mpi_pi just puts the angle in the range -180 to 180

    // If this is the smallest Delta Phi we've seen so far, save it!
    if (dPhi < min_delta_phi)
    {
      min_delta_phi = dPhi;
      best_t = t_test;
      best_pcz_intersect = test_pt;
    }
  }
  return best_pcz_intersect;
}

inline std::vector<std::vector<std::tuple<int, double, double>>>
PW::Make_Clusters(std::unordered_map<int, std::tuple<int, double, double>> wireEvents)
{
  std::vector<std::vector<std::tuple<int, double, double>>> wireClusters;
  std::vector<std::tuple<int, double, double>> wireCluster;
  // TODO: Write a macro once, call it twice
  int wirecount = 0;
  while (wirecount < 24)
  {
    if (wireEvents.find(wirecount) == wireEvents.end())
    {
      wirecount++;
      continue;
    }
    wireCluster.clear();
    int ctr2 = wirecount;
    do
    {
      wireCluster.emplace_back(wireEvents[ctr2]);
      ctr2 += 1;
      if (ctr2 == 24 || ctr2 - wirecount == 7)
        break; // loose logic, needs to be looked at.
    } while (wireEvents.find(ctr2) != wireEvents.end());
    wireClusters.push_back(std::move(wireCluster));
    wirecount = ctr2; // we already dealt with wires until the last value of ctr2
  }

  if (wireClusters.size() > 1)
  {                                            // Deal with wraparound if required
    auto first_cluster = wireClusters.front(); // front and back provide references to the elements themselves. less copy, can modify etc
    auto last_cluster = wireClusters.back();
    if (std::get<0>(last_cluster.back()) == 23 && std::get<0>(first_cluster.front()) == 0)
    {
      last_cluster.insert(last_cluster.end(), first_cluster.begin(), first_cluster.end());
    }
    wireClusters.erase(wireClusters.begin()); // canonically, erase() needs an iterator, hence begin() not front()
    // TODO: Can also deal with 'gaps' of missing wires similarly. end of one segment and beginning of another segment will be separated by missing wire --> combine the two
    // TODO: Also needs some development regarding the time-correlation. Don't put wires in the same cluster if they aren't time coincident
  }
  return wireClusters;

  /*if(aClusters.size()>1 || cClusters.size() > 1) {
    std::cout << " ============== " << std::endl;
  }
  if(aClusters.size()>1 && cClusters.size() >=1) {
    std::cout << aClusters.size() << " new anode clusters ----> " << std::endl;
    int cc=1;
    for(auto ac : aClusters) {
      std::cout << "  Cluster " << cc << std::endl;
      double first_ts = std::get<2>(ac.at(0));
      for(auto item : ac) {
        std::cout << "  \t" << std::get<0>(item) << " " << std::get<1>(item) << " " << std::get<2>(item)-first_ts << std::endl;
      }
      std::cout << "  ------" << std::endl;
      cc++;
    }
  }

  if(cClusters.size()>=1 ) {
    std::cout << cClusters.size() << " new cathode clusters ----> " << std::endl;
    int cc=1;
    for(auto ac : cClusters) {
      std::cout << "  Cluster " << cc << std::endl;
      double first_ts = std::get<2>(ac.at(0));
      for(auto item : ac) {
        std::cout << "  \t" << std::get<0>(item) << " " << std::get<1>(item) << " " << std::get<2>(item)-first_ts << std::endl;
      }
      std::cout << "  ------" << std::endl;
      cc++;
    }
  } 	*/
}

inline std::tuple<std::pair<TVector3, TVector3>, double, double, double>
PW::GetPseudoWire(const std::vector<std::tuple<int, double, double>> &cluster, std::string type)
{
  std::pair<TVector3, TVector3> avgvec = std::pair(TVector3(0, 0, 0), TVector3(0, 0, 0));
  double sumEnergy = 0;
  double maxEnergy = 0;
  double tsMaxEnergy = 0;
  if (type == "ANODE")
  {
    // if(cluster.size()>1) std::cout << " -------anodes" << std::endl;
    for (auto wire : cluster)
    {
      avgvec.first += std::get<1>(wire) * TVector3(An.at(std::get<0>(wire)).first.X(), An.at(std::get<0>(wire)).first.Y(), 0);
      avgvec.second += std::get<1>(wire) * TVector3(An.at(std::get<0>(wire)).second.X(), An.at(std::get<0>(wire)).second.Y(), 0);
      sumEnergy += std::get<1>(wire);
      if (std::get<1>(wire) > maxEnergy)
      {
        maxEnergy = std::get<1>(wire);
        tsMaxEnergy = std::get<2>(wire);
      }
      /*if(cluster.size()>1) {
        std::cout << "\t\t ch:" << std::get<0>(wire) << " " << std::get<1>(wire) << " " << std::get<2>(wire) << std::endl;
        std::cout << "\t\t w1(r,phi,z):" << An.at(std::get<0>(wire)).first.Perp() << " " << An.at(std::get<0>(wire)).first.Phi()*180/M_PI << " " << An.at(std::get<0>(wire)).first.Z() << std::endl;
        std::cout << "\t\t w2(r,phi,z):" << An.at(std::get<0>(wire)).second.Perp() << " " << An.at(std::get<0>(wire)).second.Phi()*180/M_PI << " " << An.at(std::get<0>(wire)).second.Z() << std::endl;
      }*/
    }
    avgvec.first = avgvec.first * (1.0 / sumEnergy);
    avgvec.second = avgvec.second * (1.0 / sumEnergy);
    double phi1 = avgvec.first.Phi();
    double phi2 = avgvec.second.Phi();
    avgvec.first.SetXYZ(radiusA * TMath::Cos(phi1), radiusA * TMath::Sin(phi1), -zLen / 2);
    avgvec.second.SetXYZ(radiusA * TMath::Cos(phi2), radiusA * TMath::Sin(phi2), zLen / 2);
    /*if(cluster.size()>1) {
      std::cout << "\t\t avg1(r,phi,z):" << avgvec.first.Perp() << " " << avgvec.first.Phi()*180/M_PI << " " << avgvec.first.Z() << std::endl;
      std::cout << "\t\t avg2(r,phi,z):" << avgvec.second.Perp() << " " << avgvec.second.Phi()*180/M_PI << " " << avgvec.second.Z() << std::endl;
    }*/
  }
  else if (type == "CATHODE")
  {
    for (auto wire : cluster)
    {
      avgvec.first += std::get<1>(wire) * TVector3(Ca.at(std::get<0>(wire)).first.X(), Ca.at(std::get<0>(wire)).first.Y(), 0);
      avgvec.second += std::get<1>(wire) * TVector3(Ca.at(std::get<0>(wire)).second.X(), Ca.at(std::get<0>(wire)).second.Y(), 0);
      sumEnergy += std::get<1>(wire);
      if (std::get<1>(wire) > maxEnergy)
      {
        maxEnergy = std::get<1>(wire);
        tsMaxEnergy = std::get<2>(wire);
      }
    }
    avgvec.first = avgvec.first * (1.0 / sumEnergy);
    avgvec.second = avgvec.second * (1.0 / sumEnergy);
    double phi1 = avgvec.first.Phi();
    double phi2 = avgvec.second.Phi();
    avgvec.first.SetXYZ(radiusC * TMath::Cos(phi1), radiusC * TMath::Sin(phi1), -zLen / 2);
    avgvec.second.SetXYZ(radiusC * TMath::Cos(phi2), radiusC * TMath::Sin(phi2), zLen / 2);
  }
  return std::tuple(avgvec, sumEnergy, maxEnergy, tsMaxEnergy);
}

inline std::tuple<TVector3, double, double, double, double, double, double, double> PW::FindCrossoverProperties(const std::vector<std::tuple<int, double, double>> &a_cluster,
                                                                                                                const std::vector<std::tuple<int, double, double>> &c_cluster)
{
  // std::pair<TVector3, TVector3> apwire = GetPseudoWire(a_cluster,"ANODE",anodeSumE);
  // std::pair<TVector3, TVector3> cpwire = GetPseudoWire(c_cluster,"CATHODE",cathodeSumE);
  auto [apwire, apSumE, apMaxE, apTSMaxE] = GetPseudoWire(a_cluster, "ANODE");
  auto [cpwire, cpSumE, cpMaxE, cpTSMaxE] = GetPseudoWire(c_cluster, "CATHODE");

  TVector3 crossover;
  crossover.Clear();
  TVector3 a, c, diff;
  double a2, ac, c2, adiff, cdiff, denom, alpha = 0;

  if (apSumE && cpSumE)
  {
    a = apwire.first - apwire.second;
    c = cpwire.first - cpwire.second;
    diff = apwire.first - cpwire.first;
    a2 = a.Dot(a);
    c2 = c.Dot(c);
    ac = a.Dot(c);
    adiff = a.Dot(diff);
    cdiff = c.Dot(diff);
    denom = a2 * c2 - ac * ac;
    alpha = (ac * cdiff - c2 * adiff) / denom;
    crossover = apwire.first + alpha * a;
    if (crossover.z() < -190 || crossover.Z() > 190)
    {
      alpha = 9999999;
      apSumE = -1;
      cpSumE = -1;
      apMaxE = -1;
      cpMaxE = -1;
      apTSMaxE = -1;
      cpTSMaxE = -1;
    }
  }
  // std::cout << apSumE << " " << cpSumE << " " << " " <<  crossover.Perp() << std::endl;
  return std::tuple(crossover, alpha, apSumE, cpSumE, apMaxE, cpMaxE, apTSMaxE, cpTSMaxE);
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

inline void PW::CalTrack2(TVector3 siPos, TVector3 anodeInt, bool verbose)
{

  double mx, my;
  double z;
  mx = siPos.X() / (siPos.X() - anodeInt.X());
  my = siPos.Y() / (siPos.Y() - anodeInt.Y());
  z = siPos.Z() + mx * (anodeInt.Z() - siPos.Z());
  // if (mx == my)
  {
    trackVec = TVector3(0, 0, z);
  }

  if (verbose)
    printf("X slope = %f and Y slope = %f \n", mx, my);
}

/*inline TVector3 PW::CalTrack3(TVector3 siPos, TVector3 anodeInt, bool verbose)
{

  TVector3 v = anodeInt-siPos;
  double t_minimum = -1.0*(siPos.X()*v.X()+siPos.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
  TVector3 vector_closest_to_z = siPos + t_minimum*v;

  return vector_closest_to_z;
  if (verbose)
    printf("X slope = %f and Y slope = %f \n", mx, my);
}*/

inline double PW::GetZ0()
{

  double x = trackPos.X();
  double y = trackPos.Y();
  double rho = TMath::Sqrt(x * x + y * y);
  double theta = trackVec.Theta();

  return trackVec.Z();
}

#endif
