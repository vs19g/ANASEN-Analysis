#define Analyzer_cxx

#include "Analyzer.h"
#include "Armory/ClassSX3.h"
#include "Armory/ClassPW.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include "TVector3.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>

TH2F *hsx3IndexVE;
TH2F *hqqqIndexVE;
TH2F *hpcIndexVE;

TH2F *hpcIndexVE_GM;
TH2F *hsx3Coin;
TH2F *hqqqCoin;
TH2F *hpcCoin;
TH2F *hAVCcoin;

TH2F *hqqqPolar;
TH2F *hsx3VpcIndex;
TH2F *hqqqVpcIndex;
TH2F *hqqqVpcE;
TH2F *hsx3VpcE;
TH2F *hanVScatsum;
TH2F *hanVScatsum_a[24];
TH1F *hPC_E[48];
TH1F *hCat4An;
TH1F *hCat0An;
TH1F *hAnodehits;
TH2F *hNosvAe;

int padID = 0;

SX3 sx3_contr;
PW pw_contr;
PW pwinstance;
TVector3 hitPos;
// TVector3 anodeIntersection;
std::map<int, std::pair<double, double>> slopeInterceptMap;

const int MAX_DET = 24;
const int MAX_UP = 4;
const int MAX_DOWN = 4;
const int MAX_BK = 4;
double backGain[MAX_DET][MAX_BK] = {{0}};
bool backGainValid[MAX_DET][MAX_BK] = {{false}};
double frontGain[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool frontGainValid[MAX_DET][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

TH1F *hZProj;
TH1F *hPCZProj;

void Analyzer::Begin(TTree * /*tree*/)
{
  TString option = GetOption();

  hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
  hsx3IndexVE->SetNdivisions(-612, "x");
  hqqqIndexVE = new TH2F("hqqqIndexVE", "QQQ index vs Energy; QQQ index ; Energy", 4 * 2 * 16, 0, 4 * 2 * 16, 400, 0, 5000);
  hqqqIndexVE->SetNdivisions(-1204, "x");
  hpcIndexVE = new TH2F("hpcIndexVE", "PC index vs Energy; PC index ; Energy", 2 * 24, 0, 2 * 24, 400, 0, 16000);
  hpcIndexVE->SetNdivisions(-1204, "x");
  hpcIndexVE_GM = new TH2F("hpcIndexVE_GM", "PC index vs Energy; PC index ; Energy", 2 * 24, 0, 2 * 24, 400, 0, 16000);
  hpcIndexVE_GM->SetNdivisions(-1204, "x");

  hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);
  hqqqCoin = new TH2F("hqqqCoin", "QQQ Coincident", 4 * 2 * 16, 0, 4 * 2 * 16, 4 * 2 * 16, 0, 4 * 2 * 16);
  hpcCoin = new TH2F("hpcCoin", "PC Coincident", 2 * 24, 0, 2 * 24, 2 * 24, 0, 2 * 24);
  hAVCcoin = new TH2F("hAVCcoin", "Anode vs Cathode Coincident", 24, 0, 24, 24, 0, 24);

  hqqqPolar = new TH2F("hqqqPolar", "QQQ Polar ID", 16 * 4, -TMath::Pi(), TMath::Pi(), 16, 10, 50);

  hsx3VpcIndex = new TH2F("hsx3Vpcindex", "sx3 vs pc; sx3 index; pc index", 24 * 12, 0, 24 * 12, 48, 0, 48);
  hsx3VpcIndex->SetNdivisions(-612, "x");
  hsx3VpcIndex->SetNdivisions(-12, "y");
  hqqqVpcIndex = new TH2F("hqqqVpcindex", "qqq vs pc; qqq index; pc index", 4 * 2 * 16, 0, 4 * 2 * 16, 48, 0, 48);
  hqqqVpcIndex->SetNdivisions(-612, "x");
  hqqqVpcIndex->SetNdivisions(-12, "y");

  hqqqVpcE = new TH2F("hqqqVpcEnergy", "qqq vs pc; qqq energy; pc energy", 400, 0, 5000, 800, 0, 16000);
  hqqqVpcE->SetNdivisions(-612, "x");
  hqqqVpcE->SetNdivisions(-12, "y");

  hsx3VpcE = new TH2F("hsx3VpcEnergy", "sx3 vs pc; sx3 energy; pc energy", 400, 0, 5000, 800, 0, 16000);
  hsx3VpcE->SetNdivisions(-612, "x");
  hsx3VpcE->SetNdivisions(-12, "y");

  hZProj = new TH1F("hZProj", "Z Projection", 1200, -600, 600);
  hPCZProj = new TH1F("hPCZProj", "PC Z Projection", 600, -300, 300);

  hanVScatsum = new TH2F("hanVScatsum", "Anode vs Cathode Sum; Anode E; Cathode E", 400, 0, 16000, 400, 0, 20000);
  hCat4An = new TH1F("hCat4An", "Number of Cathodes/Anode", 24, 0, 24);
  hCat0An = new TH1F("hCat0An", "Number of Cathodes without Anode", 24, 0, 24);
  hAnodehits = new TH1F("hAnodehits", "Number of Anode hits", 24, 0, 24);
  hNosvAe = new TH2F("hnosvAe", "Number of Cathodes/Anode vs Anode Energy", 20, 0, 20, 400, 0, 16000);
  // for (int i = 0; i < 24; i++)
  // {
  //   TString histName = Form("hAnodeVsCathode_%d", i);
  //   TString histTitle = Form("Anode %d vs Cathode Sum; Anode E; Cathode Sum E", i);
  //   hanVScatsum_a[i] = new TH2F(histName, histTitle, 400, 0, 16000, 400, 0, 20000);
  // }
  // for (int i = 0; i < 48; i++)
  // {
  //   TString histName = Form("hCathode_%d", i);
  //   TString histTitle = Form("Cathode_E_%d;", i);
  //   hPC_E[i] = new TH1F(histName, histTitle, 3200, 0, 32000);
  // }
  sx3_contr.ConstructGeo();
  pw_contr.ConstructGeo();

  std::ifstream inputFile("slope_intercept_results.txt");

  if (inputFile.is_open())
  {
    std::string line;
    int index;
    double slope, intercept;
    while (std::getline(inputFile, line))
    {
      std::stringstream ss(line);
      ss >> index >> slope >> intercept;
      // wires 37, 39, 44 have fit data that is incorrect or not present, they have thus been set to 1,0 (slope, intercept) for convenience
      // wire 19 the 4th point was genereated using the slope of the line produced uising the other 3 points from the wire 1 vs wire 19 plot
      if (index >= 0 && index <= 47)
      {
        slopeInterceptMap[index] = std::make_pair(slope, intercept);
      }
    }
    inputFile.close();
  }
  else
  {
    std::cerr << "Error opening slope_intercept.txt" << std::endl;
  }

  std::string filename = "sx3_GainMatchback.txt";

  std::ifstream infile(filename);
  if (!infile.is_open())
  {
    std::cerr << "Error opening " << filename << "!" << std::endl;
    return;
  }

  int id, bk;
  double gain;
  while (infile >> id >> bk >> gain)
  {
    backGain[id][bk] = gain;
    if (backGain[id][bk] > 0)
      backGainValid[id][bk] = true;
    else
      backGainValid[id][bk] = false;
  }

  infile.close();
  std::cout << "Loaded back gains from " << filename << std::endl;

  std::string filename1 = "sx3_GainMatchfront.txt";

  std::ifstream infile1(filename1);
  if (!infile1.is_open())
  {
    std::cerr << "Error opening " << filename1 << "!" << std::endl;
    return;
  }

  int idf, bkf, uf, df;
  double fgain;
  while (infile1 >> idf >> bkf >> uf >> df >> fgain)
  {
    frontGain[idf][bkf][uf][df] = fgain;
    frontGainValid[idf][bkf][uf][df] = true;
  }
}

Bool_t Analyzer::Process(Long64_t entry)
{

  // if ( entry > 100 ) return kTRUE;

  hitPos.Clear();
  HitNonZero = false;

  // if( entry > 1) return kTRUE;
  // printf("################### ev : %llu \n", entry);

  b_sx3Multi->GetEntry(entry);
  b_sx3ID->GetEntry(entry);
  b_sx3Ch->GetEntry(entry);
  b_sx3E->GetEntry(entry);
  b_sx3T->GetEntry(entry);
  b_qqqMulti->GetEntry(entry);
  b_qqqID->GetEntry(entry);
  b_qqqCh->GetEntry(entry);
  b_qqqE->GetEntry(entry);
  b_qqqT->GetEntry(entry);
  b_pcMulti->GetEntry(entry);
  b_pcID->GetEntry(entry);
  b_pcCh->GetEntry(entry);
  b_pcE->GetEntry(entry);
  b_pcT->GetEntry(entry);

  sx3.CalIndex();
  qqq.CalIndex();
  pc.CalIndex();

  // sx3.Print();

  // ########################################################### Raw data
  //  //======================= SX3
  sx3ecut = false;
  std::vector<std::pair<int, int>> ID; // first = id, 2nd = index
  for (int i = 0; i < sx3.multi; i++)
  {
    ID.push_back(std::pair<int, int>(sx3.id[i], i));
    hsx3IndexVE->Fill(sx3.index[i], sx3.e[i]);

    if (sx3.e[i] > 100)
    {
      sx3ecut = true;
    }

    for (int j = i + 1; j < sx3.multi; j++)
    {
      hsx3Coin->Fill(sx3.index[i], sx3.index[j]);
    }

    for (int j = 0; j < pc.multi; j++)
    {
      if (pc.index[j] < 24 && pc.e[j] > 100)
      {
        hsx3VpcIndex->Fill(sx3.index[i], pc.index[j]);
        //  if( sx3.ch[index] > 8 ){
        //    hsx3VpcE->Fill( sx3.e[i], pc.e[j] );
        //   }
      }
    }
  }

  if (ID.size() > 0)
  {
    std::sort(ID.begin(), ID.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b)
              { return a.first < b.first; });
    // printf("##############################\n");
    // for( size_t i = 0; i < ID.size(); i++) printf("%zu | %d %d \n", i, ID[i].first, ID[i].second );

    std::vector<std::pair<int, int>> sx3ID;
    sx3ID.push_back(ID[0]);
    bool found = false;
    for (size_t i = 1; i < ID.size(); i++)
    {
      if (ID[i].first == sx3ID.back().first)
      {
        sx3ID.push_back(ID[i]);
        if (sx3ID.size() >= 3)
        {
          found = true;
        }
      }
      else
      {
        if (!found)
        {
          sx3ID.clear();
          sx3ID.push_back(ID[i]);
        }
      }
    }

    // printf("---------- sx3ID Multi : %zu \n", sx3ID.size());

    if (found)
    {
      int sx3ChUp, sx3ChDn, sx3ChBk;
      float sx3EUp, sx3EDn;
      // printf("------ sx3 ID : %d, multi: %zu\n", sx3ID[0].first, sx3ID.size());
      for (size_t i = 0; i < sx3ID.size(); i++)
      {
        int index = sx3ID[i].second;
        // printf(" %zu | index %d | ch : %d, energy : %d \n", i, index, sx3.ch[index], sx3.e[index]);

        if (sx3.ch[index] < 8)
        {
          if (sx3.ch[index] % 2 == 0)
          {
            sx3ChDn = sx3.ch[index];
            sx3EDn = sx3.e[index];
          }
          else
          {
            sx3ChUp = sx3.ch[index];
            sx3EUp = sx3.e[index];
          }
        }
        else
        {
          sx3ChBk = sx3.ch[index];
        }
        for (int j = 0; j < pc.multi; j++)
        {
          // hsx3VpcIndex->Fill( sx3.index[i], pc.index[j] );
          if (sx3.ch[index] > 8)
          {
            hsx3VpcE->Fill(sx3.e[i], pc.e[j]);
            //  hpcIndexVE->Fill( pc.index[i], pc.e[i] );
          }
        }
      }

      sx3_contr.CalSX3Pos(sx3ID[0].first, sx3ChUp, sx3ChDn, sx3ChBk, sx3EUp, sx3EDn);
      hitPos = sx3_contr.GetHitPos();
      HitNonZero = true;
      // hitPos.Print();
    }
  }

  // //======================= QQQ

  qqqEcut = false;
  for (int i = 0; i < qqq.multi; i++)
  {
    // for( int j = 0; j < pc.multi; j++){
    // if(pc.index[j]==4){
    hqqqIndexVE->Fill(qqq.index[i], qqq.e[i]);
    // }
    // printf("QQQ ID : %d, ch : %d, e : %d \n", qqq.id[i], qqq.ch[i], qqq.e[i]);
    if (qqq.e[i] > 100)
    {
      qqqEcut = true;
    }
    // }
    for (int j = 0; j < qqq.multi; j++)
    {
      if (j == i)
        continue;
      hqqqCoin->Fill(qqq.index[i], qqq.index[j]);
    }

    for (int k = 0; k < pc.multi; k++)
    {
      if (pc.index[k] < 24 && pc.e[k] > 50)
      {
        hqqqVpcE->Fill(qqq.e[i], pc.e[k]);
        //  hpcIndexVE->Fill( pc.index[i], pc.e[i] );
        hqqqVpcIndex->Fill(qqq.index[i], pc.index[k]);
      }
    }

    // }

    for (int j = i + 1; j < qqq.multi; j++)
    {
      // if( qqq.used[i] == true ) continue;

      // if( qqq.id[i] == qqq.id[j] && (16 - qqq.ch[i]) * (16 - qqq.ch[j]) < 0  ){ // must be same detector and wedge and ring
      if (qqq.id[i] == qqq.id[j])
      { // must be same detector

        int chWedge = -1;
        int chRing = -1;
        if (qqq.ch[i] < qqq.ch[j])
        {
          chRing = qqq.ch[j] - 16;
          chWedge = qqq.ch[i];
        }
        else
        {
          chRing = qqq.ch[i];
          chWedge = qqq.ch[j] - 16;
        }
        // printf(" ID : %d , chWedge : %d, chRing : %d \n", qqq.id[i], chWedge, chRing);

        double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
        double rho = 50. + 40. / 16. * (chRing + 0.5);
        // if(qqq.e[i]>50){
        hqqqPolar->Fill(theta, rho);
        // }
        // qqq.used[i] = true;
        // qqq.used[j] = true;

        if (!HitNonZero)
        {
          double x = rho * TMath::Cos(theta);
          double y = rho * TMath::Sin(theta);
          hitPos.SetXYZ(x, y, 23 + 75 + 30);
          HitNonZero = true;
        }
      }
    }
  }
  // //======================= PC

  // Calculate the crossover points and put them into an array

  pwinstance.ConstructGeo();
  Coord Crossover[24][24][2];
  TVector3 a, c, diff;
  double a2, ac, c2, adiff, cdiff, denom, alpha, beta;
  int index = 0;
  for (int i = 0; i < pwinstance.An.size(); i++)
  {
    a = pwinstance.An[i].first - pwinstance.An[i].second;

    for (int j = 0; j < pwinstance.Ca.size(); j++)
    {
      // Ok so this method uses what is essentially the solution of 2 equations to find the point of intersection between the anode and cathode wires
      // here a and c are the vectors of the anode and cathode wires respectively
      // diff is the perpendicular vector between the anode and cathode wires
      // The idea behind this is to then find the scalars alpha and beta that give a ratio between 0 and -1,

      c = pwinstance.Ca[j].first - pwinstance.Ca[j].second;
      diff = pwinstance.An[i].first - pwinstance.Ca[j].first;
      a2 = a.Dot(a);
      c2 = c.Dot(c);
      ac = a.Dot(c);
      adiff = a.Dot(diff);
      cdiff = c.Dot(diff);
      denom = a2 * c2 - ac * ac;
      alpha = (ac * cdiff - c2 * adiff) / denom;
      beta = (a2 * cdiff - ac * adiff) / denom;

      Crossover[i][j][0].x = pwinstance.An[i].first.X() + alpha * a.X();
      Crossover[i][j][0].y = pwinstance.An[i].first.Y() + alpha * a.Y();
      Crossover[i][j][0].z = pwinstance.An[i].first.Z() + alpha * a.Z();
      if (Crossover[i][j][0].z < -190 || Crossover[i][j][0].z > 190)
      {
        Crossover[i][j][0].z = 9999999;
      }
      // placeholder variable Crossover[i][j][1].x has nothing to do with the geometry of the crossover and is being used to store the alpha value-
      //-so that it can be used to sort "good" hits later
      Crossover[i][j][1].x = alpha;
      Crossover[i][j][1].y = 0;
      // if(i==0){
      // printf("CID, Crossover z and alpha are :  %d %f %f \n",  j, Crossover[i][j][0].z, Crossover[i][j][1].x /*this is alpha*/);
      // }
      //   }
      // }
    }
  }
  // printf("Anode and cathode indices, alpha, denom, andiff, cndiff : %d %d %f %f %f %f\n", i, j, alpha, denom, adiff, cdiff);

  // anodeIntersection.Clear();
  for (int i = 0; i < pc.multi; i++)
  {

    if (pc.e[i] > 100)
    {
      hpcIndexVE->Fill(pc.index[i], pc.e[i]); // non gain matched energy
    }

    // Gain Matching of PC wires
    if (pc.index[i] >= 0 && pc.index[i] < 48)
    {
      // printf("index: %d, Old cathode energy: %d \n", pc.index[i],pc.e[i]);
      auto it = slopeInterceptMap.find(pc.index[i]);
      if (it != slopeInterceptMap.end())
      {
        double slope = it->second.first;
        double intercept = it->second.second;
        // printf("slope: %f, intercept:%f\n" ,slope, intercept);
        pc.e[i] = slope * pc.e[i] + intercept;
        // printf("index: %d, New cathode energy: %d \n",pc.index[i], pc.e[i]);
      }
      hpcIndexVE_GM->Fill(pc.index[i], pc.e[i]);
      // hPC_E[pc.index[i]]->Fill(pc.e[i]); // gain matched energy per channel
    }
  }

  std::vector<std::pair<int, double>> anodeHits = {};
  std::vector<std::pair<int, double>> cathodeHits = {};
  std::vector<std::pair<int, double>> corrcatMax = {};
  std::vector<std::pair<int, double>> corrcatnextMax = {};
  std::vector<std::pair<int, double>> commcat = {};
  int aID = 0;
  int cID = 0;
  float aE = 0;
  float cE = 0;
  float aESum = 0;
  float cESum = 0;
  float aEMax = 0;
  float cEMax = 0;
  float aEnextMax = 0;
  float cEnextMax = 0;
  int aIDMax = 0;
  int cIDMax = 0;
  int aIDnextMax = 0;
  int cIDnextMax = 0;

  // Define the excluded SX3 and QQQ channels
  // std::unordered_set<int> excludeSX3 = {34, 35, 36, 37, 61, 62, 67, 73, 74, 75, 76, 77, 78, 79, 80, 93, 97, 100, 103, 108, 109, 110, 111, 112};
  // std::unordered_set<int> excludeQQQ = {0, 17, 109, 110, 111, 112, 113, 119, 127, 128};
  // inCuth=false;
  // inCutl=false;
  // inPCCut=false;
  for (int i = 0; i < pc.multi; i++)
  {
    if (pc.e[i] > 100 /*&& pc.multi < 7*/)
    {
      // creating a vector of pairs of anode and cathode hits
      if (pc.index[i] < 24)
      {
        anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
      }
      else if (pc.index[i] >= 24)
      {
        cathodeHits.push_back(std::pair<int, double>(pc.index[i] - 24, pc.e[i]));
      }

      for (int j = i + 1; j < pc.multi; j++)
      {
        // if(PCCoinc_cut1->IsInside(pc.index[i], pc.index[j]) || PCCoinc_cut2->IsInside(pc.index[i], pc.index[j])){
        //   // hpcCoin->Fill(pc.index[i], pc.index[j]);
        //   inPCCut = true;
        // }
        hpcCoin->Fill(pc.index[i], pc.index[j]);
      }
    }
  }
  // sorting the anode and cathode hits in descending order of energy
  std::sort(anodeHits.begin(), anodeHits.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
            { return a.second > b.second; });
  std::sort(cathodeHits.begin(), cathodeHits.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
            { return a.second > b.second; });

  bool SiPCflag;

  corrcatMax.clear();
  if (anodeHits.size() >= 1 && cathodeHits.size() > 1)
  {
    if (((TMath::TanH(hitPos.Y() / hitPos.X())) > (TMath::TanH(a.Y() / a.X()) - TMath::PiOver4())) || ((TMath::TanH(hitPos.Y() / hitPos.X())) < (TMath::TanH(a.Y() / a.X()) + TMath::PiOver4())))
    {

      for (const auto &anode : anodeHits)
      {
        aID = anode.first;
        aE = anode.second;
        aESum += aE;
        if (aE > aEMax)
        {
          aEMax = aE;
          aIDMax = aID;
        }
        if (aE > aEnextMax && aE < aEMax)
        {
          aEnextMax = aE;
          aIDnextMax = aID;
        }
        // for(const auto &cat : cathodeHits){
        //   hAVCcoin->Fill(aID, cat.first);
        // }
      }

      // std::cout << " Anode iD : " << aIDMax << " Energy : " << aEMax << std::endl;

      // printf("aID : %d, aE : %f, cE : %f\n", aID, aE, cE);
      for (const auto &cathode : cathodeHits)
      {
        cID = cathode.first;
        cE = cathode.second;
        // std::cout << "Cathode ID : " << cID << " Energy : " << cE << std::endl;

        hAVCcoin->Fill(aIDMax, cID);

        // This section of code is used to find the cathodes are correlated with the max and next max anodes, as well as to figure out if there are any common cathodes
        // the anodes are correlated with the cathodes +/-3 from the anode number in the reverse order

        for (int j = -4; j < 3; j++)
        {
          if ((aIDMax + 24 + j) % 24 == 23 - cID)
          /* the 23-cID is used to accomodate for the fact that the order of the cathodes was reversed relative top the physical geometry */
          // if (Crossover[aIDMax][cID][0].z != 9999999)
          {
            corrcatMax.push_back(std::pair<int, double>(cID, cE));
            cESum += cE;
            // printf("Max Anode : %d Correlated Cathode : %d Anode Energy : %f z value : %f \n", aIDMax, cID, cESum, Crossover[aIDMax][cID][1].z /*prints alpha*/);
            // std::cout << " Cathode iD : " << cID << " Energy : " << cE << std::endl;
          }
        }
      }
    }
  }

  TVector3 anodeIntersection;
  anodeIntersection.Clear();
  // Implementing a method for PC reconstruction using a single Anode event
  // if (anodeHits.size() == 1)
  {
    float x, y, z = 0;
    for (const auto &corr : corrcatMax)
    {
      if (cESum > 0)
      {
        x += (corr.second) / cESum * Crossover[aIDMax][corr.first][0].x;
        y += (corr.second) / cESum * Crossover[aIDMax][corr.first][0].y;
        z += (corr.second) / cESum * Crossover[aIDMax][corr.first][0].z;
        // printf("Max Anode : %d Correlated Cathode : %d cathode Energy : %f cESum Energy : %f z value : %f \n", aIDMax, corr.first, corr.second, cESum, Crossover[aIDMax][corr.first][1].z /*prints alpha*/);
      }
      else
      {
        printf("Warning: No valid cathode hits to correlate with anode %d! \n", aIDMax);
      }
      // printf("EventID : %llu, Max Anode : %d Cathode: %d PC X and Y : (%f, %f) \n", entry, aIDMax, cID, Crossover[aIDMax][cID][0].x, Crossover[aIDMax][cID][0].y);
      // for (int i = 0; i < sx3.multi; i++)
      // {
      //   printf("EventID : %llu, HitPos X, Y, Z:  %f %f %f SX3ID : %d %d \n", entry, hitPos.X(), hitPos.Y(), hitPos.Z(), sx3.id[i], sx3.ch[i]);
      // }

      // for (int i = 0; i < qqq.multi; i++)
      // {
      //   printf("Max Anode : %d Cathode: %d PC X and Y : %f %f \n", aIDMax, cID, Crossover[aIDMax][cID][0].x, Crossover[aIDMax][cID][0].y);
      //   printf("HitPos X, Y, Z, QQQID : %f %f %f %d \n", hitPos.X(), hitPos.Y(), hitPos.Z(), qqq.id[i]);
      // }
    }
    anodeIntersection = TVector3(x, y, z);
    // std::cout << "Anode Intersection " << anodeIntersection.Z() << " " << x << " " << y << " " << z << std::endl;
  }

  if (anodeIntersection.Z() != 0)
  {
    hPCZProj->Fill(anodeIntersection.Z());
  }
  // Filling the PC Z projection histogram
  // std::cout << anodeIntersection.Z() << std::endl;
  // hPCZProj->Fill(anodeIntersection.Z());

  // }

  // inCuth = false;
  // inCutl = false;
  // inPCCut = false;
  // for(int j=i+1;j<pc.multi;j++){
  //   if(PCCoinc_cut1->IsInside(pc.index[i], pc.index[j]) || PCCoinc_cut2->IsInside(pc.index[i], pc.index[j])){
  //     // hpcCoin->Fill(pc.index[i], pc.index[j]);
  //     inPCCut = true;
  //   }
  //   hpcCoin->Fill(pc.index[i], pc.index[j]);
  // }

  // Check if the accumulated energies are within the defined ranges
  // if (AnCatSum_high && AnCatSum_high->IsInside(aESum, cESum)) {
  //     inCuth = true;
  // }
  // if (AnCatSum_low && AnCatSum_low->IsInside(aESum, cESum)) {
  //     inCutl = true;
  // }

  // Fill histograms based on the cut conditions
  // if (inCuth && inPCCut) {
  //     hanVScatsum_hcut->Fill(aESum, cESum);
  // }
  // if (inCutl && inPCCut) {
  //     hanVScatsum_lcut->Fill(aESum, cESum);
  // }
  // for(auto anode : anodeHits){

  // float aE = anode.second;
  // aESum += aE;
  // if(inPCCut){
  hanVScatsum->Fill(aEMax, cESum);
  // }

  // if (sx3ecut)
  // {
  hCat4An->Fill(corrcatMax.size());
  hNosvAe->Fill(corrcatMax.size(), aEMax);
  hAnodehits->Fill(anodeHits.size());
  // }

  // }
  if (anodeHits.size() < 1)
  {
    hCat0An->Fill(cathodeHits.size());
  }

  if (HitNonZero && anodeIntersection.Z() != 0)
  {
    pw_contr.CalTrack2(hitPos, anodeIntersection);
    hZProj->Fill(pw_contr.GetZ0());
  }

  // ########################################################### Track constrcution

  // ############################## DO THE KINEMATICS

  return kTRUE;
}

void Analyzer::Terminate()
{

  // gStyle->SetOptStat("neiou");
  // TCanvas *canvas = new TCanvas("cANASEN", "ANASEN", 2000, 2000);
  // canvas->Divide(3, 3);

  // // hsx3VpcIndex->Draw("colz");

  // //=============================================== pad-1
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // hsx3IndexVE->Draw("colz");

  // //=============================================== pad-2
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // hqqqIndexVE->Draw("colz");

  // //=============================================== pad-3
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // hpcIndexVE->Draw("colz");

  // //=============================================== pad-4
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // hsx3Coin->Draw("colz");

  // //=============================================== pad-5
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // canvas->cd(padID)->SetLogz(true);

  // hqqqCoin->Draw("colz");

  // //=============================================== pad-6
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // hpcCoin->Draw("colz");

  // //=============================================== pad-7
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // // hsx3VpcIndex ->Draw("colz");
  // hsx3VpcE->Draw("colz");

  // //=============================================== pad-8
  // padID++;
  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);

  // // hqqqVpcIndex ->Draw("colz");

  // hqqqVpcE->Draw("colz");
  // //=============================================== pad-9
  // padID++;

  // // canvas->cd(padID)->DrawFrame(-50, -50, 50, 50);
  // // hqqqPolar->Draw("same colz pol");

  // canvas->cd(padID);
  // canvas->cd(padID)->SetGrid(1);
  // //  hZProj->Draw();
  // hanVScatsum->Draw("colz");

  // // TFile *outRoot = new TFile("Histograms.root", "RECREATE");

  // // if (!outRoot->IsOpen())
  // // {
  // //   std::cerr << "Error opening file for writing!" << std::endl;
  // //   return;
  // // }

  // // // Loop through histograms and write them to the ROOT file
  // // for (int i = 0; i < 48; i++)
  // // {
  // //   if (hPC_E[i] != nullptr)
  // //   {
  // //     hPC_E[i]->Write(); // Write histogram to file
  // //   }
  // // }

  // // outRoot->Close();
}
