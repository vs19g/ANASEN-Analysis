#define Calibration_cxx

#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <fstream>
#include <utility>
#include <algorithm>
#include <TProfile.h>
#include <TVector3.h>
#include "Armory/ClassSX3.h"
#include "Armory/ClassPW.h"
#include "TGraphErrors.h"
#include "Calibration.h"

int padID = 0;

SX3 sx3_contr;
PW pw_contr;
PW pwinstance;
TVector3 hitPos;
// TVector3 anodeIntersection;
std::map<int, std::pair<double, double>> slopeInterceptMap;

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

TH2F *hSX3FvsB;
TH2F *hSX3FvsB_g;
TH2F *hSX3;
TH1F *hZProj;
TH2F *hsx3IndexVE;
TH2F *hsx3IndexVE_gm;
TH2F *hqqqIndexVE;
TH2F *hqqqIndexVE_gm;
TH2F *hsx3Coin;
TH2F *hqqqCoin;
TH2F *hqqqPolar;

TCutG *cut;
TCutG *cut1;

// Gain arrays

const int MAX_SX3 = 24;
const int MAX_UP = 4;
const int MAX_DOWN = 4;
const int MAX_BK = 4;
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
double backGain[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool backGainValid[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};
double frontGain[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
bool frontGainValid[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{false}}}};
double uvdslope[MAX_SX3][MAX_BK][MAX_UP][MAX_DOWN] = {{{{0}}}};
double qqqGain[MAX_QQQ][MAX_BK][MAX_UP] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_BK][MAX_UP] = {{{false}}};

void Calibration::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hSX3FvsB_g = new TH2F("hSX3FvsB_g", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);
    hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hSX3 = new TH2F("hSX3", "SX3 Front v Back; Fronts; Backs", 8, 0, 8, 4, 0, 4);
    hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);
    hsx3IndexVE = new TH2F("hsx3IndexVE", "SX3 index vs Energy; sx3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hsx3IndexVE_gm = new TH2F("hsx3IndexVE_cal", "SX3 index vs Energy (calibrated); SX3 index ; Energy", 24 * 12, 0, 24 * 12, 400, 0, 5000);
    hqqqIndexVE = new TH2F("hqqqIndexVE", "QQQ index vs Energy; QQQ index ; Energy", 4 * 2 * 16, 0, 4 * 2 * 16, 400, 0, 5000);
    hqqqIndexVE_gm = new TH2F("hqqqIndexVE_cal", "QQQ index vs Energy (calibrated); QQQ index ; Energy", 4 * 2 * 16, 0, 4 * 2 * 16, 400, 0, 5000);
    hsx3Coin = new TH2F("hsx3Coin", "SX3 Coincident", 24 * 12, 0, 24 * 12, 24 * 12, 0, 24 * 12);
    hqqqCoin = new TH2F("hqqqCoin", "QQQ Coincident", 4 * 2 * 16, 0, 4 * 2 * 16, 4 * 2 * 16, 0, 4 * 2 * 16);

    hqqqPolar = new TH2F("hqqqPolar", "QQQ Polar ID", 16 * 4, -TMath::Pi(), TMath::Pi(), 16, 10, 50);

    sx3_contr.ConstructGeo();
    pw_contr.ConstructGeo();
    // ----------------------- Load Back Gains
    {
        std::string filename = "sx3_GainMatchback.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int id, bk, u, d;
            double gain;
            while (infile >> id >> bk >> u >> d >> gain)
            {
                backGain[id][bk][u][d] = gain;
                backGainValid[id][bk][u][d] = (gain > 0);
            }
            infile.close();
            std::cout << "Loaded back gains from " << filename << std::endl;
        }
    }

    // ----------------------- Load Front Gains
    {
        std::string filename = "sx3_GainMatchfront.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int id, bk, u, d;
            double gain;
            while (infile >> id >> bk >> u >> d >> gain)
            {
                frontGain[id][bk][u][d] = gain;
                frontGainValid[id][bk][u][d] = (gain > 0);
            }
            infile.close();
            std::cout << "Loaded front gains from " << filename << std::endl;
        }
    }

    // ----------------------- Load QQQ Gains
    {
        std::string filename = "qqq_GainMatch.txt";
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Error opening " << filename << "!" << std::endl;
        }
        else
        {
            int det, ring, wedge;
            double gain;
            while (infile >> det >> ring >> wedge >> gain)
            {
                qqqGain[det][ring][wedge] = gain;
                qqqGainValid[det][ring][wedge] = (gain > 0);
            }
            infile.close();
            std::cout << "Loaded QQQ gains from " << filename << std::endl;
        }
    }

    SX3 sx3_contr;
}
Bool_t Calibration::Process(Long64_t entry)
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
    }

    // --- safe SX3 handling (replace your existing block that builds sx3ID) ---
    if (ID.size() > 0)
    {
        std::sort(ID.begin(), ID.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                  { return a.first < b.first; });

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
                    found = false;
                }
            }
        }

        if (found)
        {
            // initialize to sentinel values
            int sx3ChUp = -1;
            int sx3ChDn = -1;
            int sx3ChBk = -1;
            float sx3EUp = 0.0f;
            float sx3EDn = 0.0f;
            float sx3EBk = 0.0f;

            // collect channels/energies
            for (size_t i = 0; i < sx3ID.size(); i++)
            {
                int index = sx3ID[i].second;
                int ch = sx3.ch[index];
                float e = sx3.e[index];

                if (ch < 8) // front channels
                {
                    // you used even/odd to denote down/up — keep that convention
                    if ((ch % 2) == 0) // down
                    {
                        sx3ChDn = ch;
                        sx3EDn = e;
                    }
                    else // up
                    {
                        sx3ChUp = ch;
                        sx3EUp = e;
                    }
                }
                else // back channels (assuming back channels are 8..11 or so)
                {
                    sx3ChBk = ch; // store as raw channel number; adapt if you index bk differently
                    sx3EBk = e;   // if you want to track back energy too
                }
            }

            // Basic sanity checks before using indices:
            bool haveFrontPair = (sx3ChUp >= 0 && sx3ChDn >= 0);
            bool haveBack = (sx3ChBk >= 0);

            // convert raw channel numbers to array indices if needed:
            int bk_index = (haveBack ? (sx3ChBk - 8) : -1);
            int up_index = (haveFrontPair ? sx3ChUp : -1);
            int dn_index = (haveFrontPair ? sx3ChDn : -1);
            auto sx3Id = sx3ID[0].first;

            double calibEUp, calibEDn, calibEBack = 0.0;

            if (haveFrontPair && haveBack)
            {
                // If you stored front gains indexed by [id][bk][up][down]
                if (frontGainValid[sx3Id][bk_index][up_index][dn_index])
                {
                    calibEUp = frontGain[sx3Id][bk_index][up_index][dn_index] * sx3EUp;
                    // calibEDn = frontGain[sx3Id][bk_index][up_index][dn_index] * sx3EDn;
                }
                if (backGainValid[sx3Id][bk_index][up_index][dn_index])
                {
                    calibEBack = backGain[sx3Id][bk_index][up_index][dn_index] * sx3EBk;
                }
            }

            // Only call CalSX3Pos if we have reasonable energies (avoid calling with zeros/uninitialized)
            if (haveFrontPair && (calibEUp > 50.0) && haveBack && (calibEBack > 50.0))
            {
                // find exact back energy value from sx3 entries if you tracked it above
                float backEnergyRaw = 0.0f;
                // locate the back index in sx3ID if needed
                for (size_t k = 0; k < sx3ID.size(); ++k)
                {
                    int idx = sx3ID[k].second;
                    if (sx3.ch[idx] >= 8)
                    {
                        backEnergyRaw = sx3.e[idx];
                        break;
                    }
                }

                hsx3IndexVE_gm->Fill(sx3.index[sx3ID[0].second], calibEUp);
                hSX3->Fill(sx3ChDn + 4, sx3ChBk);
                hSX3->Fill(sx3ChUp, sx3ChBk);

                // Fill the histogram for the front vs back
                hSX3FvsB->Fill(sx3EUp + sx3EDn, calibEBack);

                sx3_contr.CalSX3Pos(sx3Id, sx3ChUp, sx3ChDn, sx3ChBk, static_cast<float>(calibEUp), static_cast<float>(calibEDn));
                hitPos = sx3_contr.GetHitPos();
                HitNonZero = true;
            }
        } // found
    }

    // //======================= QQQ
    for (int i = 0; i < qqq.multi; i++)
    {

        int det = qqq.id[i]; // detector ID (0–3)
        int ch = qqq.ch[i];  // raw channel (0–31)

        // Separate ring vs wedge channel
        int ring = -1;
        int wedge = -1;
        if (ch < 16)
        { // wedge
            wedge = ch;
        }
        else
        { // ring
            ring = ch - 16;
        }

        double Ecal = qqq.e[i]; // default = raw
        if (ring >= 0 && wedge >= 0 && qqqGainValid[det][ring][wedge])
        {
            Ecal *= qqqGain[det][ring][wedge];
        }
        // for( int j = 0; j < pc.multi; j++){
        // if(pc.index[j]==4){
        hqqqIndexVE_gm->Fill(qqq.index[i], Ecal);
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

    return kTRUE;
}
void Calibration::Terminate()
{
}