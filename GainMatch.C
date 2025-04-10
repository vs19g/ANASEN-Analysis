#define GainMatch_cxx

#include "GainMatch.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>

#include <utility>
#include <algorithm>

#include "Armory/ClassSX3.h"
#include "Armory/ClassPW.h"

#include "TVector3.h"

TH2F *hSX3FvsB;
TH2F *hQQQFVB;

int padID = 0;

SX3 sx3_contr;
PW pw_contr;
TVector3 hitPos;
bool HitNonZero;

void GainMatch::Begin(TTree * /*tree*/)
{
    TString option = GetOption();

    hSX3FvsB = new TH2F("hSX3FvsB", "SX3 Front vs Back; Front E; Back E", 400, 0, 16000, 400, 0, 16000);

    hQQQFVB = new TH2F("hQQQFVB", "number of good QQQ vs QQQ id", 400, 0, 16000, 400, 0, 16000);
    sx3_contr.ConstructGeo();
    pw_contr.ConstructGeo();
}

Bool_t GainMatch::Process(Long64_t entry)
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
    // Initialize vector to store pairs of ID and index
    std::vector<std::pair<int, int>> ID; // first = id, second = index

    // Loop through each entry in sx3.multi
    for (int i = 0; i < sx3.multi; i++)
    {
        // Store ID and index pair in ID vector
        ID.push_back(std::pair<int, int>(sx3.id[i], i));
    }

    // Check if the ID vector is not empty
    if (ID.size() > 0)
    {
        // Sort the ID vector by the first element (ID)
        std::sort(ID.begin(), ID.end(), [](const std::pair<int, int> &a, const std::pair<int, int> &b)
                  { return a.first < b.first; });

        // Create a new vector to store IDs that are the same
        std::vector<std::pair<int, int>> sx3ID;
        sx3ID.push_back(ID[0]); // Start with the first ID
        bool found = false;

        // Loop through the sorted IDs and group by same IDs
        for (size_t i = 1; i < ID.size(); i++)
        {
            if (ID[i].first == sx3ID.back().first)
            {
                // If the ID matches the last one in sx3ID, add to the group
                sx3ID.push_back(ID[i]);

                // If 3 or more IDs are grouped, set found to true
                if (sx3ID.size() >= 3)
                {
                    found = true;
                }
            }
            else
            {
                // If a new ID is encountered and no group is found, reset the group
                if (!found)
                {
                    sx3ID.clear();
                    sx3ID.push_back(ID[i]);
                }
            }
        }

        // If a group of 3 or more IDs is found, process the channels and energies
        if (found)
        {
            int sx3ChUp = -1, sx3ChDn = -1, sx3ChBk = -1;
            float sx3EUp = 0.0, sx3EDn = 0.0, sx3EBk = 0.0;

            // Loop through the grouped IDs
            for (size_t i = 0; i < sx3ID.size(); i++)
            {
                int index = sx3ID[i].second; // Get the index from the pair

                // Categorize channel and energy
                if (sx3.ch[index] < 8)
                {
                    // If channel is less than 8, assign it to up or down
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
                    // If channel is >= 8, assign it to back
                    sx3ChBk = sx3.ch[index];
                    sx3EBk = sx3.e[index];
                }
            }
            // make a histogram for Sx3 back energy vs sum of sx3 front energies
            hSX3FvsB->Fill(sx3EUp + sx3EDn, sx3EBk);
            // Further energy matching calculations can be added here...
        }
    }
    // //======================= QQQ
    TH2F *hist1 = NULL;
    // Loop through each entry in qqq.multi
    for (int i = 0; i < qqq.multi; i++)
    {
        for (int j = i + 1; j < qqq.multi; j++)
        {
            // if( qqq.used[i] == true ) continue;

            // if( qqq.id[i] == qqq.id[j] && (16 - qqq.ch[i]) * (16 - qqq.ch[j]) < 0  ){ // must be same detector and wedge and ring
            if (qqq.id[i] == qqq.id[j])
            { // must be same detector

                int chWedge = -1;
                int chRing = -1;
                float eWedge = 0.0;
                float eRing = 0.0;
                if (qqq.ch[i] < qqq.ch[j])
                {
                    chRing = qqq.ch[j] - 16;
                    eRing = qqq.e[j];
                    chWedge = qqq.ch[i];
                    eWedge = qqq.e[i];
                }
                else
                {
                    chRing = qqq.ch[i];
                    eRing = qqq.e[i];
                    chWedge = qqq.ch[j] - 16;
                    eWedge = qqq.e[j];
                }
                // printf(" ID : %d , chWedge : %d, chRing : %d \n", qqq.id[i], chWedge, chRing);
                hQQQFVB->Fill(qqq.e[i], qqq.e[j]);
                // 1. Create/get individual 2D histogram
                TString histName = Form("hQQQFVB_r%d_w%d_id%d", chRing, chWedge, qqq.id[i]);
                TH2F *hist2d = (TH2F *)gDirectory->Get(histName);
                if (!hist2d)
                {
                    hist2d = new TH2F(histName, Form("QQQ GainMatch ID%d R%d W%d;Wedge E;Ring E", qqq.id[i], chRing, chWedge),
                                      400, 0, 16000, 400, 0, 16000);
                }

                hist2d->Fill(eWedge, eRing);
            }
        }
    }
    return kTRUE;
}
void GainMatch::Terminate()
{
    TIter next(gDirectory->GetList());
    TObject *obj;

    while ((obj = next()))
    {
        if (!obj->InheritsFrom("TH2F"))
            continue;

        TH2F *hist2d = (TH2F *)obj;
        TString name = hist2d->GetName();
        if (!name.BeginsWith("hQQQFVB_r"))
            continue;

        if (hist2d->GetEntries() < 1000)
            continue;

        std::vector<double> wedge_vals;
        std::vector<double> peak_vals;

        for (int binX = 1; binX <= hist2d->GetNbinsX(); ++binX)
        {
            TH1D *projY = hist2d->ProjectionY("_py", binX, binX);
            if (projY->GetEntries() < 30)
            {
                delete projY;
                continue;
            }

            double mean = projY->GetMean();
            TF1 *fit = new TF1("fit", "gaus", mean - 100, mean + 100);
            projY->Fit(fit, "QNR");

            double wedgeE = hist2d->GetXaxis()->GetBinCenter(binX);
            double ringPeak = fit->GetParameter(1);

            wedge_vals.push_back(wedgeE);
            peak_vals.push_back(ringPeak);

            delete projY;
            delete fit;
        }

        if (wedge_vals.size() >= 5)
        {
            TGraph *gr = new TGraph(wedge_vals.size(), &wedge_vals[0], &peak_vals[0]);
            gr->SetName(name + "_fit");
            TF1 *line = new TF1("line", "pol1", 0, 16000);
            gr->Fit(line, "Q");

            double gain = line->GetParameter(1);
            double offset = line->GetParameter(0);
            printf("Gain match %s → Gain = %.4f, Offset = %.2f\n", name.Data(), gain, offset);

            TCanvas *c1 = new TCanvas(name + "_c", name + "_c");
            gr->SetTitle(Form("Gain Match: %s", name.Data()));
            gr->GetXaxis()->SetTitle("Wedge Energy");
            gr->GetYaxis()->SetTitle("Ring Energy");
            gr->Draw("AP");
            line->Draw("same");
            c1->SaveAs(Form("%s_fit.png", name.Data()));
            delete c1;
        }
    }
}
