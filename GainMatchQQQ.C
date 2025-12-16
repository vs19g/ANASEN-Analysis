#define GainMatchQQQ_cxx

#include "GainMatchQQQ.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCutG.h>
#include <fstream>
#include <utility>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "Armory/HistPlotter.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <cmath>

TH2F *hQQQFVB;
HistPlotter *plotter;

int padID = 0;

TCutG *cut;
std::map<std::tuple<int, int, int>, std::vector<std::pair<double, double>>> dataPoints;

void GainMatchQQQ::Begin(TTree * /*tree*/)
{
    plotter = new HistPlotter("GainQQQ.root", "TFILE");
    TString option = GetOption();

    hQQQFVB = new TH2F("hQQQFVB", "QQQ Front vs Back; Front E; Back E", 800, 0, 16000, 800, 0, 16000);

    // Load the TCutG object
    TFile *cutFile = TFile::Open("qqqcorr.root");
    if (!cutFile || cutFile->IsZombie())
    {
        std::cerr << "Error: Could not open qqqcorr.root" << std::endl;
        return;
    }
    cut = dynamic_cast<TCutG *>(cutFile->Get("qqqcorr"));
    if (!cut)
    {
        std::cerr << "Error: Could not find TCutG named 'qqqcorr' in qqqcorr.root" << std::endl;
        return;
    }
    cut->SetName("qqqcorr"); // Ensure the cut has the correct name
}

Bool_t GainMatchQQQ::Process(Long64_t entry)
{

    int ringMults[16] = {0};
    int wedgeMults[16] = {0};
    std::vector<std::tuple<int, int, int, double, double>> events;

    b_qqqMulti->GetEntry(entry);
    b_qqqID->GetEntry(entry);
    b_qqqCh->GetEntry(entry);
    b_qqqE->GetEntry(entry);
    b_qqqT->GetEntry(entry);

    qqq.CalIndex();

    for (int i = 0; i < qqq.multi; i++)
    {
        for (int j = i + 1; j < qqq.multi; j++)
        {
            if (qqq.id[i] == qqq.id[j])
            {
                int chWedge = -1;
                int chRing = -1;
                float eWedge = 0.0;
                float eRing = 0.0;
                if (qqq.ch[i] < 16 && qqq.ch[j] >= 16)
                {
                    chWedge = qqq.ch[i];
                    eWedge = qqq.e[i];
                    chRing = qqq.ch[j] - 16;
                    eRing = qqq.e[j];
                }
                else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16)
                {
                    chWedge = qqq.ch[j];
                    eWedge = qqq.e[j];
                    chRing = qqq.ch[i] - 16;
                    eRing = qqq.e[i];
                }
                else
                    continue;
                ringMults[chRing]++;
                wedgeMults[chWedge]++;
                hQQQFVB->Fill(eWedge, eRing);
                events.emplace_back(qqq.id[i], chRing, chWedge, eRing, eWedge);
                plotter->Fill2D(Form("hRaw_qqq%d_ring%d_wedge%d", qqq.id[i], chRing, chWedge), 800, 0, 3000, 800, 0, 3000, eWedge, eRing, "hRawQQQ");
                // double ratio = (eWedge > 0.0) ? (eRing / eWedge) : 0.0;
                // double maxslope = 1.5;

                // bool validPoint = false;
                // if (ratio < maxslope && ratio > 1. / maxslope)
                // {
                //     // Accumulate data for gain matching
                //     dataPoints[{qqq.id[i], chRing, chWedge}].emplace_back(eWedge, eRing);
                //     plotter->Fill2D("hAll_in", 4000, 0, 16000, 4000, 0, 16000, eWedge, eRing);
                //     validPoint = true;
                // }

                // if (!validPoint)
                // {
                //     plotter->Fill2D("hAll_out", 4000, 0, 16000, 4000, 0, 16000, eWedge, eRing);
                // }
            }
        }
    }

    for (auto tuple : events)
    {
        auto [id, chr, chw, er, ew] = tuple;
        if (ringMults[chr] > 1 || wedgeMults[chw] > 1)
            continue; // ignore multiplicity > 1 events
        double ratio = (ew > 0.0) ? (er / ew) : 0.0;
        double maxslope = 1.5;

        bool validPoint = false;
        if (ratio < maxslope && ratio > 1. / maxslope)
        {
                // Accumulate data for gain matching
                dataPoints[{id, chr, chw}].emplace_back(ew, er);
                plotter->Fill2D("hAll_in", 4000, 0, 16000, 4000, 0, 16000, ew, er);
                validPoint = true;
                    }
        if (!validPoint)
        {
            plotter->Fill2D("hAll_out", 4000, 0, 16000, 4000, 0, 16000, ew, er);
        }
    }

    return kTRUE;
}


void GainMatchQQQ::Terminate()
{
    const int MAX_DET = 4;
    const int MAX_RING = 16;
    const int MAX_WEDGE = 16;

    // We store gains locally just for the "corrected" plot, 
    // but the file will output Slopes for the global minimizer.
    double gainW[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    double gainR[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    bool gainValid[MAX_DET][MAX_RING][MAX_WEDGE] = {{{false}}};

    // Output file for the Minimizer
    std::ofstream outFile("qqq_GainMatch.txt");
    
    // Benchmark/Debug file
    std::ofstream benchFile("benchmark_diff.txt"); 
    benchFile << "ID Wedge Ring Chi2NDF Slope SlopeErr" << std::endl;

    if (!outFile.is_open()) { std::cerr << "Error opening output file!" << std::endl; return; }

    const int MIN_POINTS = 50;
    const int MAX_ITER = 3;        // Outlier rejection passes
    const double CLIP_SIGMA = 2.5; // Sigma threshold for outliers

    for (const auto &kv : dataPoints)
    {
        auto key = kv.first;
        auto [id, ring, wedge] = key;
        const auto &pts = kv.second;

        if (pts.size() < (size_t)MIN_POINTS) continue;

        std::vector<std::pair<double, double>> current_pts = pts;
        
        double finalSlope = 0.0;
        double finalSlopeErr = 0.0;
        bool converged = false;

        // --- Iterative Fitting ---
        for (int iter = 0; iter < MAX_ITER; ++iter)
        {
            if (current_pts.size() < (size_t)MIN_POINTS) break;

            std::vector<double> x, y, ex, ey;

            for (const auto &p : current_pts)
            {
                x.push_back(p.first);  // Wedge E
                y.push_back(p.second); // Ring E
                ex.push_back(std::sqrt(std::abs(p.first))); // Error in X (Poisson)
                ey.push_back(std::sqrt(std::abs(p.second))); // Error in Y (Poisson)
                
                // Sanity check to avoid 0 error
                if(ex.back() < 1.0) ex.back() = 1.0;
                if(ey.back() < 1.0) ey.back() = 1.0;
            }

            // 2. Create Graph
            TGraphErrors *gr = new TGraphErrors(current_pts.size(), x.data(), y.data(), ex.data(), ey.data());
            
            // 3. Fit Linear Function through Origin
            TF1 *f1= new TF1("calibFit", "[0]*x", 0, 16000);
            f1->SetParameter(0, 1.0);   

            // "Q"=Quiet, "N"=NoDraw, "S"=ResultPtr
            // We do NOT use "W" (Ignore weights), we want to use the errors we set.
            int fitStatus = gr->Fit(f1, "QNS");

            if (fitStatus != 0) {
                delete gr; delete f1;
                break; 
            }

            finalSlope = f1->GetParameter(0);
            double chi2 = f1->GetChisquare();
            double ndf = f1->GetNDF();
            
            // Get the statistical error on the slope
            double rawErr = f1->GetParError(0);
            
            // SCALING ERROR:
            // If Chi2/NDF > 1, the data scatters more than Poisson stats predict.
            // // We inflate the error by sqrt(Chi2/NDF) to be conservative for the Minimizer.
            // double redChi2 = (ndf > 0) ? (chi2 / ndf) : 1.0;
            // double inflation = (redChi2 > 1.0) ? std::sqrt(redChi2) : 1.0;
            
            // finalSlopeErr = rawErr * inflation;

            // 4. Outlier Rejection
            if (iter == MAX_ITER - 1) {
                converged = true;
                delete gr; delete f1;
                break;
            }

            // Calculate Residuals
            std::vector<double> residuals;
            double sumSqResid = 0.0;
            for(size_t k=0; k<current_pts.size(); ++k) {
                double val = f1->Eval(current_pts[k].first);
                double res = current_pts[k].second - val;
                residuals.push_back(res);
                sumSqResid += res*res;
            }
        //     double sigma = std::sqrt(sumSqResid / current_pts.size());
            
        //     // Filter
        //     std::vector<std::pair<double, double>> next_pts;
        //     for(size_t k=0; k<current_pts.size(); ++k) {
        //         if(std::abs(residuals[k]) < CLIP_SIGMA * sigma) {
        //             next_pts.push_back(current_pts[k]);
        //         }
        //     }

        //     if (next_pts.size() == current_pts.size()) {
        //         converged = true; 
        //         delete gr; delete f1;
        //         break;
        //     }
        //     current_pts = next_pts;
        //     delete gr; delete f1;
        }

        if (!converged || finalSlope <= 0) continue;

        // --- Store/Output ---
        
        // 1. Save locally for the verification plot (hAll)
        // Approximate local gain for plotting purposes only
        double gW_local = std::sqrt(finalSlope);
        double gR_local = 1.0 / gW_local;
        gainW[id][ring][wedge] = gW_local;
        gainR[id][ring][wedge] = gR_local;
        gainValid[id][ring][wedge] = true;

        // 2. Write to File for Minimizer
        // Format: ID Wedge Ring Slope Error
        outFile << id << " " << wedge << " " << ring << " " << finalSlope << " " << finalSlopeErr << std::endl;

        // 3. Benchmark Info
        benchFile << id << " " << wedge << " " << ring << " " 
                  << finalSlope << " " << finalSlopeErr << std::endl;
    }

    outFile.close();
    benchFile.close();
    std::cout << "Gain matching with Errors complete." << std::endl;

    // Plotting the corrected data (Visual check using local approx gains)
    for (auto &kv : dataPoints)
    {
        int id, ring, wedge;
        std::tie(id, ring, wedge) = kv.first;
        if (!gainValid[id][ring][wedge]) continue;
        auto &pts = kv.second;
        for (auto &pr : pts)
        {
            double corrWedge = pr.first * gainW[id][ring][wedge];
            double corrRing = pr.second * gainR[id][ring][wedge];
            plotter->Fill2D("hAll", 4000, 0, 16000, 4000, 0, 16000, corrWedge, corrRing);
        }
    }

    plotter->FlushToDisk();
}