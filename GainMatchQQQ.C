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

    double gainW[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    double gainR[MAX_DET][MAX_RING][MAX_WEDGE] = {{{0}}};
    bool gainValid[MAX_DET][MAX_RING][MAX_WEDGE] = {{{false}}};

    std::ofstream outFile("qqq_GainMatch.txt");
    if (!outFile.is_open())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }

    // Parameters for sigma-clipping
    const int MIN_POINTS = 5;
    const int MAX_ITER = 5;
    const double CLIP_SIGMA = 3.0;

    for (const auto &kv : dataPoints)
    {
        auto key = kv.first;
        auto [id, ring, wedge] = key;
        const auto &pts_in = kv.second;
        if (pts_in.size() < (size_t)MIN_POINTS)
            continue;

        // Make a working copy of the points for clipping
        std::vector<std::pair<double, double>> pts = pts_in;

        bool solved = false;
        double final_gW = 0.0;
        double final_gR = 0.0;

        // Iterative sigma-clipping loop
        for (int iter = 0; iter < MAX_ITER; ++iter)
        {
            // Compute sums
            double sum_w2 = 0.0;
            double sum_wr = 0.0;
            double sum_r2 = 0.0;
            for (const auto &p : pts)
            {
                double w = p.first;
                double r = p.second;
                sum_w2 += w * w;
                sum_wr += w * r;
                sum_r2 += r * r;
            }

            // Guard against degenerate cases
            if (sum_w2 <= 0.0 || sum_wr <= 0.0)
            {
                // // fallback to single-parameter linear fit (original method)
                // // Use ROOT TGraph fitting as fallback
                // if (pts.size() >= 2)
                // {
                //     std::vector<double> wE, rE;
                //     wE.reserve(pts.size());
                //     rE.reserve(pts.size());
                //     for (const auto &pr : pts)
                //     {
                //         wE.push_back(pr.first);
                //         rE.push_back(pr.second);
                //     }
                //     TGraph g(static_cast<int>(wE.size()), wE.data(), rE.data());
                //     TF1 f("f_fallback", "[0]*x", 0, 16000);
                //     g.Fit(&f, "QNR");
                //     double slope = f.GetParameter(0);
                //     if (slope > 0)
                //     {
                //         // distribute correction symmetrically:
                //         double alpha = slope; // r ≈ slope * w => alpha = slope
                //         double gW_try = std::sqrt(alpha);
                //         double gR_try = 1.0 / gW_try;
                //         final_gW = gW_try;
                //         final_gR = gR_try;
                //         solved = true;
                //     }
                // }
                break;
            }

            // alpha = sum(w*r) / sum(w^2)
            double alpha = sum_wr / sum_w2;

            if (!(alpha > 0.0) || !std::isfinite(alpha))
            {
                // // degenerate; fallback to TF1 fit as above
                // if (pts.size() >= 2)
                // {
                //     std::vector<double> wE, rE;
                //     wE.reserve(pts.size());
                //     rE.reserve(pts.size());
                //     for (const auto &pr : pts)
                //     {
                //         wE.push_back(pr.first);
                //         rE.push_back(pr.second);
                //     }
                //     TGraph g(static_cast<int>(wE.size()), wE.data(), rE.data());
                //     TF1 f("f_fallback2", "[0]*x", 0, 16000);
                //     g.Fit(&f, "QNR");
                //     double slope = f.GetParameter(0);
                //     if (slope > 0)
                //     {
                //         double gW_try = std::sqrt(slope);
                //         double gR_try = 1.0 / gW_try;
                //         final_gW = gW_try;
                //         final_gR = gR_try;
                //         solved = true;
                //     }
                // }
                break;
            }

            // distribute correction between W and R
            double gW = std::sqrt(alpha);
            double gR = 1.0 / gW;

            // compute residuals and sigma
            std::vector<double> residuals;
            residuals.reserve(pts.size());
            for (const auto &p : pts)
            {
                double w = p.first;
                double r = p.second;
                double res = gW * w - gR * r;
                residuals.push_back(res);
            }

            // mean and stddev (use population sigma here)
            double mean = 0.0;
            for (double v : residuals)
                mean += v;
            mean /= residuals.size();

            double var = 0.0;
            for (double v : residuals)
                var += (v - mean) * (v - mean);
            var /= residuals.size();
            double sigma = std::sqrt(var);

            // If sigma is NaN or zero, accept and break
            if (!std::isfinite(sigma) || sigma == 0.0)
            {
                final_gW = gW;
                final_gR = gR;
                solved = true;
                break;
            }

            // Clip > CLIP_SIGMA and build new pts
            size_t before = pts.size();
            std::vector<std::pair<double, double>> new_pts;
            new_pts.reserve(pts.size());
            for (size_t k = 0; k < pts.size(); ++k)
            {
                if (std::fabs(residuals[k] - mean) <= CLIP_SIGMA * sigma)
                    new_pts.push_back(pts[k]);
            }
            pts.swap(new_pts);
            size_t after = pts.size();

            // If no points removed or too few remain, accept current solution
            final_gW = gW;
            final_gR = gR;
            solved = true;
            if (before == after || pts.size() < (size_t)MIN_POINTS)
                break;
            // otherwise iterate again with clipped pts
        } // end iter loop

        if (!solved)
            continue;

        // sanity checks: avoid ridiculous gains
        if (!(final_gW > 0.0) || !(final_gR > 0.0) || !std::isfinite(final_gW) || !std::isfinite(final_gR))
            continue;

        // store gains
        gainW[id][ring][wedge] = final_gW;
        gainR[id][ring][wedge] = final_gR;
        gainValid[id][ring][wedge] = true;

        // write out both gains: id wedge ring gW gR
        outFile << id << " " << wedge << " " << ring << " " << final_gW << " " << final_gR << std::endl;
        printf("Gain match Det%d Ring%d Wedge%d → gW=%.6f gR=%.6f\n", id, ring, wedge, final_gW, final_gR);
    }

    outFile.close();

    std::cout << "Gain matching complete." << std::endl;

    // === Plot all gain-matched QQQ points together with a 2D histogram ===

    // Fill the combined TH2F with corrected data
    for (auto &kv : dataPoints)
    {
        int id, ring, wedge;
        std::tie(id, ring, wedge) = kv.first;
        if (!gainValid[id][ring][wedge])
            continue;
        auto &pts = kv.second;
        for (auto &pr : pts)
        {
            double corrWedge = pr.first * gainW[id][ring][wedge];
            double corrRing = pr.second * gainR[id][ring][wedge];
            // hAll->Fill(corrWedge, corrRing);
            plotter->Fill2D("hAll", 4000, 0, 16000, 4000, 0, 16000, corrWedge, corrRing);
            plotter->Fill2D(Form("hGMQQQ_id%d_ring%d_wedge%d", id, ring, wedge), 400, 0, 16000, 400, 0, 16000, corrWedge, corrRing, "GainMatched");
        }
    }

    // Optionally keep previous global histos too
    plotter->FlushToDisk();
}
