#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TApplication.h>

void plot_slope_scan()
{
    // =========================================================================
    // --- Configuration ---
    // =========================================================================
    bool checkFits = false;
    gStyle->SetOptStat(0);

    // double slopes[] = {
    //     0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60,
    //     0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80,
    //     0.82, 0.84, 0.86, 0.88, 0.90};
    // const int N = 25;

    // double slopes[] = {0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64};
    // const int N = 12; // Updated to 14 to match the size of your slopes array
    const int N = 1; // Updated to 14 to match the size of your slopes array
    double slopes[] = {0.52};
    
    int runs[] = {9, 12, 18, 19, 20, 21};
    const int NRUNS = 6;
    // int runs[] = {9, 12};
    // const int NRUNS = 2;
    // int runs[] = {18, 19, 20, 21};
    // const int NRUNS = 4;
    Int_t runColors[6] = {kRed + 1, kBlue + 2, kGreen + 2, kOrange + 1, kViolet + 2, kCyan + 2};

    // --- Initialize Graphs for both SX3 and QQQ ---
    TGraph *gStdDevRun_sx3[NRUNS], *gMeanRun_sx3[NRUNS];
    TGraph *gStdDevRun_qqq[NRUNS], *gMeanRun_qqq[NRUNS];

    for (int ir = 0; ir < NRUNS; ir++)
    {
        gStdDevRun_sx3[ir] = new TGraph(N);
        gMeanRun_sx3[ir] = new TGraph(N);
        gStdDevRun_qqq[ir] = new TGraph(N);
        gMeanRun_qqq[ir] = new TGraph(N);
    }

    TGraph *gStdDevSum_sx3 = new TGraph(N), *gMeanSum_sx3 = new TGraph(N);
    TGraph *gStdDevSum_qqq = new TGraph(N), *gMeanSum_qqq = new TGraph(N);

    // Arrays to hold histograms
    TH1 *h1_sx3[NRUNS][N];
    TH1 *h1_qqq[NRUNS][N];
    TH1 *h_vtx_sx3[NRUNS][N];
    TH1 *h_vtx_qqq[NRUNS][N];

    // --- Load Data ---
    for (int ir = 0; ir < NRUNS; ir++)
    {
        for (int i = 0; i < N; i++)
        {
            h1_sx3[ir][i] = nullptr;
            h1_qqq[ir][i] = nullptr;
            h_vtx_sx3[ir][i] = nullptr;
            h_vtx_qqq[ir][i] = nullptr;

            TString path = TString::Format("slope_scan/run%03d/slope_%.2f.root", runs[ir], slopes[i]);
            TFile *f = TFile::Open(path);
            if (!f || f->IsZombie())
                continue;

            TH1 *t1_s = (TH1 *)f->Get("pczfix-sx3pczguess_A1C2");
            if (t1_s)
            {
                t1_s->SetDirectory(nullptr);
                h1_sx3[ir][i] = t1_s;
            }

            TH1 *t1_q = (TH1 *)f->Get("pczfix-qqqpczguess_A1C2");
            if (t1_q)
            {
                t1_q->SetDirectory(nullptr);
                h1_qqq[ir][i] = t1_q;
            }

            TH1 *t_sx3 = (TH1 *)f->Get("VertexRecon_pczfix_sx3");
            if (t_sx3)
            {
                t_sx3->SetDirectory(nullptr);
                h_vtx_sx3[ir][i] = t_sx3;
            }

            TH1 *t_qqq = (TH1 *)f->Get("VertexRecon_pczfix_qqq");
            if (t_qqq)
            {
                t_qqq->SetDirectory(nullptr);
                h_vtx_qqq[ir][i] = t_qqq;
            }

            f->Close();
        }
    }

    // =========================================================================
    // --- Helper Lambda for Gaussian Fitting ---
    // =========================================================================
    auto fitPeak = [](TH1 *h, double fitWindow, bool debug, TCanvas *c, double &outMean, double &outStdDev)
    {
        if (!h || h->GetEntries() <= 0)
            return false;

        int maxBin = h->GetMaximumBin();
        double peakX = h->GetXaxis()->GetBinCenter(maxBin);
        TF1 *gfit = new TF1("gfit", "gaus", peakX - fitWindow, peakX + fitWindow);
        TString fitOpt = debug ? "RQS" : "RQS0";

        if (debug && c)
            c->cd();
        TFitResultPtr fitRes = h->Fit(gfit, fitOpt);

        if (debug && c)
        {
            h->GetXaxis()->SetRangeUser(peakX - 50, peakX + 50);
            h->SetMaximum(h->GetMaximum() * 1.2);
            h->Draw();
            c->Modified();
            c->Update();
            while (c->WaitPrimitive())
                ;
        }

        if (fitRes >= 0)
        {
            outMean = gfit->GetParameter(1);
            outStdDev = gfit->GetParameter(2);
        }
        else
        {
            outMean = h->GetMean();
            outStdDev = h->GetStdDev();
        }
        delete gfit;
        return true;
    };

    // =========================================================================
    // --- Calculate Mean & StdDev for SX3 & QQQ ---
    // =========================================================================
    std::vector<int> nPtsRun_sx3(NRUNS, 0), nPtsRun_qqq(NRUNS, 0);
    int nPtsSum_sx3 = 0, nPtsSum_qqq = 0;
    double fitWin = 15.0;

    TCanvas *cDebug = checkFits ? new TCanvas("cDebug", "Debugger", 0, 0, 800, 600) : nullptr;

    for (int i = 0; i < N; i++)
    {
        TH1 *hsum_sx3 = nullptr, *hsum_qqq = nullptr;

        for (int ir = 0; ir < NRUNS; ir++)
        {
            double mean, stddev;

            // Fit SX3
            if (fitPeak(h1_sx3[ir][i], fitWin, checkFits, cDebug, mean, stddev))
            {
                gStdDevRun_sx3[ir]->SetPoint(nPtsRun_sx3[ir], slopes[i], stddev);
                gMeanRun_sx3[ir]->SetPoint(nPtsRun_sx3[ir]++, slopes[i], mean);
                if (!hsum_sx3)
                    hsum_sx3 = (TH1 *)h1_sx3[ir][i]->Clone(TString::Format("hsum_sx3_%.2f", slopes[i]));
                else
                    hsum_sx3->Add(h1_sx3[ir][i]);
            }
            if (runs[ir] == 20 || runs[ir] == 21)
                continue;
            // Fit QQQ
            if (fitPeak(h1_qqq[ir][i], fitWin, checkFits, cDebug, mean, stddev))
            {
                gStdDevRun_qqq[ir]->SetPoint(nPtsRun_qqq[ir], slopes[i], stddev);
                gMeanRun_qqq[ir]->SetPoint(nPtsRun_qqq[ir]++, slopes[i], mean);
                if (!hsum_qqq)
                    hsum_qqq = (TH1 *)h1_qqq[ir][i]->Clone(TString::Format("hsum_qqq_%.2f", slopes[i]));
                else
                    hsum_qqq->Add(h1_qqq[ir][i]);
            }
        }

        // Fit Sums
        double sMean, sStdDev;
        if (fitPeak(hsum_sx3, fitWin, checkFits, cDebug, sMean, sStdDev))
        {
            gStdDevSum_sx3->SetPoint(nPtsSum_sx3, slopes[i], sStdDev);
            gMeanSum_sx3->SetPoint(nPtsSum_sx3++, slopes[i], sMean);
        }
        if (fitPeak(hsum_qqq, fitWin, checkFits, cDebug, sMean, sStdDev))
        {
            gStdDevSum_qqq->SetPoint(nPtsSum_qqq, slopes[i], sStdDev);
            gMeanSum_qqq->SetPoint(nPtsSum_qqq++, slopes[i], sMean);
        }
        delete hsum_sx3;
        delete hsum_qqq;
    }

    if (cDebug)
        delete cDebug;

    // =========================================================================
    // ---- Canvas 1: StdDev and Mean vs Slope (Both SX3 and QQQ) ----
    // =========================================================================
    TCanvas *c1 = new TCanvas("c_stats", "1D Residual Stats vs Slope", 0, 0, 1800, 800);

    // 1. Create a Master Legend in the top 8% of the canvas
    c1->cd();
    TLegend *leg1_master = new TLegend(0.1, 0.92, 0.9, 1.0);
    leg1_master->SetNColumns(6);
    leg1_master->SetBorderSize(0);
    leg1_master->SetFillStyle(0);
    leg1_master->SetTextSize(0.03);

    for (int ir = 0; ir < NRUNS; ir++)
    {
        TLine *dSX3 = new TLine();
        dSX3->SetLineColor(runColors[ir]);
        dSX3->SetLineWidth(2);
        dSX3->SetLineStyle(kSolid);
        leg1_master->AddEntry(dSX3, TString::Format("Run %d (SX3)", runs[ir]), "l");
    }
    for (int ir = 0; ir < NRUNS; ir++)
    {
        TLine *dQQQ = new TLine();
        dQQQ->SetLineColor(runColors[ir]);
        dQQQ->SetLineWidth(2);
        dQQQ->SetLineStyle(kDashed);
        leg1_master->AddEntry(dQQQ, TString::Format("Run %d (QQQ)", runs[ir]), "l");
    }
    leg1_master->Draw();

    // 2. Create the Grid Pad for Canvas 1
    TPad *pad1_grid = new TPad("pad1_grid", "Grid", 0.0, 0.0, 1.0, 0.92);
    pad1_grid->Draw();
    pad1_grid->Divide(2, 1);

    // --- Pad 1: StdDev ---
    pad1_grid->cd(1);
    gPad->SetGrid(1, 1);
    TMultiGraph *mg1 = new TMultiGraph("mg_stddev", "StdDev of Residuals vs Slope;Slope;StdDev (mm)");

    for (int ir = 0; ir < NRUNS; ir++)
    {
        if (nPtsRun_sx3[ir] > 0)
        {
            gStdDevRun_sx3[ir]->SetLineColor(runColors[ir]);
            gStdDevRun_sx3[ir]->SetLineWidth(2);
            gStdDevRun_sx3[ir]->SetMarkerStyle(20);
            gStdDevRun_sx3[ir]->SetMarkerColor(runColors[ir]);
            mg1->Add(gStdDevRun_sx3[ir], "PL");
        }
        if (runs[ir] == 20 || runs[ir] == 21)
            continue;
        if (nPtsRun_qqq[ir] > 0)
        {
            gStdDevRun_qqq[ir]->SetLineColor(runColors[ir]);
            gStdDevRun_qqq[ir]->SetLineWidth(2);
            gStdDevRun_qqq[ir]->SetLineStyle(kDashed);
            gStdDevRun_qqq[ir]->SetMarkerStyle(24);
            gStdDevRun_qqq[ir]->SetMarkerColor(runColors[ir]);
            mg1->Add(gStdDevRun_qqq[ir], "PL");
        }
    }
    mg1->Draw("A");

    double bestSlope = 0, bestStdDev = 1e9;
    if (gStdDevRun_sx3[3]->GetN() > 0)
    {
        for (int i = 0; i < gStdDevRun_sx3[3]->GetN(); i++)
        {
            double x, y;
            gStdDevRun_sx3[3]->GetPoint(i, x, y);
            if (y > 0 && y < bestStdDev)
            {
                bestStdDev = y;
                bestSlope = x;
            }
        }
        TLine *lb = new TLine(bestSlope, mg1->GetHistogram()->GetMinimum(), bestSlope, mg1->GetHistogram()->GetMaximum());
        lb->SetLineColor(kRed);
        lb->SetLineStyle(kDashed);
        lb->Draw();
        printf(">>> Best slope (Run %d) = %.2f  StdDev=%.2f mm\n", runs[3], bestSlope, bestStdDev);
    }

    // --- Pad 2: Mean ---
    pad1_grid->cd(2);
    // gPad->SetGrid(1, 1);
    TMultiGraph *mg2 = new TMultiGraph("mg_mean", "Mean of Residuals vs Slope;Slope;Mean Offset (mm)");

    for (int ir = 0; ir < NRUNS; ir++)
    {
        if (nPtsRun_sx3[ir] > 0)
        {
            gMeanRun_sx3[ir]->SetLineColor(runColors[ir]);
            gMeanRun_sx3[ir]->SetLineWidth(2);
            gMeanRun_sx3[ir]->SetMarkerStyle(20);
            gMeanRun_sx3[ir]->SetMarkerColor(runColors[ir]);
            mg2->Add(gMeanRun_sx3[ir], "PL");
        }
        if (runs[ir] == 20 || runs[ir] == 21)
            continue;
        if (nPtsRun_qqq[ir] > 0)
        {
            gMeanRun_qqq[ir]->SetLineColor(runColors[ir]);
            gMeanRun_qqq[ir]->SetLineWidth(2);
            gMeanRun_qqq[ir]->SetLineStyle(kDashed);
            gMeanRun_qqq[ir]->SetMarkerStyle(24);
            gMeanRun_qqq[ir]->SetMarkerColor(runColors[ir]);
            mg2->Add(gMeanRun_qqq[ir], "PL");
        }
    }
    mg2->Draw("A");
    TLine *lz = new TLine(slopes[0], 0, slopes[N - 1], 0);
    lz->SetLineColor(kGray + 2);
    lz->SetLineStyle(kDotted);
    lz->Draw();

    c1->SaveAs("slope_scan_1d_metric.png");

    // =========================================================================
    // ---- Canvas 2: Per-slope pads, 1D Residuals (SX3 + QQQ) ----
    // =========================================================================
    TCanvas *c2 = new TCanvas("c_perslope_1d", "Per-slope 1D Residuals (SX3 & QQQ)", 0, 50, 2000, 1600);

    c2->cd();
    TLegend *leg2_master = new TLegend(0.1, 0.92, 0.9, 1.0);
    leg2_master->SetNColumns(6);
    leg2_master->SetBorderSize(0);
    leg2_master->SetFillStyle(0);
    leg2_master->SetTextSize(0.025);

    for (int ir = 0; ir < NRUNS; ir++)
    {
        TLine *dSX3 = new TLine();
        dSX3->SetLineColor(runColors[ir]);
        dSX3->SetLineWidth(2);
        dSX3->SetLineStyle(kSolid);
        leg2_master->AddEntry(dSX3, TString::Format("Run %d (SX3)", runs[ir]), "l");
    }
    for (int ir = 0; ir < NRUNS; ir++)
    {
        TLine *dQQQ = new TLine();
        dQQQ->SetLineColor(runColors[ir]);
        dQQQ->SetLineWidth(2);
        dQQQ->SetLineStyle(kDashed);
        leg2_master->AddEntry(dQQQ, TString::Format("Run %d (QQQ)", runs[ir]), "l");
    }
    leg2_master->Draw();

    TPad *pad2_grid = new TPad("pad2_grid", "Grid", 0.0, 0.0, 1.0, 0.92);
    pad2_grid->Draw();
    pad2_grid->Divide(3,4);

    for (int i = 0; i < N; i++)
    {
        pad2_grid->cd(i + 1);
        gPad->SetGrid(1, 1);

        double maxY = 0;
        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (h1_sx3[ir][i] && h1_sx3[ir][i]->GetEntries() > 0)
            {
                h1_sx3[ir][i]->GetXaxis()->SetRangeUser(-50, 50);
                if (h1_sx3[ir][i]->GetMaximum() > maxY)
                    maxY = h1_sx3[ir][i]->GetMaximum();
            }
            if (runs[ir] == 20 || runs[ir] == 21)
                continue;
            if (h1_qqq[ir][i] && h1_qqq[ir][i]->GetEntries() > 0)
            {
                h1_qqq[ir][i]->GetXaxis()->SetRangeUser(-50, 50);
                if (h1_qqq[ir][i]->GetMaximum() > maxY)
                    maxY = h1_qqq[ir][i]->GetMaximum();
            }
        }

        bool first = true;
        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (h1_sx3[ir][i] && h1_sx3[ir][i]->GetEntries() > 0)
            {
                h1_sx3[ir][i]->SetLineColor(runColors[ir]);
                h1_sx3[ir][i]->SetLineWidth(2);
                h1_sx3[ir][i]->SetLineStyle(kSolid);
                h1_sx3[ir][i]->SetTitle(TString::Format("slope=%.2f;Residual (mm);Counts", slopes[i]));
                if (first)
                {
                    h1_sx3[ir][i]->SetMaximum(maxY * 1.15);
                    h1_sx3[ir][i]->Draw("HIST");
                    first = false;
                }
                else
                    h1_sx3[ir][i]->Draw("HIST SAME");
            }
            if (h1_qqq[ir][i] && h1_qqq[ir][i]->GetEntries() > 0)
            {
                h1_qqq[ir][i]->SetLineColor(runColors[ir]);
                h1_qqq[ir][i]->SetLineWidth(2);
                h1_qqq[ir][i]->SetLineStyle(kDashed);
                if (first)
                {
                    h1_qqq[ir][i]->SetTitle(TString::Format("slope=%.2f;Residual (mm);Counts", slopes[i]));
                    h1_qqq[ir][i]->SetMaximum(maxY * 1.15);
                    h1_qqq[ir][i]->Draw("HIST");
                    first = false;
                }
                else
                    h1_qqq[ir][i]->Draw("HIST SAME");
            }
        }
    }
    c2->SaveAs("slope_scan_perslope_1d.png");

    // =========================================================================
    // ---- Canvas 3: Vertex Recon (SX3 & QQQ Overlaid) ----
    // =========================================================================
    TCanvas *c3 = new TCanvas("c_vtx_recon", "Vertex Recon: SX3 vs QQQ", 0, 100, 2000, 1600);

    c3->cd();
    TLegend *leg3_master = new TLegend(0.1, 0.92, 0.9, 1.0);
    leg3_master->SetNColumns(6);
    leg3_master->SetBorderSize(0);
    leg3_master->SetFillStyle(0);
    leg3_master->SetTextSize(0.025);

    for (int ir = 0; ir < NRUNS; ir++)
    {
        TLine *dSX3 = new TLine();
        dSX3->SetLineColor(runColors[ir]);
        dSX3->SetLineWidth(2);
        dSX3->SetLineStyle(kSolid);
        leg3_master->AddEntry(dSX3, TString::Format("Run %d (SX3)", runs[ir]), "l");
    }
    for (int ir = 0; ir < NRUNS; ir++)
    {
        TLine *dQQQ = new TLine();
        dQQQ->SetLineColor(runColors[ir]);
        dQQQ->SetLineWidth(2);
        dQQQ->SetLineStyle(kDashed);
        leg3_master->AddEntry(dQQQ, TString::Format("Run %d (QQQ)", runs[ir]), "l");
    }
    leg3_master->Draw();

    TPad *pad3_grid = new TPad("pad3_grid", "Grid", 0.0, 0.0, 1.0, 0.92);
    pad3_grid->Draw();
    pad3_grid->Divide(3,4);

    for (int i = 0; i < N; i++)
    {
        pad3_grid->cd(i + 1);
        gPad->SetGrid(1, 1);

        double maxY3 = 0;
        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (h_vtx_sx3[ir][i] && h_vtx_sx3[ir][i]->GetEntries() > 0)
            {
                h_vtx_sx3[ir][i]->GetXaxis()->SetRangeUser(-200, 200);
                if (h_vtx_sx3[ir][i]->GetMaximum() > maxY3)
                    maxY3 = h_vtx_sx3[ir][i]->GetMaximum();
            }
            if (h_vtx_qqq[ir][i] && h_vtx_qqq[ir][i]->GetEntries() > 0)
            {
                h_vtx_qqq[ir][i]->GetXaxis()->SetRangeUser(-200, 200);
                if (h_vtx_qqq[ir][i]->GetMaximum() > maxY3)
                    maxY3 = h_vtx_qqq[ir][i]->GetMaximum();
            }
        }

        bool first3 = true;
        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (h_vtx_sx3[ir][i] && h_vtx_sx3[ir][i]->GetEntries() > 0)
            {
                h_vtx_sx3[ir][i]->SetLineColor(runColors[ir]);
                h_vtx_sx3[ir][i]->SetLineWidth(2);
                h_vtx_sx3[ir][i]->SetLineStyle(kSolid);
                h_vtx_sx3[ir][i]->SetTitle(TString::Format("slope=%.2f;Z_{Vertex} (mm);Counts", slopes[i]));
                if (first3)
                {
                    h_vtx_sx3[ir][i]->SetMaximum(maxY3 * 1.15);
                    h_vtx_sx3[ir][i]->Draw("HIST");
                    first3 = false;
                }
                else
                    h_vtx_sx3[ir][i]->Draw("HIST SAME");
            }
            if (h_vtx_qqq[ir][i] && h_vtx_qqq[ir][i]->GetEntries() > 0)
            {
                h_vtx_qqq[ir][i]->SetLineColor(runColors[ir]);
                h_vtx_qqq[ir][i]->SetLineWidth(2);
                h_vtx_qqq[ir][i]->SetLineStyle(kDashed);
                if (first3)
                {
                    h_vtx_qqq[ir][i]->SetTitle(TString::Format("slope=%.2f;Z_{Vertex} (mm);Counts", slopes[i]));
                    h_vtx_qqq[ir][i]->SetMaximum(maxY3 * 1.15);
                    h_vtx_qqq[ir][i]->Draw("HIST");
                    first3 = false;
                }
                else
                    h_vtx_qqq[ir][i]->Draw("HIST SAME");
            }
        }
    }
    c3->SaveAs("slope_scan_vtx_recon.png");

    c1->Modified();
    c1->Update();
    c2->Modified();
    c2->Update();
    c3->Modified();
    c3->Update();
}