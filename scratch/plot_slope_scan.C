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

void plot_slope_scan()
{
    // =========================================================================
    // --- Configuration ---
    // =========================================================================

    // Set to true to inspect each Gaussian fit. Double-click the canvas to advance.
    bool checkFits = false;
    gStyle->SetOptStat(0);
    double slopes[] = {0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66};
    const int N = 14; // Updated to 14 to match the size of your slopes array

    // int runs[] = {9, 12, 18, 19, 20, 21};
    // const int NRUNS = 6;
    // int runs[] = {9, 12};
    // const int NRUNS = 2;
    int runs[] = {18, 19, 20, 21};
    const int NRUNS = 4;
    Int_t runColors[6] = {kRed + 1, kBlue + 2, kGreen + 2, kOrange + 1, kViolet + 2, kCyan + 2};

    // --- Initialize Graphs ---
    TGraph *gStdDevRun[NRUNS];
    TGraph *gMeanRun[NRUNS];
    for (int ir = 0; ir < NRUNS; ir++)
    {
        gStdDevRun[ir] = new TGraph(N);
        gMeanRun[ir] = new TGraph(N);
    }
    TGraph *gStdDevSum = new TGraph(N);
    TGraph *gMeanSum = new TGraph(N);

    TH1 *h1[NRUNS][N];

    // --- Load Data ---
    for (int ir = 0; ir < NRUNS; ir++)
    {
        for (int i = 0; i < N; i++)
        {
            h1[ir][i] = nullptr;
            TString path = TString::Format("slope_scan/run%03d/slope_%.2f.root", runs[ir], slopes[i]);
            TFile *f = TFile::Open(path);

            if (!f || f->IsZombie())
                continue;

            TH1 *t1 = (TH1 *)f->Get("pczfix-sx3pczguess_A1C2");
            if (t1)
            {
                t1->SetDirectory(nullptr);
                h1[ir][i] = t1;
            }
            f->Close();
        }
    }

    // =========================================================================
    // --- Calculate Mean & StdDev with Bounded Gaussian Fit ---
    // =========================================================================
    std::vector<int> nPtsRun(NRUNS, 0);
    int nPtsSum = 0;
    double fitWindow = 15.0; // Window around the peak to fit (mm)

    TCanvas *cDebug = nullptr;
    if (checkFits)
    {
        cDebug = new TCanvas("cDebug", "Fit Diagnostic Debugger", 0, 0, 800, 600);
    }

    for (int i = 0; i < N; i++)
    {
        TH1 *hsum = nullptr;

        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (!h1[ir][i] || h1[ir][i]->GetEntries() <= 0)
                continue;

            // 1. Find center of the highest peak
            int maxBin = h1[ir][i]->GetMaximumBin();
            double peakX = h1[ir][i]->GetXaxis()->GetBinCenter(maxBin);

            // 2. Define the fit
            TF1 *gfit = new TF1("gfit", "gaus", peakX - fitWindow, peakX + fitWindow);

            // 3. Fit Conditionally (Draw if checking, otherwise run quietly in background)
            TString fitOpt = checkFits ? "RQS" : "RQS0";

            if (checkFits)
                cDebug->cd(); // Ensure we draw on the debug canvas
            TFitResultPtr fitRes = h1[ir][i]->Fit(gfit, fitOpt);

            // --- DIAGNOSTIC PAUSE ---
            if (checkFits)
            {
                h1[ir][i]->SetTitle(TString::Format("Run %d | Slope %.2f", runs[ir], slopes[i]));
                h1[ir][i]->GetXaxis()->SetRangeUser(peakX - 50, peakX + 50); // Zoom in X

                // Scale Y-axis to the maximum of this specific zoomed window + 20% headroom
                h1[ir][i]->SetMaximum(h1[ir][i]->GetMaximum() * 1.2);

                h1[ir][i]->Draw();
                cDebug->Modified();
                cDebug->Update();
                printf("Visual Check: Run %d | Slope %.2f. Double-click canvas to continue...\n", runs[ir], slopes[i]);
                while (cDebug->WaitPrimitive())
                    ;
            }

            double mean, stddev;
            if (fitRes >= 0)
            { // If fit succeeds
                mean = gfit->GetParameter(1);
                stddev = gfit->GetParameter(2);
            }
            else
            { // Fallback to raw stats if fit fails
                mean = h1[ir][i]->GetMean();
                stddev = h1[ir][i]->GetStdDev();
            }

            gStdDevRun[ir]->SetPoint(nPtsRun[ir], slopes[i], stddev);
            gMeanRun[ir]->SetPoint(nPtsRun[ir]++, slopes[i], mean);

            delete gfit;

            // Build the cumulative sum histogram
            if (!hsum)
                hsum = (TH1 *)h1[ir][i]->Clone(TString::Format("hsum_%.2f", slopes[i]));
            else
                hsum->Add(h1[ir][i]);
        }

        // --- Process Summed Histogram ---
        if (hsum && hsum->GetEntries() > 0)
        {
            int maxBinSum = hsum->GetMaximumBin();
            double peakXSum = hsum->GetXaxis()->GetBinCenter(maxBinSum);
            TF1 *gfitSum = new TF1("gfitSum", "gaus", peakXSum - fitWindow, peakXSum + fitWindow);

            TString fitOptSum = checkFits ? "RQS" : "RQS0";
            if (checkFits)
                cDebug->cd();

            TFitResultPtr fitResSum = hsum->Fit(gfitSum, fitOptSum);

            // --- DIAGNOSTIC PAUSE FOR SUM ---
            if (checkFits)
            {
                hsum->SetTitle(TString::Format("SUMMED RUNS | Slope %.2f", slopes[i]));
                hsum->GetXaxis()->SetRangeUser(peakXSum - 50, peakXSum + 50);

                // Scale Y-axis for the sum + 20% headroom
                hsum->SetMaximum(hsum->GetMaximum() * 1.2);

                hsum->Draw();
                cDebug->Modified();
                cDebug->Update();
                printf("Visual Check: SUMMED | Slope %.2f. Double-click canvas to continue...\n", slopes[i]);
                while (cDebug->WaitPrimitive())
                    ;
            }

            // if (fitResSum >= 0)
            // {
            //     gStdDevSum->SetPoint(nPtsSum, slopes[i], gfitSum->GetParameter(2));
            //     gMeanSum->SetPoint(nPtsSum++, slopes[i], gfitSum->GetParameter(1));
            // }
            // else
            // {
            //     gStdDevSum->SetPoint(nPtsSum, slopes[i], hsum->GetStdDev());
            //     gMeanSum->SetPoint(nPtsSum++, slopes[i], hsum->GetMean());
            // }

            // delete gfitSum;
            // delete hsum;
        }
    }

    if (checkFits && cDebug)
    {
        delete cDebug; // Clean up debug canvas before plotting final results
    }

    // =========================================================================
    // ---- Canvas 1: StdDev and Mean vs Slope
    // =========================================================================
    TCanvas *c1 = new TCanvas("c_stats", "1D Residual Stats vs Slope", 0, 0, 1600, 600);
    c1->Divide(2, 1);

    c1->cd(1);
    gPad->SetGrid(1, 1);
    TLegend *leg1 = new TLegend(0.78, 0.45, 0.97, 0.88);
    leg1->SetBorderSize(1);
    TMultiGraph *mg1 = new TMultiGraph("mg_stddev", "StdDev(pczfix - sx3pczguess) vs Slope;Slope;StdDev (mm)");

    for (int ir = 0; ir < NRUNS; ir++)
    {
        if (nPtsRun[ir] == 0)
            continue;
        gStdDevRun[ir]->SetLineColor(runColors[ir]);
        gStdDevRun[ir]->SetMarkerColor(runColors[ir]);
        gStdDevRun[ir]->SetLineWidth(2);
        gStdDevRun[ir]->SetMarkerStyle(20 + ir);
        mg1->Add(gStdDevRun[ir], "PL");
        leg1->AddEntry(gStdDevRun[ir], TString::Format("run %d", runs[ir]), "lp");
    }
    if (nPtsSum > 0)
    {
        gStdDevSum->SetLineColor(kBlack);
        gStdDevSum->SetLineWidth(3);
        gStdDevSum->SetLineStyle(kDashed);
        gStdDevSum->SetMarkerStyle(29);
        mg1->Add(gStdDevSum, "PL");
        leg1->AddEntry(gStdDevSum, "sum", "lp");
    }
    mg1->Draw("A");
    leg1->Draw();

    double bestSlope = 0, bestStdDev = 1e9;
    for (int i = 0; i < gStdDevRun[0]->GetN(); i++)
    {
        double x, y;
        gStdDevRun[0]->GetPoint(i, x, y);
        if (y > 0 && y < bestStdDev)
        {
            bestStdDev = y;
            bestSlope = x;
        }
    }

    TLine *lb = new TLine(bestSlope, mg1->GetHistogram()->GetMinimum(),
                          bestSlope, mg1->GetHistogram()->GetMaximum());
    lb->SetLineColor(kRed);
    lb->SetLineStyle(kDashed);
    lb->Draw();
    printf(">>> Best slope (Run %d) = %.2f  StdDev=%.2f mm\n", runs[0], bestSlope, bestStdDev);

    c1->cd(2);
    gPad->SetGrid(1, 1);
    TLegend *leg2 = new TLegend(0.78, 0.45, 0.97, 0.88);
    leg2->SetBorderSize(1);
    TMultiGraph *mg2 = new TMultiGraph("mg_mean", "Mean(pczfix - sx3pczguess) vs Slope;Slope;Mean Offset (mm)");

    for (int ir = 0; ir < NRUNS; ir++)
    {
        if (nPtsRun[ir] == 0)
            continue;
        gMeanRun[ir]->SetLineColor(runColors[ir]);
        gMeanRun[ir]->SetMarkerColor(runColors[ir]);
        gMeanRun[ir]->SetLineWidth(2);
        gMeanRun[ir]->SetMarkerStyle(20 + ir);
        mg2->Add(gMeanRun[ir], "PL");
        leg2->AddEntry(gMeanRun[ir], TString::Format("run %d", runs[ir]), "lp");
    }
    if (nPtsSum > 0)
    {
        gMeanSum->SetLineColor(kBlack);
        gMeanSum->SetLineWidth(3);
        gMeanSum->SetLineStyle(kDashed);
        gMeanSum->SetMarkerStyle(29);
        mg2->Add(gMeanSum, "PL");
        leg2->AddEntry(gMeanSum, "sum", "lp");
    }
    mg2->Draw("A");

    TLine *lz = new TLine(slopes[0], 0, slopes[N - 1], 0);
    lz->SetLineColor(kGray + 2);
    lz->SetLineStyle(kDotted);
    lz->Draw();
    leg2->Draw();

    c1->SaveAs("slope_scan_1d_metric.png");

    // =========================================================================
    // ---- Canvas 2: Per-slope pads, overlaying the 1D Residual distributions
    // =========================================================================
    TCanvas *c2 = new TCanvas("c_perslope_1d", "Per-slope 1D Residuals: all runs overlaid", 0, 100, 1800, 1200);
    c2->Divide(4, 4);

    for (int i = 0; i < N; i++)
    {
        c2->cd(i + 1);
        gPad->SetGrid(1, 1);
        TLegend *legS = new TLegend(0.62, 0.55, 0.88, 0.88);
        legS->SetBorderSize(0);

        // --- STEP A: Pre-scan to find the highest Y value across all runs for this slope ---
        double maxY = 0;
        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (!h1[ir][i] || h1[ir][i]->GetEntries() <= 0)
                continue;

            // Set X range first so GetMaximum() doesn't accidentally grab a tail far away
            h1[ir][i]->GetXaxis()->SetRangeUser(-50, 50);
            if (h1[ir][i]->GetMaximum() > maxY)
            {
                maxY = h1[ir][i]->GetMaximum();
            }
        }

        // --- STEP B: Draw with the dynamically scaled Y-axis ---
        bool first = true;
        for (int ir = 0; ir < NRUNS; ir++)
        {
            if (!h1[ir][i] || h1[ir][i]->GetEntries() <= 0)
                continue;

            h1[ir][i]->SetLineColor(runColors[ir]);
            h1[ir][i]->SetLineWidth(2);
            h1[ir][i]->SetTitle(TString::Format("slope=%.2f;Residual (mm);Counts", slopes[i]));

            if (first)
            {
                // Set the pad's Y-axis using the global max + 15% buffer so the legend fits cleanly
                h1[ir][i]->SetMaximum(maxY * 1.15);
                h1[ir][i]->Draw("HIST");
                first = false;
            }
            else
            {
                h1[ir][i]->Draw("HIST SAME");
            }
            legS->AddEntry(h1[ir][i], TString::Format("run %d", runs[ir]), "l");
        }
        if (!first)
            legS->Draw();
    }

    c2->SaveAs("slope_scan_perslope_1d.png");
    c2->Modified();
    c2->Update();
    while (gPad->WaitPrimitive())
        ;
}