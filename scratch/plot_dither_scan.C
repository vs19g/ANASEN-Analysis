#include <TFile.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>

void plot_dither_scan()
{
    gStyle->SetOptStat(0);

    // =========================================================================
    // --- Configuration ---
    // =========================================================================
    double dithers[] = {0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    const int N_DITHER = sizeof(dithers) / sizeof(dithers[0]);

    int runs[] = {12, 18, 19, 20, 21};
    const int N_RUNS = sizeof(runs) / sizeof(runs[0]);

    int runColors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta};
    TString runNames[] = {"27Al Run 12", "17F Run 18", "17F Run 19", "17F Run 20", "17F Run 21"};

    // Graph arrays
    TGraphErrors *g_rms_A1C0[N_RUNS];
    TGraphErrors *g_rms_A1C1[N_RUNS];
    TGraph *g_kurt_A1C0[N_RUNS];
    TGraph *g_kurt_A1C1[N_RUNS];

    // Initialize graphs
    for (int ir = 0; ir < N_RUNS; ir++)
    {
        // RMS (Spread) Graphs
        g_rms_A1C0[ir] = new TGraphErrors(N_DITHER);
        g_rms_A1C0[ir]->SetTitle(runNames[ir]);
        g_rms_A1C0[ir]->SetMarkerStyle(20);
        g_rms_A1C0[ir]->SetMarkerColor(runColors[ir]);
        g_rms_A1C0[ir]->SetLineColor(runColors[ir]);

        g_rms_A1C1[ir] = new TGraphErrors(N_DITHER);
        g_rms_A1C1[ir]->SetTitle(runNames[ir]);
        g_rms_A1C1[ir]->SetMarkerStyle(21);
        g_rms_A1C1[ir]->SetMarkerColor(runColors[ir]);
        g_rms_A1C1[ir]->SetLineColor(runColors[ir]);

        // Kurtosis Graphs
        g_kurt_A1C0[ir] = new TGraph(N_DITHER);
        g_kurt_A1C0[ir]->SetTitle(runNames[ir]);
        g_kurt_A1C0[ir]->SetMarkerStyle(24);
        g_kurt_A1C0[ir]->SetMarkerColor(runColors[ir]);
        g_kurt_A1C0[ir]->SetLineColor(runColors[ir]);

        g_kurt_A1C1[ir] = new TGraph(N_DITHER);
        g_kurt_A1C1[ir]->SetTitle(runNames[ir]);
        g_kurt_A1C1[ir]->SetMarkerStyle(25);
        g_kurt_A1C1[ir]->SetMarkerColor(runColors[ir]);
        g_kurt_A1C1[ir]->SetLineColor(runColors[ir]);
    }

    // =========================================================================
    // --- Data Extraction Loop ---
    // =========================================================================
    for (int ir = 0; ir < N_RUNS; ir++)
    {
        for (int id = 0; id < N_DITHER; id++)
        {
            TString fileName = TString::Format("dither_scan/run%03d/dither_%.1f.root", runs[ir], dithers[id]);
            TFile *f = TFile::Open(fileName);

            if (!f || f->IsZombie()) {
                std::cout << "Skipping " << fileName << " (Not found)" << std::endl;
                continue;
            }

            auto fetchStats = [&](const char *tag, TGraphErrors *graph_rms, TGraph *graph_kurt)
            {
                TH1D *h_sx3 = (TH1D *)f->Get(TString::Format("Benchmark_SX3/Benchmark_SX3_PCZ-sx3pczguess_%s", tag));
                TH1D *h_qqq = (TH1D *)f->Get(TString::Format("Benchmark_QQQ/Benchmark_QQQ_PCZ-qqqpczguess_%s", tag));

                TH1D *h_combo = nullptr;
                if (h_sx3 && h_qqq) {
                    h_combo = (TH1D *)h_sx3->Clone(TString::Format("combo_%s", tag));
                    h_combo->Add(h_qqq);
                } else if (h_sx3) {
                    h_combo = (TH1D *)h_sx3->Clone(TString::Format("combo_%s", tag));
                } else if (h_qqq) {
                    h_combo = (TH1D *)h_qqq->Clone(TString::Format("combo_%s", tag));
                }

                if (h_combo && h_combo->GetEntries() > 50)
                {
                    // Use raw statistical properties instead of a Gaussian fit
                    double rms = h_combo->GetStdDev();
                    double rmsErr = h_combo->GetStdDevError();
                    double kurtosis = h_combo->GetKurtosis(); // < 0 means bifurcated, ~0 means Gaussian

                    graph_rms->SetPoint(id, dithers[id], rms);
                    graph_rms->SetPointError(id, 0.0, rmsErr);
                    graph_kurt->SetPoint(id, dithers[id], kurtosis);
                    
                    delete h_combo;
                }
            };

            fetchStats("A1C0", g_rms_A1C0[ir], g_kurt_A1C0[ir]);
            fetchStats("A1C1", g_rms_A1C1[ir], g_kurt_A1C1[ir]);

            f->Close();
            delete f;
        }
    }

    // =========================================================================
    // --- Plotting ---
    // =========================================================================
    TCanvas *c1 = new TCanvas("c1", "Dither Optimization", 1400, 1000);
    c1->Divide(2, 2);

    // -- Pad 1: A1C0 RMS
    c1->cd(1);
    gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_rms = new TMultiGraph();
    mg_A1C0_rms->SetTitle("A1C0 Emulation Resolution; Dither Sigma (mm); Raw RMS Spread (mm)");
    TLegend *leg0 = new TLegend(0.15, 0.65, 0.45, 0.85);
    leg0->SetBorderSize(1);
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_rms_A1C0[ir]->GetN() > 0) {
            mg_A1C0_rms->Add(g_rms_A1C0[ir], "PL");
            leg0->AddEntry(g_rms_A1C0[ir], runNames[ir], "pl");
        }
    }
    mg_A1C0_rms->Draw("A");
    leg0->Draw();

    // -- Pad 2: A1C1 RMS
    c1->cd(2);
    gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_rms = new TMultiGraph();
    mg_A1C1_rms->SetTitle("A1C1 Emulation Resolution; Dither Sigma (mm); Raw RMS Spread (mm)");
    TLegend *leg1 = new TLegend(0.15, 0.65, 0.45, 0.85);
    leg1->SetBorderSize(1);
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_rms_A1C1[ir]->GetN() > 0) {
            mg_A1C1_rms->Add(g_rms_A1C1[ir], "PL");
            leg1->AddEntry(g_rms_A1C1[ir], runNames[ir], "pl");
        }
    }
    mg_A1C1_rms->Draw("A");
    leg1->Draw();

    // -- Pad 3: A1C0 Kurtosis
    c1->cd(3);
    gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_kurt = new TMultiGraph();
    mg_A1C0_kurt->SetTitle("A1C0 Bifurcation Metric; Dither Sigma (mm); Excess Kurtosis (0 = Gaussian)");
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_kurt_A1C0[ir]->GetN() > 0) mg_A1C0_kurt->Add(g_kurt_A1C0[ir], "PL");
    }
    mg_A1C0_kurt->Draw("A");
    
    // Draw a horizontal line at y=0 to easily spot where it becomes Gaussian
    c1->Update();
    TLine *line0 = new TLine(c1->cd(3)->GetUxmin(), 0, c1->cd(3)->GetUxmax(), 0);
    line0->SetLineColor(kBlack); line0->SetLineStyle(2); line0->SetLineWidth(2);
    line0->Draw();

    // -- Pad 4: A1C1 Kurtosis
    c1->cd(4);
    gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_kurt = new TMultiGraph();
    mg_A1C1_kurt->SetTitle("A1C1 Bifurcation Metric; Dither Sigma (mm); Excess Kurtosis (0 = Gaussian)");
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_kurt_A1C1[ir]->GetN() > 0) mg_A1C1_kurt->Add(g_kurt_A1C1[ir], "PL");
    }
    mg_A1C1_kurt->Draw("A");
    
    c1->Update();
    TLine *line1 = new TLine(c1->cd(4)->GetUxmin(), 0, c1->cd(4)->GetUxmax(), 0);
    line1->SetLineColor(kBlack); line1->SetLineStyle(2); line1->SetLineWidth(2);
    line1->Draw();

    c1->SaveAs("dither_optimization_results.png");
}