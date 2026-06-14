#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
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
    // double dithers[] = { 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
    // 19.0, 20.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0};
    double dithers[] = { 8.0 ,12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0};
    const int N_DITHER = sizeof(dithers) / sizeof(dithers[0]);

    int runs[] = {9, 12, 18, 19, 20, 21};
    const int N_RUNS = sizeof(runs) / sizeof(runs[0]);

    // Safety: 6 runs defined, 6 colors and 6 names provided.
    int runColors[] = {kOrange + 1, kBlack, kRed, kBlue, kGreen + 2, kMagenta};
    TString runNames[] = {"27Al Run 9", "27Al Run 12", "17F Run 18", "17F Run 19", "17F Run 20", "17F Run 21"};

    // =========================================================================
    // --- Graph Arrays ---
    // =========================================================================
    // A1C0 Arrays
    TGraphErrors *g_rms_A1C0_orig[N_RUNS];
    TGraphErrors *g_rms_A1C0_hyb[N_RUNS];
    TGraphErrors *g_fit_A1C0_orig[N_RUNS];
    TGraphErrors *g_fit_A1C0_hyb[N_RUNS];
    TGraph *g_kurt_A1C0_orig[N_RUNS];
    TGraph *g_kurt_A1C0_hyb[N_RUNS];

    // A1C1 Arrays
    TGraphErrors *g_rms_A1C1_orig[N_RUNS];
    TGraphErrors *g_rms_A1C1_hyb[N_RUNS];
    TGraphErrors *g_fit_A1C1_orig[N_RUNS];
    TGraphErrors *g_fit_A1C1_hyb[N_RUNS];
    TGraph *g_kurt_A1C1_orig[N_RUNS];
    TGraph *g_kurt_A1C1_hyb[N_RUNS];

    // Initialize all graphs
    for (int ir = 0; ir < N_RUNS; ir++)
    {
        // ------------------ A1C0 Initialization ------------------
        g_rms_A1C0_orig[ir] = new TGraphErrors(N_DITHER);
        g_rms_A1C0_orig[ir]->SetTitle(runNames[ir]); g_rms_A1C0_orig[ir]->SetMarkerStyle(20); g_rms_A1C0_orig[ir]->SetMarkerColor(runColors[ir]); g_rms_A1C0_orig[ir]->SetLineColor(runColors[ir]);
        
        g_rms_A1C0_hyb[ir] = new TGraphErrors(N_DITHER);
        g_rms_A1C0_hyb[ir]->SetTitle(runNames[ir]); g_rms_A1C0_hyb[ir]->SetMarkerStyle(21); g_rms_A1C0_hyb[ir]->SetMarkerColor(runColors[ir]); g_rms_A1C0_hyb[ir]->SetLineColor(runColors[ir]);
        
        g_fit_A1C0_orig[ir] = new TGraphErrors(N_DITHER);
        g_fit_A1C0_orig[ir]->SetTitle(runNames[ir]); g_fit_A1C0_orig[ir]->SetMarkerStyle(22); g_fit_A1C0_orig[ir]->SetMarkerColor(runColors[ir]); g_fit_A1C0_orig[ir]->SetLineColor(runColors[ir]);
        
        g_fit_A1C0_hyb[ir] = new TGraphErrors(N_DITHER);
        g_fit_A1C0_hyb[ir]->SetTitle(runNames[ir]); g_fit_A1C0_hyb[ir]->SetMarkerStyle(23); g_fit_A1C0_hyb[ir]->SetMarkerColor(runColors[ir]); g_fit_A1C0_hyb[ir]->SetLineColor(runColors[ir]);
        
        g_kurt_A1C0_orig[ir] = new TGraph(N_DITHER);
        g_kurt_A1C0_orig[ir]->SetTitle(runNames[ir]); g_kurt_A1C0_orig[ir]->SetMarkerStyle(24); g_kurt_A1C0_orig[ir]->SetMarkerColor(runColors[ir]); g_kurt_A1C0_orig[ir]->SetLineColor(runColors[ir]);
        
        g_kurt_A1C0_hyb[ir] = new TGraph(N_DITHER);
        g_kurt_A1C0_hyb[ir]->SetTitle(runNames[ir]); g_kurt_A1C0_hyb[ir]->SetMarkerStyle(25); g_kurt_A1C0_hyb[ir]->SetMarkerColor(runColors[ir]); g_kurt_A1C0_hyb[ir]->SetLineColor(runColors[ir]);

        // ------------------ A1C1 Initialization ------------------
        g_rms_A1C1_orig[ir] = new TGraphErrors(N_DITHER);
        g_rms_A1C1_orig[ir]->SetTitle(runNames[ir]); g_rms_A1C1_orig[ir]->SetMarkerStyle(20); g_rms_A1C1_orig[ir]->SetMarkerColor(runColors[ir]); g_rms_A1C1_orig[ir]->SetLineColor(runColors[ir]);
        
        g_rms_A1C1_hyb[ir] = new TGraphErrors(N_DITHER);
        g_rms_A1C1_hyb[ir]->SetTitle(runNames[ir]); g_rms_A1C1_hyb[ir]->SetMarkerStyle(21); g_rms_A1C1_hyb[ir]->SetMarkerColor(runColors[ir]); g_rms_A1C1_hyb[ir]->SetLineColor(runColors[ir]);
        
        g_fit_A1C1_orig[ir] = new TGraphErrors(N_DITHER);
        g_fit_A1C1_orig[ir]->SetTitle(runNames[ir]); g_fit_A1C1_orig[ir]->SetMarkerStyle(22); g_fit_A1C1_orig[ir]->SetMarkerColor(runColors[ir]); g_fit_A1C1_orig[ir]->SetLineColor(runColors[ir]);
        
        g_fit_A1C1_hyb[ir] = new TGraphErrors(N_DITHER);
        g_fit_A1C1_hyb[ir]->SetTitle(runNames[ir]); g_fit_A1C1_hyb[ir]->SetMarkerStyle(23); g_fit_A1C1_hyb[ir]->SetMarkerColor(runColors[ir]); g_fit_A1C1_hyb[ir]->SetLineColor(runColors[ir]);
        
        g_kurt_A1C1_orig[ir] = new TGraph(N_DITHER);
        g_kurt_A1C1_orig[ir]->SetTitle(runNames[ir]); g_kurt_A1C1_orig[ir]->SetMarkerStyle(24); g_kurt_A1C1_orig[ir]->SetMarkerColor(runColors[ir]); g_kurt_A1C1_orig[ir]->SetLineColor(runColors[ir]);
        
        g_kurt_A1C1_hyb[ir] = new TGraph(N_DITHER);
        g_kurt_A1C1_hyb[ir]->SetTitle(runNames[ir]); g_kurt_A1C1_hyb[ir]->SetMarkerStyle(25); g_kurt_A1C1_hyb[ir]->SetMarkerColor(runColors[ir]); g_kurt_A1C1_hyb[ir]->SetLineColor(runColors[ir]);
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

            auto fetchStats = [&](const char *tag, TGraphErrors *graph_rms, TGraphErrors *graph_fit, TGraph *graph_kurt)
            {
                TString sx3_path = TString::Format("Benchmark_SX3_ref/Benchmark_SX3_PCZ_%s_minus_ref", tag);
                TString qqq_path = TString::Format("Benchmark_QQQ_ref/Benchmark_QQQ_PCZ_%s_minus_ref", tag);
                
                TH1D *h_sx3 = (TH1D *)f->Get(sx3_path);
                TH1D *h_qqq = (TH1D *)f->Get(qqq_path);

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
                    // 1. Raw Stats
                    double rms = h_combo->GetStdDev();
                    double rmsErr = h_combo->GetStdDevError();
                    double kurtosis = h_combo->GetKurtosis(); 

                    // 2. Gaussian Fit Stats
                    h_combo->Fit("gaus", "Q0"); // Q = Quiet, 0 = Do not draw
                    TF1 *gausFit = h_combo->GetFunction("gaus");
                    double fitSigma = 0.0;
                    double fitSigmaErr = 0.0;

                    if (gausFit) {
                        fitSigma = gausFit->GetParameter(2); // Param 2 is Sigma
                        fitSigmaErr = gausFit->GetParError(2);
                    }

                    // 3. Fill Graphs
                    graph_rms->SetPoint(id, dithers[id], rms);
                    graph_rms->SetPointError(id, 0.0, rmsErr);
                    
                    graph_fit->SetPoint(id, dithers[id], fitSigma);
                    graph_fit->SetPointError(id, 0.0, fitSigmaErr);
                    
                    graph_kurt->SetPoint(id, dithers[id], kurtosis);
                    
                    delete h_combo;
                }
            };

            // Fetch Normal A1C0 vs Hybrid A1C0
            fetchStats("A1C0", g_rms_A1C0_orig[ir], g_fit_A1C0_orig[ir], g_kurt_A1C0_orig[ir]);
            fetchStats("A1C0_Hyb", g_rms_A1C0_hyb[ir], g_fit_A1C0_hyb[ir], g_kurt_A1C0_hyb[ir]);

            // Fetch Normal A1C1 vs Hybrid A1C1
            fetchStats("A1C1", g_rms_A1C1_orig[ir], g_fit_A1C1_orig[ir], g_kurt_A1C1_orig[ir]);
            fetchStats("A1C1_Hyb", g_rms_A1C1_hyb[ir], g_fit_A1C1_hyb[ir], g_kurt_A1C1_hyb[ir]);

            f->Close();
            delete f;
        }
    }

    // =========================================================================
    // --- CANVAS 1: A1C0 Plotting ---
    // =========================================================================
    TCanvas *c1 = new TCanvas("c1", "A1C0 Normal vs Hybrid", 1400, 1400);
    c1->Divide(2, 3);

    // ROW 1: RAW RMS SPREAD
    c1->cd(1); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_orig_rms = new TMultiGraph();
    mg_A1C0_orig_rms->SetTitle("Original A1C0 (Pure Z-Dither); Dither Sigma (mm); PC Z Raw RMS (mm)");
    TLegend *leg0 = new TLegend(0.15, 0.65, 0.45, 0.85); leg0->SetBorderSize(1);
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_rms_A1C0_orig[ir]->GetN() > 0) {
            mg_A1C0_orig_rms->Add(g_rms_A1C0_orig[ir], "PL");
            leg0->AddEntry(g_rms_A1C0_orig[ir], runNames[ir], "pl");
        }
    }
    mg_A1C0_orig_rms->Draw("A"); leg0->Draw();

    c1->cd(2); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_hyb_rms = new TMultiGraph();
    mg_A1C0_hyb_rms->SetTitle("Hybrid A1C0 (Si Smear + Z-Dither); Dither Sigma (mm); PC Z Raw RMS (mm)");
    TLegend *leg1 = new TLegend(0.15, 0.65, 0.45, 0.85); leg1->SetBorderSize(1);
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_rms_A1C0_hyb[ir]->GetN() > 0) {
            mg_A1C0_hyb_rms->Add(g_rms_A1C0_hyb[ir], "PL");
            leg1->AddEntry(g_rms_A1C0_hyb[ir], runNames[ir], "pl");
        }
    }
    mg_A1C0_hyb_rms->Draw("A"); leg1->Draw();

    // ROW 2: GAUSSIAN FIT SIGMA
    c1->cd(3); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_orig_fit = new TMultiGraph();
    mg_A1C0_orig_fit->SetTitle("Original A1C0 Fit; Dither Sigma (mm); Gaussian Fit #sigma (mm)");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_fit_A1C0_orig[ir]->GetN() > 0) mg_A1C0_orig_fit->Add(g_fit_A1C0_orig[ir], "PL");
    mg_A1C0_orig_fit->Draw("A");

    c1->cd(4); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_hyb_fit = new TMultiGraph();
    mg_A1C0_hyb_fit->SetTitle("Hybrid A1C0 Fit; Dither Sigma (mm); Gaussian Fit #sigma (mm)");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_fit_A1C0_hyb[ir]->GetN() > 0) mg_A1C0_hyb_fit->Add(g_fit_A1C0_hyb[ir], "PL");
    mg_A1C0_hyb_fit->Draw("A");

    // ROW 3: KURTOSIS
    c1->cd(5); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_orig_kurt = new TMultiGraph();
    mg_A1C0_orig_kurt->SetTitle("Original A1C0 Bifurcation Metric; Dither Sigma (mm); Excess Kurtosis");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_kurt_A1C0_orig[ir]->GetN() > 0) mg_A1C0_orig_kurt->Add(g_kurt_A1C0_orig[ir], "PL");
    mg_A1C0_orig_kurt->Draw("A");
    
    c1->Update(); TLine *line0 = new TLine(c1->cd(5)->GetUxmin(), 0, c1->cd(5)->GetUxmax(), 0);
    line0->SetLineColor(kBlack); line0->SetLineStyle(2); line0->SetLineWidth(2); line0->Draw();

    c1->cd(6); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C0_hyb_kurt = new TMultiGraph();
    mg_A1C0_hyb_kurt->SetTitle("Hybrid A1C0 Bifurcation Metric; Dither Sigma (mm); Excess Kurtosis");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_kurt_A1C0_hyb[ir]->GetN() > 0) mg_A1C0_hyb_kurt->Add(g_kurt_A1C0_hyb[ir], "PL");
    mg_A1C0_hyb_kurt->Draw("A");
    
    c1->Update(); TLine *line1 = new TLine(c1->cd(6)->GetUxmin(), 0, c1->cd(6)->GetUxmax(), 0);
    line1->SetLineColor(kBlack); line1->SetLineStyle(2); line1->SetLineWidth(2); line1->Draw();

    c1->SaveAs("a1c0_normal_vs_hybrid_results.png");

    // =========================================================================
    // --- CANVAS 2: A1C1 Plotting ---
    // =========================================================================
    TCanvas *c2 = new TCanvas("c2", "A1C1 Normal vs Hybrid", 1400, 1400);
    c2->Divide(2, 3);

    // ROW 1: RAW RMS SPREAD
    c2->cd(1); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_orig_rms = new TMultiGraph();
    mg_A1C1_orig_rms->SetTitle("Original A1C1 (Pure Z-Dither); Dither Sigma (mm); PC Z Raw RMS (mm)");
    TLegend *leg2 = new TLegend(0.15, 0.65, 0.45, 0.85); leg2->SetBorderSize(1);
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_rms_A1C1_orig[ir]->GetN() > 0) {
            mg_A1C1_orig_rms->Add(g_rms_A1C1_orig[ir], "PL");
            leg2->AddEntry(g_rms_A1C1_orig[ir], runNames[ir], "pl");
        }
    }
    mg_A1C1_orig_rms->Draw("A"); leg2->Draw();

    c2->cd(2); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_hyb_rms = new TMultiGraph();
    mg_A1C1_hyb_rms->SetTitle("Hybrid A1C1 (Si Smear + Z-Dither); Dither Sigma (mm); PC Z Raw RMS (mm)");
    TLegend *leg3 = new TLegend(0.15, 0.65, 0.45, 0.85); leg3->SetBorderSize(1);
    for (int ir = 0; ir < N_RUNS; ir++) {
        if (g_rms_A1C1_hyb[ir]->GetN() > 0) {
            mg_A1C1_hyb_rms->Add(g_rms_A1C1_hyb[ir], "PL");
            leg3->AddEntry(g_rms_A1C1_hyb[ir], runNames[ir], "pl");
        }
    }
    mg_A1C1_hyb_rms->Draw("A"); leg3->Draw();

    // ROW 2: GAUSSIAN FIT SIGMA
    c2->cd(3); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_orig_fit = new TMultiGraph();
    mg_A1C1_orig_fit->SetTitle("Original A1C1 Fit; Dither Sigma (mm); Gaussian Fit #sigma (mm)");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_fit_A1C1_orig[ir]->GetN() > 0) mg_A1C1_orig_fit->Add(g_fit_A1C1_orig[ir], "PL");
    mg_A1C1_orig_fit->Draw("A");

    c2->cd(4); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_hyb_fit = new TMultiGraph();
    mg_A1C1_hyb_fit->SetTitle("Hybrid A1C1 Fit; Dither Sigma (mm); Gaussian Fit #sigma (mm)");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_fit_A1C1_hyb[ir]->GetN() > 0) mg_A1C1_hyb_fit->Add(g_fit_A1C1_hyb[ir], "PL");
    mg_A1C1_hyb_fit->Draw("A");

    // ROW 3: KURTOSIS
    c2->cd(5); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_orig_kurt = new TMultiGraph();
    mg_A1C1_orig_kurt->SetTitle("Original A1C1 Bifurcation Metric; Dither Sigma (mm); Excess Kurtosis");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_kurt_A1C1_orig[ir]->GetN() > 0) mg_A1C1_orig_kurt->Add(g_kurt_A1C1_orig[ir], "PL");
    mg_A1C1_orig_kurt->Draw("A");
    
    c2->Update(); TLine *line2 = new TLine(c2->cd(5)->GetUxmin(), 0, c2->cd(5)->GetUxmax(), 0);
    line2->SetLineColor(kBlack); line2->SetLineStyle(2); line2->SetLineWidth(2); line2->Draw();

    c2->cd(6); gPad->SetGrid(1, 1);
    TMultiGraph *mg_A1C1_hyb_kurt = new TMultiGraph();
    mg_A1C1_hyb_kurt->SetTitle("Hybrid A1C1 Bifurcation Metric; Dither Sigma (mm); Excess Kurtosis");
    for (int ir = 0; ir < N_RUNS; ir++) if (g_kurt_A1C1_hyb[ir]->GetN() > 0) mg_A1C1_hyb_kurt->Add(g_kurt_A1C1_hyb[ir], "PL");
    mg_A1C1_hyb_kurt->Draw("A");
    
    c2->Update(); TLine *line3 = new TLine(c2->cd(6)->GetUxmin(), 0, c2->cd(6)->GetUxmax(), 0);
    line3->SetLineColor(kBlack); line3->SetLineStyle(2); line3->SetLineWidth(2); line3->Draw();

    c2->SaveAs("a1c1_normal_vs_hybrid_results.png");
}