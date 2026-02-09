#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <iostream>

void RunTimeSummary(int startRun, int endRun)
{
    TString fileDir = "/mnt/d/Remapped_files/17F_data/root_data/";
    TString histName = "AnodeQQQ_Time";
    TString filePattern = "Run_%03d_mapped_histograms.root";
    TString filePatternAlt = "ProtonRun_%d_mapped_histograms.root";
    TString filePatternAlt2 = "Source_%d_mapped_histograms.root";

    int nBinsTime = 0;
    double timeMin = 0, timeMax = 0;
    bool foundRef = false;

    for (int r = startRun; r <= endRun; r++)
    {
        TString tempName;

        // 1. Try Pattern 1: Run_XXX...
        tempName = fileDir + Form(filePattern, r);
        if (gSystem->AccessPathName(tempName))
        { // Returns true if MISSING

            // 2. Try Pattern 2: ProtonRun_X...
            tempName = fileDir + Form(filePatternAlt, r);
            if (gSystem->AccessPathName(tempName))
            {

                // 3. Try Pattern 3: Source_X...
                tempName = fileDir + Form(filePatternAlt2, r);
                if (gSystem->AccessPathName(tempName))
                {
                    // All 3 patterns failed. Skip this run.
                    continue;
                }
            }
        }

        // If we get here, 'tempName' holds the valid filename that was found
        TFile *fTemp = TFile::Open(tempName);
        if (!fTemp || fTemp->IsZombie())
        {
            if (fTemp)
                delete fTemp;
            continue;
        }

        TH1F *hRef = (TH1F *)fTemp->Get(histName);
        if (hRef)
        {

            nBinsTime = hRef->GetNbinsX();
            timeMin = hRef->GetXaxis()->GetXmin();
            timeMax = hRef->GetXaxis()->GetXmax();
            foundRef = true;

            delete hRef;
            fTemp->Close();
            delete fTemp;
            printf("Reference found in Run %d: %d bins, Range [%.1f, %.1f]\n", r, nBinsTime, timeMin, timeMax);
            break;
        }
        fTemp->Close();
        delete fTemp;
    }

    if (!foundRef)
    {
        printf("Error: No valid histograms found in the entire range. Exiting.\n");
        return;
    }

    int nRuns = endRun - startRun + 1;
    TH2F *hSummary = new TH2F("hSummary",
                              Form("Timing Summary (Runs %d-%d);Timing;Run Number", startRun, endRun),
                              nBinsTime, timeMin, timeMax,
                              nRuns, startRun, endRun + 1);

    for (int run = startRun; run <= endRun; run++)
    {

        TString filename = fileDir + Form(filePattern, run);

        if (gSystem->AccessPathName(filename))
            continue;

        TFile *fin = TFile::Open(filename);
        if (!fin || fin->IsZombie())
        {
            if (fin)
                delete fin;
            continue;
        }

        TH1F *hin = (TH1F *)fin->Get(histName);

        if (hin)
        {
            // Determine which ROW (Y-bin) corresponds to this Run
            // Note: ROOT bins start at 1.
            // If startRun=10 and run=10 -> binY=1.
            int binY_Run = run - startRun + 1;

            // Loop through the Time bins (X-bins in the 1D hist)
            for (int binX_Time = 1; binX_Time <= hin->GetNbinsX(); binX_Time++)
            {

                double content = hin->GetBinContent(binX_Time);

                // Copy content to: (Time, Run)
                if (content > 0)
                {
                    hSummary->SetBinContent(binX_Time, binY_Run, content);
                }
            }
            delete hin;
        }

        fin->Close();
        delete fin;

        if ((run - startRun) % 10 == 0)
            printf("Stitched Run %d...\n", run);
    }

    TFile *fOut = new TFile("SummaryPlot.root", "RECREATE");
    hSummary->Write();

    TCanvas *c1 = new TCanvas("c1", "Time Summary Plot", 1000, 800);
    hSummary->SetStats(0);
    hSummary->Draw("COLZ");

    printf("Done! Saved to SummaryPlot.root\n");
}