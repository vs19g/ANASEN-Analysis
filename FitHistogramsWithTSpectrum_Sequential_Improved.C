#include <TFile.h>
#include <TH1.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <TText.h>

void FitHistogramsWithTSpectrum_Sequential_Improved() {
    TFile *inputFile = new TFile("../Histograms_anodes.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening the input file!" << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "Histogram Viewer", 800, 600);

    // Open the output ASCII file to save the centroids
    std::ofstream outFile("centroids.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        return;
    }
    outFile << "HistogramIndex\tPeakNumber\tCentroid\tAmplitude\tSigma" << std::endl;

    for (int i = 0; i < 24; ++i) {
        TH1 *histogram = dynamic_cast<TH1*>(inputFile->Get(Form("hCathode_%d", i)));
        if (!histogram) {
            std::cerr << "Failed to retrieve histogram_" << i << " from the file." << std::endl;
            continue;
        }

        // Set range for peak search
        double minX = 700;
        double maxX = 25000;
        histogram->GetXaxis()->SetRangeUser(minX, maxX);

        // Draw the histogram
        c1->cd();
        histogram->Draw();

        // Peak search using TSpectrum
        const int maxPeaks = 5;
        TSpectrum spectrumFinder(maxPeaks);
        int nFound = spectrumFinder.Search(histogram, 2, "", 0.01);

        if (nFound <= 0) {
            std::cerr << "No peaks found for histogram " << i << std::endl;
            continue;
        }

        Double_t *xPositions = spectrumFinder.GetPositionX();
        Double_t *yPositions = spectrumFinder.GetPositionY();
        std::vector<std::pair<Double_t, Double_t>> peaks;

        // Collect and sort peaks by X position
        for (int j = 0; j < nFound; ++j) {
            peaks.emplace_back(xPositions[j], yPositions[j]);
        }
        std::sort(peaks.begin(), peaks.end());

        // Fit each peak with a Gaussian
        for (int j = 0; j < peaks.size(); ++j) {
            Double_t peakX = peaks[j].first;
            Double_t peakY = peaks[j].second;
            Double_t initialAmplitude = peakY;  // Better initial guess
            Double_t initialCentroid = peakX;   // Centroid based on peak position
            Double_t initialSigma = 60.0;
            // Define Gaussian with initial parameters
            TF1 *gaussFit = new TF1(Form("gauss_%d", j), "gaus", peakX - 200, peakX + 200);
            //gaussFit->SetParameters(peakY, peakX, 25.0); // Initial guesses for amplitude, mean, sigma
            gaussFit->SetParameters(initialAmplitude, initialCentroid, initialSigma);
            // Perform fit
            int fitStatus = histogram->Fit(gaussFit, "RQ+");
            if (fitStatus != 0) {
                std::cerr << "Fit failed for peak " << j + 1 << " in histogram " << i << std::endl;
                delete gaussFit;
                continue;
            }

            // Retrieve fit parameters
            double amplitude = gaussFit->GetParameter(0);
            double centroid = gaussFit->GetParameter(1);
            double sigma = gaussFit->GetParameter(2);
            double amplitudeError = gaussFit->GetParError(0);
            double centroidError = gaussFit->GetParError(1);
            double sigmaError = gaussFit->GetParError(2);

            // Chi-squared value
            double chi2 = gaussFit->GetChisquare();
            int ndf = gaussFit->GetNDF();
            outFile << i << "\t" << j + 1 << "\t" << centroid << std::endl;
            gaussFit->SetLineColor(kRed);
            gaussFit->Draw("SAME");
            TText *text = new TText();
            text->SetNDC();
            text->SetTextSize(0.03);
            text->SetTextColor(kRed);
            //text->DrawText(0.15, 0.8 - j * 0.05, Form("Peak %d: Amp=%.2f, Mean=%.2f, Sigma=%.2f", j + 1, amplitude, centroid, sigma));
            text->DrawText(0.15, 0.8 - j * 0.05, 
                Form("Peak %d: Amp=%.2f±%.2f, Mean=%.2f±%.2f, Sigma=%.2f±%.2f, Chi2/NDF=%.2f",
                     j + 1, amplitude, amplitudeError, centroid, centroidError, sigma, sigmaError, chi2 / ndf));
            // Save results
            
            
            // Clean up
            delete gaussFit;
        }

        // Update canvas for visualization
        c1->Update();
        std::cout << "Press Enter to view the next histogram..." << std::endl;
        c1->WaitPrimitive(); // Wait until Enter is pressed in the ROOT console
    }

    // Close resources
    inputFile->Close();
    outFile.close();
    delete c1;
}

