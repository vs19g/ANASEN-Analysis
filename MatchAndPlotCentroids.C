#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1.h>



void MatchAndPlotCentroids() {
    // Open the centroid data file
    std::ifstream inputFile("centroids.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open Centroids.txt" << std::endl;
        return;
    }

    // Data structure to store centroids by histogram and peak number
    std::map<int, std::map<int, double>> centroidData;

    // Read data from the file
    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        int histogramIndex, peakNumber;
        double centroid;
        if (iss >> histogramIndex >> peakNumber >> centroid) {
            centroidData[histogramIndex][peakNumber] = centroid;
        }
    }

    inputFile.close();

    // Ensure histogram 24 exists and has data
    if (centroidData.find(1) == centroidData.end()) {
        std::cerr << "Error: Histogram 0 not found in the data!" << std::endl;
        return;
    }

    // Reference centroids from histogram 24
    const auto& referenceCentroids = centroidData[1];
    std::ofstream outputFile("slope_intercept_results.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the output file for writing!" << std::endl;
        return;
    }
    outputFile << "Histogram Number\tSlope\tIntercept\n";
    // Loop through histograms 25 to 47
    for (int targetHist = 0; targetHist <= 23; targetHist++) {
        // Ensure the target histogram exists and matches in peak numbers
        if (centroidData.find(targetHist) == centroidData.end() || centroidData[targetHist].size() != referenceCentroids.size()) {
            //4th cnetroid data point for 19 was generated using the 3 datqa points for the slope of wires 0 and 19
            std::cout << "Skipping Histogram " << targetHist << " due to mismatched or missing data." << std::endl;
            continue;
        }

        // Prepare x and y values for TGraph
        std::vector<double> xValues, yValues;
        for (const auto& [peakNumber, refCentroid] : referenceCentroids) {
            if (centroidData[targetHist].find(peakNumber) != centroidData[targetHist].end()) {
                yValues.push_back(refCentroid);
                xValues.push_back(centroidData[targetHist][peakNumber]);
            } else {
                std::cerr << "Warning: Peak " << peakNumber << " missing in histogram " << targetHist << std::endl;
            }
        }

        if (xValues.size() < 3) {
            std::cout << "Skipping Histogram " << targetHist << " as it has less than 3 matching centroids." << std::endl;
            continue;
        }

        // Create a TGraph
        TCanvas *c1 = new TCanvas(Form("c_centroid_1_vs_%d", targetHist), Form("Centroid 1 vs %d", targetHist), 800, 600);
        TGraph *graph = new TGraph(xValues.size(), &xValues[0], &yValues[0]);
        graph->SetTitle(Form("Centroid of Histogram  %d vs 1", targetHist));
        graph->GetYaxis()->SetTitle("Centroid of Histogram 1");
        graph->GetXaxis()->SetTitle(Form("Centroid of Histogram %d", targetHist));
        graph->SetMarkerStyle(20); // Full circle marker
        graph->SetMarkerSize(1.0);
        graph->SetMarkerColor(kBlue);
        // Draw the graph
        graph->Draw("AP");
        double minX = *std::min_element(xValues.begin(), xValues.end());
        double maxX = *std::max_element(xValues.begin(), xValues.end());
        // Fit the data with a linear function
        TF1 *fitLine = new TF1("fitLine", "pol1", minX, maxX);  // Adjust range as needed
	fitLine->SetLineColor(kRed); // Set the line color to distinguish it
        fitLine->SetLineWidth(2);    // Thicker line for visibility
        graph->Fit(fitLine, "M");
        fitLine->Draw("same");
        fitLine->SetParLimits(0, -10, 10); // Limit intercept between -10 and 10
        fitLine->SetParLimits(1, 0, 2); 
        // Extract slope and intercept
        double slope = fitLine->GetParameter(1);
        double intercept = fitLine->GetParameter(0);
        outputFile << targetHist << "\t" << slope << "\t" << intercept << "\n";
        std::cout << "Histogram 24 vs " << targetHist << ": Slope = " << slope << ", Intercept = " << intercept << std::endl;
        std::vector<double> residuals;
    for (size_t i = 0; i < xValues.size(); ++i) {
        double fittedY = fitLine->Eval(xValues[i]); // Evaluate fitted function at x
        double residual = yValues[i] - fittedY;    // Residual = observed - fitted
        residuals.push_back(residual);
    }

    // Create a graph for the residuals
    /*TGraph *residualGraph = new TGraph(residuals.size(), &xValues[0], &residuals[0]);
    residualGraph->SetTitle(Form("Residuals for Histogram 24 vs %d", targetHist));
    residualGraph->GetYaxis()->SetTitle("Residuals");
    residualGraph->GetXaxis()->SetTitle(Form("Centroid of Histogram %d", targetHist));
    residualGraph->SetMarkerStyle(20);
    residualGraph->SetMarkerSize(1.0);
    residualGraph->SetMarkerColor(kGreen);

    // Draw the residuals plot below the original plot (can be on a new canvas if preferred)
    TCanvas *c2 = new TCanvas(Form("c_residuals_24_vs_%d", targetHist), Form("Residuals for Centroid 24 vs %d", targetHist), 800, 400);
    residualGraph->Draw("AP");*/
        c1->Update();
        //c2->Update();
        std::cout << "Press Enter to continue..." << std::endl;
        
        //std::cin.get();
        c1->WaitPrimitive();
        //c2->WaitPrimitive();
        //std::cin.get();
        //std::cin.get();
    }
    outputFile.close();
    std::cout << "Results written to slope_intercept_results.txt" << std::endl;
}
