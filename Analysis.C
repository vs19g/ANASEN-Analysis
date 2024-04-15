#include "TChain.h"
#include <iostream>



void Analysis(int start, int end) {
    // Create a TChain
    TChain *chain = new TChain("tree");
    for(int i = start; i < end+1; i++) {
        chain->Add(Form("data/root_data/Run_%03d_mapped.root", i));
    }

    // Process the chain using Analyzer.C+
    chain->Process("Analyzer.C+");
}

// Define a macro with the same name as the script
void Analysis() {
    Analysis(72, 194); // Adjust the range if needed
}