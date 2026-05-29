#include "../Armory/PC_StepLadder_Correction.h"
int quit = 0;
void handler(int) { quit = 1; }

int colors[] = {kSpring + 3, kRed, kGreen + 3, kBlue + 3, kViolet, kOrange, kSpring - 7, kAzure - 5};
void PhiCorr()
{
    signal(SIGINT, handler);
    TCanvas c("c1", "c1", 0, 0, 1600, 800);

    std::vector<TFile *> files;

    for (int i = 12; i <= 300; i++)
    {
        TFile *f = new TFile(Form("results_run%d.root", i));
        if (!f || f->IsZombie() || !f->IsOpen()) { delete f; continue; }

        TH1F *h1 = (TH1F *)(f->Get("pcphi_minus_sx3phi_TC"));
        TH1F *h2 = (TH1F *)(f->Get("phiPC_minus_phiQQQ_TC"));
        if (!h1 || !h2) { f->Close(); delete f; continue; }

        files.push_back(f);

        c.Clear();
        auto c1 = c.cd(1);
        c1->SetGrid(1, 1);

        h1->SetLineColorAlpha(kRed, 0.75);
        h2->GetXaxis()->SetRangeUser(-100, 100);
        h2->Draw("");

        h2->SetTitle(Form("run%d", i));
        h2->SetLineColorAlpha(kBlue, 0.75);
        h1->Draw("SAME");

        c1->Modified();
        c1->Update();
        c.SaveAs(Form("scratch/PhiCorr_run%d.png", i));
        while (gPad->WaitPrimitive());
        if (quit) break;
    }

    for (auto file : files) file->Close();
}