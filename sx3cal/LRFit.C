{
    TFile *f = new TFile("../results_run19.root");
    f->cd("l_vs_r");
    gDirectory->ls();
    int clkpos = 13;
    std::ofstream ofile(Form("rightgains%d.dat",clkpos));
    for(int i=1; i<4; i++) {
        TH2F h2(*(TH2F*)(f->Get(Form("l_vs_r/l_vs_r_sx3_id_%d_f%d",clkpos,i))));
        h2.Draw();

        TH1F hproj(*(TH1F*)(h2.ProjectionX("_px")));
        /*hproj.Draw("SAME");
        gPad->Modified();
        gPad->Update();
        while(gPad->WaitPrimitive());*/

        int leftbin = hproj.FindFirstBinAbove(hproj.GetMaximum()*0.4);
        int rightbin = hproj.FindLastBinAbove(hproj.GetMaximum()*0.1);

        TH1F h1(*(TH1F*)(h2.ProfileX("_pfx",leftbin,rightbin)));
        h1.Draw("histo same");
        TLine L1(h1.GetBinCenter(leftbin),0,h1.GetBinCenter(leftbin),1000); L1.SetLineColor(kRed); L1.Draw("SAME");
        TLine L2(h1.GetBinCenter(rightbin),0,h1.GetBinCenter(rightbin),1000); L2.SetLineColor(kRed); L2.Draw("SAME");
        //h2.GetYaxis()->SetRangeUser(0,2000);
        //h2.GetXaxis()->SetRangeUser(hproj.GetBinCenter(leftbin),hproj.GetBinCenter(rightbin));
        h2.Fit("pol1","","SAME",h1.GetBinCenter(leftbin),h1.GetBinCenter(rightbin));

        TF1 *f1 = (TF1*)h2.GetFunction("pol1");
        f1->Draw("SAME");
        ofile << clkpos << " " << i << " " << f1->GetParameter(0) << " " << TMath::Abs(f1->GetParameter(1)) << std::endl;
        gPad->Modified();
        gPad->Update();
        while(gPad->WaitPrimitive());
    }
    ofile.close();
    f->Close();
}
