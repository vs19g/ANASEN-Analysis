void anode_gainmatch(){
    TFile *f = new TFile("../results_run16.root");

    TH2F *pc_index_h2d = (TH2F*)(f->Get("hRawPC/PC_Index_Vs_Energy"));
    std::cout << pc_index_h2d << std::endl;
    TCanvas c("c1","c1",0,0,1600,800);
    //TCanvas c_g("cg","cg",0,900,400,400);
    c.Divide(2,1);
    auto c1=c.cd(1);
    pc_index_h2d->Draw("COLZ");
    pc_index_h2d->GetYaxis()->SetRangeUser(240,5000);
    auto c2=c.cd(2);
    c2->SetLogy();
    TH1F *h_1d=NULL;
    int bin_index=1;
    std::vector<std::vector<double>> all_peaks;
    std::vector<int> found_wire_list;
    while(bin_index<=24) {
        h_1d=(TH1F*)(pc_index_h2d->ProjectionY("_py",bin_index,bin_index));
        auto c1 = c.cd(1);
        TBox box(pc_index_h2d->GetXaxis()->GetBinLowEdge(bin_index),0,pc_index_h2d->GetXaxis()->GetBinUpEdge(bin_index),pc_index_h2d->GetYaxis()->GetXmax());
        box.SetFillColorAlpha(kYellow+3,0.3);
        box.Draw("SAME");
        c1->Modified(); c1->Update();
        //while(c1->WaitPrimitive());

        TSpectrum s;
        auto c2 = c.cd(2);
        h_1d->Draw();
        c2->Modified(); c2->Update();
        int npeaks = s.Search(h_1d,8,"",0.02); std::cout << npeaks << std::endl;
        if(npeaks>=3) {
            std::vector<double> xpeaks(s.GetPositionX(),s.GetPositionX()+npeaks);
            std::sort(xpeaks.begin(),xpeaks.end(),std::greater());
            found_wire_list.push_back((int)pc_index_h2d->GetXaxis()->GetBinCenter(bin_index));
            all_peaks.push_back(xpeaks);
        }
        while(c2->WaitPrimitive());
        bin_index++;
    }
    c.cd(2)->SetLogy(kFALSE);
    gStyle->SetOptFit(1111);

    std::ofstream outfile("anode_gm_coeffs.dat");
    outfile <<  found_wire_list.at(0) << " "
            << 1.0  << " "
            << 0.0  << std::endl;

    for(int i=0; i<all_peaks.size(); i++){
        if(i==1) continue;
        TGraph g(all_peaks.at(i).size(), all_peaks.at(i).data(), all_peaks.at(1).data());
        auto c2 = c.cd(2);
        g.SetMarkerStyle(20);
        //g.Print();
        g.Draw("AP");
        g.Fit("pol1");
        outfile <<  found_wire_list.at(i) << " "
                << ((TF1*)g.FindObject("pol1"))->GetParameter(1) << " "
                << ((TF1*)g.FindObject("pol1"))->GetParameter(0) << std::endl;
        c2->Modified();
        c2->Update();
        while(c2->WaitPrimitive());
    }
    outfile.close();
    f->Close();
    return;
}
