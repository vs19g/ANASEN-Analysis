{
    int index = 1;
    TFile *f = new TFile("../results_run19.root");
    TH2F *h2=NULL;
    TH1F *h1x=NULL, *h1y=NULL;
    //f->cd("evsx");
    //f->ls();

    double known_xpos[2][2] = {{0,18.75},{-18.75,0}};
    std::vector<double> xpos, xposkn; //first = x = known position, second = y = unknown position

    std::ofstream ofb(Form("backgains%d.dat",index));
    std::ofstream off(Form("frontgains%d.dat",index));
    for(int i=1; i<4; i++) {
        //do it for pad#2
        int backnum=2;
        h2 = (TH2F*)(f->Get(Form("evsx/be_vs_x_sx3_id_%d_f%d_b%d",index,i,backnum)));
        auto macro = [&]() {
            h1x = (TH1F*)(h2->ProjectionX("_px"));
            double xleft = h1x->GetBinCenter(h1x->FindFirstBinAbove(h1x->GetMaximum()*0.4));
            double xright = h1x->GetBinCenter(h1x->FindLastBinAbove(h1x->GetMaximum()*0.4));
            //h1x->GetXaxis()->SetRangeUser(4*xleft, xright*4);
            h1x->Draw();
            TLine L1(xleft,0,xleft,h1x->GetMaximum()); L1.SetLineColor(kRed); L1.Draw("SAME");
            TLine L2(xright,0,xright,h1x->GetMaximum()); L2.SetLineColor(kRed); L2.Draw("SAME");
            gPad->Modified();
            gPad->Update();
            xpos.push_back(xleft); xposkn.push_back(known_xpos[backnum-1][0]);
            xpos.push_back(xright); xposkn.push_back(known_xpos[backnum-1][1]);
            while(gPad->WaitPrimitive());

            h1y = (TH1F*)(h2->ProjectionY("_py"));
            double ycenter = h1y->GetBinCenter(h1y->GetMaximumBin());
//            std::cout << "front " << i << " back " << backnum << " "  << xleft << " " << xright << " " << ycenter << " " << 5486/ycenter << std::endl;
            ofb << index <<" front " << i << " back " << backnum << " "  << 5486/ycenter << std::endl;
            h1y->GetXaxis()->SetRangeUser(ycenter-200,ycenter+200);
            h1y->Draw();
            TLine L3(ycenter,0,ycenter,h1y->GetMaximum()*1.1); L3.SetLineColor(kRed); L3.Draw("SAME");

            gPad->Modified();
            gPad->Update();
            while(gPad->WaitPrimitive());
        };
        if(h2)
            macro();

        //repeat for pad#1
        backnum=1;
        h2 = (TH2F*)(f->Get(Form("evsx/be_vs_x_sx3_id_%d_f%d_b%d",index,i,backnum)));
        if(h2)
            macro();

        double xtofit[] = {xpos[0],xpos[3]};
        double xktofit[] = {xposkn[0],xposkn[3]};
        TGraph G1(xpos.size(),xpos.data(),xposkn.data());
        G1.Draw("APL*");
        G1.Fit("pol1","Q");
        off << index<<" lengthcal front " << i << " " << G1.GetFunction("pol1")->GetParameter(0) << " " << G1.GetFunction("pol1")->GetParameter(1) << std::endl;
        gPad->Modified(); gPad->Update();
        while(gPad->WaitPrimitive());
        xpos.clear();
        xposkn.clear();
    }
    ofb.close();
    off.close();
    f->Close();
}
