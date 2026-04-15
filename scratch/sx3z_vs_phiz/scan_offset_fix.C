#include "testmodel.h"

int quit=0;
void handler(int){quit=1;}

int colors[] = {kSpring+3, kRed, kGreen+3, kBlue+3, kViolet, kOrange, kSpring-7, kAzure-5};
void scan_offset_fix(){
   signal(SIGINT,handler);
   TCanvas c("c1","c1",0,0,1600,800);
   c.Divide(2,1);

    TF1 f1("model",model,-200,200,2);
    f1.SetNpx(10000);
    std::vector<double> pars = {0.0,1.};
    f1.SetParameters(pars.data());
    f1.SetLineColor(kGreen+2);
    f1.SetLineStyle(kLine);




   TFile* f=NULL;
   std::vector<TFile*> files;
   int ctr=0;
    for(int i=12; i<=21; i++) {
        auto c1=c.cd(1);
        c1->SetGrid(1,1);
        f = new TFile(Form("../../results_run%d.root",i));

        if(i==12) {
            //TH2F *h2 = (TH2F*)(f->Get("phicut/pczguess_vs_pc_int"));
            TH2F *h23 = (TH2F*)(f->Get("pczfix_vs_qqqpczguess_A1C2"));
            h23->SetLineColorAlpha(kOrange,0.75);
            h23->Draw("SAME");

        }  else {
            //TH2F *h2 = (TH2F*)(f->Get("phicut/pczguess_vs_pc_int"));
            //TH2F *h2 = (TH2F*)(f->Get("pcz_vs_sx3pczguess_A1C2_strip12"));
            TH2F *h2 = (TH2F*)(f->Get("pczfix_vs_sx3pczguess_A1C2"));
            //TH2F *h2 = (TH2F*)(f->Get("hPCQQQ/PC_XY_Projection_QQQ2"));
            if(!h2) continue;
            h2->SetTitle(Form("case%d",i));
            //h2->Draw("colz same");
            h2->SetLineColorAlpha(colors[ctr],0.75);
            h2->Draw("box same");
            //f1.Draw("same");
        }
        TF1 eqline("x","x",-200,200);
        eqline.Draw("SAME");
        c1->Modified();
        c1->Update();
        ctr+=1;


        auto c2=c.cd(2);
        c2->SetGrid(1,1);

        TH2F *h3 = (TH2F*)(f->Get("sx3phi_vs_pcphi1"));
//        TH2F *h2 = (TH2F*)(f->Get("hPCQQQ/PC_XY_Projection_QQQ2"));
        if(!h3) continue;
        h3->SetTitle(Form("case%d",i));
        h3->Draw("colz");
        eqline.Draw("SAME");
        c2->Modified();
        c2->Update();

        while(gPad->WaitPrimitive());

        files.emplace_back(f);
        if(i==21) {
            i=11;
            c.Clear();
            c.Divide(2,1);
            ctr=0;
        }
        if(quit) break;
    }
    for(auto file : files) {
        file->Close();
    }
}
