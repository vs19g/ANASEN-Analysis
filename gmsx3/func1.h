#include "TF1.h"

double model2(double *x, double *par) {
    /*  'Potential Well' of width 2a from from xx-a to xx+a
        xx is coordinate about the point of origin, set at x=center
        v0 is the y-offset of the potential
        k is the 'steepness' of the potential

        continuous across xx-a and xx+a, and differentiable
    */

    double center= par[3];
    double xx = x[0]-center;
    double a = TMath::Abs(par[0]);
    double k = TMath::Abs(par[1]);
    double v0 = par[2];

    if(xx < -a)
        return k*(xx+a)*(xx+a) + v0;
    else if(xx > a)
        return k*(xx-a)*(xx-a) + v0;
    else
        return v0;
}

void func1() {
    //TF1 f1("bowl",model,-2.,2.,2);
    TCanvas c("c1","c1",800,600);
    TF1 f1("bowl",model2,-10.,10.,4);
    f1.SetMaximum(10);

    for(int i=-4; i<4; i++) {
        f1.SetParameters(.4,100,2,i); //a, k, v0, center
        f1.SetNpx(100000);
        if(i==-4) f1.Draw("L");
        f1.DrawCopy("L SAME");
        c.Modified(); c.Update();
        //c.SaveAs(Form("%d.png",out));
        while(c.WaitPrimitive());
   }
}
