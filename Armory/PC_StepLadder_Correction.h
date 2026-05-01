#include <TF1.h>
double model(double* x, double* p) {
    double result = x[0];
    double factor = 29.0;
    double slope = 0.7;
    if(TMath::Abs(x[0]) < 16.2) result=x[0]*slope;
    else if(TMath::Abs(x[0]) < 49.8 ) result=x[0]*slope+TMath::Sign(1.0,x[0])*factor;
    else if(TMath::Abs(x[0]) < 85.2 ) result=x[0]*slope+TMath::Sign(1.0,x[0])*factor*2;
    else result=x[0]*slope+TMath::Sign(1.0,x[0])*factor*2;
    return result;
}

double model_invert(double *y, double *q) {
    double result=y[0];
    double slope = 0.7;
    double factor = 40.0;
    if(TMath::Abs(y[0]) < 16.2/slope) result = y[0]/slope;
    else if(TMath::Abs(y[0]) < 49.8/slope ) result=y[0]/slope-TMath::Sign(1.0,y[0])*factor;
    else if(TMath::Abs(y[0]) < 85.2/slope ) result=y[0]/slope-TMath::Sign(1.0,y[0])*factor*2;
    else result=y[0]/slope-TMath::Sign(1.0,y[0])*factor*2;
    return result+40;
}

/*void testmodel() {
    TF1 eqline("x","x",-200,200);
    eqline.Draw("");
    eqline.SetLineStyle(kDashed);

    //TF1 f1("model",model,-200,200,2);
    TF1 f1("model_inv",model_invert,-200,200,2);
    eqline.SetNpx(10000);
    f1.SetNpx(10000);
    std::vector<double> pars = {0.0,1.};
    f1.SetParameters(pars.data());
    f1.SetLineColor(kGreen+2);
    f1.SetLineStyle(kLine);
    f1.Draw("L SAME");

    gPad->Modified(); gPad->Update();
    while(gPad->WaitPrimitive());

}*/
