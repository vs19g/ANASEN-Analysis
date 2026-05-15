#include <TF1.h>
/*double model(double* x, double* p) {
    double result = x[0];
    double factor = 29.0;
    double slope = 0.7;
    if(TMath::Abs(x[0]) < 16.2) result=x[0]*slope;
    else if(TMath::Abs(x[0]) < 49.8 ) result=x[0]*slope+TMath::Sign(1.0,x[0])*factor;
    else if(TMath::Abs(x[0]) < 85.2 ) result=x[0]*slope+TMath::Sign(1.0,x[0])*factor*2;
    else result=x[0]*slope+TMath::Sign(1.0,x[0])*factor*3;
    return result;
}

double model_invert(double *y, double *q) {
    double result=y[0];
    double slope = 0.7;
    double factor = 0.0;
    if(TMath::Abs(y[0]) < 16.2/slope) result = y[0]/slope;
    else if(TMath::Abs(y[0]) < 49.8/slope ) result=y[0]/slope-TMath::Sign(1.0,y[0])*factor;
    else if(TMath::Abs(y[0]) < 85.2/slope ) result=y[0]/slope-TMath::Sign(1.0,y[0])*factor*2;
    else result=y[0]/slope-TMath::Sign(1.0,y[0])*factor*3;
    return result;
}*/

double model_invert(double* y, double* p) {
    double result = y[0];
    double slope = 0.6;
	double z_grid[8] = {147.998,101.946,59.7634,19.6965,-19.6965,-59.7634,-101.946,-147.998};
 	for(int i=0;i<7;i++) {
 		if(y[0] <= z_grid[i] && y[0] > z_grid[i+1]) {
 			double zavg = (z_grid[i] + z_grid[i+1])*0.5; //midpoint about which we pivot
 			result = (y[0]-zavg)/slope + zavg;
 			break;
 		}
 	}
    return result+80;
}

double model_a1c1(double* x, double* p) {
    double result = x[0];
    double factor = 29.0;
    double slope = 0.0;
    if(TMath::Abs(x[0]) < 16.2) result=x[0]*slope;
    else if(TMath::Abs(x[0]) < 49.8 ) result=x[0]*slope+TMath::Sign(1.0,x[0])*factor;
    else if(TMath::Abs(x[0]) < 85.2 ) result=x[0]*slope+TMath::Sign(1.0,x[0])*factor*2;
    else result=x[0]*slope+TMath::Sign(1.0,x[0])*factor*3;
    return result;
}

double model_invert_a1c1(double *y, double *q) {
    double result=y[0];
/*    double slope = 1.0;
    double factor = 5.0;
    if(TMath::Abs(y[0]) < 16.2/slope) result = y[0]/slope;
    else if(TMath::Abs(y[0]) < 49.8/slope ) result=y[0]/slope-TMath::Sign(1.0,y[0])*factor;
    else if(TMath::Abs(y[0]) < 85.2/slope ) result=y[0]/slope-TMath::Sign(1.0,y[0])*factor*2;
    else result=y[0]/slope-TMath::Sign(1.0,y[0])*factor**/;
    return result+40;
}


/*void testmodel() {
    TF1 eqline("x","x",-200,200);
    eqline.Draw("");
    eqline.SetLineStyle(kDashed);

 //TF1 f1("model",model,-200,200,2);
    TF1 f1a("model_inv",model_a1c1,-200,200,2);
    eqline.SetNpx(10000);
    f1a.SetNpx(10000);
    std::vector<double> pars = {0.0,1.};
    f1a.SetParameters(pars.data());
    f1a.SetLineColor(kGreen+2);
    f1a.SetLineStyle(kLine);
    f1a.Draw("L SAME");

    TF1 f1("model",model,-200,200,2);
    //TF1 f1("model_inv",model_invert,-200,200,2);
    eqline.SetNpx(10000);
    f1.SetNpx(10000);
    //std::vector<double> pars = {0.0,1.};
    f1.SetParameters(pars.data());
    f1.SetLineColor(kGreen+2);
    f1.SetLineStyle(kLine);
    f1.Draw("L SAME");

    gPad->Modified(); gPad->Update();
    while(gPad->WaitPrimitive());
}*/
