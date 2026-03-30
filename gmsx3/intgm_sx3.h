#include "../Armory/HistPlotter.h"
#include <Minuit2/FCNBase.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TMath.h>
#include <TPad.h>
#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <TF1.h>
#include "func1.h"
static long iters=0;

//class intgm_sx3 : public ROOT::Minuit2::FCNBase {
class intgm_sx3 {
        int N;

        //L.at(0).at(3).at(n) is front strip = 0, back pad = 3, nth datapoint
        std::array<std::array<std::vector<double>,4>,4> L,R,B;
        //std::array<std::array<double,5>,4> stripedge; //stripedge.at(i).at(j) is the jth edge of the ith strip. there are five edges for the four strips 'i'=0 to 3, (0,1) (1,2) (2,3) (3,4) for each
        //the edges are at  -2a, -a, 0, a, 2a respectively so we enforce four ratios in the chi2 value - 'a' can be held constant, no need to fit it.

        //assume z = M*(aL-bR)
        //stripedge[i][1] = max(z) when pad==0 = min(z) when pad==1  this should be -1
        //stripedge[i][2] = max(z) when pad==1 = min(z) when pad==2. this should be 0
        //stripedge[i][3] = max(z) when pad==2 = min(z) when pad==3. this should be 1

        //i.e. stripedge[i][j] = max(z) when pad == j-1, min(z) when pad==j, for i= 1,2,3

        //ncounts.at(frontch).at(backch) is the number of (L,R,B) tuples we've filled (frontch,backch) coordinates in the detector
        std::array<std::array<long,4>,4> ncounts;
        TH1F *localhists[4][4]; //one histogram for each fc, bc combination
        HistPlotter *plotter;
        TF1 *pos_weight[4];
        TF1 *energywell;
    public:
        intgm_sx3() {
            for(int bc=0; bc<4; bc++) {
                for(int fc=0; fc<4; fc++) {
                    L[fc][bc].reserve(1000);
                    R[fc][bc].reserve(1000);
                    B[fc][bc].reserve(1000);
                    //localhists[fc][bc] = new TH1F(Form("h_%d_%d",fc,bc),Form("h_%d_%d",fc,bc),1000,-4.,4.);
                    ncounts[fc][bc] = 0;
                }
                pos_weight[bc] = new TF1(Form("b_strip%d",bc),model2,-10,10,4); //from -10, to 10, 4 parameters
                pos_weight[bc]->SetParameters(1.0,10,1.,3-2*bc); //centers at 1, 3.,5,7 Width 2a with a=1.0
                pos_weight[bc]->SetNpx(1'000'000);
            }
            energywell = new TF1("ewell",model2,0,2000,4); //0 to 2000 channels, 4 params
            energywell->SetParameters(1000,20,1,1500); //center the back E values at 1500 +/- 500
            energywell->SetNpx(1'000'000);
            N=0;
        }
        void set_plotter(HistPlotter *p) {plotter=p;}
        void set_iters(long i) { iters=i;}
        intgm_sx3(HistPlotter *p) : plotter(p) {
            for(int bc=0; bc<4; bc++) {
                for(int fc=0; fc<4; fc++) {
                    L[fc][bc].reserve(1000);
                    R[fc][bc].reserve(1000);
                    B[fc][bc].reserve(1000);
                    //localhists[fc][bc] = new TH1F(Form("h_%d_%d",fc,bc),Form("h_%d_%d",fc,bc),1000,-4.,4.);
                    ncounts[fc][bc] = 0;
                }
                pos_weight[bc] = new TF1(Form("b_strip%d",bc),model2,-10,10,4); //from -10, to 10, 4 parameters
                                                //a/2, k, v0, center
                pos_weight[bc]->SetParameters(0.92,10,1.,-1.*(3-2*bc)); //centers at 7, 5.,3,1 Width 2a with a=1.0
                pos_weight[bc]->SetNpx(1'000'000);
            }
            energywell = new TF1("ewell",model2,0,8000,4); //0 to 2000 channels, 4 params
//            energywell->SetParameters(60,10,0,1430); //center the back E values at 1430 +/- 60
            energywell->SetParameters(400,10,0,5246); //center the back E values at 5486 +/- 600
            energywell->SetNpx(1'000'000);
            N=0;
        }

        inline void fill(int fc, int bc, double leftE, double rightE, double backE) {
            /*
             *
             */
            assert(fc>=0 && fc<=3 && "Front channels should fit the range 0 to 3 inclusive!");
            assert(bc>=0 && bc<=3 && "Back channels should fit the range 0 to 3 inclusive!");
            if(leftE>0 && rightE >0 && backE>0) {
                L[fc][bc].emplace_back(leftE);
                R[fc][bc].emplace_back(rightE);
                B[fc][bc].emplace_back(backE);
                ncounts[fc][bc]+=1;
                N+=1;
            }
        }

        inline void print() {
            for(int i=0; i<16; i++) {
                std::cout << ncounts[i%4][i/4] << std::endl;
            }
        }
        inline void plot(std::string comment, const double* params) {
            std::array<double,4> l,r,b,bo,ro,lo,offset,stretch; //aliases to help with book-keeping
            std::array<std::array<double,4>,4> back_gains;// back_gains[fc][bc] are for fc,bc firing in combo
            for(int ctr=0; ctr<4; ctr++) {
                r[ctr] = params[ctr];
            }
            for(int ctr=4; ctr<20; ctr++) {
                int bch = (ctr-4)%4;
                int fch = (ctr-4)/4;
                back_gains[bch][fch] = params[ctr];
            }
            for(int ctr=20; ctr<24; ctr++) {
                stretch[ctr-20] = params[ctr];
            }
            for(int ctr=24; ctr<28; ctr++) {
                l[ctr-24] = params[ctr];
            }
            for(int fc=0; fc<4; fc++) {
                for(int bc=0; bc<4; bc++) {
                    for(int n=0; n<ncounts[fc][bc]; n++) {
                        if(plotter) {
                            double left = l[fc]*L[fc][bc].at(n);
                            double right = r[fc]*R[fc][bc].at(n);
                            double back = back_gains[bc][fc]*B[fc][bc].at(n);
                            //double zpos = (left - right)/(left+right);
                            double zpos = stretch[fc]*(left - right)/(left+right);// + offset[fc]; //back;
                            plotter->Fill2D(Form("normlf_fc%d_%d_%s",fc,bc,comment.c_str()),800,0,1.,800,0,1.,left/back, right/back,"l_vs_r");
                            plotter->Fill2D(Form("normlf_all_%s",comment.c_str()),800, 0, 1., 800, 0, 1.,left/back, right/back);
                            plotter->Fill2D(Form("case_f%d_b%d_%s",fc,bc,comment.c_str()),800,0,8192,800,0,8192,left+right,back,"l_vs_r");
                            plotter->Fill2D(Form("case_all_%s",comment.c_str()),800,0,8192,800,0,8192,left+right,back);
                            //plotter->Fill2D(Form("z_vs_backe_f%d_b%d_%s",fc,bc,comment.c_str()),800,-10,10,800,0,8192,zpos,back,"z_vs_be");
                            plotter->Fill2D(Form("z_vs_backe_all_%s",comment.c_str()),800,-10,10,800,0,8192,zpos,back);
                        } //end if plotter
                    }// end for-n
                }//end for-bc
            }//end for-fc
        }//end plot()

    //    double operator()(const std::vector<double>& params) const override{
        double eval(const double* params) const {
            iters+=1;

            std::array<double,4> l,r,b,bo,ro,lo, offset, stretch; //aliases to help with book-keeping
            std::array<std::array<double,4>,4> back_gains;// back_gains[fc][bc] are for fc,bc firing in combo
            for(int ctr=0; ctr<16; ctr++) {
                int bch = (ctr)%4;
                int fch = (ctr)/4;
                back_gains[bch][fch] = params[ctr];
            }
            for(int ctr=16; ctr<20; ctr++) {
                r[ctr-16] = params[ctr];
                l[ctr-16] = 1.0;
            }
            for(int ctr=20; ctr<24; ctr++) {
            	stretch[ctr-20] = params[ctr];
            }
            double result=0, sumcount=0;
            for(int fc=0; fc<4; fc++) {
                for(int bc=0; bc<4; bc++) {
                    //if(bc >= 1 || fc >= 1 ) continue;
                    if(ncounts[fc][bc] == 0 && iters ==0) {
                        std::cout << "Missing any data in front:" << fc << " back:" << bc << " combination." << std::endl;
                    }

                    for(int n=0; n<ncounts[fc][bc] ; n++) {
                        //double left = l[fc]*L[fc][bc].at(n) + lo[fc];
                        //double right = r[fc]*R[fc][bc].at(n) + ro[fc];
                        //double back = b[bc]*B[fc][bc].at(n) + bo[bc];
                        //double add = TMath::Power(left + right - back,2);
                        if(B[fc][bc].at(n)<100) continue;//ignore events too close to noise threshold

                        double left = l[fc]*L[fc][bc].at(n);
                        double right = r[fc]*R[fc][bc].at(n);
                        double back = back_gains[bc][fc]*B[fc][bc].at(n);
                        double lnorm = left/B[fc][bc].at(n);
                        double rnorm = right/B[fc][bc].at(n);

                        //double add = TMath::Power(left/back + right/back - 1.0,2);
                        double add = TMath::Power(left + right - back,2);
                        double zpos = stretch[fc]*(left - right)/(left+right); //back;
                        std::cout << zpos << " " << pos_weight[bc]->Eval(zpos) << " " << bc << std::endl;
                        double add_position = pos_weight[bc]->Eval(zpos);
                        double eback_align_penalty = energywell->Eval(back);
/*                      if(back>1000) zmid[fc][bc] += zpos;
                        if(back>1000 && zpos < zmin[fc][bc]) zmin[fc][bc] = zpos;
                        if(back> 1000 && zpos > zmax[fc][bc]) zmax[fc][bc] = zpos;

                        if(back>1000) {
                            localhists[fc][bc]->Fill(zpos);
                        }*/

                        result += add_position;
                        //result += add;
                        result += eback_align_penalty;
                        sumcount+=1;
                        //if(bc==0) std::cout << add << " " << add_position << " " <<  zpos << std::endl;
                        //To avoid drift towards (0,0,0) trivial solution. This value ~1 close to (1,1,1)
                        //result+=(1e-3/(TMath::Power(l[fc],2)+TMath::Power(r[fc],2)+TMath::Power(b[bc],2)+1e-9));
						//result+=(1e-3/(TMath::Power(l[fc],2)+TMath::Power(r[fc],2)+TMath::Power(b[bc],2)+1e-9));
                    } //end for-n
                } //end for-bc
            } //end for-fc
            result/=sumcount; //normalize, so the value doesn't scream
            if(iters%1'000==0) {
                std::cout << "iters : " << iters << " params: " << std::endl ;
                for(int i=0 ; i< 10; i++) std::cout << params[i] << " " << std::flush;
                std::cout<< std::endl;
                for(int i=10 ; i< 20; i++) std::cout << params[i] << " " << std::flush;
                std::cout << std::endl <<  " result: " << result << std::endl;
            } //end if
            return result;
        } //end eval()

        //double Up() const override { return 1.0; } // Required by minuit2 FCBase

};

