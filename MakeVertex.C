#define MakeVertex_cxx
//#define ADD_NEW_BRANCHES 0

Int_t colors[40] = {
	kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange,
	kSpring, kTeal, kAzure, kViolet, kPink, kGray, kWhite,
	kRed+2, kGreen+2, kBlue+2, kYellow+2, kMagenta+2, kCyan+2, kOrange+2,
	kSpring+2, kTeal+2, kAzure+2, kViolet+2, kPink+2,
	kRed-7, kGreen-7, kBlue-7, kYellow-7, kMagenta-7, kCyan-7, kOrange-7,
	kSpring-7, kTeal-7, kAzure-7, kViolet-7, kPink-7, kGray+2
};

#include "MakeVertex.h"
#include "Armory/ClassPW.h"
#include "Armory/HistPlotter.h"
#include "Armory/SX3Geom.h"
#include "Armory/PC_StepLadder_Correction.h"
#include "Armory/Kinematics.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TBranch.h>
#include <TGraph2D.h>
#include <TView.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TH3D.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <algorithm>

bool realtime = true;
bool process_alpha_proton_scattering = true;
const double qqq_z = 100.0;
const double z_entrance = -174.3 - 9.7 - 100.0;
const double anode_gain = 1.5146e-5; //channels --> MeV
bool Seven_MeV_Cut=false;

TF1 pcfix_func("func",model_invert,-200,200);
TF1 pcfix_func_a1c1("func_a1c1",model_invert_a1c1,-200,200);
TGraph *MeV_to_cm=NULL,*cm_to_MeV=NULL;
TGraph *cm_to_MeV_27Al=NULL,  *MeV_to_cm_27Al=NULL;
TGraph *MeV_to_cm_p=NULL,*cm_to_MeVp=NULL;
TGraph *cm_to_MeV_17F=NULL, *MeV_to_cm_17F=NULL;
TGraph *MeV_to_cm_d=NULL,*cm_to_MeVd=NULL;
TGraph *MeV_to_cm_t=NULL,*cm_to_MeVt=NULL;
TGraph *MeV_to_cm_3he=NULL, *cm_to_MeV3he=NULL;
TApplication *app=NULL;
TH1F *hha=NULL,*hhc=NULL;
TH3D *frame=NULL;
TCanvas *can1=NULL,*can2=NULL;

TPolyLine3D *pla[24]={NULL};
TPolyLine3D *plc[24]={NULL};
TPolyLine3D *qqqw[16][4]={NULL};
TPolyLine3D *trajectory=NULL;
TGraph2D *qqqg=NULL, *crossoverg=NULL, *guessg=NULL;


double z_to_crossover_rho(double z) {
//    return 9.20645e-5*z*z + 34.1973;
	return 0.000165896*z*z + 4.61626e-08*z + 32.067;

}

// Global instances
PW pwinstance;
TRandom3 rnd_qqq, rnd_sx3;
TVector3 hitPos;
double qqqenergy, qqqtimestamp;
class Event {
public:
	Event(TVector3 p, double e1, double e2, double t1, double t2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2) {}
	Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2) {}
	//Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2, int m1, int m2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2), multi1(m1), multi2(m2) {}

	TVector3 pos;
	int ch1=-1; //int(ch1/16) gives qqq id, ch1%16 gives ring#
	int ch2=-1; //int(ch2/16) gives qqq id, ch2%16 gives wedge#
	double Energy1=-1; //Front for QQQ, Anode for PC
	double Energy2=-1; //Back for QQQ, Cathode for PC
	double Time1=-1;
	double Time2=-1;

	//misc elements;
	int multi1=-1, multi2=-1;
};

// Calibration globals
const int MAX_QQQ = 4;
const int MAX_RING = 16;
const int MAX_WEDGE = 16;
double qqqGain[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqGainValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};
double qqqCalib[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{0}}};
bool qqqCalibValid[MAX_QQQ][MAX_RING][MAX_WEDGE] = {{{false}}};

double sx3BackGain[24][4][4] = {{{1.}}};
double sx3FrontGain[24][4] = {{1.}};
double sx3FrontOffset[24][4] = {{0.}};
double sx3RightGain[24][4] = {{1.}};

// PC Arrays
double pcSlope[48];
double pcIntercept[48];

HistPlotter *plotter;

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

void protonAlphaHistograms(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);

void paMiscHistograms(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);
void paMiscHistograms_oneWire(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters);

void ppMiscHistograms(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);
void ppMiscHistograms_oneWire(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters);
void ppMiscHistograms_sx3(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);

void miscHistograms_27Alaa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);
void miscHistograms_oneWire_27Alaa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters);
void miscHistograms_27Alaa_sx3(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events, std::string globaltag);

void miscHistograms_17Faa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events, std::string);
void miscHistograms_oneWire_17Faa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters);
void miscHistograms_17Faa_sx3(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events);

void MakeVertex::Begin(TTree * /*tree*/)
{
	pcfix_func.SetNpx(100000);
	TString option = GetOption();
	if(option!="")
	  plotter = new HistPlotter(option.Data(),"TFILE");
	else
	  plotter = new HistPlotter("Analyzer_SX3.root", "TFILE");
	plotter->ReadCuts("cutlist.txt");
	// ---------------------------------------------------------
	// 1. CRITICAL FIX: Initialize PC Arrays to Default (Raw)
	// ---------------------------------------------------------
	for (int i = 0; i < 48; i++)
	{
		pcSlope[i] = 1.0;     // Default slope = 1 (preserves Raw energy)
		pcIntercept[i] = 0.0; // Default intercept = 0
	}
	rnd_qqq.SetSeed(0);
	rnd_sx3.SetSeed(0);
	if(getenv("DATASET")) 
		dataset = std::string(getenv("DATASET"));
	if(getenv("source_vertex"))
		source_vertex = (double)std::atof(std::string(getenv("source_vertex")).c_str());
	std::cout << "Dataset set to  " << dataset << std::endl;
	std::cout << "source_vertex set to  " << source_vertex << std::endl;
	if(getenv("reactiondata")) {
		reactiondata = std::atoi(getenv("reactiondata"));
		std::cout << "Analyzing dataset as reactiondata" << std::endl;
	}

	int aoffset = 0;
	int coffset = 0;
	if(getenv("anode_offset")) {
		aoffset = std::atoi(getenv("anode_offset"));
		std::cout << "Offseting anodes by " << aoffset << " wires." << std::endl;
	}
	if(getenv("cathode_offset")) {
		coffset = std::atoi(getenv("cathode_offset"));
		std::cout << "Offseting cathodes by " << coffset << " wires." << std::endl;
	}
	pwinstance.ConstructGeo(aoffset,coffset);
	//pwinstance.PrintGeometry();


	fflush(stdout);
	//usleep(4e5);
	// Load PC Calibrations
	std::ifstream inputFile("slope_intercept_results_"+dataset+".dat");
	if (inputFile.is_open())
	{
		std::string line;
		int index;
		double slope, intercept;
		while (std::getline(inputFile, line))
		{
			std::stringstream ss(line);
			ss >> index >> slope >> intercept;
			if (index >= 0 && index <= 47)
			{
				pcSlope[index] = slope;
				pcIntercept[index] = intercept;
			}
		}
		inputFile.close();
	}
	else
	{
		std::cerr << "Error opening slope_intercept.dat" << std::endl;
	}

	// ... (Load QQQ Gains and Calibs - same as before) ...
	{
		std::string filename = "qqq_GainMatch.dat";
		std::ifstream infile(filename);
		if (infile.is_open())
		{
			int det, ring, wedge;
			double gainw, gainr;
			while (infile >> det >> wedge >> ring >> gainw >> gainr)
			{
				qqqGain[det][wedge][ring] = gainw;
				qqqGainValid[det][wedge][ring] = (gainw > 0);
				// std::cout << "QQQ Gain Loaded: Det " << det << " Ring " << ring << " Wedge " << wedge << " GainW " << gainw << " GainR " << gainr << std::endl;
			}
			infile.close();
		}
	}

	{
		std::string filename = "qqq_Calib.dat";
		std::ifstream infile(filename);
		if (infile.is_open())
		{
			int det, ring, wedge;
			double slope;
			while (infile >> det >> wedge >> ring >> slope)
			{
				qqqCalib[det][wedge][ring] = slope;
				qqqCalibValid[det][wedge][ring] = (slope > 0);
				// std::cout << "QQQ Calib Loaded: Det " << det << " Ring " << ring << " Wedge " << wedge << " Slope " << slope << std::endl;
			}
			infile.close();
		}
	}

	{
		std::ifstream infile("sx3cal/"+dataset+"/backgains.dat");
		std::string temp;
		int backpos, frontpos, clkpos;
		if (infile.is_open())
			while(infile>>clkpos>>temp>>frontpos>>temp>>backpos>>sx3BackGain[clkpos][frontpos][backpos])
				;//std::cout << sx3BackGain[clkpos][frontpos][backpos] << std::endl;
		infile.close();

		infile.open("sx3cal/"+dataset+"/frontgains.dat");
		if (infile.is_open())
			while(infile>>clkpos>>temp>>temp>>frontpos>>sx3FrontOffset[clkpos][frontpos]>>sx3FrontGain[clkpos][frontpos])
				;//std::cout << sx3FrontOffset[clkpos][frontpos] << " " << sx3FrontGain[clkpos][frontpos] << std::endl;
		infile.close();

		infile.open("sx3cal/"+dataset+"/rightgains.dat");
		if (infile.is_open())
			while(infile>>clkpos>>frontpos>>temp>>sx3RightGain[clkpos][frontpos]) {
				sx3RightGain[clkpos][frontpos]=TMath::Abs(sx3RightGain[clkpos][frontpos]);
			}
		infile.close();
	}
//    MeV_to_cm = new TGraph("eloss_calculations/alphas_in_250torr_mix_filtered_6MeV.txt","%lf %*lf %lf");
	MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_40MeV_3percent.dat","%lf %*lf %lf");
	cm_to_MeV=  new TGraph(MeV_to_cm->GetN(), MeV_to_cm->GetY(), MeV_to_cm->GetX());

	MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV_3percent.dat","%lf %*lf %lf");
	cm_to_MeVp=  new TGraph(MeV_to_cm_p->GetN(), MeV_to_cm_p->GetY(), MeV_to_cm_p->GetX());

	MeV_to_cm_d = new TGraph("eloss_calculations/deuteron_lookup_20MeV_3percent.dat","%lf %*lf %lf");
	cm_to_MeVd=  new TGraph(MeV_to_cm_d->GetN(), MeV_to_cm_d->GetY(), MeV_to_cm_d->GetX());

	MeV_to_cm_t = new TGraph("eloss_calculations/triton_lookup_20MeV_3percent.dat","%lf %*lf %lf");
	cm_to_MeVt=  new TGraph(MeV_to_cm_t->GetN(), MeV_to_cm_t->GetY(), MeV_to_cm_t->GetX());

	MeV_to_cm_3he = new TGraph("eloss_calculations/3he_lookup_40MeV_3percent.dat","%lf %*lf %lf");
	cm_to_MeV3he=  new TGraph(MeV_to_cm_3he->GetN(), MeV_to_cm_3he->GetY(), MeV_to_cm_3he->GetX());

	MeV_to_cm_27Al = new TGraph("eloss_calculations/aluminum_lookup_80MeV.dat", "%lf %*lf %lf");
	cm_to_MeV_27Al = new TGraph(MeV_to_cm_27Al->GetN(), MeV_to_cm_27Al->GetY(), MeV_to_cm_27Al->GetX());

	MeV_to_cm_17F = new TGraph("eloss_calculations/fluorine_lookup_70MeV.dat", "%lf %*lf %lf");
	cm_to_MeV_17F = new TGraph(MeV_to_cm_17F->GetN(), MeV_to_cm_17F->GetY(), MeV_to_cm_17F->GetX());

	//cm_to_MeV.Eval(MeV_to_cm.Eval(detectedE)-PathLength) gives energy of particle before it traversed 'path length'

}

Bool_t MakeVertex::Process(Long64_t entry)
{
	hitPos.Clear();
	qqqenergy = -1;
	qqqtimestamp=-1;
	HitNonZero = false;
	bool qqq1000cut = false;
	b_sx3Multi->GetEntry(entry);
	b_sx3ID->GetEntry(entry);
	b_sx3Ch->GetEntry(entry);
	b_sx3E->GetEntry(entry);
	b_sx3T->GetEntry(entry);
	b_qqqMulti->GetEntry(entry);
	b_qqqID->GetEntry(entry);
	b_qqqCh->GetEntry(entry);
	b_qqqE->GetEntry(entry);
	b_qqqT->GetEntry(entry);
	b_pcMulti->GetEntry(entry);
	b_pcID->GetEntry(entry);
	b_pcCh->GetEntry(entry);
	b_pcE->GetEntry(entry);
	b_pcT->GetEntry(entry);
#ifdef ADD_NEW_BRANCHES
		//std:cout << " aaa" << std::endl;
		b_miscMulti->GetEntry(entry);
		b_miscID->GetEntry(entry);
		b_miscCh->GetEntry(entry);
		b_miscE->GetEntry(entry);
		b_miscT->GetEntry(entry);
		b_miscTf->GetEntry(entry);
#endif
	double timecut_low = getenv("timecut_low")?std::atof(getenv("timecut_low")):0;
	double timecut_high = getenv("timecut_high")?std::atof(getenv("timecut_high")):1e15  ;

	if( pc.multi>0) {
		for(int i=0; i<pc.multi; i++) {
			if(pc.t[i]*1e-9 < timecut_high && pc.t[i]*1e-9 >= timecut_low) {
				//good, keep it moving
			} else {
				return kTRUE;
			}
		}
	}

	sx3.CalIndex();
	qqq.CalIndex();
	pc.CalIndex();

	std::vector<Event> SX3_Events;
	if(sx3.multi>1) {
		std::array<sx3det,24> Fsx3;
		//std::cout << "-----" << std::endl;
		bool found_upstream_sx3=0;
		for(int i=0; i<sx3.multi; i++) {
			int id = sx3.id[i];
			if(id>=12) continue;
			if(sx3.ch[i]>=8) {
				int sx3ch=sx3.ch[i]-8;
				sx3ch=(sx3ch+3)%4;
				if(id>=12) {
					found_upstream_sx3=1;
					//std::cout << Form("f%d(",id) << sx3ch << "," << sx3.e[i] << ") " << std::flush;
				}
				//if(sx3ch==0 || sx3ch==3) continue;
				double value=sx3.e[i];
				int gch = sx3.id[i]*4+(sx3.ch[i]-8);
				if(id<12)  Fsx3.at(id).fillevent("BACK",sx3ch,value);
				Fsx3.at(id).ts = static_cast<double>(sx3.t[i])+(rnd_sx3.Uniform(16.0)-8.0);
				plotter->Fill2D("sx3backs_all_raw",100,0,100,800,0,4096,gch,sx3.e[i]);
			} else {
				int sx3ch=sx3.ch[i]/2;
				double value=sx3.e[i];
				if(id>=12) {
					found_upstream_sx3=1;
					//std::cout << Form("b%d(",id) << sx3ch << "," << value << ") " << std::flush;
				}
				if(sx3.ch[i]%2==0) {
					Fsx3.at(id).fillevent("FRONT_L",sx3ch,value*sx3RightGain[id][sx3ch]);
				} else {
					Fsx3.at(id).fillevent("FRONT_R",sx3ch,value);
				}
			}
		} //end for (i in sx3.multi)
		//if(found_upstream_sx3) std::cout << std::endl;

		for(int id=0; id<24; id++) {
			//std::cout << id << " " << Fsx3.at(id).valid_front_chans.size() << " " << Fsx3.at(id).valid_back_chans.size() << std::endl;;
			try {
				Fsx3.at(id).validate();
			} catch(std::exception exc) {
				std::cout << "oops! anyway " << std::endl;
				continue;
			}
			auto det = Fsx3.at(id);
			bool no_charge_sharing_strict = det.valid_front_chans.size()==1 && det.valid_back_chans.size()==1;
			if(det.valid) {
				//std::cout << det.frontEL << " " << det.frontEL*sx3RightGain[id][det.stripF] << std::endl;
				//plotter->Fill2D("be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_b"+std::to_string(det.stripB),200,-1,1,800,0,8192,det.frontX,det.backE,"evsx");

				plotter->Fill2D("matched_be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF),200,-60,60,800,0,8192,
								det.frontX*sx3FrontGain[id][det.stripF]+sx3FrontOffset[id][det.stripF],det.backE*sx3BackGain[id][det.stripF][det.stripB],"evsx_matched");
				//plotter->Fill2D("fe_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_"+std::to_string(det.stripB),200,-1,1,800,0,4096,det.frontX,det.backE,"evsx");
				//plotter->Fill2D("l_vs_r_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF),800,0,4096,800,0,4096,det.frontEL,det.frontER,"l_vs_r");
			}
			if(det.valid && (id ==9 || id==7 || id == 1 || id==3) && det.stripF!=DEFAULT_NULL && det.stripB!=DEFAULT_NULL) {
				double z = det.frontX*sx3FrontGain[id][det.stripF]+sx3FrontOffset[id][det.stripF];
				z = z+(75.0/2.0)-3.0;//convert local sx3z to detector global coordinate system as indicated by measurements. 
				// Note that this will be different for the upstream barrel, when it gets implemented
				double backE = det.backE*sx3BackGain[id][det.stripF][det.stripB]*5.2/5.0;
				//if(backE<2000) continue;
				//det.stripF=3-det.stripF;
				if(id == 9 && backE< 2000) continue;
				if(Seven_MeV_Cut && backE<5000) continue;
				double alpha_n = TMath::ATan2((2*(3-det.stripF)-3)*40.30, 8.0*88.0*TMath::Cos(15.0*M_PI/180.0))*180./M_PI; //angle subtended w.r.t the radial perpendicular bisector of each sx3
				double beta_n = 15.0+alpha_n; //how much to add per strip to the starting position? this is the angle w.r.t an edge of the sx3, the above values run as (-10.08deg, -3.39deg, 3.39deg, 10.08deg)
				double phi_n = ((-id+0.5)*30+beta_n);
				phi_n+=45;
				double rho_at_strip = 88.0/TMath::Cos(alpha_n*M_PI/180.0); //TMath::Cos(15.0*M_PI/180.0) if the edge-length is 88mm
				phi_n*=M_PI/180.; //starting-position phi + strip contribution
				//Event sx3ev(TVector3(88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n),z),backE*0.001,-1,det.ts,-1,det.stripB+4*id,det.stripF+4*id);
				Event sx3ev(TVector3(rho_at_strip*TMath::Cos(phi_n),rho_at_strip*TMath::Sin(phi_n),z),backE*0.001,-1,det.ts,-1,det.stripB+4*id,det.stripF+4*id);
				SX3_Events.push_back(sx3ev);
				plotter->Fill2D("sx3backs_gm",100,0,100,800,0,8192,det.stripB+4*id,backE);
				plotter->Fill2D("SX3CartesianPlot" + std::to_string(id), 200, -100, 100, 200, -100, 100, 88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n), "hCalSX3");
			}
			if(det.valid && det.stripF!=DEFAULT_NULL && det.stripB!=DEFAULT_NULL) {
				plotter->Fill2D("sx3backs_raw",100,0,100,800,0,8192,det.stripB+4*id,det.backE);
			}

		}
	}
	//return kTRUE;
	// QQQ Processing

	int qqqCount = 0;
	int qqqAdjCh = 0;
	// REMOVE WHEN RERUNNING USING THE NEW CALIBRATION FILE
	std::vector<Event> QQQ_Events, PC_Events,PC_Events_OnlyAnode, PC_Events_OnlyCathode;
	std::vector<Event> QQQ_Events_Raw, PC_Events_Raw;
	std::vector<Event> QQQ_Events2; //clustering done

	std::unordered_map<int,std::tuple<int,int,double,double>> qvecr[4], qvecw[4];
	if(qqq.multi>1) {
		//if(qqq.multi>=3) std::cout << "-----" << std::endl;
		for(int i=0; i<qqq.multi; i++) {
			if(qqq.ch[i]/16) {
				if(qvecr[qqq.id[i]].find(qqq.ch[i])!=qvecr[qqq.id[i]].end()) std::cout << "mayday!" << std::endl;
				qvecr[qqq.id[i]][qqq.ch[i]] = std::tuple(qqq.id[i],qqq.ch[i],qqq.e[i],qqq.t[i]);
			} else {
				if(qvecw[qqq.id[i]].find(qqq.ch[i])!=qvecw[qqq.id[i]].end()) std::cout << "mayday!" << std::endl;
				qvecw[qqq.id[i]][qqq.ch[i]] = std::tuple(qqq.id[i],qqq.ch[i],qqq.e[i],qqq.t[i]);
			}
		}
	}

	bool PCSX3TimeCut = false;
	bool PCASX3TimeCut = false;
	bool PCCSX3TimeCut = false;
	
	bool PCQQQTimeCut = false;
	bool PCAQQQTimeCut = false;
	bool PCCQQQTimeCut = false;
	for (int i = 0; i < qqq.multi; i++) {
		plotter->Fill2D("QQQ_Index_Vs_Energy", 16 * 8, 0, 16 * 8, 2000, 0, 16000, qqq.index[i], qqq.e[i], "hRawQQQ");

		for (int j = 0; j < qqq.multi; j++) {
			if (j == i)
				continue;
			plotter->Fill2D("QQQ_Coincidence_Matrix", 16 * 8, 0, 16 * 8, 16 * 8, 0, 16 * 8, qqq.index[i], qqq.index[j], "hRawQQQ");
		}

		for (int k = 0; k < pc.multi; k++) {
			if (pc.index[k] < 24 && pc.e[k] > 10) {
				plotter->Fill2D("QQQ_Vs_Anode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
				plotter->Fill2D("QQQ_Vs_PC_Index", 16 * 8, 0, 16 * 8, 24, 0, 24, qqq.index[i], pc.index[k], "hRawQQQ");
			}
			else if (pc.index[k] >= 24 && pc.e[k] > 10) {
				plotter->Fill2D("QQQ_Vs_Cathode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
			}
		}

		for (int j = i + 1; j < qqq.multi; j++) {
			if (qqq.id[i] == qqq.id[j]) {
				qqqCount++;
				int chWedge = -1;
				int chRing = -1;
				double eWedge = 0.0;
				double eWedgeMeV = 0.0;
				double eRing = 0.0;
				double eRingMeV = 0.0;
				double tRing = 0.0;
				double tWedge = 0.0;

				if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16]) {
					chWedge = qqq.ch[i];
					eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
					chRing = qqq.ch[j] - 16;
					eRing = qqq.e[j];
					tRing = static_cast<double>(qqq.t[j])+(rnd_qqq.Uniform(16.0)-8.0);
					tWedge = static_cast<double>(qqq.t[i])+(rnd_qqq.Uniform(16.0)-8.0);
				}
				else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16]) {
					chWedge = qqq.ch[j];
					eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
					chRing = qqq.ch[i] - 16;
					eRing = qqq.e[i];
					tRing = static_cast<double>(qqq.t[i])+(rnd_qqq.Uniform(16.0)-8.0);
					tWedge = static_cast<double>(qqq.t[j])+(rnd_qqq.Uniform(16.0)-8.0);
				}
				else
					continue;

				plotter->Fill1D("Wedgetime_Vs_Ringtime", 100, -1000, 1000, tWedge - tRing, "hTiming");
				plotter->Fill2D("RingE_vs_Index", 16 * 4, 0, 16 * 4, 1000, 0, 16000, chRing + qqq.id[i] * 16, eRing, "hRawQQQ");
				plotter->Fill2D("WedgeE_vs_Index", 16 * 4, 0, 16 * 4, 1000, 0, 16000, chWedge + qqq.id[i] * 16, eWedge, "hRawQQQ");
				plotter->Fill2D("WedgeE_Vs_RingECal", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");
				if(chWedge + qqq.id[i] * 16 == 49 || chWedge + qqq.id[i] * 16 == 48 ) continue;
				if(chRing + qqq.id[i] * 16 == 63) continue;

				if (qqqCalibValid[qqq.id[i]][chWedge][chRing]) {
					eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;
					eRingMeV = eRing * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;

					if(eRingMeV/eWedgeMeV > 3.0 || eRingMeV/eWedgeMeV<1.0/3.0) continue;
					//if(eRingMeV<1.2 || eWedgeMeV<1.2) continue;

					//double theta = 2 * TMath::Pi() * (-qqq.id[i] * 16 + (15-chWedge) + 0.5)/(16*4);
					double theta = (M_PI/180.)*(-90*qqq.id[i]+(87./16.)*((15-chWedge)+0.5)+3.0);
					double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
					
					
					//z used to be 75+30+23=128
					//we found a 12mm shift towards the vertex later --> 116
					Event qqqevent(TVector3(rho*TMath::Cos(theta),rho*TMath::Sin(theta),qqq_z), eRingMeV, eWedgeMeV, tRing, tWedge,chRing+qqq.id[i]*16, chWedge+qqq.id[i]*16);
					Event qqqeventr(TVector3(rho*TMath::Cos(theta),rho*TMath::Sin(theta),qqq_z), eRing, eWedge, tRing, tWedge,chRing+qqq.id[i]*16, chWedge+qqq.id[i]*16);
					assert(qqq.id[i]>=0);
					if(Seven_MeV_Cut &&(eRingMeV<6.6 || eWedgeMeV< 6.6)) continue;                   
					QQQ_Events.push_back(qqqevent);
					QQQ_Events_Raw.push_back(qqqeventr);
					plotter->Fill2D("WedgeE_Vs_RingECal_selected", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");                    
					plotter->Fill2D("QQQCartesianPlot", 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
					plotter->Fill2D("QQQCartesianPlot" + std::to_string(qqq.id[i]), 200, -100, 100, 200, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hCalQQQ");
					plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(theta), rho * TMath::Sin(theta), "hPCQQQ");
				}
				else
					continue;


				for (int k = 0; k < pc.multi; k++)
				{
					plotter->Fill2D("RingCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
					plotter->Fill2D("WedgeCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
					plotter->Fill2D("WedgeCh_vs_Anode_Index" + std::to_string(qqq.id[i]), 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k]);
					plotter->Fill2D("RingCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
					plotter->Fill2D("WedgeCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");

					if (pc.index[k] < 24 && pc.e[k] > 10)
					{
						plotter->Fill2D("Timing_Difference_QQQ_PC", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
						plotter->Fill2D("DelT_Vs_QQQRingECal", 500, -2000, 2000, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hTiming");
						//if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
						if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
						{
							PCAQQQTimeCut = true;
							plotter->Fill2D("CalibratedQQQEvsPCE_R", 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
							plotter->Fill2D("CalibratedQQQEvsPCE_W", 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
						}
					}

					if (pc.index[k] >= 24 && pc.e[k] > 10) {
						if (tRing - static_cast<double>(pc.t[k]) < -200) PCCQQQTimeCut = true;
						//if (tRing - static_cast<double>(pc.t[k]) > 200) PCCQQQTimeCut = true;
						plotter->Fill2D("Timing_Difference_QQQ_PC_Cathode", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
					}
				} //end of pc k loop

				if (!HitNonZero) {
					//double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
					//double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
					//double theta = 2 * TMath::Pi() * (-qqq.id[i] * 16 + (15-chWedge) + 0.5)/(16*4);
					double theta = -90.+90*qqq.id[i]+(87./16.)*((15-chWedge)+0.5)+3.0;
					theta *= (M_PI/180.);
					double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
					double x = rho * TMath::Cos(theta);
					double y = rho * TMath::Sin(theta);
					hitPos.SetXYZ(x, y, qqq_z);
					qqqenergy = eRingMeV;
					qqqtimestamp = tRing;
					HitNonZero = true;
				}
			} // if j==i
		} //j loop end
	} //i loop end

	PCQQQTimeCut = PCAQQQTimeCut && PCCQQQTimeCut;
	plotter->Fill1D("QQQ_Multiplicity", 10, 0, 10, qqqCount, "hRawQQQ");

	typedef std::unordered_map<int,std::tuple<int,double,double>> WireEvent; //this stores nearest neighbour wire events, or a 'cluster'
	WireEvent aWireEvents, cWireEvents; //naming for book keeping

	aWireEvents.clear();
	aWireEvents.reserve(24);
	cWireEvents.clear();
	cWireEvents.reserve(24);

	// PC Gain Matching and Filling
	double anodeT = -99999;
	double cathodeT = 99999;
	int anodeIndex = -1;
	int cathodeIndex = -1;
	for (int i = 0; i < pc.multi; i++)
	{
		//std::cout << pc.index[i] << " " << pc.e[i] << " " << std::endl;
		if (pc.e[i] > 25)
		{
			plotter->Fill2D("PC_Index_Vs_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], static_cast<double>(pc.e[i]), "hRawPC");
		} else {
			continue;
		}

		if (pc.index[i] < 48)
		{
			pc.e[i] = pcSlope[pc.index[i]] * pc.e[i] + pcIntercept[pc.index[i]];
			plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], pc.e[i], "hGMPC");
		}

		if (pc.index[i] < 24)
		{
			anodeT = static_cast<double>(pc.t[i]);
			anodeIndex = pc.index[i];
			aWireEvents[pc.index[i]] = std::tuple(pc.index[i],pc.e[i],static_cast<double>(pc.t[i]));
		}
		else
		{
			cathodeT = static_cast<double>(pc.t[i]);
			cathodeIndex = pc.index[i] - 24;
			cWireEvents[pc.index[i]-24] = std::tuple(pc.index[i]-24,pc.e[i],static_cast<double>(pc.t[i]));
		}

		if (anodeT != -99999 && cathodeT != 99999)
		{
			for (int j = 0; j < qqq.multi; j++)
			{
				plotter->Fill1D("PC_Time_qqq", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
				plotter->Fill2D("PC_Time_Vs_QQQ_ch", 200, -2000, 2000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, qqq.ch[j], "hTiming");
				plotter->Fill2D("PC_Time_vs_AIndex", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
				plotter->Fill2D("PC_Time_vs_CIndex", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
				// plotter->Fill1D("PC_Time_A" + std::to_string(anodeIndex) + "_C" + std::to_string(cathodeIndex), 200, -1000, 1000, anodeT - cathodeT, "TimingPC");
			}

			for (int j = 0; j < sx3.multi; j++)
			{
				plotter->Fill1D("PC_Time_sx3", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
			}
			for(auto sx3event : SX3_Events) {
					bool TCC = sx3event.Time1 - cathodeT   < 0;
					bool TCA = sx3event.Time1 - anodeT   < 0;
					//plotter->Fill2D("sx3_z_phi_awire"+std::to_string(anodeIndex)+"_TC"+std::to_string(TCA), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
					//plotter->Fill2D("sx3_z_phi_cwire"+std::to_string(cathodeIndex)+"_TC"+std::to_string(TCC), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );	    		
			}

			plotter->Fill1D("PC_Time", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
		}

		for (int j = i + 1; j < pc.multi; j++)
		{
			plotter->Fill2D("PC_Coincidence_Matrix", 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
			plotter->Fill2D("PC_Coincidence_Matrix_anodeMinusCathode_lt_-200_" + std::to_string(anodeT - cathodeT < -200), 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
			plotter->Fill2D("Anode_V_Anode", 24, 0, 24, 24, 0, 24, pc.index[i], pc.index[j], "hGMPC");
		}
	}
	anodeHits.clear();
	cathodeHits.clear();
	corrcatMax.clear();

	int aID = 0;
	int cID = 0;
	double aE = 0;
	double cE = 0;
	double aESum = 0;
	double cESum = 0;
	double aEMax = 0;
	double cEMax = 0;
	int aIDMax = 0;
	int cIDMax = 0;

	for (int i = 0; i < pc.multi; i++) {
		// if (pc.e[i] > 100)
		{
			if (pc.index[i] < 24) {
				anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
			}
			else if (pc.index[i] >= 24) {
				cathodeHits.push_back(std::pair<int, double>(pc.index[i] - 24, pc.e[i]));
			}
		}
	}

	std::sort(anodeHits.begin(),anodeHits.end(),[](std::pair<int,double> a, std::pair<int,double> b) {
		return a.first < b.first;
	});

	std::sort(cathodeHits.begin(),cathodeHits.end(),[](std::pair<int,double> a, std::pair<int,double> b) {
		return a.first < b.first;
	});

	//clusters = collection of (collection of wires) where each wire is (index, energy, timestamp)
	std::vector<std::vector<std::tuple<int,double,double>>> aClusters = pwinstance.Make_Clusters(aWireEvents);
	std::vector<std::vector<std::tuple<int,double,double>>> cClusters = pwinstance.Make_Clusters(cWireEvents);

	std::vector<std::pair<double,double>> sumE_AC;
	int actr = 0;
	for(auto aCluster: aClusters) {
		int cctr = 0;
		for(auto cCluster: cClusters) {
			//both have at least 1, here. Keep the a1, c1 events
			auto [crossover,alpha,apSumE,cpSumE,apMaxE,cpMaxE,apTSMaxE,cpTSMaxE] = pwinstance.FindCrossoverProperties(aCluster, cCluster);
			if(alpha!=9999999 && apSumE!=-1 && aCluster.size() && cCluster.size()) { //needs both anodes and cathodes, AND for the crossover to fall in [-173.6, 173.6]
				//Event PCEvent(crossover,apMaxE,cpMaxE,apTSMaxE,cpTSMaxE);
				//Event PCEvent(crossover,apSumE,cpSumE,apTSMaxE,cpTSMaxE);
				Event PCEvent(crossover,apSumE,cpMaxE,apTSMaxE,cpTSMaxE); //run12 shows cathode-max and anode-sum provide best dE signals.
				//std::cout << apSumE << " " << crossover.Perp() << " " << apMaxE << " " << apTSMaxE << std::endl;
				PCEvent.multi1=aCluster.size();
				PCEvent.multi2=cCluster.size();
				PC_Events.push_back(PCEvent);
				sumE_AC.push_back(std::pair(apSumE,cpSumE));
			} 
			if(cCluster.size()!=0 && actr==0) {
					Event PCEvent_OnlyC(crossover,-1,cpMaxE,-1,cpTSMaxE); //run12 shows cathode-max and anode-sum provide best dE signals.
					PCEvent_OnlyC.multi1=0;
					PCEvent_OnlyC.multi2=cCluster.size();
					PC_Events_OnlyCathode.push_back(PCEvent_OnlyC);            	
			}
			
			if(aCluster.size()!=0 && cctr==0) { //avoid 'libeio' alternativesdouble-counting
				Event PCEvent_OnlyA(crossover,apSumE,-1,apTSMaxE,-1); //run12 shows cathode-max and anode-sum provide best dE signals.
				PCEvent_OnlyA.multi1=aCluster.size();
				PCEvent_OnlyA.multi2=0;
				PC_Events_OnlyAnode.push_back(PCEvent_OnlyA);
			}
		 	cctr++;   
		}
		actr++;
	}
	bool is_fluorine=false;
#ifdef ADD_NEW_BRANCHES
	TRandom3 rnd;
	rnd.SetSeed();//random seed set
	if(dataset=="17F" && reactiondata) {
		int ctr=0;
		for(auto qqqevent: QQQ_Events) {
			double ts_rf = -987654321;
			double ts_needle = -987654321;
			double ts_mcp =-987654321;
			double needle_energy = -987654321;
			double ts_qqq = static_cast<double>(qqqevent.Time1) + (rnd.Uniform(16.0)-8.0);
			bool found_rf=false;
			bool found_mcp=false;
			bool found_needle=false;
			for(int j=0; j< misc.multi; j++) {
				plotter->Fill1D("channels_misc",20,0,20,misc.ch[j],"misc");
				if(misc.ch[j]==2) { //Needle 
					plotter->Fill2D("needle_vs_qqqE",800,0,16384,800,0,10,misc.e[j],qqqevent.Energy1,"misc");
					ts_needle = static_cast<double>(misc.t[j])+static_cast<double>(misc.tf[j]);
					needle_energy = static_cast<double>(misc.e[j]);
					found_needle=1;
					plotter->Fill1D("dt_qqq_needle",800,-2000,2000,ts_qqq-ts_needle,"misc");
				}
				if(misc.ch[j]==3) { //RF 
					ts_rf = static_cast<double>(misc.t[j])+static_cast<double>(misc.tf[j]);
					found_rf=1;
					plotter->Fill1D("dt_qqq_rf",800,-2000,2000,ts_qqq-ts_rf,"misc");
				}
				if(misc.ch[j]==4) { //mcp
					ts_mcp = static_cast<double>(misc.t[j])+static_cast<double>(misc.tf[j]);
					found_mcp=1;
					plotter->Fill1D("dt_qqq_mcp",800,-2000,2000,ts_qqq-ts_mcp,"misc");
				}
			}
			if(found_rf && found_mcp) {
				if(ctr==0) plotter->Fill1D("dt_rf_mcp",500,-1000,1000,ts_rf-ts_mcp,"misc");
				double dt_rf_mcp = ts_rf - ts_mcp;
				double dt_qqq_rf = ts_qqq-ts_rf;
				double dt_qqq_mcp = ts_qqq-ts_mcp;
				plotter->Fill2D("dt(qqq,rf)_vs_(rf,mcp)",640,-2000,2000,640,-2000,2000,dt_qqq_rf,dt_rf_mcp,"misc");
				plotter->Fill2D("dt_(qqq,mcp)_vs_(qqq,rf)",640,-1400,2000,640,-2000,2000,dt_qqq_mcp,dt_qqq_rf,"misc");
				plotter->Fill2D("dt_(qqq,mcp)_vs_(rf,mcp)",640,-1400,-600,640,-2000,2000,dt_qqq_mcp,dt_rf_mcp,"misc");
				if(dt_qqq_mcp > -1160 && dt_qqq_mcp < -960 && dt_rf_mcp > 60 && dt_rf_mcp < 160) {
					plotter->Fill2D("dt_(qqq,mcp)_vs_(rf,mcp)_beam1",640,-1400,-600,640,-2000,2000,dt_qqq_mcp,dt_rf_mcp,"misc");
					is_fluorine=true;
				}
				if(found_needle) {
					plotter->Fill2D("dt_(rf,mcp)_vs_needleE_beam"+std::to_string(is_fluorine),640,-2000,2000,800,0,16384,dt_rf_mcp,needle_energy,"misc");
					plotter->Fill2D("dt_(rf,mcp)_vs_dt(qqq,needle)_beam"+std::to_string(is_fluorine),640,-2000,2000,800,-2000,2000,dt_rf_mcp,ts_qqq-ts_needle,"misc");
				}
			}
			ctr+=1;
		}
		
		for(auto sx3event: SX3_Events) {
			double ts_rf = -987654321;
			double ts_needle = -987654321;
			double ts_mcp =-987654321;
			double needle_energy = -987654321;
			double ts_sx3 = static_cast<double>(sx3event.Time1) + (rnd.Uniform(16.0)-8.0);
			bool found_rf=false;
			bool found_mcp=false;
			bool found_needle=false;
			for(int j=0; j< misc.multi; j++) {
				plotter->Fill1D("channels_misc",20,0,20,misc.ch[j],"misc");
				if(misc.ch[j]==2) { //Needle 
					plotter->Fill2D("needle_vs_sx3E",800,0,16384,800,0,10,misc.e[j],sx3event.Energy1,"misc");
					ts_needle = static_cast<double>(misc'libeio' alternatives.t[j])+static_cast<double>(misc.tf[j]);
					needle_energy = static_cast<double>(misc.e[j]);
					found_needle=1;
					plotter->Fill1D("dt_sx3_needle",800,-2000,2000,ts_sx3-ts_needle,"misc");
				}
				if(misc.ch[j]==3) { //RF 
					ts_rf = static_cast<double>(misc.t[j])+static_cast<double>(misc.tf[j]);
					found_rf=1;
					plotter->Fill1D("dt_sx3_rf",800,-2000,2000,ts_sx3-ts_rf,"misc");
				}
				if(misc.ch[j]==4) { //mcp
					ts_mcp = static_cast<double>(misc.t[j])+static_cast<double>(misc.tf[j]);
					found_mcp=1;
					plotter->Fill1D("dt_sx3_mcp",800,-2000,2000,ts_sx3-ts_mcp,"misc");
				}
			}
			if(found_rf && found_mcp) {
				if(ctr==0) plotter->Fill1D("dt_rf_mcp_sx3",500,-1000,1000,ts_rf-ts_mcp,"misc");
				double dt_rf_mcp = ts_rf - ts_mcp;
				double dt_sx3_rf = ts_sx3-ts_rf;
				double dt_sx3_mcp = ts_sx3-ts_mcp;
				plotter->Fill2D("dt(sx3,rf)_vs_(rf,mcp)",640,-2000,2000,640,-2000,2000,dt_sx3_rf,dt_rf_mcp,"misc");
				plotter->Fill2D("dt_(sx3,mcp)_vs_(sx3,rf)",640,-1400,2000,640,-2000,2000,dt_sx3_mcp,dt_sx3_rf,"misc");
				plotter->Fill2D("dt_(sx3,mcp)_vs_(rf,mcp)",640,-1400,-600,640,-2000,2000,dt_sx3_mcp,dt_rf_mcp,"misc");
				if(dt_sx3_mcp > -1160 && dt_sx3_mcp < -960 && dt_rf_mcp > 60 && dt_rf_mcp < 160) {
					plotter->Fill2D("dt_(sx3,mcp)_vs_(rf,mcp)_beam1",640,-1400,-600,640,-2000,2000,dt_sx3_mcp,dt_rf_mcp,"misc");
					is_fluorine=true;
				}
				if(found_needle) {
					plotter->Fill2D("dt_(rf,mcp)_vs_needleE_sx3_beam"+std::to_string(is_fluorine),640,-2000,2000,800,0,16384,dt_rf_mcp,needle_energy,"misc");
					plotter->Fill2D("dt_(rf,mcp)_vs_dt(sx3,needle)_beam"+std::to_string(is_fluorine),640,-2000,2000,800,-2000,2000,dt_rf_mcp,ts_sx3-ts_needle,"misc");
				}
			}
			ctr+=1;
		}//end sx3 loop
	}
#endif
	if(process_alpha_proton_scattering) {
		//protonAlphaHistograms(plotter,QQQ_Events,SX3_Events,PC_Events);
		//return kTRUE;
	}//end if(process_alpha_proton_scattering)
	
	paMiscHistograms(plotter,QQQ_Events,SX3_Events,PC_Events);
	paMiscHistograms_oneWire(plotter, QQQ_Events, aClusters);
	
	//ppMiscHistograms(plotter,QQQ_Events,SX3_Events,PC_Events);
	//ppMiscHistograms_oneWire(plotter, QQQ_Events, aClusters);

	//miscHistograms_oneWire_27Alaa(plotter, QQQ_Events, SX3_Events, aClusters);
	//miscHistograms_27Alaa(plotter,QQQ_Events,SX3_Events,PC_Events);
	//miscHistograms_27Alaa_sx3(plotter,QQQ_Events,SX3_Events,PC_Events,"");


	//miscHistograms_17Faa(plotter,QQQ_Events,SX3_Events,PC_Events,"beam"+std::to_string(is_fluorine));
	//miscHistograms_oneWire_17Faa(plotter, QQQ_Events, aClusters);

	//ppMiscHistograms_sx3(plotter,QQQ_Events,SX3_Events,PC_Events);
	
	return kTRUE;
	
	if(QQQ_Events.size() && PC_Events.size())
		plotter->Fill2D("PCEv_vs_QQQEv",20,0,20,20,0,20,QQQ_Events.size(),PC_Events.size());

	plotter->Fill2D("ac_vs_cc",20,0,20,20,0,20,aClusters.size(),cClusters.size(),"wiremult");
	for(auto cluster: aClusters) {
		plotter->Fill1D("aClusters"+std::to_string(aClusters.size()),20,-5,15,cluster.size(),"wiremult");
	}
	for(auto cluster: cClusters) {
		plotter->Fill1D("cClusters"+std::to_string(cClusters.size()),20,-5,15,cluster.size(),"wiremult");
	}

	if(cClusters.size() && aClusters.size()) {
		plotter->Fill2D("ac_vs_cc_ign0",20,0,20,20,0,20,aClusters.size(),cClusters.size(),"wiremult");
	}
	
	for(auto sx3event: SX3_Events) {
	
		/*for(const auto pcevent : PC_Events_OnlyAnode) {
			plotter->Fill1D("dt_oa_sx3",800,-2000,2000,sx3event.Time1 - apTSMaxE,"pcev_onlyanode");
			if(sx3event.Time1 - apTSMaxE < 150) {
				plotter->Fill1D("dt_oa_sx3_gated",800,-2000,2000,sx3event.Time1 - pcevent.Energy1,"pcev_onlyanode");
				plotter->Fill2D("dt_oa_sx3_gated_vs_sx3E",800,-2000,2000,800,0,10,sx3event.Time1 - pcevent.Time1,sx3event.Energy1,"pcev_onlyanode");
				plotter->Fill2D("dEa_oa_Esx3_TC1_ignC"+std::to_string(acluster.size()),400,0,10,800,0,40000,sx3event.Energy1,pcevent.Energy1,"pcev_onlyanode");
				plotter->Fill2D("pcPhi_oa_sx3Phi_TC1_ignC"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,sx3event.pos.Phi()*180./M_PI,"pcev_onlyanode");
				plotter->Fill2D("pcZ_oa_sx3Z_TC1_ignC"+std::to_string(acluster.size()),300,-100,200,400,-200,200,sx3event.pos.Z(),pc_closest.Z(),"pcev_onlyanode");

				double sx3theta = TMath::ATan2(88.0,sx3event.pos.Z()-source_vertex);
				double pczguess = z_to_crossover_rho(pc_closest.Z())/TMath::Tan(sx3theta)  + source_vertex;
				double sinTheta = TMath::Sin(sx3theta);

				plotter->Fill2D("pcZ_oa_sx3pczguess_TC1_ignC"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"pcev_onlyanode");
				TVector3 x2(pc_closest), x1(sx3event.pos);
				TVector3 v = x2-x1;
				double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
				TVector3 vector_closest_to_axis = x1 + t_minimum*v;

				double sinTheta2 = TMath::Sin(TMath::ATan2(88.0,sx3event.pos.Z()-vector_closest_to_axis.Z()));
				plotter->Fill2D("dEa3_oa_Esx3_TC1_ignC"+std::to_string(acluster.size()),400,0,10,800,0,30000,sx3event.Energy1,apSumE*sinTheta2,"pcev_
				
				
				");

				plotter->Fill2D("vertexZ_oa_sx3Z_TC1_ignC"+std::to_string(acluster.size()),300,-100,200,400,-200,200,sx3event.pos.Z(),vector_closest_to_axis.Z(),"pcev_onlyanode");
				plotter->Fill2D("vertexXY_oa_TC1_ignC"+std::to_string(acluster.size()),200,-100,100,200,-100,100,vector_closest_to_axis.X(),vector_closest_to_axis.Y(),"pcev_onlyanode");
			}
		}*/
		if(sx3event.Energy1<1.2) continue;
		for(const auto acluster: aClusters) {
			auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster,"ANODE");
			int a_number = acluster.size();
			TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire,sx3event.pos.Phi());
			plotter->Fill1D("dt_anode_interp_sx3",800,-2000,2000,sx3event.Time1 - apTSMaxE,"ainterp_noc");
			if(sx3event.Time1 - apTSMaxE < 150) {
				bool phicut = sx3event.pos.Phi() <= pc_closest.Phi()+TMath::Pi()/4. && sx3event.pos.Phi() >= pc_closest.Phi()-TMath::Pi()/4.;
				double sx3theta = TMath::ATan2(88.0,sx3event.pos.Z()-source_vertex);
				double pczguess = z_to_crossover_rho(pc_closest.Z())/TMath::Tan(sx3theta)  + source_vertex;
				double sinTheta = TMath::Sin(sx3theta);

				plotter->Fill1D("dt_anode_ainterp_sx3_gated",800,-2000,2000,sx3event.Time1 - apTSMaxE,"ainterp_noc");
				plotter->Fill2D("dt_anode_ainterp_sx3_gated_vs_sx3E",800,-2000,2000,800,0,10,sx3event.Time1 - apTSMaxE,sx3event.Energy1,"ainterp_noc");
				plotter->Fill2D("dEa_ainterp_Esx3_TC1_ignC_a"+std::to_string(acluster.size()),400,0,10,800,0,40000,sx3event.Energy1,apSumE,"ainterp_noc");
				plotter->Fill2D("pcPhi_ainterp_sx3Phi_TC1_ignC_a"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,sx3event.pos.Phi()*180./M_PI,"ainterp_noc");
				plotter->Fill2D("pcZ_ainterp_sx3Z_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),300,-100,200,400,-200,200,sx3event.pos.Z(),pc_closest.Z(),"ainterp_noc");

				plotter->Fill2D("pcZ_ainterp_sx3pczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc");
				TVector3 x2(pc_closest), x1(sx3event.pos);
				TVector3 v = x2-x1;
				double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
				TVector3 vector_closest_to_axis = x1 + t_minimum*v;
				double sinTheta2 =  TMath::Sin((sx3event.pos - vector_closest_to_axis).Theta());;
				plotter->Fill2D("dEa3_ainterp_Esx3_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),1200,0,30,800,0,30000,sx3event.Energy1,apSumE*sinTheta2*3.,"ainterp_noc");
				plotter->Fill2D("vertexZ_ainterp_sx3Z_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,800,-400,400,sx3event.pos.Z(),vector_closest_to_axis.Z(),"ainterp_noc");
				plotter->Fill1D("vertexZ1d_ainterp_sx3Z_TC1_ignC_a"+std::to_string(acluster.size()),800,-400,400,vector_closest_to_axis.Z(),"ainterp_noc");
				plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),200,-100,100,200,-100,100,vector_closest_to_axis.X(),vector_closest_to_axis.Y(),"ainterp_noc");
			}
		}

		for(int i=0; i<24; i++) {
			if(aWireEvents.find(i) != aWireEvents.end()) {
				auto awire = aWireEvents[i];
				if(sx3event.Time1 -(double)std::get<2>(awire)< 150) {
					//plotter->Fill2D("sx3_z_phi2_awire"+std::to_string(std::get<0>(awire)), 400,-100,100, 100, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );	    		
					//plotter->Fill2D("sx3_z_strip#_awire"+std::to_string(std::get<0>(awire)), 400,-100,100, 100, -50,50,sx3event.pos.Z(), sx3event.ch2);	    		
					plotter->Fill2D("anodeNum_vs_stripNum",64,0,64,24,0,24,sx3event.ch2,i,"onewire");
					bool sx3diagonalphi = (!plotter->FindCut("anode_sx3_diag1")->IsInside(sx3event.ch2,i)) || (!plotter->FindCut("anode_sx3_diag2")->IsInside(sx3event.ch2,i));
					plotter->Fill2D("anodeNum_vs_stripNum_diag"+std::to_string(sx3diagonalphi),64,0,64,24,0,24,sx3event.ch2,i,"onewire");
					plotter->Fill2D("onewire_dEa_Esx3_TC1_fullev"+std::to_string(PC_Events.size()>0)+"_PC"+std::to_string(sx3diagonalphi),400,0,10,800,0,40000,sx3event.Energy1,std::get<1>(awire),"onewire");
					//plotter->Fill2D("onewire_aNum_sx3Phi_TC1_fullev"+std::to_string(PC_Events.size()>0),24,0,24,120,-360,360,i,sx3event.pos.Phi()*180./M_PI,"onewire");
					//plotter->Fill2D("sx3_z_phi_ow_awire"+std::to_string(anodeIndex)+"_sx3strip"+std::to_string(sx3event.ch2), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
				}
			}

			if(cWireEvents.find(i) != cWireEvents.end()) {
				auto cwire = cWireEvents[i];
				if(sx3event.Time1 -(double)std::get<2>(cwire) < 150) {
					//plotter->Fill2D("sx3_z_phi2_cwire"+std::to_string(std::get<0>(cwire)),400,-100,100, 100, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );	    	
					//plotter->Fill2D("sx3_z_strip#_cwire"+std::to_string(std::get<0>(cwire)),400,-100,100, 100, -50,50,sx3event.pos.Z(), sx3event.ch2 );	    	
					plotter->Fill2D("onewire_dEc_Esx3_fullev"+std::to_string(PC_Events.size()>0),400,0,10,800,0,40000,sx3event.Energy1,std::get<1>(cwire),"onewire");
					plotter->Fill2D("onewire_cNum_sx3Phi_TC1_fullev"+std::to_string(PC_Events.size()>0),24,0,24,120,-360,360,i,sx3event.pos.Phi()*180./M_PI,"onewire");
				}
			}
	   }//for 'i' loop
	}

	for(auto sx3event:SX3_Events) {
		PCSX3TimeCut=false;
		for(auto pcevent:PC_Events) {
			plotter->Fill1D("dt_pcA_sx3B"+std::to_string(sx3event.ch2),640,-2000,2000,sx3event.Time1 - pcevent.Time1,"hTiming");
			plotter->Fill1D("dt_pcC_sx3B"+std::to_string(sx3event.ch2),640,-2000,2000,sx3event.Time1 - pcevent.Time2,"hTiming");
			if(sx3event.Time1 - pcevent.Time1 < 200)//-150 for alphas
				PCASX3TimeCut = 1;
			if(sx3event.Time1 - pcevent.Time2 < 200)//-200 for alphas
				PCCSX3TimeCut = 1;
			PCSX3TimeCut = PCASX3TimeCut && PCCSX3TimeCut;

			bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && sx3event.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4.;

			plotter->Fill1D("dt_pcA_sx3B",640,-2000,2000,sx3event.Time1 - pcevent.Time1);
			plotter->Fill1D("dt_pcC_sx3B",640,-2000,2000,sx3event.Time1 - pcevent.Time2);
			plotter->Fill2D("dt_pcA_vs_sx3RE",640,-2000,2000,400,0,30,sx3event.Time1-pcevent.Time1, sx3event.Energy1);
			plotter->Fill2D("dE_E_Anodesx3B",400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1);
			plotter->Fill2D("dE_E_Cathodesx3B",400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2);
			if(pcevent.multi1==1 && pcevent.multi2==2) plotter->Fill2D("dE_E_Anodesx3B_a1c2",400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1);
			if(pcevent.multi1==1 && pcevent.multi2==2) plotter->Fill2D("dE_E_Cathodesx3B_a1c2",400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2);
			if(pcevent.multi1==2 && pcevent.multi2==1) plotter->Fill2D("dE_E_Anodesx3B_a2c1",400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1);
			if(pcevent.multi1==2 && pcevent.multi2==1) plotter->Fill2D("dE_E_Cathodesx3B_a2c1",400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2);
			plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<150),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
			plotter->Fill1D("d_sx3phi_minus_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<150),180,-360,360,sx3event.pos.Phi()*180/M_PI-pcevent.pos.Phi()*180/M_PI);
			if(PCSX3TimeCut) {
				plotter->Fill1D("dt_pcA_sx3B_timecut",640,-2000,2000,sx3event.Time1 - pcevent.Time1);
				plotter->Fill1D("dt_pcC_sx3B_timecut",640,-2000,2000,sx3event.Time1 - pcevent.Time2);	
				plotter->Fill2D("xyplot_sx3"+std::to_string(sx3event.ch2/4),100,-100,100,100,-100,100,sx3event.pos.X(),sx3event.pos.Y());
				plotter->Fill2D("xyplot_sx3"+std::to_string(sx3event.ch2/4),100,-100,100,100,-100,100,pcevent.pos.X(),pcevent.pos.Y());
				plotter->Fill2D("pcz_vs_pcphi_TimeCut",600,-200,200,120,-360,360,pcevent.pos.Z(),pcevent.pos.Phi()*180/M_PI); //x-axis is all Si det, y-axis is PC anode+cathode only
			}

			double sx3rho = 88.0;//approximate barrel radius
			double sx3z = sx3event.pos.Z(); //w.r.t target origin at 90 for run12
			double pcz = pcevent.pos.Z();
			double calcsx3theta = TMath::ATan2(sx3rho-z_to_crossover_rho(pcz),sx3z-pcz);
			plotter->Fill2D("dE4_E_Anodesx3B",400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1*TMath::Sin(calcsx3theta));
			plotter->Fill2D("dE4_E_Cathodesx3B",400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2*TMath::Sin(calcsx3theta));

			double sx3theta = TMath::ATan2(sx3rho,sx3z-source_vertex);
			double pczguess_int = z_to_crossover_rho(pcevent.pos.Z())/TMath::Tan(sx3theta)  + source_vertex;
			double pcz_guess_int_self = pcevent.pos.Perp()/TMath::Tan(sx3theta)  + source_vertex;
			plotter->Fill1D("d_guess_self_vs_int_sx3",4000,-200,200,pczguess_int-pcz_guess_int_self);
			plotter->Fill2D("guess_self_vs_int_sx3",400,-200,200,400,-200,200,pczguess_int,pcz_guess_int_self);

			double sinTheta = TMath::Sin(sx3theta);
			plotter->Fill2D("dE2_E_Anodesx3B",400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1*TMath::Sin(sx3theta));
			plotter->Fill2D("dE2_E_Cathodesx3B",400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2*TMath::Sin(sx3theta));

			TVector3 x2(pcevent.pos), x1(sx3event.pos);
			TVector3 v = x2-x1;
			double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
			TVector3 vector_closest_to_z_sx3 = x1 + t_minimum*v;
			if(vector_closest_to_z_sx3.Perp()>20.0) continue;
			double path_length_s = (sx3event.pos-vector_closest_to_z_sx3).Mag()*0.1;
			double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1)-path_length_s);
			double sinTheta2 = TMath::Sin((sx3event.pos - vector_closest_to_z_sx3).Theta());///TMath::Sin((TVector3(51.5,0,128.) - TVector3(0,0,85)).Theta());
			plotter->Fill2D(std::string("dE3_Ef_Anodesx3B")+"_PC"+std::to_string(phicut),400,0,30,800,0,40000,sx3Efix,pcevent.Energy1*sinTheta2*3.);
			plotter->Fill2D(std::string("dE3_E_Anodesx3B")+"_PC"+std::to_string(phicut),400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1*sinTheta2*3.);
			plotter->Fill2D(std::string("dE3_Ef_Cathodesx3B")+"_PC"+std::to_string(phicut),400,0,30,800,0,10000,sx3Efix,pcevent.Energy2*sinTheta2*3.);
			
			plotter->Fill1D("VertexReconZ_SX3"+std::to_string(PCSX3TimeCut),600,-1300,1300,vector_closest_to_z_sx3.Z(),"hPCZSX3");
			plotter->Fill2D("VertexReconXY_SX3"+std::to_string(PCSX3TimeCut),100,-100,100,100,-100,100,vector_closest_to_z_sx3.X(),vector_closest_to_z_sx3.Y(),"hPCZSX3");

			plotter->Fill2D("pcz_vs_time",2000,0,2000,600,-200,200,pcevent.Time1*1e-9,pcevent.pos.Z()); //x-axis is all Si det, y-axis is PC anode+cathode only
			plotter->Fill2D("pcphi_vs_time",2000,0,2000,180,-360,360,pcevent.Time1*1e-9,pcevent.pos.Phi()*180./M_PI); //x-axis is all Si det, y-axis is PC anode+cathode only
			//plotter->Fill2D("pcz_vs_time_strip"+std::to_string(sx3event.ch2),2000,0,2000,600,-200,200,pcevent.Time1*1e-9,pcevent.pos.Z()); //x-axis is all Si det, y-axis is PC anode+cathode only
			plotter->Fill2D("sx3phi_vs_time",2000,0,2000,180,-360,360,pcevent.Time1*1e-9,sx3event.pos.Phi()*180./M_PI); //x-axis is all Si det, y-axis is PC anode+cathode only
			plotter->Fill2D("pcz_vs_sx3pczguess",400,-200,200,400,-200,200,pczguess_int,pcevent.pos.Z()); //x-axis is all Si det, y-axis is PC anode+cathode only
			plotter->Fill2D("pcz_vs_sx3pczguess_self",400,-200,200,400,-200,200,pcz_guess_int_self,pcevent.pos.Z());
			plotter->Fill1D("d_pcz_vs_sx3pczguess",400,-200,200,pczguess_int-pcevent.pos.Z());
			plotter->Fill1D("d_pcz_vs_sx3pczguess_self",400,-200,200,pcz_guess_int_self-pcevent.pos.Z());

			if(pcevent.multi1>0) { //any anodes at all
				//plotter->Fill2D("anodeNum_vs_WedgeNum",64,0,64,24,0,24,qqqevent.ch2,i,"onewire");
			}

			if(pcevent.multi1 == 1 && pcevent.multi2==1) {
				plotter->Fill2D("pcz_vs_sx3pczguess_A1C1",600,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
				double pcz_fix_a1c1 = pcfix_func_a1c1.Eval(pcevent.pos.Z());

				TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix_a1c1);
				TVector3 v = x2f-x1;
				double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
				TVector3 r_rhoMin_fix = x1 + t_minimum*v;
				plotter->Fill1D("VertexRecon_pczfix_sx3_a1c1",800,-300,300,r_rhoMin_fix.Z());
				plotter->Fill1D("pczfix_A1C1_1d_sx3",600,-200,200,pcz_fix_a1c1);
				plotter->Fill2D("pczfix_vs_sx3pczguess_A1C1",600,-200,200,600,-200,200,pczguess_int,pcz_fix_a1c1);
				plotter->Fill2D("pcz_vs_sx3pczguess_A1C1_strip"+std::to_string(sx3event.ch2),300,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
			}
			if(pcevent.multi1==1 && pcevent.multi2==2) { //a1c2
				bool TCC = sx3event.Time1 - cathodeT   < 0;
				bool TCA = sx3event.Time1 - anodeT   < 0;
				plotter->Fill2D("pcz_vs_sx3pczguess_A1C2",600,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
				plotter->Fill2D("pcz_vs_sx3pczguess_self_A1C2",400,-200,200,400,-200,200,pcz_guess_int_self,pcevent.pos.Z());
				plotter->Fill1D("d_sx3pczguess_minus_pcz_a1c2",400,-200,200,pczguess_int-pcevent.pos.Z());
				plotter->Fill1D("d_sx3pczguess_self_minus_pcz_a1c2",400,-200,200,pcz_guess_int_self-pcevent.pos.Z());
				double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
				plotter->Fill1D("d_sx3pczguess_minus_pczfix_a1c2",400,-200,200,pczguess_int-pcz_fix);
				TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
				TVector3 v = x2f-x1;
				double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
				TVector3 r_rhoMin_fix = x1 + t_minimum*v;
				plotter->Fill1D("VertexRecon_pczfix_sx3",800,-300,300,r_rhoMin_fix.Z());
				plotter->Fill1D("pczfix_A1C2_1d_sx3",600,-200,200,pcz_fix);
				plotter->Fill2D("pczfix_vs_sx3pczguess_A1C2",600,-200,200,600,-200,200,pczguess_int,pcz_fix);
				plotter->Fill2D("pcz_vs_sx3pczguess_A1C2_strip"+std::to_string(sx3event.ch2),300,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
				double sinTheta_customV = TMath::Sin((sx3event.pos - r_rhoMin_fix).Theta());
				plotter->Fill2D("dE3_E_CathodeSX3_A1C2_TC"+std::to_string(PCSX3TimeCut)+"_PC"+std::to_string(phicut),400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2*sinTheta_customV);
				plotter->Fill2D("dE3_E_AnodeSX3_A1C2_TC"+std::to_string(PCSX3TimeCut)+"_PC"+std::to_string(phicut),400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1*sinTheta_customV);
				/*if(TMath::Abs(r_rhoMin_fix.Z())<200.0) {
					plotter->Fill2D("dE3_E_AnodeSX3B_A1C2_(vertex_fix_z/100)="+std::to_string(floor(r_rhoMin_fix.Z()/100.0)),400,0,30,800,0,40000,sx3event.Energy1,pcevent.Energy1*sinTheta_customV);
					plotter->Fill2D("dE3_E_CathodeSX3B_A1C2_(vertex_fix_z/100)="+std::to_string(floor(r_rhoMin_fix.Z()/100.0)),400,0,30,800,0,10000,sx3event.Energy1,pcevent.Energy2*sinTheta_customV);
				}*/
			}
			if(pcevent.multi1==1 && pcevent.multi2==3) {
				plotter->Fill2D("pcz_vs_sx3pczguess_A1C3",600,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
				plotter->Fill2D("pcz_vs_sx3pczguess_A1C3_strip"+std::to_string(sx3event.ch2),300,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
			}

			if(pcevent.multi1==2 && pcevent.multi2==1) {
				plotter->Fill2D("pcz_vs_sx3pczguess_A2C1",600,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
				double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
				TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
				TVector3 v = x2f-x1;
				double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
				TVector3 r_rhoMin_fix = x1 + t_minimum*v;
				plotter->Fill1D("VertexRecon_pczfix_sx3",800,-300,300,r_rhoMin_fix.Z());
				plotter->Fill1D("pczfix_A2C1_1d_sx3",600,-200,200,pcz_fix);
				plotter->Fill2D("pczfix_vs_sx3pczguess_A2C1",600,-200,200,600,-200,200,pczguess_int,pcz_fix);
				//plotter->Fill2D("pcz_vs_sx3pczguess_A2C1_strip"+std::to_string(sx3event.ch2),300,-200,200,600,-200,200,pczguess,pcevent.pos.Z());
			}
			plotter->Fill2D("pcz_vs_sx3pczguess",600,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z()); //x-axis is all Si det, y-axis is PC anode+cathode only
			plotter->Fill2D("pcz_vs_sx3pczguess_strip"+std::to_string(sx3event.ch2),300,-200,200,600,-200,200,pczguess_int,pcevent.pos.Z());
			//plotter->Fill2D("pcz_vs_sx3pczguess_phi"+std::to_string(sx3event.pos.Phi()*180/M_PI),300,0,200,600,-200,200,pczguess,pcevent.pos.Z());
			/*plotter->Fill2D("pcz_vs_sx3z_strip="+std::to_string(sx3event.ch2),300,0,100,600,-200,200,sx3z,pcevent.pos.Z(),"sx3_vs_pc_zcorr");
			plotter->Fill2D("pcz_vs_sx3z_strip="+std::to_string(sx3event.ch2)+"_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),300,0,100,600,-200,200,sx3z,pcevent.pos.Z(),"sx3_vs_pc_zcorr");

			plotter->Fill2D("pcdEC_vs_sx3z_strip="+std::to_string(sx3event.ch2)+"_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),800,0,20000,600,-200,200,pcevent.Energy2,sx3z,"sx3_vs_pc_zcorr");
			plotter->Fill2D("pcdEA_vs_sx3z_strip="+std::to_string(sx3event.ch2)+"_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),800,0,20000,600,-200,200,pcevent.Energy1,sx3z,"sx3_vs_pc_zcorr");*/

			/*for(auto cc: cClusters)
			 for(auto ac: aClusters) {
				plotter->Fill2D("pcz_sx3_phicut_a"+std::to_string(ac.size())+"_c"+std::to_string(cc.size())+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
				if(ac.size()==2 && cc.size()==1) {
					plotter->Fill2D("pcz_sx3_phicut_a("+std::to_string(std::get<0>(ac.at(0)))+","+std::to_string(std::get<0>(ac.at(1)))+")_c"+std::to_string(std::get<0>(cc.at(0)))+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
					plotter->Fill2D("a2c1_vs_sx3_strip",24,0,24,64,0,64,0.5*(std::get<0>(ac.at(0))+std::get<0>(ac.at(1))),sx3event.ch2,"hPCZSX3");
					//plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<-150)+"_a("+std::to_string(std::get<0>(ac.at(0)))+","+std::to_string(std::get<0>(ac.at(1)))+")_c"+std::to_string(std::get<0>(cc.at(0))),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
				}
				if(cc.size()==2 && ac.size()==1) {
					plotter->Fill2D("pcz_sx3_phicut_c("+std::to_string(std::get<0>(cc.at(0)))+","+std::to_string(std::get<0>(cc.at(1)))+")_a"+std::to_string(std::get<0>(ac.at(0)))+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
					plotter->Fill2D("c2a1_vs_sx3_strip",24,0,24,64,0,64,0.5*(std::get<0>(cc.at(0))+std::get<0>(cc.at(1))),sx3event.ch2,"hPCZSX3");
					plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<-150)+"_c("+std::to_string(std::get<0>(cc.at(0)))+","+std::to_string(std::get<0>(cc.at(1)))+")_a"+std::to_string(std::get<0>(ac.at(0))),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
					//plotter->Fill2D("pcz_vs_sx3z_2C1A_phiCut_TC"+std::to_string(PCSX3TimeCut),300,0,200,600,-400,400,sx3z,pcevent.pos.Z());
				}

				if(ac.size()==1 && cc.size()==1) {
					plotter->Fill2D("pcz_sx3_phicut_a("+std::to_string(std::get<0>(ac.at(0)))+")_c"+std::to_string(std::get<0>(cc.at(0)))+"_sx3guess",300,0,200,600,-200,200,sx3z,pcevent.pos.Z(),"hPCZSX3");
					//plotter->Fill2D("a2c1_vs_sx3_strip",24,0,24,64,0,64,0.5*(std::get<0>(ac.at(0))+std::get<0>(ac.at(1))),sx3event.ch2,"hPCZSX3");
					//plotter->Fill2D("sx3phi_vs_pcphi"+std::to_string(sx3event.Time1 - pcevent.Time1<-150)+"_a("+std::to_string(std::get<0>(ac.at(0)))+")_c"+std::to_string(std::get<0>(cc.at(0))),100,-360,360,100,-360,360,sx3event.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
				}
			}*/  //end for

			bool sx3PhiCut = (TMath::Abs(sx3event.pos.Phi()-pcevent.pos.Phi()) < 45.0*M_PI/180.);
			plotter->Fill1D("pcz_sx3Coinc_phiCut"+std::to_string(sx3PhiCut)+"_TC"+std::to_string(PCSX3TimeCut),300,0,200,sx3z);
			plotter->Fill2D("pcz_vs_sx3z_phiCut"+std::to_string(sx3PhiCut)+"_TC"+std::to_string(PCSX3TimeCut),300,0,200,600,-400,400,sx3z,pcevent.pos.Z());

			//plotter->Fill2D("sx3E_vs_sx3z"+std::to_string(sx3event.ch2),400,0,30,300,0,200,sx3event.Energy1,sx3z);
			plotter->Fill2D("sx3E_vs_sx3z",400,0,30,300,0,200,sx3event.Energy1,sx3z);

			plotter->Fill2D("pcdEA_vs_sx3z",800,0,20000,300,0,200,pcevent.Energy1,sx3z);
			plotter->Fill2D("pcdEC_vs_sx3z",800,0,20000,300,0,200,pcevent.Energy2,sx3z);

			plotter->Fill2D("pcdEA_vs_sx3z"+std::to_string(sx3event.ch2),800,0,20000,300,0,200,pcevent.Energy1,sx3z,"pcE_vs_sx3pos");
			plotter->Fill2D("pcdEC_vs_sx3z"+std::to_string(sx3event.ch2),800,0,20000,300,0,200,pcevent.Energy2,sx3z,"pcE_vs_sx3pos");

			plotter->Fill2D("pcdE2A_vs_sx3z",800,0,20000,300,0,200,pcevent.Energy1*sinTheta,sx3z);
			plotter->Fill2D("pcdE2C_vs_sx3z",800,0,20000,300,0,200,pcevent.Energy2*sinTheta,sx3z);
			plotter->Fill2D("phi_vs_stripnum",180,-180,180,48,0,48,pcevent.pos.Phi()*180./M_PI,sx3event.ch2);
			plotter->Fill2D("E_theta_AnodeSX3",400,-20,180,300,0,15,sx3theta*180/M_PI,sx3event.Energy1);
		}
	}//end PC-SX3 coincidence

	for(auto qqqevent: QQQ_Events) {
		/* Events with QQQ present, but PC events don't have a reliable cathode signal, so we scan just the anode wire clusters */
		for(auto pcevent: PC_Events) {
			bool phicut = qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4.;
			plotter->Fill1D("dt_pcA_qqqR",640,-2000,2000,qqqevent.Time1 - pcevent.Time1);
			plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE_PC"+std::to_string(phicut),640,-2000,2000,400,0,30,qqqevent.Time1-pcevent.Time1, qqqevent.Energy1);
			plotter->Fill1D("dt_pcC_qqqW",640,-2000,2000,qqqevent.Time2 - pcevent.Time2);
			plotter->Fill2D("phiPC_vs_phiQQQ",180,-360,360,180,-360,360,qqqevent.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
			double sinTheta = TMath::Sin((qqqevent.pos - TVector3(0,0,source_vertex)).Theta());///TMath::Sin((TVector3(51.5,0,128.) - TVector3(0,0,85)).Theta());

			TVector3 x2(pcevent.pos);
			TVector3 x1(qqqevent.pos);
			TVector3 v = x2-x1;
			double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
			TVector3 r_rhoMin = x1 + t_minimum*v;
			if(r_rhoMin.Perp()>20.0) continue;
			double path_length_q = (qqqevent.pos-r_rhoMin).Mag()*0.1;
			double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
			double sinTheta2 = TMath::Sin((qqqevent.pos - r_rhoMin).Theta());///TMath::Sin((TVector3(51.5,0,128.) - TVector3(0,0,85)).Theta());

			//bool timecut = (qqqevent.Time1 - pcevent.Time1 < -150);
			bool timecut = (qqqevent.Time1 - pcevent.Time1 < 150);
			bool lowercut_cath = pcevent.Energy2*sinTheta < 250 && (qqqevent.Energy2 < 5.0 || qqqevent.Energy1 < 5.0) ;
			
			if(lowercut_cath && phicut) {
				plotter->Fill1D("dt_pcA_qqqR_pidlow_PC1",640,-2000,2000,qqqevent.Time1 - pcevent.Time1);
				plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE_pidlow_PC1",640,-2000,2000,400,0,30,qqqevent.Time1-pcevent.Time1, qqqevent.Energy1);
				plotter->Fill1D("dt_pcC_qqqW_pidlow_PC1",640,-2000,2000,qqqevent.Time2 - pcevent.Time2);
			}
			if(timecut) {// && qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4. ) {

				plotter->Fill2D("dE_E_AnodeQQQR",400,0,30,800,0,40000,qqqevent.Energy1,pcevent.Energy1);
				plotter->Fill2D("dE_E_CathodeQQQR",400,0,30,800,0,10000,qqqevent.Energy2,pcevent.Energy2);
				plotter->Fill2D(std::string("dE3_Ef_AnodeQQQR")+"_PC"+std::to_string(phicut),400,0,30,800,0,40000,qqqEfix,pcevent.Energy1*sinTheta2*3.);
				plotter->Fill2D(std::string("dE3_Ef_CathodeQQQR")+"_PC"+std::to_string(phicut),400,0,30,800,0,10000,qqqEfix,pcevent.Energy2*sinTheta2*3.);
				plotter->Fill2D(std::string("dE3_E_AnodeQQQR")+"_PC"+std::to_string(phicut),400,0,30,800,0,40000,qqqevent.Energy1,pcevent.Energy1*sinTheta2*3.);
				plotter->Fill2D(std::string("dE3_E_CathodeQQQR")+"_PC"+std::to_string(phicut),400,0,30,800,0,10000,qqqevent.Energy1,pcevent.Energy2*sinTheta2*3.);

				
				if(pcevent.multi1==1 && pcevent.multi2==2) plotter->Fill2D("dE_E_AnodeQQQR_a1c2",400,0,30,800,0,40000,qqqevent.Energy1,pcevent.Energy1);
				if(pcevent.multi1==1 && pcevent.multi2==2) plotter->Fill2D("dE_E_CathodeQQQR_a1c2",400,0,30,800,0,10000,qqqevent.Energy1,pcevent.Energy2);
				if(pcevent.multi1==2 && pcevent.multi2==1) plotter->Fill2D("dE_E_AnodeQQQR_a2c1",400,0,30,800,0,40000,qqqevent.Energy1,pcevent.Energy1);
				if(pcevent.multi1==2 && pcevent.multi2==1) plotter->Fill2D("dE_E_CathodeQQQR_a2c1",400,0,30,800,0,10000,qqqevent.Energy1,pcevent.Energy2);

				/*if(phicut) {
					plotter->Fill2D("dE2_E_AnodeQQQR_TC1PC1_pidlow"+std::to_string(lowercut_cath),400,0,30,800,0,4000,qqqevent.Energy1,pcevent.Energy1*sinTheta);
					plotter->Fill2D("dE2_E_CathodeQQQW_TC1PC1_pidlow"+std::to_string(lowercut_cath),400,0,30,800,0,1000,qqqevent.Energy2,pcevent.Energy2*sinTheta);
					//plotter->Fill2D("E_theta_AnodeQQQR_TC1PC1_pidlow"+std::to_string(lowercut_cath),75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
					plotter->Fill2D("E_theta_zoomin_AnodeQQQR_TC1PC1_pidlow"+std::to_string(lowercut_cath),60,0,30,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);

				}*/

				plotter->Fill2D("dE2_E_AnodeQQQR_TC1_PC"+std::to_string(phicut),400,0,30,800,0,4000,qqqevent.Energy1,pcevent.Energy1*sinTheta);
				plotter->Fill2D("dE2_E_CathodeQQQR_TC1_PC"+std::to_string(phicut),400,0,30,800,0,1000,qqqevent.Energy2,pcevent.Energy2*sinTheta);
				plotter->Fill2D("dEC_vs_dEA_TC1_PC"+std::to_string(phicut),800,0,40000,800,0,10000,pcevent.Energy1,pcevent.Energy2);
				plotter->Fill2D("qqqphi_vs_time",2000,0,2000,180,-360,360,pcevent.Time1*1e-9,qqqevent.pos.Phi()*180./M_PI); //x-axis is all Si det, y-axis is PC anode+cathode only

				plotter->Fill1D("dt_pcA_qqqR_timecut",640,-2000,2000,qqqevent.Time1 - pcevent.Time1);
				plotter->Fill1D("dt_pcC_qqqW_timecut",640,-2000,2000,qqqevent.Time2 - pcevent.Time2);
				plotter->Fill2D("dE_theta_AnodeQQQR",90,0,90,400,0,20000,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,pcevent.Energy1);
				
				plotter->Fill2D("dE2_theta_AnodeQQQR",90,0,90,400,0,20000,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,pcevent.Energy1*sinTheta);
				plotter->Fill2D("phiPC_vs_phiQQQ_TimeCut",180,-360,360,180,-360,360,qqqevent.pos.Phi()*180/M_PI,pcevent.pos.Phi()*180/M_PI);
				plotter->Fill1D("d_phiPC_phiQQQ_TimeCut",180,-360,360,qqqevent.pos.Phi()*180/M_PI-pcevent.pos.Phi()*180/M_PI);
				//plotter->Fill2D("E_theta_AnodeQQQR_TC1_PC"+std::to_string(phicut),75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
				//plotter->Fill2D("E_theta_zoomin_AnodeQQQR_TC1_PC"+std::to_string(phicut),60,0,30,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
				//plotter->Fill2D("E2_theta_AnodeQQQR",75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1);
				plotter->Fill2D("Etot2_theta_AnodeQQQR",75,0,90,300,0,15,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,qqqevent.Energy1+pcevent.Energy1*anode_gain*sinTheta);

				//plotter->Fill2D("dE_theta_CathodeQQQR",75,0,90,800,0,10000,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,pcevent.Energy2);
				//plotter->Fill2D("dE2_theta_CathodeQQQR",75,0,90,800,0,10000,(qqqevent.pos - TVector3(0,0,source_vertex)).Theta()*180/M_PI,pcevent.Energy2*sinTheta);
				

				//plotter->Fill2D("dE_phi_AnodeQQQR",100,-180,180,800,0,40000,(qqqevent.pos - TVector3(0,0,source_vertex)).Phi()*180/M_PI,pcevent.Energy1);
				//plotter->Fill2D("dE_phi_CathodeQQQR",100,-180,180,800,0,10000,(qqqevent.pos - TVector3(0,0,source_vertex)).Phi()*180/M_PI,pcevent.Energy2);
				plotter->Fill1D("PCZ",800,-200,200,pcevent.pos.Z(),"phicut");
				//plotter->Fill1D("PCZ_phicut_a"+std::to_string(aClusters.at(0).size())+"_c"+std::to_string(cClusters.at(0).size()),800,-200,200,pcevent.pos.Z(),"wiremult");

				double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z())/TMath::Tan((qqqevent.pos-TVector3(0,0,source_vertex)).Theta())  + source_vertex;
				double pcz_guess_int_self = pcevent.pos.Perp()/TMath::Tan((qqqevent.pos-TVector3(0,0,source_vertex)).Theta())  + source_vertex;
				plotter->Fill1D("d_guess_self_vs_int_qqq",400,-200,200,pcz_guess_int-pcz_guess_int_self);
				plotter->Fill2D("guess_self_vs_int_qqq",400,-200,200,400,-200,200,pcz_guess_int,pcz_guess_int_self);
				


				//plotter->Fill2D("pczguess_vs_pc_int",180,0,200,150,0,200,pcz_guess_int,pcevent.pos.Z(),"phicut");
				plotter->Fill2D("pczguess_vs_pc_int",400,-200,200,600,-400,400,pcz_guess_int,pcevent.pos.Z(),"phicut");
				plotter->Fill2D("pczguess_vs_pc_int_self",400,-200,200,600,-200,200,pcz_guess_int_self,pcevent.pos.Z(),"phicut");
				plotter->Fill1D("d_pczqqq_vs_pc_int",400,-200,200,pcz_guess_int-pcevent.pos.Z());
				plotter->Fill1D("d_pczqqq_vs_pc_int_self",400,-200,200,pcz_guess_int_self-pcevent.pos.Z());
				if(pcevent.multi1==1 && pcevent.multi2==1) {
					double pcz_fix_a1c1 = pcfix_func_a1c1.Eval(pcevent.pos.Z());
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix_a1c1);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					plotter->Fill1D("pczfix_A1C1_1d_qqq",600,-200,200,pcz_fix_a1c1);
					plotter->Fill2D("pczfix_vs_qqqpczguess_A1C1",600,-200,200,600,-200,200,pcz_guess_int,pcz_fix_a1c1);
					plotter->Fill2D("pczguess_vs_pc_int_A1C1",400,-200,200,600,-400,400,pcz_guess_int,pcevent.pos.Z(),"phicut");
				}
				if(pcevent.multi1==1 && pcevent.multi2==2) { //a1c2
					double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;

 					double sinTheta_customV = TMath::Sin((qqqevent.pos - r_rhoMin_fix).Theta());            	
					plotter->Fill2D("dE3_E_CathodeQQQW_A1C2_TC1_PC"+std::to_string(phicut),400,0,30,800,0,10000,qqqevent.Energy2,pcevent.Energy2*sinTheta_customV);
					plotter->Fill2D("dE3_E_AnodeQQQR_A1C2_TC1_PC"+std::to_string(phicut),400,0,30,800,0,10000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV);

					plotter->Fill1D("VertexRecon_pczfix_qqq",800,-400,400,r_rhoMin_fix.Z());      
					plotter->Fill1D("VertexRecon_pczfix_qqq_PC"+std::to_string(phicut)+"_pidlow"+std::to_string(lowercut_cath),800,-400,400,r_rhoMin_fix.Z());      
					plotter->Fill1D("d_qqqpczguess_minus_pcz_a1c2",400,-200,200,pcz_guess_int-pcevent.pos.Z());
					plotter->Fill1D("d_qqqpczguess_self_minus_pcz_a1c2",400,-200,200,pcz_guess_int_self-pcevent.pos.Z());
					/*if(TMath::Abs(r_rhoMin_fix.Z())<200.0) {
						plotter->Fill2D("dE3_E_AnodeQQQR_A1C2_(vertex_fix_z/100)="+std::to_string(floor(r_rhoMin_fix.Z()/100.0)),400,0,30,800,0,40000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV);                        					
						plotter->Fill2D("dE3_E_CathodeQQQR_A1C2_(vertex_fix_z/100)="+std::to_string(floor(r_rhoMin_fix.Z()/100.0)),400,0,30,800,0,10000,qqqevent.Energy1,pcevent.Energy2*sinTheta_customV);                        			
					} */     	

					plotter->Fill1D("pczfix_A1C2_1d_qqq",600,-200,200,pcz_fix);     
					plotter->Fill1D("d_qqqpczguess_minus_pczfix_a1c2",400,-200,200,pcz_guess_int-pcz_fix);
					plotter->Fill1D("d_qqqpczguess_self_minus_pczfix_a1c2",400,-200,200,pcz_guess_int_self-pcz_fix);
					plotter->Fill2D("pczfix_vs_qqqpczguess_A1C2",600,-200,200,600,-200,200,pcz_guess_int,pcz_fix);
					plotter->Fill2D("pczfix_vs_qqqpczguess_self_A1C2",600,-200,200,600,-200,200,pcz_guess_int_self,pcz_fix);
					plotter->Fill2D("pczguess_vs_pc_int_A1C2",400,-200,200,600,-400,400,pcz_guess_int,pcevent.pos.Z(),"phicut");
					plotter->Fill2D("pczguess_self_vs_pc_int_A1C2",400,-200,200,600,-400,400,pcz_guess_int_self,pcevent.pos.Z(),"phicut");

					double path_length = (qqqevent.pos-TVector3(0,0,r_rhoMin_fix.Z())).Mag()*0.1;
					//std::cout << path_length << std::endl;
					double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length);
					double qqqEfix_p = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1)-path_length);

					plotter->Fill2D("E_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut),180,0,180,600,0,15,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqevent.Energy1);
					if(lowercut_cath) 
						plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut)+"_pidlow"+std::to_string(lowercut_cath),180,0,180,600,0,15,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqEfix_p);
					else {
						std::string zcut = "_"+std::to_string((TMath::Abs(r_rhoMin_fix.Z())<180));
						plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut)+"_pidlow"+std::to_string(lowercut_cath)+zcut,180,0,180,600,0,15,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqEfix);					
					}
					std::string morecuts = "_pidlow"+std::to_string(lowercut_cath)+"_vertexfix="+std::to_string(floor(r_rhoMin_fix.Z()/20)*20+10);
					//plotter->Fill2D("E_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut)+morecuts,180,0,180,800,0,8,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqevent.Energy1,"morecuts");
					//plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC"+std::to_string(phicut)+morecuts,180,0,180,800,0,8,(qqqevent.pos - TVector3(0,0,r_rhoMin_fix.Z())).Theta()*180/M_PI,qqqEfix,"morecuts");
					plotter->Fill2D("dE3_Ef_AnodeQQQR_TC1"+std::to_string(phicut)+"_pidlow"+std::to_string(lowercut_cath),600,0,15,800,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV);
					plotter->Fill2D("dE3_Ef_CathodeQQQR_TC1PC"+std::to_string(phicut)+"_pidlow"+std::to_string(lowercut_cath),600,0,15,800,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV);

				}
				if(pcevent.multi1==2 && pcevent.multi2==1) {
					plotter->Fill2D("pcz_int_vs_qqqpczguess_A2C1",600,-200,200,600,-200,200,pcz_guess_int,pcevent.pos.Z());
					plotter->Fill2D("pczguess_vs_pc_int_self_A2C1",400,-200,200,600,-200,200,pcz_guess_int_self,pcevent.pos.Z(),"phicut");
					double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					plotter->Fill1D("VertexRecon_a2c1_pczfix_qqq",800,-300,300,r_rhoMin_fix.Z());
					plotter->Fill1D("pczfix_A2C1_1d_qqqs",600,-200,200,pcz_fix);
					plotter->Fill2D("pczfix_vs_qqqpczguess_A2C1",600,-200,200,600,-200,200,pcz_guess_int,pcz_fix);
					//plotter->Fill2D("pcz_vs_sx3pczguess_A2C1_strip"+std::to_string(sx3event.ch2),300,-200,200,600,-200,200,pczguess,pcevent.pos.Z());
				}

				double qqqrho = qqqevent.pos.Perp();
				double qqqz = (qqqevent.pos - TVector3(0,0,source_vertex)).Z();
				double tan_theta = qqqrho/qqqz;
				double pcz_guess_int2 = z_to_crossover_rho(pcevent.pos.Z())/tan_theta + source_vertex;
				plotter->Fill2D("pczguess_vs_pc_int2",180,0,200,150,0,200,pcz_guess_int2,pcevent.pos.Z(),"phicut");

				double qqqz2 = (qqqevent.pos - r_rhoMin).Z();
				double tan_theta2 = qqqrho/qqqz2;
				double pcz_guess_int3 = z_to_crossover_rho(pcevent.pos.Z())/tan_theta2 + r_rhoMin.Z();
				plotter->Fill2D("pczguess_vs_pc_int3",180,0,200,150,0,200,pcz_guess_int3,pcevent.pos.Z(),"phicut");
				//plotter->Fill2D("pczguess_vs_pc_int2_a"+std::to_string(pcevent.multi1)+"_c"+std::to_string(pcevent.multi2),180,0,200,150,0,200,pcz_guess_int2,pcevent.pos.Z(),"phicut");

				plotter->Fill2D("pctheta_vs_qqqtheta_sv",180,-360,360,180,-360,360,(qqqevent.pos-TVector3(0,0,source_vertex)).Theta()*180/M_PI,(pcevent.pos-TVector3(0,0,source_vertex)).Theta()*180/M_PI,"phicut");
				plotter->Fill2D("pctheta_vs_qqqtheta_rmz",180,-360,360,180,-360,360,(qqqevent.pos-TVector3(0,0,r_rhoMin.Z())).Theta()*180/M_PI,(pcevent.pos-TVector3(0,0,r_rhoMin.Z())).Theta()*180/M_PI,"phicut");
				plotter->Fill2D("pctheta_vs_qqqtheta_rm",180,-360,360,180,-360,360,(qqqevent.pos-r_rhoMin).Theta()*180/M_PI,(pcevent.pos-r_rhoMin).Theta()*180/M_PI,"phicut");
				plotter->Fill2D("pczguess_vs_pc_phi="+std::to_string(qqqevent.pos.Phi()*180./M_PI),300,0,200,150,0,200,pcz_guess_int,pcevent.pos.Z(),"phicut");
			}
		}
	}//end PC QQQ coincidence

	return kTRUE;
}

void MakeVertex::Terminate()
{
	plotter->FlushToDisk(10);
}


void miscHistograms_oneWire_17Faa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters) {
		//consider the 'proton-like' QQQ branch seen in a,p data
		TRandom3 rand;
		rand.SetSeed();//random seed set
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 0.6) continue; //coarse gating
				//if(qqqevent.Energy1 > 5.0) continue; //coarse gating
					for(const auto acluster: aClusters) {
						auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster,"ANODE");
						//if(apSumE<6000) continue;
						int a_number = acluster.size();
						TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire,qqqevent.pos.Phi());
						plotter->Fill1D("dt_anode_interp_qqq",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
						if(qqqevent.Time1 - apTSMaxE < 150) {
							bool phicut = qqqevent.pos.Phi() <= pc_closest.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pc_closest.Phi()-TMath::Pi()/4.;
							
							TVector3 x2(pc_closest), x1(qqqevent.pos);
							TVector3 v = x2-x1;
							double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
							TVector3 r_rhoMin_fix = x1 + t_minimum*v;
													
							double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
		 					double sinTheta2 = TMath::Sin(theta_q);

							if(r_rhoMin_fix.Perp()>6.0) continue;
							//if(r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100) continue;
							if(!phicut) continue;
							plotter->Fill1D("dt_anode_ainterp_qqq_gated",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
							plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_qqqE",800,-2000,2000,800,0,30,qqqevent.Time1 - apTSMaxE,qqqevent.Energy1,"ainterp_noc");
							plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size()),400,0,30,800,0,40000,qqqevent.Energy1,apSumE,"ainterp_noc");
							plotter->Fill2D("pcPhi_ainterp_qqqPhi_TC1_ignC_a"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,qqqevent.pos.Phi()*180./M_PI,"ainterp_noc");
							plotter->Fill2D("pcZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),300,-100,200,400,-200,200,qqqevent.pos.Z(),pc_closest.Z(),"ainterp_noc");
			
			
							//plotter->Fill2D("pcZ_ainterp_qqqpczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc");
							plotter->Fill2D("dEa3_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),1200,0,30,800,0,30000,qqqevent.Energy1,apSumE*sinTheta2*3.,"ainterp_noc");
			
							plotter->Fill2D("vertexZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,800,-400,400,qqqevent.pos.Z(),r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill1D("vertexZ1d_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),800,-400,400,r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),"ainterp_noc");

							double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
							double qqqEfix;
							//qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
							qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1)-path_length_q);
							
	  						double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
		  					double beam_energy_at_vertex = cm_to_MeV_17F->Eval(MeV_to_cm_17F->Eval(65.0) - beam_path_length);
							Kinematics aakin_17F(17.0020952,4.002603254,4.002603254,17.0020952,beam_energy_at_vertex/17.0020952); //m3 is alpha
	  						//Kinematics apkin_17F(17.0020952,4.002603254,1.008664916,19.9924401753,beam_energy_at_vertex/17.0020952); //m3 is alpha
		
	  						
							auto kin_17F = aakin_17F;
							plotter->Fill1D("pmisc_ow_Ex_from_alpha",400,-20,20,kin_17F.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
							plotter->Fill2D("pmisc_ow_Ef_vs_theta_qqq",100,0,180,800,0,30,theta_q*180/M_PI,qqqEfix,"ainterp_noc");
							plotter->Fill2D("pmisc_ow_VertexReconZ_vs_Ef",800,-400,400,800,0,30,r_rhoMin_fix.Z(),qqqEfix,"ainterp_noc");
						}	
					}
		}//end QQQEvents loop
}


void miscHistograms_17Faa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events, std::string globaltag="") {
		TRandom3 rand;
		rand.SetSeed();//random seed set
		//Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,7.0); //m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.2, 250 torr)
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 0.6) continue; //coarse gating
				//if(qqqevent.Energy1 > 5.0) continue; //coarse gating
				for(auto pcevent: PC_Events) {
					if(!(pcevent.multi1==1 && pcevent.multi2<=2)) continue;
					//if(pcevent.Energy1 > 11000) continue; //coarse gating
					bool phicut = qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4.;
					if(!phicut) continue;
					//if(pcevent.Time1-qqqevent.Time1<-150 || pcevent.Time1-qqqevent.Time1 >850) continue;
					
					double pcz_fix, pcz_dith=pcevent.pos.Z();
					if(pcevent.multi2==2) 
						pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					else {
						pcz_fix = rand.Gaus(pcevent.pos.Z(),8.0);//dither for a1c1 events
						pcz_dith = pcz_fix;
					}
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(qqqevent.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					if(r_rhoMin_fix.Perp()>6) continue;
					//if(vertex_z > 90) continue;

					//double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
					double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_q);
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);
					
					double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
		  			double beam_energy_at_vertex = cm_to_MeV_17F->Eval(MeV_to_cm_17F->Eval(65.0) - beam_path_length);
					Kinematics aakin_17F(17.0020952,4.002603254,4.002603254,17.0020952,beam_energy_at_vertex/17.0020952); //m3 is alpha
	  				Kinematics apkin_17F(17.0020952,4.002603254,1.008664916,19.9924401753,beam_energy_at_vertex/17.0020952); //m3 is alpha
	  				auto kin_27Al = aakin_17F;
					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = globaltag+"_proton+misc"+tag;
						plotter->Fill2D("pmisc_dE_E_AnodeQQQ"+globaltag+tag,400,0,30,800,0,120000,qqqevent.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE_E_CathodeQQQ"+globaltag+tag,400,0,30,800,0,40000,qqqevent.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_AnodeQQQ"+globaltag+tag,400,0,30,400,0,120000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV*3.,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_CathodeQQQ"+globaltag+tag,400,0,30,400,0,40000,qqqevent.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmisc_dPhi_QQQ_PC"+globaltag+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,qqqevent.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmisc_dt_Anode_QQQ_PC"+std::to_string(phicut)+globaltag+tag,600,-2000,2000,pcevent.Time1-qqqevent.Time1,pmlabel);
						plotter->Fill1D("pmisc_dt_Cathode_QQQ"+globaltag+tag,600,-2000,2000,pcevent.Time2-qqqevent.Time1,pmlabel);
						plotter->Fill2D("pmisc_dt_Anode_E_QQQ_PC"+std::to_string(phicut)+globaltag+tag,600,-2000,2000,400,0,30,pcevent.Time1-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_AnodeQQQ_vsPCPhi"+globaltag+tag,600,-2000,2000,180,-360,360,pcevent.Time1-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill2D("pmisc_dt_Cathode_E_QQQ"+globaltag+tag,600,-2000,2000,400,0,30,pcevent.Time2-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_CathodeQQQ_vsPCPhi"+globaltag+tag,600,-2000,2000,180,-360,360,pcevent.Time2-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmisc_pczfix"+globaltag+tag,600,-300,300,pcz_fix,pmlabel);
						if(pcevent.multi2==2) {
							plotter->Fill1D("pmisc_pcz"+globaltag+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz2"+globaltag+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
						}
						if(pcevent.multi2==1) {
							plotter->Fill1D("pmisc_pcz"+globaltag+tag,600,-300,300,pcz_fix,pmlabel);
							plotter->Fill1D("pmisc_pcz1"+globaltag+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz_dith"+globaltag+tag,600,-300,300,pcz_dith,pmlabel);
						}

						//double path_length_q = (qqqevent.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
						double qqqEfix;
						if(tag == "_cathode_alphas") {//satisfied when find succeeds 
							qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
						} else {
							qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1)-path_length_q);
						}
						plotter->Fill1D("pmisc_Ex_"+globaltag+tag,400,-20,20,kin_27Al.getExc(qqqEfix,theta_q*180/M_PI),pmlabel);
						//plotter->Fill2D("qqqEf_sx3E_matrix_all"+tag,400,0,10,400,0,10,qqqEfix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_AnodeQQQ"+globaltag+tag,400,0,30,400,0,120000,qqqEfix,pcevent.Energy1*sinTheta_customV*3,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_CathodeQQQ"+globaltag+tag,400,0,30,400,0,40000,qqqEfix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill1D("pmisc_VertexReconZ"+globaltag+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmisc_VertexReconXY"+globaltag+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);						
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+globaltag+tag,800,-400,400,800,0,30,vertex_z,qqqEfix,pmlabel);
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+globaltag+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						
						plotter->Fill2D("pmisc_Ef_vs_theta_qqq"+globaltag+tag,100,0,180,800,0,30,theta_q*180/M_PI,qqqEfix,pmlabel);
						if(pcevent.multi2==1) {
							plotter->Fill2D("pmisc_Ef_vs_theta_qqq_a1c1"+globaltag+tag,100,0,180,800,0,30,theta_q*180/M_PI,qqqEfix,pmlabel);						
							plotter->Fill2D("pmisc_VertexReconZ_vs_Ef_a1c1"+globaltag+tag,800,-400,400,800,0,30,vertex_z,qqqEfix,pmlabel);
						}
					};
					
					if(cathode_alpha_select) 
						plot_with_tag("_cathode_alphas");
					else 
						plot_with_tag("_cathode_protons");
					//plot_with_tag();
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(qqqEfix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end QQQEvents loop
}

void miscHistograms_17Faa_sx3(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events, std::string globaltag="") {
		//consider the 'proton-like' QQQ branch seen in a,p data
		for(auto sx3event: SX3_Events) {
				//if(sx3event.Energy1 < 1.2) continue; //coarse gating
				//if(sx3event.Energy1 > 5.0) continue; //coarse gating
				for(auto pcevent: PC_Events) {
					if(!(pcevent.multi1==1 && pcevent.multi2==2)) continue;
					//if(pcevent.Energy1 > 11000) continue; //coarse gating
					bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/3. && sx3event.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/3.;
					if(!phicut) continue;
					//if(pcevent.Time1-sx3event.Time1<-150 || pcevent.Time1-sx3event.Time1 >850) continue;
					
					double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(sx3event.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (sx3event.pos - TVector3(0,0,vertex_z)).Theta();

					if(r_rhoMin_fix.Perp()>10.0) continue;
					double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_s);
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);

		

					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = "proton+miscsx3"+tag;
						plotter->Fill2D("pmiscs_dE_E_Anodesx3"+tag,400,0,10,800,0,40000,sx3event.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dE_E_Cathodesx3"+tag,400,0,10,800,0,10000,sx3event.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmiscs_dE3_E_Anodesx3"+tag,400,0,10,400,0,40000,sx3event.Energy1,pcevent.Energy1*sinTheta_customV*3.,pmlabel);
						plotter->Fill2D("pmiscs_dE3_E_Cathodesx3"+tag,400,0,10,400,0,10000,sx3event.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmiscs_dPhi_sx3_PC"+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,sx3event.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmiscs_dt_Anode_sx3_PC"+std::to_string(phicut)+tag,600,-2000,2000,pcevent.Time1-sx3event.Time1,pmlabel);
						plotter->Fill1D("pmiscs_dt_Cathode_sx3"+tag,600,-2000,2000,pcevent.Time2-sx3event.Time1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Anode_E_sx3_PC"+std::to_string(phicut)+tag,600,-2000,2000,400,0,10,pcevent.Time1-sx3event.Time1,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Cathode_E_sx3"+tag,600,-2000,2000,400,0,10,pcevent.Time2-sx3event.Time1,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Cathodesx3_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time2-sx3event.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmiscs_pczfix"+tag,600,-300,300,pcz_fix,pmlabel);
						plotter->Fill1D("pmiscs_pcz"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
	
						//double path_length_q = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_s = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;
						double sx3Efix = cm_to_MeV->Eval(MeV_to_cm->Eval(sx3event.Energy1)-path_length_s);
	
						//plotter->Fill2D("sx3Ef_sx3E_matrix_all"+tag,400,0,10,400,0,10,sx3Efix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dE3_Ef_Anodesx3"+tag,400,0,10,400,0,40000,sx3Efix,pcevent.Energy1*sinTheta_customV*3,pmlabel);
						plotter->Fill2D("pmiscs_dE3_Ef_Cathodesx3"+tag,400,0,10,400,0,10000,sx3Efix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill2D("pmiscs_Ef_vs_theta_sx3"+tag,100,0,180,800,0,20,theta_s*180/M_PI,sx3Efix,pmlabel);
						plotter->Fill1D("pmiscs_VertexReconZ"+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmiscs_VertexReconXY"+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);						
						plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef"+tag,800,-400,400,800,0,20,vertex_z,sx3Efix,pmlabel);
						plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef"+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,sx3Efix,pmlabel);
					};
					
					plot_with_tag();
					if(cathode_alpha_select) 
						plot_with_tag(globaltag+"_cathode_alphas");
					else 
						plot_with_tag(globaltag+"_cathode_protons");
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end sx3Events loop
}


void miscHistograms_27Alaa_sx3(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events, std::string globaltag="") {
		//consider the 'proton-like' QQQ branch seen in a,p data
		for(auto sx3event: SX3_Events) {
				//if(sx3event.Energy1 < 1.2) continue; //coarse gating
				//if(sx3event.Energy1 > 5.0) continue; //coarse gating
				for(auto pcevent: PC_Events) {
					if(!(pcevent.multi1==1 && pcevent.multi2==2)) continue;
					//if(pcevent.Energy1 > 11000) continue; //coarse gating
					bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/3. && sx3event.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/3.;
					if(!phicut) continue;
					//if(pcevent.Time1-sx3event.Time1<-150 || pcevent.Time1-sx3event.Time1 >850) continue;
					
					double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(sx3event.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (sx3event.pos - TVector3(0,0,vertex_z)).Theta();

					if(r_rhoMin_fix.Perp()>10.0) continue;
					double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_s);
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);

					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = "27Alaa+miscsx3"+tag;
						plotter->Fill2D("pmiscs_dE_E_Anodesx3"+tag,400,0,10,800,0,40000,sx3event.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dE_E_Cathodesx3"+tag,400,0,10,800,0,10000,sx3event.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmiscs_dE3_E_Anodesx3"+tag,400,0,10,400,0,40000,sx3event.Energy1,pcevent.Energy1*sinTheta_customV*3.,pmlabel);
						plotter->Fill2D("pmiscs_dE3_E_Cathodesx3"+tag,400,0,10,400,0,10000,sx3event.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmiscs_dPhi_sx3_PC"+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,sx3event.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmiscs_dt_Anode_sx3_PC"+std::to_string(phicut)+tag,600,-2000,2000,pcevent.Time1-sx3event.Time1,pmlabel);
						plotter->Fill1D("pmiscs_dt_Cathode_sx3"+tag,600,-2000,2000,pcevent.Time2-sx3event.Time1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Anode_E_sx3_PC"+std::to_string(phicut)+tag,600,-2000,2000,400,0,10,pcevent.Time1-sx3event.Time1,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Cathode_E_sx3"+tag,600,-2000,2000,400,0,10,pcevent.Time2-sx3event.Time1,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Cathodesx3_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time2-sx3event.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmiscs_pczfix"+tag,600,-300,300,pcz_fix,pmlabel);
						plotter->Fill1D("pmiscs_pcz"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
	
						//double path_length_q = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_s = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;
						double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1)-path_length_s);
						double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
						double beam_energy_at_vertex = cm_to_MeV_27Al->Eval(MeV_to_cm_27Al->Eval(72.0) - beam_path_length);
						//Kinematics aakin_27Al(26.9815384,4.002603254,4.002603254,26.9815384,beam_energy_at_vertex/26.9815384); //m3 is alpha
						Kinematics apkin_27Al(26.9815384,4.002603254,1.00727647,29.9661,beam_energy_at_vertex/26.97441); //m3 is alpha
						auto kin_27Al = apkin_27Al;

						//plotter->Fill2D("sx3Ef_sx3E_matrix_all"+tag,400,0,10,400,0,10,sx3Efix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dE3_Ef_Anodesx3"+tag,400,0,10,400,0,40000,sx3Efix,pcevent.Energy1*sinTheta_customV*3,pmlabel);
						plotter->Fill2D("pmiscs_dE3_Ef_Cathodesx3"+tag,400,0,10,400,0,10000,sx3Efix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill2D("pmiscs_Ef_vs_theta_sx3"+tag,100,0,180,800,0,20,theta_s*180/M_PI,sx3Efix,pmlabel);
						plotter->Fill1D("pmisc_sx3_Ex_from_alpha"+tag,400,-20,20,kin_27Al.getExc(sx3Efix,theta_s*180/M_PI),pmlabel);

						plotter->Fill1D("pmiscs_VertexReconZ"+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmiscs_VertexReconXY"+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);						
						plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef"+tag,800,-400,400,800,0,20,vertex_z,sx3Efix,pmlabel);
						plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef"+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,sx3Efix,pmlabel);
					};
					
					plot_with_tag();
					if(cathode_alpha_select) 
						plot_with_tag(globaltag+"_cathode_alphas");
					else 
						plot_with_tag(globaltag+"_cathode_protons");
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end sx3Events loop
}

void miscHistograms_oneWire_27Alaa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters) {
		//consider the 'proton-like' QQQ branch seen in a,p data
		TRandom3 rand;
		rand.SetSeed();//random seed set
		bool sort_alphas=false;
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 1.0) continue; //coarse gating
				//if(qqqevent.Energy1 < 9.0) continue; //coarse gating
					for(const auto acluster: aClusters) {
						auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster,"ANODE");
						//if(apSumE<9000) continue;
						int a_number = acluster.size();
						if(a_number != 2) continue;
						TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire,qqqevent.pos.Phi());
						pc_closest.SetZ(rand.Gaus(pc_closest.Z(),8.0)); //dither
						plotter->Fill1D("dt_anode_interp_qqq",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
						if(qqqevent.Time1 - apTSMaxE < 150) {
							bool phicut = qqqevent.pos.Phi() <= pc_closest.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pc_closest.Phi()-TMath::Pi()/4.;
							TVector3 x2(pc_closest), x1(qqqevent.pos);
							TVector3 v = x2-x1;
							double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
							TVector3 r_rhoMin_fix = x1 + t_minimum*v;

							double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
		 					double sinTheta2 = TMath::Sin(theta_q);
		 					//if(theta_q*180/M_PI < 25 || theta_q*180/M_PI > 90.0 ) continue;
							if(r_rhoMin_fix.Perp()>10.0) continue;
							//if(r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100) continue;
							//if(r_rhoMin_fix.Z() < -50) continue;
							if(qqqevent.pos.Phi()*180./M_PI > 55 && qqqevent.pos.Phi()*180./M_PI < 72) continue; 
							//bool mystery = (plotter->FindCut("mystery_gate")->IsInside(qqqevent.Energy1,apSumE));
							//if(!mystery) continue;
							
							if(!phicut) continue;
							if(r_rhoMin_fix.Z()<-200) 
								sort_alphas=false;
							else
								sort_alphas=true;
								
							int zbin = std::floor(r_rhoMin_fix.Z()/50.0);
							//if(zbin!=-4) continue;
							
							plotter->Fill1D("dt_anode_ainterp_qqq_gated",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
							plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_qqqE",800,-2000,2000,800,0,10,qqqevent.Time1 - apTSMaxE,qqqevent.Energy1,"ainterp_noc");
							plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size()),400,0,10,800,0,40000,qqqevent.Energy1,apSumE,"ainterp_noc");
							plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size())+"_zbin50="+std::to_string(zbin),400,0,10,800,0,40000,qqqevent.Energy1,apSumE,"ainterp_noc");
							plotter->Fill2D("pcPhi_ainterp_qqqPhi_TC1_ignC_a"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,qqqevent.pos.Phi()*180./M_PI,"ainterp_noc");
							plotter->Fill2D("pcZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),300,-100,200,400,-200,200,qqqevent.pos.Z(),pc_closest.Z(),"ainterp_noc");

							//plotter->Fill2D("pcZ_ainterp_qqqpczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc");
							plotter->Fill2D("dEa3_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),400,0,10,800,0,40000,qqqevent.Energy1,apSumE*sinTheta2*3.,"ainterp_noc");
							plotter->Fill2D("vertexZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,800,-400,400,qqqevent.pos.Z(),r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill1D("vertexZ1d_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),800,-400,400,r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),"ainterp_noc");

							double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
							double qqqEfix;
							if(sort_alphas)
								qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
							else
								qqqEfix = cm_to_MeVd->Eval(MeV_to_cm_d->Eval(qqqevent.Energy1)-path_length_q);
							
							double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
							double beam_energy_at_vertex = cm_to_MeV_27Al->Eval(MeV_to_cm_27Al->Eval(56.155950) - beam_path_length);
							Kinematics aakin_27Al(26.9815384,4.002603254,4.002603254,26.9815384,beam_energy_at_vertex/26.9815384); //m3 is alpha
							Kinematics apkin_27Al(26.9815384,4.002603254,1.00727647,29.9661,beam_energy_at_vertex/26.9815384); //m3 is proton
							Kinematics adkin_27Al(26.9815384,4.002603254, 2.014101777844, 28.97649466434,beam_energy_at_vertex/26.9815384);
							Kinematics atkin_27Al(26.9815384,4.002603254, 3.016049, 27.97692653442,beam_energy_at_vertex/26.9815384);
							//Kinematics a3hekin_27Al(26.9815384,4.002603254, 3.016029321967, 26.981538408,beam_energy_at_vertex/26.9815384);

							if(sort_alphas) {
								auto kin_27Al = aakin_27Al;
								plotter->Fill1D("pmisc_ow_Ex_from_alpha",700,-20,50,kin_27Al.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
								plotter->Fill2D("pmisc_ow_Ex_from_alpha_vs_VertexZs",800,-400,400,700,-20,50,r_rhoMin_fix.Z(),kin_27Al.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
							} else {
								auto kin_27Al = adkin_27Al;
								plotter->Fill1D("pmisc_ow_Ex_from_d",700,-20,50,kin_27Al.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
								plotter->Fill2D("pmisc_ow_Ex_from_d_vs_VertexZs",800,-400,400,700,-20,50,r_rhoMin_fix.Z(),kin_27Al.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
							}
							plotter->Fill2D("pmisc_ow_Ef_vs_theta_qqq",100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,"ainterp_noc");
							plotter->Fill2D("pmisc_ow_E_vs_theta_qqq",100,0,180,800,0,20,theta_q*180/M_PI,qqqevent.Energy1,"ainterp_noc");
							plotter->Fill2D("pmisc_ow_VertexReconZ_vs_Ef",800,-400,400,800,0,20,r_rhoMin_fix.Z(),qqqEfix,"ainterp_noc");
							plotter->Fill2D("pmisc_ow_VertexReconZ_vs_E",800,-400,400,800,0,20,r_rhoMin_fix.Z(),qqqevent.Energy1,"ainterp_noc");
						}	
					}
		}//end QQQEvents loop


		for(auto sx3event: SX3_Events) {
				if(sx3event.Energy1 < 1.0) continue; //coarse gating
				if(sx3event.pos.Phi()*180/M_PI>-50 && sx3event.pos.Phi()*180/M_PI < 0 ) continue;
				//if(sx3event.Energy1 < 8.0) continue; //coarse gating
					for(const auto acluster: aClusters) {
						auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster,"ANODE");
						//if(apSumE<9000) continue;
						int a_number = acluster.size();
						//if(a_number != 2) continue;
						TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire,sx3event.pos.Phi());
						pc_closest.SetZ(rand.Gaus(pc_closest.Z(),8.0)); //dither
						plotter->Fill1D("pmiscsx3_ow_dt_anode_interp_sx3",800,-2000,2000,sx3event.Time1 - apTSMaxE,"ainterp_noc_sx3");
						if(sx3event.Time1 - apTSMaxE < 150) {
							bool phicut = sx3event.pos.Phi() <= pc_closest.Phi()+TMath::Pi()/4. && sx3event.pos.Phi() >= pc_closest.Phi()-TMath::Pi()/4.;
							TVector3 x2(pc_closest), x1(sx3event.pos);
							TVector3 v = x2-x1;
							double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
							TVector3 r_rhoMin_fix = x1 + t_minimum*v;
													
							double theta_q = (sx3event.pos - r_rhoMin_fix).Theta();
		 					double sinTheta2 = TMath::Sin(theta_q);
							//if(theta_q*180/M_PI < 25 || theta_q*180/M_PI > 90.0 ) continue;

							if(r_rhoMin_fix.Perp()>10.0) continue;
							//if(r_rhoMin_fix.Z() < -50) continue;
							//if(r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100) continue;
							if(!phicut) continue;
							
							//if(r_rhoMin_fix.Z()<-200) 
							if(apSumE<3000)
								sort_alphas=false;
							else
								sort_alphas=true;
							int zbin = std::floor(r_rhoMin_fix.Z()/50.0);
							//if(zbin!=-4) continue;

							plotter->Fill1D("pmiscsx3_ow_dt_anode_ainterp_sx3_gated",800,-2000,2000,sx3event.Time1 - apTSMaxE,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_dt_anode_ainterp_sx3_gated_vs_apSumE",800,-2000,2000,800,0,40000,sx3event.Time1 - apTSMaxE,apSumE,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_dt_anode_ainterp_sx3_gated_vs_sx3E",800,-2000,2000,800,0,10,sx3event.Time1 - apTSMaxE,sx3event.Energy1,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_dEa_ainterp_Esx3_TC1_ignC_a"+std::to_string(acluster.size()),400,0,10,800,0,40000,sx3event.Energy1,apSumE,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_dEa_ainterp_Esx3_TC1_ignC_a"+std::to_string(acluster.size())+"_zbin50="+std::to_string(zbin),400,0,10,800,0,40000,sx3event.Energy1,apSumE,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_pcPhi_ainterp_sx3Phi_TC1_ignC_a"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,sx3event.pos.Phi()*180./M_PI,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_pcZ_ainterp_sx3Z_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),300,-100,200,400,-200,200,sx3event.pos.Z(),pc_closest.Z(),"ainterp_noc_sx3");
							//plotter->Fill2D("pmiscsx3_ow_pcZ_ainterp_sx3pczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_dEa3_ainterp_Esx3_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),400,0,10,800,0,40000,sx3event.Energy1,apSumE*sinTheta2*3.,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_dEa3_ainterp_Esx3_TC1_ignC_a"+std::to_string(acluster.size())+"_PC1_"+std::string(sort_alphas?"a":"d"),400,0,10,800,0,40000,sx3event.Energy1,apSumE*sinTheta2*3.,"ainterp_noc_sx3");

							plotter->Fill2D("pmiscsx3_ow_vertexZ_ainterp_sx3Z_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,800,-400,400,sx3event.pos.Z(),r_rhoMin_fix.Z(),"ainterp_noc_sx3");
							plotter->Fill1D("pmiscsx3_ow_vertexZ1d_ainterp_sx3Z_TC1_ignC_a"+std::to_string(acluster.size()),800,-400,400,r_rhoMin_fix.Z(),"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_vertexXY_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),"ainterp_noc_sx3");

							double path_length_q = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;
							double sx3Efix;
							if(sort_alphas)
								sx3Efix = cm_to_MeV->Eval(MeV_to_cm->Eval(sx3event.Energy1)-path_length_q);
							else
								sx3Efix = cm_to_MeVd->Eval(MeV_to_cm_d->Eval(sx3event.Energy1)-path_length_q);
							
							double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
		  					double beam_energy_at_vertex = cm_to_MeV_27Al->Eval(MeV_to_cm_27Al->Eval(56.155950) - beam_path_length);
							Kinematics aakin_27Al(26.9815384,4.002603254,4.002603254,26.9815384,beam_energy_at_vertex/26.9815384); //m3 is alpha
	  						Kinematics apkin_27Al(26.9815384,4.002603254,1.00727647,29.9661,beam_energy_at_vertex/26.9815384); //m3 is alpha
							Kinematics adkin_27Al(26.9815384,4.002603254, 2.014101777844, 28.97649466434,beam_energy_at_vertex/26.9815384);
							Kinematics atkin_27Al(26.9815384,4.002603254, 3.016049, 27.97692653442,beam_energy_at_vertex/26.9815384);

							if(sort_alphas) {
								auto kin_27Al = aakin_27Al;
								plotter->Fill1D("pmiscsx3_ow_Ex_from_alpha",700,-20,50,kin_27Al.getExc(sx3Efix,theta_q*180/M_PI),"ainterp_noc_sx3");
								plotter->Fill2D("pmiscsx3_ow_Ex_from_alpha_vs_VertexZs",800,-400,400,700,-20,50,r_rhoMin_fix.Z(),kin_27Al.getExc(sx3Efix,theta_q*180/M_PI),"ainterp_noc_sx3");
							} else {
								auto kin_27Al = adkin_27Al;
								plotter->Fill1D("pmiscsx3_ow_Ex_from_d",700,-20,50,kin_27Al.getExc(sx3Efix,theta_q*180/M_PI),"ainterp_noc_sx3");
								plotter->Fill2D("pmiscsx3_ow_Ex_from_d_vs_VertexZs",800,-400,400,700,-20,50,r_rhoMin_fix.Z(),kin_27Al.getExc(sx3Efix,theta_q*180/M_PI),"ainterp_noc_sx3");
							}
							plotter->Fill2D("pmiscsx3_ow_Ef_vs_theta_sx3",100,0,180,800,0,20,theta_q*180/M_PI,sx3Efix,"ainterp_noc_sx3");
							plotter->Fill2D("pmiscsx3_ow_VertexReconZ_vs_Ef",800,-400,400,800,0,20,r_rhoMin_fix.Z(),sx3Efix,"ainterp_noc_sx3");
						}	
					}
		}//end sx3Events loop
}



void miscHistograms_27Alaa(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events) {
		TRandom3 rand;
		rand.SetSeed();//random seed set
		//Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,7.0); //m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.2, 250 torr)
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 0.6) continue; //coarse gating
				//if(qqqevent.Energy1 > 5.0) continue; //coarse gating
				for(auto pcevent: PC_Events) {
					if(!(pcevent.multi1==1 && pcevent.multi2<=2)) continue;
					plotter->Fill2D("pmisc_dE_E_AnodeQQQ",400,0,10,800,0,40000,qqqevent.Energy1,pcevent.Energy1,"proton+misc");
					plotter->Fill2D("pmisc_dE_E_CathodeQQQ",400,0,10,800,0,10000,qqqevent.Energy1,pcevent.Energy2,"proton+misc");

					//if(pcevent.Energy1 > 11000) continue; //coarse gating
					bool phicut = qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4.;
					if(!phicut) continue;
					//if(pcevent.Time1-qqqevent.Time1<-150 || pcevent.Time1-qqqevent.Time1 >850) continue;
					
					double pcz_fix, pcz_dith=pcevent.pos.Z();
					if(pcevent.multi2==2) 
						pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					else {
						pcz_fix = rand.Gaus(pcevent.pos.Z(),8.0);//dither for a1c1 events
						pcz_dith = pcz_fix;
					}
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(qqqevent.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
					double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_q);
					if(r_rhoMin_fix.Perp()>6) continue;
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);
					if(vertex_z < -173.6 || vertex_z > 100) continue;

					double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
		  			double beam_energy_at_vertex = cm_to_MeV_27Al->Eval(MeV_to_cm_27Al->Eval(72.0) - beam_path_length);
					Kinematics aakin_27Al(26.9815384,4.002603254,4.002603254,26.9815384,beam_energy_at_vertex/26.9815384); //m3 is alpha
	  				//Kinematics apkin_27Al(26.9815384,4.002603254,1.00727647,29.9661,beam_energy_at_vertex/26.97441); //m3 is alpha
	  				auto kin_27Al = aakin_27Al;
					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = "proton+misc"+tag;
						plotter->Fill2D("pmisc_dE_E_AnodeQQQ"+tag,400,0,10,800,0,40000,qqqevent.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE_E_CathodeQQQ"+tag,400,0,10,800,0,10000,qqqevent.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_AnodeQQQ"+tag,400,0,10,400,0,40000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV*3.,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_CathodeQQQ"+tag,400,0,10,400,0,10000,qqqevent.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmisc_dPhi_QQQ_PC"+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,qqqevent.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmisc_dt_Anode_QQQ_PC"+std::to_string(phicut)+tag,600,-2000,2000,pcevent.Time1-qqqevent.Time1,pmlabel);
						plotter->Fill1D("pmisc_dt_Cathode_QQQ"+tag,600,-2000,2000,pcevent.Time2-qqqevent.Time1,pmlabel);
						plotter->Fill2D("pmisc_dt_Anode_E_QQQ_PC"+std::to_string(phicut)+tag,600,-2000,2000,400,0,10,pcevent.Time1-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_AnodeQQQ_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time1-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill2D("pmisc_dt_Cathode_E_QQQ"+tag,600,-2000,2000,400,0,10,pcevent.Time2-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_CathodeQQQ_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time2-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmisc_pczfix"+tag,600,-300,300,pcz_fix,pmlabel);
						if(pcevent.multi2==2) {
							plotter->Fill1D("pmisc_pcz"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz2"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
						}
						if(pcevent.multi2==1) {
							plotter->Fill1D("pmisc_pcz"+tag,600,-300,300,pcz_fix,pmlabel);
							plotter->Fill1D("pmisc_pcz1"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz_dith"+tag,600,-300,300,pcz_dith,pmlabel);
						}

						//double path_length_q = (qqqevent.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
						double qqqEfix;
						qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
						plotter->Fill1D("pmisc_Ex_from_alpha",400,-20,20,kin_27Al.getExc(qqqEfix,theta_q*180/M_PI),pmlabel);
						//plotter->Fill2D("qqqEf_sx3E_matrix_all"+tag,400,0,10,400,0,10,qqqEfix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_AnodeQQQ"+tag,400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV*3,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_CathodeQQQ"+tag,400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill1D("pmisc_VertexReconZ"+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmisc_VertexReconXY"+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);						
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+tag,800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						
						plotter->Fill2D("pmisc_Ef_vs_theta_qqq"+tag,100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,pmlabel);
						if(pcevent.multi2==1) {
							plotter->Fill2D("pmisc_Ef_vs_theta_qqq_a1c1"+tag,100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,pmlabel);						
							plotter->Fill2D("pmisc_VertexReconZ_vs_Ef_a1c1"+tag,800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						}
					};
					
					if(cathode_alpha_select) 
						plot_with_tag("_cathode_alphas");
					//else 
					//	plot_with_tag("_cathode_protons");
					//plot_with_tag();
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(qqqEfix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end QQQEvents loop
}


void paMiscHistograms_oneWire(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters) {
		//consider the 'proton-like' QQQ branch seen in a,p data
		TRandom3 rand;
		rand.SetSeed();  //random seed set
		Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,6.88445); //m3 is alpha, 6.88445 MeV is 7.0 MeV proton energy after 8um kapton
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 0.6) continue; //coarse gating
				//if(qqqevent.Energy1 > 5.0) continue; //coarse gating
					for(const auto acluster: aClusters) {
						auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster,"ANODE");
						//if(apSumE<2000) continue;
						if(!plotter->FindCut("anode_qqq_alphas_27Al")->IsInside(qqqevent.Energy1,apSumE)) continue;
						int a_number = acluster.size();
						TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire,qqqevent.pos.Phi());
						plotter->Fill1D("dt_anode_interp_qqq",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
						if(qqqevent.Time1 - apTSMaxE < 150) {
							bool phicut = qqqevent.pos.Phi() <= pc_closest.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pc_closest.Phi()-TMath::Pi()/4.;
							
							TVector3 x2(pc_closest), x1(qqqevent.pos);
							TVector3 v = x2-x1;
							double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
							TVector3 r_rhoMin_fix = x1 + t_minimum*v;
							double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
		 					double sinTheta2 = TMath::Sin(theta_q)*3.0;

							if(r_rhoMin_fix.Perp()>6.0) continue;
							if(r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100) continue;
							if(!phicut) continue;
							plotter->Fill1D("dt_anode_ainterp_qqq_gated",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
							plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_qqqE",800,-2000,2000,800,0,10,qqqevent.Time1 - apTSMaxE,qqqevent.Energy1,"ainterp_noc");
							plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_anodeMaxE",800,-2000,2000,800,0,40000,qqqevent.Time1 - apTSMaxE,apMaxE,"ainterp_noc");
							plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_anodeSumE",800,-2000,2000,800,0,40000,qqqevent.Time1 - apTSMaxE,apSumE,"ainterp_noc");
							plotter->Fill2D("anodeSumE_vs_Thetaq",800,0,40000,180,0,180,apSumE,theta_q*180/M_PI,"ainterp_noc");
							plotter->Fill2D("anodeMaxE_vs_Thetaq",800,0,40000,180,0,180,apMaxE,theta_q*180/M_PI,"ainterp_noc");
							
							plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size()),400,0,10,800,0,40000,qqqevent.Energy1,apSumE,"ainterp_noc");
							plotter->Fill2D("dEaMax_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size()),400,0,10,800,0,40000,qqqevent.Energy1,apMaxE,"ainterp_noc");
							plotter->Fill2D("pcPhi_ainterp_qqqPhi_TC1_ignC_a"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,qqqevent.pos.Phi()*180./M_PI,"ainterp_noc");
							plotter->Fill1D("pcZ_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),400,-200,200,pc_closest.Z(),"ainterp_noc");

							//plotter->Fill2D("pcZ_ainterp_qqqpczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc");
							plotter->Fill2D("dEa3_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),400,0,10,800,0,40000,qqqevent.Energy1,apSumE*sinTheta2,"ainterp_noc");
			
							//plotter->Fill2D("vertexZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,800,-400,400,qqqevent.pos.Z(),r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill1D("vertexZ1d_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),800,-400,400,r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),"ainterp_noc");

							double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
							double qqqEfix;
							qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
							plotter->Fill1D("pmisc_ow_Ex_from_alpha",200,-10,10,apkin_a.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
							plotter->Fill2D("pmisc_ow_Ef_vs_theta_qqq",100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,"ainterp_noc");
							plotter->Fill2D("pmisc_ow_VertexReconZ_vs_Ef",800,-400,400,800,0,20,r_rhoMin_fix.Z(),qqqEfix,"ainterp_noc");
							//plotter->Fill2D("pmisc_ow_VertexReconRhoZ",800,-400,400,40,-20,20,r_rhoMin_fix.Z(),r_rhoMin_fix.Perp(),"ainterp_noc");

							double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
		  					double beam_energy_at_vertex = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(6.88445) - beam_path_length);
		  					//double beame_guess = (1+4.002603254/1.007825)*qqqEfix/(4.*(4.002603254/1.007825)*cos(theta_q)*cos(theta_q));
							double beame_guess = TMath::Power(1.007825+4.002603254,2)*qqqEfix/(4.*4.002603254*1.007825*cos(theta_q)*cos(theta_q));
		  					plotter->Fill2D("omisc_ow_beame_guess_vs_beame_catima",800,4,10,800,6,8,beame_guess,beam_energy_at_vertex,"ainterp_noc");
		  					Kinematics apkin_a_2(1.007825,4.002603254,4.002603254,1.007825,beam_energy_at_vertex/1.007825); //m3 is alpha
							plotter->Fill1D("pmisc_ow_Ex_from_alpha_beampathfix",200,-10,10,apkin_a_2.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
						}
					}
		}//end QQQEvents loop
}

/*
	Analyze p+a data for (p,a) alphas going into QQQ via PC, phi-phi correlation being true
*/
void paMiscHistograms(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events) {
		//consider the 'proton-like' QQQ branch seen in a,p data
		TRandom3 rand;
		rand.SetSeed();//random seed set
		double beame=7.0;
		if(dataset=="17F") 
			beame=6.6877; //havar 5um + kapton 8um and 100mm gas
		else
			beame=6.88445; //kapton 8um and gas only
		Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,beame); //m3 is alpha, 6.88445 MeV is 7.0 MeV proton energy after 8um kapton
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 0.6) continue; //coarse gating
				//if(qqqevent.Energy1 > 5.0) continue; //coarse gating
				for(auto pcevent: PC_Events) {

					if(!(pcevent.multi1==1 && pcevent.multi2<=2)) continue;
					plotter->Fill2D("pmisc_dE_E_AnodeQQQ",400,0,10,800,0,40000,qqqevent.Energy1,pcevent.Energy1,"proton+misc");
					plotter->Fill2D("pmisc_dE_E_CathodeQQQ",400,0,10,800,0,10000,qqqevent.Energy1,pcevent.Energy2,"proton+misc");

					//if(pcevent.Energy1 < 2000) continue; //coarse gating
					bool phicut = qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4.;
					if(!phicut) continue;
					if(pcevent.Time1-qqqevent.Time1<-150 || pcevent.Time1-qqqevent.Time1 >850) continue;
					
					double pcz_fix, pcz_dith=pcevent.pos.Z();
					if(pcevent.multi2==2) 
						pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					else {
						pcz_fix = rand.Gaus(pcevent.pos.Z(),8.0);//dither for a1c1 events
						pcz_dith = pcz_fix;
					}
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(qqqevent.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
					double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_q)*3.0;
					if(r_rhoMin_fix.Perp()>6) continue;
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);
					//if(vertex_z < -173.6 || vertex_z > 100) continue;


					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = "proton+misc"+tag;
						plotter->Fill2D("pmisc_dE_E_AnodeQQQ"+tag,400,0,10,800,0,40000,qqqevent.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE_Theta_anodeQQQ"+tag,180,0,180,800,0,40000,theta_q*180/M_PI,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE_E_CathodeQQQ"+tag,400,0,10,800,0,10000,qqqevent.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_AnodeQQQ"+tag,400,0,10,400,0,40000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_CathodeQQQ"+tag,400,0,10,400,0,10000,qqqevent.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmisc_dPhi_QQQ_PC"+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,qqqevent.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmisc_dt_Anode_QQQ_PC"+std::to_string(phicut)+tag,600,-2000,2000,pcevent.Time1-qqqevent.Time1,pmlabel);
						plotter->Fill1D("pmisc_dt_Cathode_QQQ"+tag,600,-2000,2000,pcevent.Time2-qqqevent.Time1,pmlabel);
						plotter->Fill2D("pmisc_dt_Anode_E_QQQ_PC"+std::to_string(phicut)+tag,600,-2000,2000,400,0,10,pcevent.Time1-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_AnodeQQQ_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time1-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill2D("pmisc_dt_Cathode_E_QQQ"+tag,600,-2000,2000,400,0,10,pcevent.Time2-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_CathodeQQQ_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time2-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmisc_pczfix"+tag,600,-300,300,pcz_fix,pmlabel);
						if(pcevent.multi2==2) {
							plotter->Fill1D("pmisc_pcz"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz2"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
						}
						if(pcevent.multi2==1) {
							plotter->Fill1D("pmisc_pcz"+tag,600,-300,300,pcz_fix,pmlabel);
							plotter->Fill1D("pmisc_pcz1"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz_dith"+tag,600,-300,300,pcz_dith,pmlabel);
						}

						//double path_length_q = (qqqevent.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
						double qqqEfix;
						if(tag == "_cathode_alphas") {//satisfied when find succeeds 
							qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
							plotter->Fill1D("pmisc_Ex_from_alpha",200,-10,10,apkin_a.getExc(qqqEfix,theta_q*180/M_PI),pmlabel);
							
							double beam_path_length = TMath::Abs(r_rhoMin_fix.Z() - z_entrance) * 0.1;
		  					double beam_energy_at_vertex = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(6.88445) - beam_path_length);
		  					Kinematics apkin_a_2(1.007825,4.002603254,4.002603254,1.007825,beam_energy_at_vertex/1.007825); //m3 is alpha
							plotter->Fill1D("pmisc_Ex_from_alpha_beampathfix",200,-10,10,apkin_a_2.getExc(qqqEfix,theta_q*180/M_PI),pmlabel);
							double beame_guess = TMath::Power(1.007825+4.002603254,2)*qqqEfix/(4.*4.002603254*1.007825*cos(theta_q)*cos(theta_q));
		  					plotter->Fill2D("pmisc_beame_guess_vs_beame_catima",800,4,10,800,6,8,beame_guess,beam_energy_at_vertex,pmlabel);
						}
						else
							qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1)-path_length_q);
						//plotter->Fill2D("qqqEf_sx3E_matrix_all"+tag,400,0,10,400,0,10,qqqEfix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_AnodeQQQ"+tag,400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_CathodeQQQ"+tag,400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill1D("pmisc_VertexReconZ"+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmisc_VertexReconXY"+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);	
						//plotter->Fill2D("pmisc_VertexReconRhoZ"+tag,800,-400,400,40,-20,20,r_rhoMin_fix.Z(),r_rhoMin_fix.Perp(),pmlabel);
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+tag,800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						
						plotter->Fill2D("pmisc_Ef_vs_theta_qqq"+tag,100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,pmlabel);
						if(pcevent.multi2==1) {
							plotter->Fill2D("pmisc_Ef_vs_theta_qqq_a1c1"+tag,100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,pmlabel);						
							plotter->Fill2D("pmisc_VertexReconZ_vs_Ef_a1c1"+tag,800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						}

					};
					
					if(cathode_alpha_select) 
						plot_with_tag("_cathode_alphas");
					//else 
					//	plot_with_tag("_cathode_protons");
					//plot_with_tag();
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(qqqEfix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end QQQEvents loop
}

/*
	Analyze p,p data from run17 for protons going into QQQ via PC, phi-phi correlation being true
*/
void ppMiscHistograms(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events) {
		//consider the 'proton-like' QQQ branch seen in a,p data
		TRandom3 rand;
		rand.SetSeed();//random seed set
		Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,7.0); //m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.2, 250 torr)
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 6.6) continue; //coarse gating
				for(auto pcevent: PC_Events) {
					if(!(pcevent.multi1==1 && pcevent.multi2<=2)) continue;
					//if(pcevent.Energy1 > 11000) continue; //coarse gating
					bool phicut = qqqevent.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/4.;
					if(!phicut) continue;
					//if(pcevent.Time1-qqqevent.Time1<-150 || pcevent.Time1-qqqevent.Time1 >850) continue;
					
					double pcz_fix, pcz_dith=pcevent.pos.Z();
					if(pcevent.multi2==2) 
						pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					else {
						pcz_fix = rand.Gaus(pcevent.pos.Z(),8.0);//dither for a1c1 events
						pcz_dith = pcz_fix;
					}
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(qqqevent.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
					double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_q)*3.0;
					//if(r_rhoMin_fix.Perp()>6) continue;
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);
					if(vertex_z < -173.6 || vertex_z > 100) continue;

					double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z())/TMath::Tan((qqqevent.pos-TVector3(0,0,source_vertex)).Theta())  + source_vertex;
					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = "proton+misc"+tag;
						plotter->Fill2D("pmisc_dE_E_AnodeQQQ"+tag,400,0,10,800,0,40000,qqqevent.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE_E_CathodeQQQ"+tag,400,0,10,800,0,10000,qqqevent.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_AnodeQQQ"+tag,400,0,10,400,0,40000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV*3.,pmlabel);
						plotter->Fill2D("pmisc_dE3_E_CathodeQQQ"+tag,400,0,10,400,0,10000,qqqevent.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmisc_dPhi_QQQ_PC"+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,qqqevent.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmisc_dt_Anode_QQQ_PC"+std::to_string(phicut)+tag,600,-2000,2000,pcevent.Time1-qqqevent.Time1,pmlabel);
						plotter->Fill1D("pmisc_dt_Cathode_QQQ"+tag,600,-2000,2000,pcevent.Time2-qqqevent.Time1,pmlabel);
						plotter->Fill2D("pmisc_dt_Anode_E_QQQ_PC"+std::to_string(phicut)+tag,600,-2000,2000,400,0,10,pcevent.Time1-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_AnodeQQQ_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time1-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill2D("pmisc_dt_Cathode_E_QQQ"+tag,600,-2000,2000,400,0,10,pcevent.Time2-qqqevent.Time1,qqqevent.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dt_CathodeQQQ_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time2-qqqevent.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmisc_pczfix"+tag,600,-300,300,pcz_fix,pmlabel);
						if(pcevent.multi2==2) {
							plotter->Fill1D("pmisc_pcz"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill2D("pmisc_pcz_vs_pczguess"+tag,600,-300,300,600,-300,300,pcz_guess_int,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz2"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
						}
						if(pcevent.multi2==1) {
							plotter->Fill1D("pmisc_pcz"+tag,600,-300,300,pcz_fix,pmlabel);
							plotter->Fill1D("pmisc_pcz1"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
							plotter->Fill2D("pmisc_pcz_vs_pczguess"+tag,600,-300,300,600,-300,300,pcz_guess_int,pcevent.pos.Z(),pmlabel);
							plotter->Fill1D("pmisc_pcz_dith"+tag,600,-300,300,pcz_dith,pmlabel);
						}

						//double path_length_q = (qqqevent.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
						double qqqEfix;
						if(tag == "_cathode_alphas") {//satisfied when find succeeds 
							qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
							plotter->Fill1D("pmisc_Ex_from_alpha",200,-10,10,apkin_a.getExc(qqqEfix,theta_q*180/M_PI),pmlabel);
						}
						else
							qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1)-path_length_q);
						//plotter->Fill2D("qqqEf_sx3E_matrix_all"+tag,400,0,10,400,0,10,qqqEfix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_AnodeQQQ"+tag,400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV*3,pmlabel);
						plotter->Fill2D("pmisc_dE3_Ef_CathodeQQQ"+tag,400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill1D("pmisc_VertexReconZ"+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmisc_VertexReconXY"+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);						
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+tag,800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						plotter->Fill2D("pmisc_VertexReconZ_vs_Ef"+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						
						plotter->Fill2D("pmisc_Ef_vs_theta_qqq"+tag,100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,pmlabel);
						if(pcevent.multi2==1) {
							plotter->Fill2D("pmisc_Ef_vs_theta_qqq_a1c1"+tag,100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,pmlabel);						
							plotter->Fill2D("pmisc_VertexReconZ_vs_Ef_a1c1"+tag,800,-400,400,800,0,20,vertex_z,qqqEfix,pmlabel);
						}
					};
					
					if(cathode_alpha_select) 
						plot_with_tag("_cathode_alphas");
					else 
						plot_with_tag("_cathode_protons");
					plot_with_tag();
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(qqqEfix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end QQQEvents loop
}

//p,p data, has correlated PC/QQQ only in run17, 27Al data
void ppMiscHistograms_oneWire(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<std::vector<std::tuple<int,double,double>>> aClusters) {
		TRandom3 rand;
		rand.SetSeed();//random seed set
		Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,7.0); //m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.2, 250 torr)
		for(auto qqqevent: QQQ_Events) {
				if(qqqevent.Energy1 < 6.6) continue; //coarse gating
					for(const auto acluster: aClusters) {
						auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster,"ANODE");
						//if(apSumE<6000) continue;
						int a_number = acluster.size();
						TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire,qqqevent.pos.Phi());
						plotter->Fill1D("dt_anode_interp_qqq",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
						if(qqqevent.Time1-apTSMaxE>-550 && qqqevent.Time1 - apTSMaxE < -150) {
							bool phicut = qqqevent.pos.Phi() <= pc_closest.Phi()+TMath::Pi()/4. && qqqevent.pos.Phi() >= pc_closest.Phi()-TMath::Pi()/4.;
							
							TVector3 x2(pc_closest), x1(qqqevent.pos);
							TVector3 v = x2-x1;
							double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
							TVector3 r_rhoMin_fix = x1 + t_minimum*v;
													
							double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
		 					double sinTheta2 = TMath::Sin(theta_q)*3;

							if(r_rhoMin_fix.Perp()>6.0) continue;
							//if(r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100) continue;
							if(!phicut) continue;
							double pcz_guess_int = z_to_crossover_rho(pc_closest.Z())/TMath::Tan((qqqevent.pos-TVector3(0,0,source_vertex)).Theta())  + source_vertex;

							plotter->Fill1D("dt_anode_ainterp_qqq_gated",800,-2000,2000,qqqevent.Time1 - apTSMaxE,"ainterp_noc");
							plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_qqqE",800,-2000,2000,800,0,10,qqqevent.Time1 - apTSMaxE,qqqevent.Energy1,"ainterp_noc");
							plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size()),400,0,10,800,0,40000,qqqevent.Energy1,apSumE,"ainterp_noc");
							plotter->Fill2D("pcPhi_ainterp_qqqPhi_TC1_ignC_a"+std::to_string(acluster.size()),120,-360,360,120,-360,360,pc_closest.Phi()*180./M_PI,qqqevent.pos.Phi()*180./M_PI,"ainterp_noc");
							plotter->Fill2D("pcZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),300,-100,200,400,-200,200,qqqevent.pos.Z(),pc_closest.Z(),"ainterp_noc");
			
							plotter->Fill2D("pmisc_ow_pcz_vs_pczguess",600,-300,300,600,-300,300,pcz_guess_int,pc_closest.Z(),"ainterp_noc");

							//plotter->Fill2D("pcZ_ainterp_qqqpczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc");
							plotter->Fill2D("dEa3_ainterp_Eqqq_TC1_ignC_a"+std::to_string(acluster.size())+"_PC"+std::to_string(phicut),1200,0,30,800,0,30000,qqqevent.Energy1,apSumE*sinTheta2,"ainterp_noc");
			
							plotter->Fill2D("vertexZ_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,800,-400,400,qqqevent.pos.Z(),r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill1D("vertexZ1d_ainterp_qqqZ_TC1_ignC_a"+std::to_string(acluster.size()),800,-400,400,r_rhoMin_fix.Z(),"ainterp_noc");
							plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a"+std::to_string(acluster.size()),200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),"ainterp_noc");

							double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
							double qqqEfix;
							qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1)-path_length_q);
							plotter->Fill1D("pmisc_ow_Ex_from_alpha",200,-10,10,apkin_a.getExc(qqqEfix,theta_q*180/M_PI),"ainterp_noc");
							plotter->Fill2D("pmisc_ow_Ef_vs_theta_qqq",100,0,180,800,0,20,theta_q*180/M_PI,qqqEfix,"ainterp_noc");
							plotter->Fill2D("pmisc_ow_VertexReconZ_vs_Ef",800,-400,400,800,0,20,r_rhoMin_fix.Z(),qqqEfix,"ainterp_noc");
						}	
					}
		}//end QQQEvents loop
}


//Half baked, doesn't work
void ppMiscHistograms_sx3(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events) {
		//consider the 'proton-like' QQQ branch seen in a,p data
		for(auto sx3event: SX3_Events) {
				if(sx3event.Energy1 < 1.2) continue; //coarse gating
				//if(sx3event.Energy1 > 5.0) continue; //coarse gating
				for(auto pcevent: PC_Events) {
					if(!(pcevent.multi1==1 && pcevent.multi2==2)) continue;
					//if(pcevent.Energy1 > 11000) continue; //coarse gating
					bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi()+TMath::Pi()/3. && sx3event.pos.Phi() >= pcevent.pos.Phi()-TMath::Pi()/3.;
					if(!phicut) continue;
					//if(pcevent.Time1-sx3event.Time1<-150 || pcevent.Time1-sx3event.Time1 >850) continue;
					
					double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(sx3event.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (sx3event.pos - TVector3(0,0,vertex_z)).Theta();

					if(r_rhoMin_fix.Perp()>10.0) continue;
					double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_s);
					bool cathode_alpha_select = (pcevent.Energy2 > 1400);

		

					//What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
					auto plot_with_tag = [&](std::string tag="") {
						std::string pmlabel = "proton+miscsx3"+tag;
						plotter->Fill2D("pmiscs_dE_E_Anodesx3"+tag,400,0,10,800,0,40000,sx3event.Energy1,pcevent.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dE_E_Cathodesx3"+tag,400,0,10,800,0,10000,sx3event.Energy1,pcevent.Energy2,pmlabel);
						plotter->Fill2D("pmiscs_dE3_E_Anodesx3"+tag,400,0,10,400,0,40000,sx3event.Energy1,pcevent.Energy1*sinTheta_customV*3.,pmlabel);
						plotter->Fill2D("pmiscs_dE3_E_Cathodesx3"+tag,400,0,10,400,0,10000,sx3event.Energy1,pcevent.Energy2*sinTheta_customV,pmlabel);
						plotter->Fill2D("pmiscs_dPhi_sx3_PC"+tag,180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,sx3event.pos.Phi()*180/M_PI,pmlabel);
						plotter->Fill1D("pmiscs_dt_Anode_sx3_PC"+std::to_string(phicut)+tag,600,-2000,2000,pcevent.Time1-sx3event.Time1,pmlabel);
						plotter->Fill1D("pmiscs_dt_Cathode_sx3"+tag,600,-2000,2000,pcevent.Time2-sx3event.Time1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Anode_E_sx3_PC"+std::to_string(phicut)+tag,600,-2000,2000,400,0,10,pcevent.Time1-sx3event.Time1,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Cathode_E_sx3"+tag,600,-2000,2000,400,0,10,pcevent.Time2-sx3event.Time1,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dt_Cathodesx3_vsPCPhi"+tag,600,-2000,2000,180,-360,360,pcevent.Time2-sx3event.Time1,pcevent.pos.Phi()*180./M_PI,pmlabel);
						plotter->Fill1D("pmiscs_pczfix"+tag,600,-300,300,pcz_fix,pmlabel);
						plotter->Fill1D("pmiscs_pcz"+tag,600,-300,300,pcevent.pos.Z(),pmlabel);
	
						//double path_length_q = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
						double path_length_s = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;
						double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1)-path_length_s);
	
						//plotter->Fill2D("sx3Ef_sx3E_matrix_all"+tag,400,0,10,400,0,10,sx3Efix,sx3event.Energy1,pmlabel);
						plotter->Fill2D("pmiscs_dE3_Ef_Anodesx3"+tag,400,0,10,400,0,40000,sx3Efix,pcevent.Energy1*sinTheta_customV*3,pmlabel);
						plotter->Fill2D("pmiscs_dE3_Ef_Cathodesx3"+tag,400,0,10,400,0,10000,sx3Efix,pcevent.Energy2*sinTheta_customV,pmlabel);
						
						plotter->Fill2D("pmiscs_Ef_vs_theta_sx3"+tag,100,0,180,800,0,20,theta_s*180/M_PI,sx3Efix,pmlabel);
						plotter->Fill1D("pmiscs_VertexReconZ"+tag,800,-400,400,vertex_z,pmlabel);
						plotter->Fill2D("pmiscs_VertexReconXY"+tag,200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),pmlabel);						
						plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef"+tag,800,-400,400,800,0,20,vertex_z,sx3Efix,pmlabel);
						plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef"+tag+"_a"+std::to_string(pcevent.multi1),800,-400,400,800,0,20,vertex_z,sx3Efix,pmlabel);
					};
					
					plot_with_tag();
					if(cathode_alpha_select) 
						plot_with_tag("_cathode_alphas");
					else 
						plot_with_tag("_cathode_protons");
					
					//plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),pmlabel);

				}//end PCEvents loop
		}//end sx3Events loop
}

/*
	p(a,a) events with a in qqq and p in sx3 along with PC dE for alphas
*/
void protonAlphaHistograms(HistPlotter* plotter, std::vector<Event> QQQ_Events, std::vector<Event> SX3_Events, std::vector<Event> PC_Events) {
		//Sidetrack for a(p,p)
		std::string aplabel = "a(p,p)";
   		Kinematics apkin_p(1.008664916,4.002603254,1.008664916,4.002603254,7.0);//m3 is proton 
		Kinematics apkin_a(1.007825,4.002603254,4.002603254,1.007825,7.0); //m3 is alpha

		for(auto qqqevent: QQQ_Events) {
			for(auto sx3event:SX3_Events) {
				plotter->Fill1D("ap_qqq_sx3_dt",800,-2000,2000,qqqevent.Time1-sx3event.Time1,aplabel);	
				if(TMath::Abs(qqqevent.Time1-sx3event.Time1)>300) continue;
				//sx3event.pos.SetZ(sx3event.pos.Z()+5.0);
				plotter->Fill1D("ap_qqq_sx3_dt_timecut",800,-2000,2000,qqqevent.Time1-sx3event.Time1,aplabel);	
				plotter->Fill1D("ap_qqq_sx3_dphi",180,-360,360,qqqevent.pos.Phi()*180/M_PI - sx3event.pos.Phi()*180/M_PI,aplabel);
				plotter->Fill2D("ap_qqq_sx3_dphi_vs_qqqphi",180,-360,360,180,-360,360,qqqevent.pos.Phi()*180/M_PI - sx3event.pos.Phi()*180/M_PI,qqqevent.pos.Phi()*180/M_PI,aplabel);
				plotter->Fill2D("ap_qqq_sx3_matrix",400,0,10,400,0,10,qqqevent.Energy1,sx3event.Energy1,aplabel);

				for(auto pcevent: PC_Events) {
					double pcz_fix = pcfix_func.Eval(pcevent.pos.Z())-5.0;
					TVector3 x2f(pcevent.pos.X(),pcevent.pos.Y(),pcz_fix);
					TVector3 x1(qqqevent.pos);
					TVector3 v = x2f-x1;
					double t_minimum = -1.0*(x1.X()*v.X()+x1.Y()*v.Y())/(v.X()*v.X()+v.Y()*v.Y());
					TVector3 r_rhoMin_fix = x1 + t_minimum*v;
					double vertex_z = r_rhoMin_fix.Z();
					//double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
					double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
 					double sinTheta_customV = TMath::Sin(theta_q);
					//double theta_s = (sx3event.pos - TVector3(0,0,vertex_z)).Theta();
					double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
					double sinTheta_s = TMath::Sin(theta_s);
					//if(vertex_z<0 || vertex_z>100) continue;

					//double sinTheta = TMath::Sin((qqqevent.pos - pcevent.pos).Theta());
					//plotter->Fill2D("sinTheta2_vs_sinTheta",80,-2,2,80,-2,2,sinTheta,sinTheta_customV,aplabel);
					plotter->Fill2D("ap_dE_E_Anodesx3B",400,0,10,800,0,40000,sx3event.Energy1,pcevent.Energy1,aplabel);
					plotter->Fill2D("ap_dE_E_Cathodesx3B",400,0,10,800,0,10000,sx3event.Energy1,pcevent.Energy2,aplabel);
					plotter->Fill2D("ap_dE_E_AnodeQQQ",400,0,10,800,0,40000,qqqevent.Energy1,pcevent.Energy1,aplabel);
					plotter->Fill2D("ap_dE_E_CathodeQQQ",400,0,10,800,0,10000,qqqevent.Energy1,pcevent.Energy2,aplabel);
					plotter->Fill2D("ap_dE3_E_AnodeQQQ",400,0,10,400,0,40000,qqqevent.Energy1,pcevent.Energy1*sinTheta_customV,aplabel);
					plotter->Fill2D("ap_dE3_E_CathodeQQQ",400,0,10,400,0,10000,qqqevent.Energy1,pcevent.Energy2*sinTheta_customV,aplabel);
					plotter->Fill2D("ap_dPhi_QQQ_PC",180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,qqqevent.pos.Phi()*180/M_PI,aplabel);
					plotter->Fill2D("ap_dPhi_SX3_PC",180,-360,360,180,-360,360,pcevent.pos.Phi()*180/M_PI,sx3event.pos.Phi()*180/M_PI,aplabel);
					plotter->Fill1D("ap_dt_Anode_QQQ",600,-2000,2000,pcevent.Time1-qqqevent.Time1,aplabel);
					plotter->Fill1D("ap_dt_Cathode_QQQ",600,-2000,2000,pcevent.Time2-qqqevent.Time1,aplabel);
					plotter->Fill1D("ap_dt_Anode_SX3",600,-2000,2000,pcevent.Time1-sx3event.Time1,aplabel);
					plotter->Fill1D("ap_dt_Cathode_SX3",600,-2000,2000,pcevent.Time2-sx3event.Time1,aplabel);
					plotter->Fill1D("ap_pczfix",600,-300,300,pcz_fix,aplabel);
					plotter->Fill1D("ap_pcz",600,-300,300,pcevent.pos.Z(),aplabel);

					//double path_length_q = (qqqevent.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
					//double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
					double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
					double path_length_s = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;
						//We know that alphas predominantly are detected in QQQs, and protons in SX3s, and that protons don't leave much of a trace in dE layer.
						//Using the estimated path lengths, we correct alpha eloss in qqq, and protons in sx3. The result should (hopefully be) vertex independent.
					double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1)-path_length_q);
					double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1)-path_length_s);
					//plotter->Fill2D("qqqEf_sx3E_matrix_all",400,0,10,400,0,10,qqqEfix,sx3event.Energy1,aplabel);
					plotter->Fill2D("ap_dE3_Ef_AnodeQQQ",400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV,aplabel);
					plotter->Fill2D("ap_dE3_Ef_CathodeQQQ",400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,aplabel);
					
					plotter->Fill2D("ap_qqqEf_sx3Ef_matrix",400,0,10,400,0,10,qqqEfix,sx3Efix,aplabel);
					plotter->Fill2D("ap_Ef_vs_theta_qqq",100,0,180,400,0,10,theta_q*180/M_PI,qqqEfix,aplabel);
					plotter->Fill2D("ap_Ef_vs_theta_sx3",100,0,180,400,0,10,theta_s*180/M_PI,sx3Efix,aplabel);
					plotter->Fill2D("ap_theta_vs_theta_qqq_sx3",100,0,180,100,0,180,theta_q*180/M_PI,theta_s*180/M_PI,aplabel);
					plotter->Fill1D("ap_VertexReconZ",400,-200,200,vertex_z,aplabel);
					plotter->Fill2D("ap_VertexReconXY",200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),aplabel);
					plotter->Fill1D("ap_Ex_from_protons",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),aplabel);
					plotter->Fill1D("ap_Ex_from_alpha",200,-10,10,apkin_a.getExc(qqqEfix,theta_q*180/M_PI),aplabel);
					if(pcevent.multi1==1 && pcevent.multi2==2) { //one-anode, two-cathode events, as originally intended
						//std::cout << "Test" << std::endl;
						plotter->Fill1D("ap_VertexReconZ_a1c2",400,-200,200,vertex_z,aplabel);
						plotter->Fill2D("ap_VertexReconXY_a1c2",200,-100,100,200,-100,100,r_rhoMin_fix.X(),r_rhoMin_fix.Y(),aplabel);
						plotter->Fill2D("ap_theta_vs_theta_qqq_sx3_a1c2",100,0,180,100,0,180,theta_q*180/M_PI,theta_s*180/M_PI,aplabel);
						plotter->Fill2D("ap_Ef_vs_theta_qqq_a1c2",100,0,180,400,0,10,theta_q*180/M_PI,qqqEfix,aplabel);
						plotter->Fill1D("ap_Ex_from_protons_a1c2",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),aplabel);
						plotter->Fill1D("ap_Ex_from_alpha_a1c2",200,-10,10,apkin_a.getExc(qqqEfix,theta_q*180/M_PI),aplabel);

						//std::cout << apkin_p.getExc(sx3Efix,theta_s*180/M_PI) << " " << apkin_a.getExc(qqqEfix,theta_q*180/M_PI)<<  std::endl;
						plotter->Fill2D("ap_Ef_vs_theta_sx3_a1c2",100,0,180,400,0,10,theta_s*180/M_PI,sx3Efix,aplabel);
						//plotter->Fill2D("qqqEf_sx3E_matrix",400,0,10,400,0,10,qqqEfix,sx3event.Energy1,aplabel);
						plotter->Fill2D("ap_qqq_sx3_matrix_a1c2",400,0,10,400,0,10,qqqevent.Energy1,sx3event.Energy1,aplabel);
						plotter->Fill2D("ap_qqqEf_sx3Ef_matrix_a1c2",400,0,10,400,0,10,qqqEfix,sx3Efix,aplabel);
						//std::cout << sx3event.Energy1 	<< " " <<  path_length_s << " " << sx3Efix << std::endl;

						//plotter->Fill2D("dE3_Ef_AnodeQQQ_a1c2",400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV,aplabel);
						//plotter->Fill2D("dE3_Ef_CathodeQQQ_a1c2",400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,aplabel);
					} //end if(a1c2) loop
				}//end PC_Events for loop
			}//end SX3_Events for loop
		} //end QQQ_Events for loop, end sidetrack a(p,p)
		return;
}