#ifndef UTILS_ORR_H
#define UTILS_ORR_H
#include "datatypes.h"
#include "HistPlotter.h"
#include "Geometry_orr.h" //contains orruba geometry constants
#include <cassert>
#include <stdio.h>
#include <cassert>
#include <cstdint>
#include <fcntl.h>
#include <unistd.h>
#include <vector>
//#include <cstring>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <set>
#include <sys/stat.h>
#include "TMath.h"
#include "counters.h"

static named_counter oc("named orruba counters");

inline float get_filesize(std::string filename) {
    struct stat st;
    stat(filename.c_str(), &st);
    return st.st_size;
}

class orruba_params {
    /*
    *   
    */
public:
    int chnum=DEFAULT_NULL; //!< global channel number
    std::string type = "-";
    int id=DEFAULT_NULL;
    int layer=DEFAULT_NULL;
    int frontback=DEFAULT_NULL;
    int updown=DEFAULT_NULL;
    int subid=DEFAULT_NULL;
    int leftright=DEFAULT_NULL;

    float ped=DEFAULT_NULL;
    float offset=DEFAULT_NULL;
    float gain=DEFAULT_NULL;
    float gain2=DEFAULT_NULL; //for use with sx3's
};

class sx3_geometry_scalefactors {
public:
    //If sx3 has L, R being the left and right extremities, we choose add, stretch here such that
    // x_mm = (x_raw+add)*stretch; so add=abs(L), stretch=75/(abs(L)+R)
    float add[4];
    float stretch[4];
};

class qqq5_finegains {
public:
    std::array<std::pair<float,float>,32> front;
    //front.at(30).first = slope at clkpos 0, ring 30 for E front layer
    //front.at(30).second = intercept for the same as above
    std::array<std::pair<float,float>,4> back;
};
class sx3_fbgains {
public:
    //Order of indices is [pad][strip]
    float padoffsets[4][4];
    float padgains[4][4];

    float stripLoffsets[4][4];
    float stripLgains[4][4];

    float stripRoffsets[4][4];
    float stripRgains[4][4];
};

//Metadata ORRUBA needs to know about itself, to be configured at the start
extern std::array<orruba_params,MAX_ORRUBA_CHANS> o_params;
extern std::array<sx3_fbgains,24> sx3_xtalk_gains; //every sx3 needs to be gainmatched as a frontL-back, frontR-back pair (pad strip pair)
extern std::array<sx3_geometry_scalefactors,24> sx3gs;
extern std::array<qqq5_finegains,4> qqq5_fg_dE, qqq5_fg_E;

class type19Raw {
public:
    long int timestamp;
    std::vector<unsigned short> ch;
    std::vector<unsigned short> val;

    type19Raw() : timestamp(0), ch(50,0), val(50,0) {} //Reserve 50 for size of these vectors, initial value of zero
    void print() const {
        /*
            print()
            Prints type19Raw's contents to stdout for monitoring
        */
        std::cout << "------" << std::endl;
        for(unsigned int ii=0; ii<ch.size(); ii++) {
            std::cout << ch.at(ii) << " " << val.at(ii) << std::endl;
        }
    }
};

class sx3 {
public:
    //TODO: Convert to std::array
    //Holds all information in an event, including ped subtraction+scaling. back[2].at(0) will have the largest energy seen in ch2, if any
    std::vector<float> back[4];
    std::vector<float> frontL[4];
    std::vector<float> frontR[4];

    double ts = DEFAULT_NULL;
    //Easy lookup of final calibrated event. Only filled for valid cases, assumed for now to be 1L, 1R, 1B
    float frontX=DEFAULT_NULL;
    float frontXmm=DEFAULT_NULL;
    float frontE=DEFAULT_NULL;
    float backE=DEFAULT_NULL;
    int stripF=DEFAULT_NULL;
    int stripB=DEFAULT_NULL;
    float frontEL=DEFAULT_NULL;
    float frontER=DEFAULT_NULL;

    float phi=DEFAULT_NULL; //

    std::set<int> valid_front_chans;
    std::set<int> valid_back_chans;
    std::set<int> unmatched_front_chans; //every front channel is unmatched and invalid at first. when it gets matched, it gets removed and sent to valid

    bool foundevent=false;
    bool valid=false;//valid will be set to false in all cases where we have ambiguity
    int flags=-1;//flags settable to different types of values to indicate different invalid situations

    void fillevent(const std::string& position, const int subchannel, const float value); //make 'const' what functions don't need to change, helps with performance
    void validate(const sx3_fbgains&, const sx3_geometry_scalefactors&);
};

class qqq5 {
public:
    //Holds all information in an event, including ped subtraction+scaling. front[2].at(0) will have the largest energy seen in ch2, if any
    //TODO: Convert to std::array
    std::vector<float> back[4];
    std::vector<float> front[32];

    double ts = DEFAULT_NULL;
    float selftheta=DEFAULT_NULL,selfrho=DEFAULT_NULL;
    //Easy lookup of the final calibrated event. Only filled for valid cases.
    double frontE=DEFAULT_NULL;
    double backE=DEFAULT_NULL;
    int frontch;
    int backch;

    std::pair<int,int> adj_front_strips = {-1,-1};
    std::pair<int,int> adj_back_strips = {-1,-1};
    std::set<int> valid_front_chans; //list of channels that fire. can inspect size() of these to see if there are front/back events
    std::set<int> valid_back_chans; // we use std::set since it makes for very readable code

    bool foundevent=false;
    bool valid=false; //valid will be set to false in all cases where we have ambiguity
    int flags=-1; //flags settable to different types of values to indicate different invalid situations

    void fillevent(const std::string& position, const int subchannel, const float value); //make 'const' what functions don't need to change, helps with performance
    void validate();
};

struct orrubaevent {
    //Every clean, valid charged-particle event will have these four parts
    float dE=DEFAULT_NULL;  //!< true energy-loss in the dE layer. Found by gainmatching ADC readout to alpha data
    float E=DEFAULT_NULL;   //!< energy deposited in the  E layer. When summed with dE, gives true energy in keV deposited by the particle in ORRUBA
    float dE_PID = DEFAULT_NULL; //!< dE scaled for dE-layer's thickness, reducing the spread due to angular straggling by explicitly accounting for it. This will give a sharper pid plot that can be gated on better
    float dE_linPID = DEFAULT_NULL; //!< dE_PID, but linearized using the prescription described in, say, PhysRevC.90.034601 (2014). dE_linPID = ((dE+E)^a-E^a)^(1/a) where a ~ 1.68 is chosen empirically
    float Theta=DEFAULT_NULL; //!< Laboratory polar angle of event in radians, deprecated
    float Phi=DEFAULT_NULL; //!< Lab azimuthal angle of event in radians, deprecated

    //Helpful indices to make dE-E plots
    std::string type; //!< "endcap" vs "barrel"
    //!< Identify the position of the detector in the barrel, usually in accordance with the channel map: say we might learn detector is at "Quad 4" or "clk_pos 10", together with 'type'. Useful with HistPlotter class
    int position=DEFAULT_NULL;
	int subchdE_1=DEFAULT_NULL, subchdE_2=DEFAULT_NULL; //!< Identify the subchannels corresponding to the two sides of the dE detector. To avoid confusion, 1=strip(sx3), ring(qqq) and 2=pad(sx3), wedge(qqq)
	int subchE_1=DEFAULT_NULL, subchE_2=DEFAULT_NULL; //!< Identify the subchannels corresponding to the two sides of the E detector. Same convention as above

    float x=DEFAULT_NULL,y=DEFAULT_NULL,z=DEFAULT_NULL; //!< Laboratory x,y,z coordinates of the event from the E layer
    float r0=DEFAULT_NULL,theta0=DEFAULT_NULL,phi0=DEFAULT_NULL; //!< vector elements from hit to origin for E layer
    float r1=DEFAULT_NULL,theta1=DEFAULT_NULL,phi1=DEFAULT_NULL; //!< vector elements from hit to origin for dE layer
};

/*TODO:
    * There will be some use for a class such that it stores 
        PhysicsEvent = <Ex, Brho, THeta4, ...orrubaevent >
    * Once the 'orrubaevent' structs are made, it should be callable. like 
    std::vector<PhysicsEvent> getPhysicsFromVertices(Kinematics dpkin, std::vector<orrubaevent> orvec);?
    * should the physics just go sit in orrubaevent? maybe orruba2024 class can have a kinematics object in it?
*/

class orruba2024 {
private:
//Class expected to be changed for each version of the analysis code
public:
    bool found_trk, found_trkpresc, found_tdcq, found_s800e1, found_s800trg,
         found_rf, found_gt, found_si, found_siup;

	bool found_de, found_e, found_qqq, found_sx3;//orruba
    long long timestamp=DEFAULT_NULL;
    std::vector<type19Raw> o_rawvec;

    std::array<qqq5,4> uendcapE;
    std::array<qqq5,4> uendcapdE;
    std::array<sx3,12> ubarrelE;
    std::array<sx3,12> ubarreldE;

    std::array<sx3,2> dbarrelE;
    //Results after post-processing, including possible multiplicities
    std::vector<orrubaevent> events;

    double target_z_offset;

    //void initialize(const std::string & filename1, const std::string& filename2, const std::string& filename3);
	float tdc_trk = DEFAULT_NULL;
	float tdc_trksi = DEFAULT_NULL; //trackerpr+si stops are combined into one channel
	float tdc_trkp = DEFAULT_NULL;
	float tdc_q = DEFAULT_NULL;
	float tdc_s800e1 = DEFAULT_NULL;
	float tdc_s800trg = DEFAULT_NULL;
	float tdc_rf = DEFAULT_NULL;
	float tdc_rf_unwrap = DEFAULT_NULL;
	float tdc_gt = DEFAULT_NULL;
	float tdc_si = DEFAULT_NULL;
	float tdc_siup = DEFAULT_NULL;
    orruba2024(const std::vector<type19Raw>& orvec);
    void postprocess(); //generates 'events', and performs validations on freshly entered data. performs fine gainmatching for front-back-pairs
    void print() const {
    	std::cout << Form("TDCs\n trk:%1.4f\ntrkp:%1.4f\nq:%1.4f\ns800e1:%1.4f\ns800trg:%1.4f\nrf:%1.4f\nrfuw:%1.4f\ngt:%1.4f\nsi:%1.4f\nsiup:%1.4f\n-----\nevents_size:%lu\ntimestamp:%lld\nfound_de:%d, found_e:%d\n-------\n-------\n",
    	tdc_trk, tdc_trkp, tdc_q, tdc_s800e1, tdc_s800trg, tdc_rf, tdc_rf_unwrap, tdc_gt, tdc_si, tdc_siup, events.size(), timestamp, found_de, found_e) << std::endl;
    }
};

class trackingdet {
public:
    double timestamp=DEFAULT_NULL;
    //TODO: Convert all to std::array
    std::vector<float> xwires[MAXNWIRES_TRACK];
    std::vector<float> ywires[MAXNWIRES_TRACK];
    std::vector<float> xwiresf[MAXNWIRES_TRACK];
    std::vector<float> ywiresf[MAXNWIRES_TRACK];
    std::vector<float> xwiresf_nn[MAXNWIRES_TRACK];
    std::vector<float> ywiresf_nn[MAXNWIRES_TRACK];

    std::vector<float> xtimes[MAXNWIRES_TRACK];
    std::vector<float> ytimes[MAXNWIRES_TRACK];
    std::vector<float> xtimesf[MAXNWIRES_TRACK];
    std::vector<float> ytimesf[MAXNWIRES_TRACK];
    std::vector<float> cathode;

    int multx; //how many x-wires fired?
    int multy; //how many y-wires fired?
    int multxf; //how many filtered x-wires fired?
    int multyf; //how many filtered y-wires fired?

    int multxt; //how many x-tdcwires fired?
    int multyt; //how many y-tdcwires fired?

    int multxtf; //how many x-tdcwires fired in window?
    int multytf; //how many y-tdcwires fired in window?

    float tot_cathode=0;
    float tot_x=0;
    float tot_y=0;
    float tot_anode=0;

    //list of x and y wires fired above energy threshold, within timing gate window
    std::set<int> list_ywires;
    std::set<int> list_xwires;

    std::set<int> list_ytwires;
    std::set<int> list_xtwires;

    //position of the vertex estimated by up to 2 neighbouring wires firing together within window  
    double xpos=DEFAULT_NULL;
    double ypos=DEFAULT_NULL;
    bool clean_event = false;
    int maxnx=-124, maxny=-124;
    int nnx = 0; //nearest neighbour x wires set this to +/- 1 if present w.r.t. wire maxnx
    int nny = 0; //nearest neighbour y wires set this to +/- 1 if present w.r.t  wire maxny
    bool clean_event_no_timing = false;
    bool clean_single_xy_event = false;
    bool clean_single_xy_event_no_timing = false;
    void Reset() {
        /***
        Resets all data members.
        **/
        cathode.clear();
        for(int i=0; i<MAXNWIRES_TRACK; i++) {
            xwires[i].clear();
            ywires[i].clear();
            xwiresf[i].clear();
            ywiresf[i].clear();
            xwiresf_nn[i].clear();
            ywiresf_nn[i].clear();
            xtimes[i].clear();
            ytimes[i].clear();
            xtimesf[i].clear();
            ytimesf[i].clear();
        }
        list_xwires.clear();
        list_ywires.clear();
        list_xtwires.clear();
        list_ytwires.clear();

        nnx = 0;
        nny = 0;
        clean_event_no_timing = false;
        clean_event = false; //x, y both fire above thresh, xt, yt present within broad coinc window
        xpos=DEFAULT_NULL;
        ypos=DEFAULT_NULL;
        maxnx=-124;
        maxny=-124;
        multx=0; //how many x-wires fired?
        multy=0; //how many y-wires fired?
        multxf=0; //how many filt x-wires fired?
        multyf=0; //how many filt y-wires fired?
        multxt=0; //how many tdc x-wires fired?
        multyt=0; //how many tdc y-wires fired?
        multxtf=0; //how many tdc x-wires fired?
        multytf=0; //how many tdc y-wires fired?
    }

    trackingdet() {
        Reset();
    };
    trackingdet(const std::vector<type19Raw>& orvec);
};

const float alpha = 0.0;
int matchchantype(unsigned short chan, const std::array<orruba_params,MAX_ORRUBA_CHANS>& index, const std::string& label);
void initialize_orruba(const std::string & filename1, const std::string& filename2, const std::string& filename3, const std::string& filename4);//,
int parse_orruba_data(const unsigned short* buffer, int32_t length, type19Raw& oraw_event);
#endif
