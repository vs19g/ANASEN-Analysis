#ifndef SX3Geom_h
#define SX3Geom_h
#include <vector>

const double DEFAULT_NULL=-987654321.;

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

std::array<sx3_fbgains,24> sx3_xtalk_gains; //every sx3 needs to be gainmatched as a frontL-back, frontR-back pair (pad strip pair)
std::array<sx3_geometry_scalefactors,24> sx3gs;

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
    void validate();
};


void sx3::fillevent(const std::string& positionstring, const int subchannel, const float value) {
    assert(subchannel>=0 && subchannel<4);
    if(positionstring=="FRONT_L") {
        frontL[subchannel].push_back(value);
        unmatched_front_chans.insert(subchannel);
    } else if(positionstring=="FRONT_R") {
        frontR[subchannel].push_back(value);
        unmatched_front_chans.insert(subchannel);
    } else if(positionstring=="BACK") {
        back[subchannel].push_back(value);
        valid_back_chans.insert(subchannel);
    } else {
        std::cout << "Unknown string "+positionstring+" encountered in sx3::fillevent \n" << std::endl;
    }
    if(frontL[subchannel].size()!=0 && frontR[subchannel].size()!=0 ) {
        unmatched_front_chans.erase(subchannel);
        valid_front_chans.insert(subchannel); //std::set, so no duplication will happen
    }
}


//void sx3::validate(const sx3_fbgains& fbgains, const sx3_geometry_scalefactors& sx3gs) {
void sx3::validate() {
    if(valid_front_chans.size()!=0 && valid_back_chans.size()!=0) {
        valid=true;
        float maxFE=0;
        float maxBE=0;
        //float zpos=0;
        int bchan=-1;
        int fchan=-1;
/*            	for(auto cc: valid_front_chans) {
            		std::cout << "fc" << cc << std::endl;// " " << frontL[cc].at(0) << " " << frontR[cc].at(0) << std::endl;
            	}
            	for(auto cc: valid_back_chans) {
            		std::cout << "bc" << cc << std::endl; //" " << back[cc].at(0) << std::endl;
            	}
*/
        for(auto chan: valid_front_chans) {
            if(frontL[chan].size()>1) {
                printf("\nmultihit sx3 at Lsubchan:%d, ts:%1.13g\n",chan,ts);
                for(const auto& e: frontL[chan]) printf("e: %f\t",e);
                std::sort(frontL[chan].begin(), frontL[chan].end(), std::greater<float>());
                flags += (-1000);
            }
            if(frontR[chan].size()>1) {
                printf("\nmultihit sx3 at Rsubchan:%d, ts:%1.13g\n",chan,ts);
                for(const auto& e: frontR[chan]) printf("e: %f\t",e);
                std::sort(frontR[chan].begin(), frontR[chan].end(), std::greater<float>());
                flags += (-2000);
            }
            //assign position using max L+R value
            /*printf("chan:%d sizeL: %d sizeR: %d\n",chan, frontL[chan].size(), frontR[chan].size()); fflush(stdout);
			printf("foo\n");
			std::cout << "\nL:" << std::endl;
			for(auto thing: frontL[chan]) std::cout << thing << " " << std::flush;
			std::cout << "\nR:" << std::endl;
			for(auto thing: frontR[chan]) std::cout << thing << " " << std::flush;*/
            if(frontL[chan].at(0) + frontR[chan].at(0)> maxFE) {
                maxFE = frontL[chan].at(0) + frontR[chan].at(0);
                //zpos = (frontL[chan].at(0)-frontR[chan].at(0))/maxFE;
                fchan = chan;
            }
        }
        for(auto chan: valid_back_chans) {
            if(back[chan].size()>1) {
                printf("\nmultihit sx3 at Bsubchan:%d, ts:%1.13g\n",chan,ts);
                for(const auto& e: back[chan]) printf("e: %f\t",e);
                std::sort(back[chan].begin(), back[chan].end(), std::greater<float>());
                flags += (-3000);
            }
            if(back[chan].size() ==0 ) {
            	printf("foo\n");
            	//continue;
            }
            if(back[chan].at(0) > maxBE) {
                maxBE = back[chan].at(0);
                bchan = chan;
            }
        }
        /*
        	Cross-talk corrections are important when evaluating 'energy' signals from strips/pads.
        	They can cause unexpected behavior when used universally for all EL, ER cases, so we split scenarios in two.
        		- Positions along each strip (frontX) *are not* corrected for crosstalk.
        		- Total F and B energies (frontE, backE) *are*.
        	Sudarsan B, 31 Oct 2024
        */
        float Eleft = frontL[fchan].at(0);
        float Eright = frontR[fchan].at(0);
        frontEL = Eleft;
        frontER = Eright;
        frontX = (Eleft-Eright)/(Eleft+Eright);
        //frontXmm = (frontX+sx3gs.add[fchan])*sx3gs.stretch[fchan]; //convert to mm

        //frontE = Eleft*fbgains.stripLgains[bchan][fchan] + fbgains.stripLoffsets[bchan][fchan]
         //        + Eright*fbgains.stripRgains[bchan][fchan]  + fbgains.stripRoffsets[bchan][fchan];
        //backE = back[bchan].at(0)*fbgains.padgains[bchan][fchan]+fbgains.padoffsets[bchan][fchan];
        frontE = Eleft+Eright;
        backE = maxBE;
        stripF=fchan;
        stripB=bchan;

        flags = 0;
    } else if(valid_front_chans.size()!=0 && valid_back_chans.size()==0) {
        flags = -10;
    } else if(valid_front_chans.size()==0 && valid_back_chans.size()!=0) {
        flags = -20;
    }
}

typedef sx3 sx3det;
#endif
