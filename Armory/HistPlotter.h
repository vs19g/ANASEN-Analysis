#ifndef HISTPLOTTER_H
#define HISTPLOTTER_H
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>
#include <TFile.h>
#include <TMemFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCutG.h>
#include <signal.h>
#include <cstdlib>
#include <utility>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <TGraphErrors.h>

class HistPlotter {
private:
    long long barrier_count, barrier_limit; //meant to keep track of how often to call FillN() on histograms
    enum {TFILE, TMEMFILE} filetype;
    std::unordered_map<std::string,TObject*> oMap; //!< Maps std::string to all TH1, TH2 objects in the class
    std::unordered_map<std::string,TObject*> cutsMap; //!< Maps std::string to TCutG objects held by the class
    std::set<std::string> folderList; //!< List of all folder names used to nest objects
    std::unordered_map<TObject*,std::string> foldersForObjects; //!< Map that returns the folder corresponding to the object whose pointer is specified
    TFile *ofile=nullptr; //!< TFile pointer for the output file
    TMemFile *omfile=nullptr; //!< TFile pointer for the output memfile

    //Caches to permit FillN() calls
    std::unordered_map<std::string, std::vector<double>> onedimcache;
    std::unordered_map<std::string, std::pair<std::vector<double>, std::vector<double>>> twodimcache;
    inline void FillN_All_Histograms();
public:
    HistPlotter(std::string outfile, std::string type);
    inline void FlushToDisk(int integral); //!< Writes all objects to file before closing, nesting objects in folders as is found necessary
    inline void PrintObjects(); //!< Dump objects to std::cout for inspection
    inline void ReadCuts(std::string); 
    inline TCutG* FindCut(std::string cut) {
        return static_cast<TCutG*>(cutsMap.at(cut));
    }
    inline void set_barrier_limit(long long limit) { barrier_limit = limit; }
    inline void barrier_increment() {
        barrier_count++;
        if(barrier_count == barrier_limit) {
            FillN_All_Histograms();
            barrier_count=0;
        }
    }
    /*! \fn void FindCut()
        \brief        
            - Searches for a cut by name 'cut' in the internal list of cuts 'cutsMap'. Ugly fails (via unresolved at()) if such a cut isn't found.            
        \param filename - name of the plainxtext file containing the cut file locations and identifiers
        \return Pointer to the TCutG object that matches the name. Very useful to use this as plotter.FindCut("protonbarrelpid")->IsInside(deltaE, E) for instance.
    */

    inline void SetNewTitle(std::string name, std::string title) {
        auto result = oMap.find(name); //result is an iterator
        if(result==oMap.end()) return; //no warnings, could be changed in future
        else
        	static_cast<TNamed*>(oMap.at(name))->SetTitle(title.c_str()); // set new title
	}

    //Smart functions that create a new histogram if it doesn't exist.
    inline void FillGraph(const std::string &name, float valuex, float valuey, float errx=0, float erry=0);
    inline void Fill1D(const std::string& name,int nbinsx, float xlow, float xhigh, float value);
    inline void Fill2D(const std::string& name,int nbinsx, float xlow, float xhigh
                ,int nbinsy, float ylow, float yhigh, float valuex, float valuey);
    inline void Fill1D(const std::string& name,int nbinsx, float xlow, float xhigh, float value, const std::string& folder);
    inline void Fill2D(const std::string& name,int nbinsx, float xlow, float xhigh
                ,int nbinsy, float ylow, float yhigh, float valuex, float valuey, const std::string& folder);
    //TObject* findObject(std::string key);
};

HistPlotter::HistPlotter(std::string outfile, std::string type="") {
    /*!
        \brief Constructor. Opens a TFile instance with the specified filename
        \param outfile : std::string that holds the desired output ROOT filename
        \return None
    */
    if(type=="" || type == "TFILE") {
        ofile = new TFile(outfile.c_str(),"recreate");
        filetype = TFILE;
    } else if(type =="TMEMFILE") {
        omfile = new TMemFile(outfile.c_str(),"recreate");
        filetype=TMEMFILE;
    } else {
        std::cout << "Unknown type "<< type << " specified for HistPlotter (use \"TFILE\" or \"TMEMFILE\"), using default \"TFILE\"  " << std::endl;
        ofile = new TFile(outfile.c_str(),"recreate");
        filetype = TFILE;
    }
    barrier_count=0;
    barrier_limit=1000;
}

void HistPlotter::FillN_All_Histograms() {
    for(auto it=oMap.begin(); it!=oMap.end(); it++ ) {
        //it->first is std::string 'name', it->second is the TObject
        if(it->second->InheritsFrom("TH1F")) {
            //FillN(size, array-of-doubles, array-of-weights); //we set array-of-weights to (1,1,1,.. (size)
            static_cast<TH1F*>(it->second)->FillN(onedimcache[it->first].size(), //size
                                                  onedimcache[it->first].data(), //array
                                                  std::vector<double>(onedimcache[it->first].size(),1.0).data()); //weight of ones
            onedimcache[it->first].clear();
        }  else if(it->second->InheritsFrom("TH2F")) {
            //FillN(size, array-of-doubles, array-of-weights); //we set array-of-weights to (1,1,1,.. (size))
            static_cast<TH2F*>(it->second)->FillN(twodimcache[it->first].first.size(), //size
                                                  twodimcache[it->first].first.data(), //x array
                                                  twodimcache[it->first].second.data(), //y array
                                                  std::vector<double>(twodimcache[it->first].first.size(),1.0).data()); //weight of ones
           twodimcache[it->first].first.clear();
           twodimcache[it->first].second.clear();
        }
    }
    std::cout << "." << std::endl;
}

void HistPlotter::FlushToDisk(int min_integral=0) {
    /*! \fn void FlushToDisk()
        \brief         Function that can be used at any point to exit smoothly by saving all ROOT objects in memory
        to the output file before closing it. Obeys the binding of histograms to separate folders, if so specified.
        \return No return -- void
    */
    if(filetype==TMEMFILE && omfile) {
        std::cout << "Not flushing a TMemfile .. exiting .." << std::endl;
        delete omfile;
        return;
    }
    if(ofile->IsZombie() || !ofile) {
        std::cerr << "Output file is zombie, finishing up without writing to disk!" << std::endl;
        return;
    }
    FillN_All_Histograms();
    for(auto it=oMap.begin(); it!=oMap.end(); it++ ) {
        //omap maps: name(first) to object address(second).
        // foldersForObjects maps: object address(first) to foldername(second)
        auto result = foldersForObjects.find(it->second); //returns <TObject* histogram,std::string foldername> pair if found
        if(result!=foldersForObjects.end()) { //we try to create folder if needed and cd to it
            ofile->mkdir(result->second.c_str(),"",kTRUE); // args: name, title, returnExistingDirectory
            ofile->cd(result->second.c_str());
        } else {
            ofile->cd(); //toplevel for all default histograms. Default setting
        }
        if(((TH1F*)it->second)->Integral()>min_integral)
	        it->second->Write();
    }

    //Create a directory for all cuts, and save all cuts in them
    ofile->mkdir("gCUTS","",kTRUE);
    ofile->cd("gCUTS");
    for(auto it=cutsMap.begin(); it!=cutsMap.end(); it++) {
        (static_cast<TNamed*>(it->second))->SetName(it->first.c_str());
        it->second->Write();
    }
    ofile->Close();
    std::cout << "Wrote " << oMap.size() << " histograms to TFile " << std::string(ofile->GetName()) << std::endl;
}

void HistPlotter::FillGraph(const std::string& name, float valuex, float valuey, float errx, float erry) {
    /*! \fn void FillGraph()
        \brief
        - Creates a TGraphError in memory with name 'name' if it doesn't exist, and fills it with valuex, valuey
        - Writes present state to disk and fails with return value -1 if the name clashes with another object that's not of type TGraph*

        \param name name of the TGraph
        \param valuex The xvalue
        \param valuey The yvalue
        \param errx The x error
        \param erry The y error
        \return No return void
    */
    auto result = oMap.find(name); 
    if(result==oMap.end()) {
        TGraphErrors *tempG = new TGraphErrors();
        tempG->SetName(name.c_str());
        oMap.insert(std::make_pair(name,static_cast<TObject*>(tempG)));
    }
    if(!oMap.at(name)->InheritsFrom("TGraphErrors")) {
        std::cerr << "Object " << name << " refers to something other than a TGraph*, not filling it hence!" << std::endl;
        std::cerr << "Abort.." << std::endl;
        FlushToDisk();
        exit(-1);
    }
    // static_cast<TGraphErrors*>(oMap.at(name))->AddPointError(valuex,valuey,errx,erry);
}

void HistPlotter::Fill1D(const std::string& name, int nbinsx, float xlow, float xhigh, float value) {
    /*! \fn void Fill1D()
        \brief  
        - Creates a TH1F in memory with name 'name' if it doesn't exist, and fills it with valuex, valuey
        - Writes present state to disk and fails with return value -1 if the name clashes with another object that's not of type TH1*

        \param name name of the TH1F histogram
        \param nbinsx Number of bins in the histogram
        \param xlow Lower limit on x-axis
        \param xhigh Upper limit on x-axis
        \param value The bin corresponding to value in (nbinsx, xlow, xhigh) is incremented by 1
        \return No return void
    */
    auto result = oMap.find(name); //result is an iterator
    if(result==oMap.end()) {
        TH1F* temp1D = new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xhigh);
        oMap.insert(std::make_pair(name,static_cast<TObject*>(temp1D)));
        onedimcache.insert(std::make_pair(name, std::vector<double>()));
        onedimcache[name].reserve(16384);
    } else if(foldersForObjects.find(oMap.at(name))!=foldersForObjects.end()) { //shouldn't have a folder associated with it
            std::cerr << "Object " << name << " already registered at " << foldersForObjects[oMap[name]] <<  ", choose a different name for the histogram to be stored in toplevel .." << std::endl;
    }

    //Check if the string 'name' maps to a 1D hist. If there's any other object by this name raise issue
    if(!oMap.at(name)->InheritsFrom("TH1F")) {
        std::cerr << "Object " << name << " refers to something other than a TH1*, not filling it hence!" << std::endl;
        std::cerr << "Abort.." << std::endl;
        FlushToDisk();
        exit(-1);
    }
    onedimcache[name].emplace_back(value);
    //static_cast<TH1F*>(oMap.at(name))->Fill(value);
}

void HistPlotter::Fill1D(const std::string& name, int nbinsx, float xlow, float xhigh, float value, const std::string& foldername) {
    /*! \fn void Fill1D()
        \brief        
      -  Creates a TH1F in memory with name 'name' if it doesn't exist, and fills it with valuex, valuey
      -  Writes present state to disk and fails with return value -1 if the name clashes with another object that's not of type TH1*
      -  Remembers the foldername this particular histogram maps to, if provided. If not, defaults to toplevel.

        \param name name of the TH1F histogram
        \param nbinsx Number of bins in the histogram
        \param xlow Lower limit on x-axis
        \param xhigh Upper limit on x-axis
        \param value The bin corresponding to value in (nbinsx, xlow, xhigh) is incremented by 1
        \param foldername Name of the folder to put this histogram into. Defaults to toplevel if left empty
        \return No return -- void
    */

    auto result = oMap.find(name); //result is an iterator
    if(result==oMap.end()) {
        TH1F* temp1D = new TH1F(name.c_str(), name.c_str(), nbinsx, xlow, xhigh);
        oMap.insert(std::make_pair(name,static_cast<TObject*>(temp1D)));
        onedimcache.insert(std::make_pair(name, std::vector<double>()));
        onedimcache[name].reserve(16384);
        if(foldername!="") {
            if(folderList.find(foldername)==folderList.end()) {
                folderList.insert(foldername);
            }
            foldersForObjects.insert(std::make_pair(static_cast<TObject*>(temp1D),foldername));
        }
    } else {
        //object is present in map, but we enforce unique names
        //it must already have a folder attached to it
        if(foldersForObjects.find(oMap.at(name))==foldersForObjects.end()) {
            std::cerr << "Object " << name << " already registered at toplevel, choose a different name for the histogram to be stored in " << foldername << " folder .." << std::endl;
        } else if(foldersForObjects[oMap[name]]!=foldername) {
            std::cerr << "Object " << name << " already registered at " << foldersForObjects[oMap[name]] <<  ", choose a different name for the histogram to be stored in " << foldername << " folder .." << std::endl;
        }
    }
    //Check if the string 'name' maps to a 1D hist. If there's any other object by this name raise issue
    if(!oMap.at(name)->InheritsFrom("TH1F")) {
        std::cerr << "Object " << name << " refers to something other than a TH1*, not filling it hence!" << std::endl;
        std::cerr << "Abort.." << std::endl;
        FlushToDisk();
        exit(-1);
    }
    onedimcache[name].emplace_back(value);
    //static_cast<TH1F*>(oMap.at(name))->Fill(value);
}

void HistPlotter::Fill2D(const std::string& name, int nbinsx, float xlow, float xhigh, int nbinsy, float ylow, float yhigh, float valuex, float valuey) {
    /*! \fn void Fill2D()
        \brief        
      -  Creates a TH2F in memory with name 'name' if it doesn't exist, and fills it with valuex, valuey
      -  Writes present state to disk and fails with return value -1 if the name clashes with another object that's not of type TH2*
        \param name name of the TH1F histogram
        \param nbinsx Number of xbins in the histogram
        \param xlow Lower limit on x-axis
        \param xhigh Upper limit on x-axis
        \param nbinsy Number of ybins in the histogram
        \param ylow Lower limit on y-axis
        \param yhigh Upper limit on y-axis
        \param valuex
        \param valuey The bin corresponding to (valuex, valuey) in (nbinsx, xlow, xhigh, ybinsx, ylow, yhigh) is incremented by 1
        \return No return -- void
    */

    auto result = oMap.find(name); //result is an iterator
    if(result==oMap.end()) {
        TH2F* temp2D = new TH2F(name.c_str(), name.c_str(), nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
        oMap.insert(std::make_pair(name,static_cast<TObject*>(temp2D)));
        twodimcache.insert(std::make_pair(name, std::make_pair(std::vector<double>(),std::vector<double>())));
        twodimcache[name].first.reserve(16384);
        twodimcache[name].second.reserve(16384);
    } else if(foldersForObjects.find(oMap.at(name))!=foldersForObjects.end()) { //shouldn't have a folder associated with it
        std::cerr << "Object " << name << " already registered at " << foldersForObjects[oMap[name]] <<  ", choose a different name for the histogram to be stored in toplevel .." << std::endl;
    }

    //Check if the string 'name' maps to a 1D hist. If there's any other object by this name raise issue
    if(!oMap.at(name)->InheritsFrom("TH2F")) {
        std::cerr << "Object " << name << " refers to something other than a TH2*, not filling it hence!" << std::endl;
        std::cerr << "Abort.." << std::endl;
        FlushToDisk();
        exit(-1);
    }
    twodimcache[name].first.emplace_back(valuex);
    twodimcache[name].second.emplace_back(valuey);
    //static_cast<TH2F*>(oMap.at(name))->Fill(valuex,valuey);
}

void HistPlotter::Fill2D(const std::string& name, int nbinsx, float xlow, float xhigh, int nbinsy, float ylow, float yhigh, float valuex, float valuey, const std::string& foldername) {
    /*! \fn void Fill2D()
        \brief        
      -  Creates a TH2F in memory with name 'name' if it doesn't exist, and fills it with valuex, valuey
      -  Writes present state to disk and fails with return value -1 if the name clashes with another object that's not of type TH2*
      -  Remembers the foldername this particular histogram maps to, if provided. If not defaults to toplevel

        \param name name of the TH1F histogram
        \param nbinsx Number of xbins in the histogram
        \param xlow Lower limit on x-axis
        \param xhigh Upper limit on x-axis
        \param nbinsy Number of ybins in the histogram
        \param ylow Lower limit on y-axis
        \param yhigh Upper limit on y-axis
        \param valuex
        \param valuey The bin corresponding to (valuex, valuey) in (nbinsx, xlow, xhigh, ybinsx, ylow, yhigh) is incremented by 1
        \param foldername Name of the folder to put this histogram into. Defaults to toplevel if left empty
        \return No return -- void
    */

    auto result = oMap.find(name); //result is an iterator
    if(result==oMap.end()) {
        TH2F* temp2D = new TH2F(name.c_str(), name.c_str(), nbinsx, xlow, xhigh, nbinsy, ylow, yhigh);
        oMap.insert(std::make_pair(name,static_cast<TObject*>(temp2D)));
        twodimcache.insert(std::make_pair(name, std::make_pair(std::vector<double>(),std::vector<double>())));
        twodimcache[name].first.reserve(16384);
        twodimcache[name].second.reserve(16384);
        if(foldername!="") {
            if(folderList.find(foldername)==folderList.end()) {
                folderList.insert(foldername);
            }
            foldersForObjects.insert(std::make_pair(static_cast<TObject*>(temp2D),foldername));
        }
    } else {
        //object is present in map, but we enforce unique names
        //it must already have a folder attached to it
        if(foldersForObjects.find(oMap.at(name))==foldersForObjects.end()) {
            std::cerr << "Object " << name << " already registered at toplevel, choose a different name for the histogram to be stored in " << foldername << " folder .." << std::endl;
        } else if(foldersForObjects[oMap.at(name)]!=foldername) {
            std::cerr << "Object " << name << " already registered at " << foldersForObjects[oMap[name]] <<  ", choose a different name for the histogram to be stored in " << foldername << " folder .." << std::endl;
        }
    }

    //Check if the string 'name' maps to a 1D hist. If there's any other object by this name raise issue
    if(!oMap.at(name)->InheritsFrom("TH2F")) {
        std::cerr << "Object " << name << " refers to something other than a TH2*, not filling it hence!" << std::endl;
        std::cerr << "Abort.." << std::endl;
        FlushToDisk();
        exit(-1);
    }
    twodimcache[name].first.emplace_back(valuex);
    twodimcache[name].second.emplace_back(valuey);
    //static_cast<TH2F*>(oMap.at(name))->Fill(valuex,valuey);
}

void HistPlotter::ReadCuts(std::string filename) {
    /*! \fn void ReadCuts()
        \brief  Reads a list of cuts from a file. The file must have the format below, two columns
                   - Column#1 - path to a file that contains a single TCutG object named "CUTG", the default name in ROOT.
                   - Column#2 - The identifier name you plan to use in the code, like 'protonbarrelpid' or something, that will be searched by FindCut()
        \param filename name of the plainxtext file containing the cut file locations and identifiers
        \return No return -- void
    */

    std::ifstream infile;
    infile.open(filename);
    std::string cutfilename, cutname;
    for(std::string line; std::getline(infile, line); ) {
        if(line.size()!=0 && line[0]=='#')
            ; //don't do anything with '#' lines
        else {
            std::stringstream ss(line);
            ss>>cutfilename>>cutname;

            TFile f(cutfilename.c_str());
            if(f.IsZombie()) {
                std::cerr << "Cannot open cutfile " << cutfilename << " .. skipping.." << std::endl;
                continue;
            }
            TCutG *cut = (TCutG*)(f.Get("CUTG"));
            cutsMap.insert(std::make_pair(cutname,static_cast<TObject*>(cut)));
            f.Close();
        } //else
    }//for loop
    infile.close();
}

void HistPlotter::PrintObjects() {
    /*
        void PrintObjects()
        Prints the contents of the unordered_maps oMap and cutsMap to facilitate debugging

    */
    std::cout << "Type | Name " << std::endl;
    std::cout << "---- | --------------------- " << std::endl;
    for(auto it=oMap.begin(); it!=oMap.end(); it++ ) {
        std::cout << it->second->ClassName() << " | "<< it->first << std::endl;
    }
    for(auto it=cutsMap.begin(); it!=cutsMap.end(); it++ ) {
        std::cout << it->second->ClassName() << " | "<< it->first << std::endl;
    }
    std::cout << "---- | --------------------- " << std::endl;
}

#endif
