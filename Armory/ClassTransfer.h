#ifndef ClassTransfer_h
#define ClassTransfer_h

#include "TBenchmark.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TMacro.h"
#include "TGraph.h"
#include <vector>
#include <fstream>

#include "Isotope.h"

class ReactionConfig{

public:

  ReactionConfig(){}
  ~ReactionConfig(){}

  int beamA, beamZ;       
  int targetA, targetZ;     
  int recoilLightA, recoilLightZ;
  int recoilHeavyA, recoilHeavyZ;

  float beamEnergy;         ///MeV/u
  float beamEnergySigma;    ///beam-energy_sigma_in_MeV/u
  float beamAngle;          ///beam-angle_in_mrad
  float beamAngleSigma;     ///beam-emittance_in_mrad
  float beamX;              ///x_offset_of_Beam_in_mm
  float beamY;              ///y_offset_of_Beam_in_mm
  int numEvents;            ///number_of_Event_being_generated
  bool isTargetScattering;  ///isTargetScattering
  float targetDensity;      ///target_density_in_g/cm3
  float targetThickness;    ///targetThickness_in_cm
  std::string beamStoppingPowerFile;         ///stopping_power_for_beam
  std::string recoilLightStoppingPowerFile;  ///stopping_power_for_light_recoil
  std::string recoilHeavyStoppingPowerFile;  ///stopping_power_for_heavy_recoil
  bool isDecay;        ///isDacay
  int heavyDecayA;     ///decayNucleus_A
  int heavyDecayZ;     ///decayNucleus_Z
  bool isRedo;         ///isReDo
  std::vector<float> beamEx;        ///excitation_energy_of_A[MeV]


  void SetReaction(int beamA, int beamZ,
                   int targetA, int targetZ,
                   int recoilA, int recoilZ, float beamEnergy_AMeV){
    this->beamA = beamA;
    this->beamZ = beamZ;
    this->targetA = targetA;
    this->targetZ = targetZ;
    this->recoilLightA = recoilA;
    this->recoilLightZ = recoilZ;

    recoilHeavyA = this->beamA + this->targetA - recoilLightA;
    recoilHeavyZ = this->beamZ + this->targetZ - recoilLightZ;

  } 

  void LoadReactionConfig(TMacro * macro){

    if( macro == NULL ) return ;

    int numLine = macro->GetListOfLines()->GetSize();

    for( int i = 0; i < numLine; i ++){

      std::vector<std::string> str = SplitStr(macro->GetListOfLines()->At(i)->GetName(), " ");

      ///printf("%d | %s\n", i,  str[0].c_str());

      if( str[0].find_first_of("#") == 0 ) break;

      if( i ==  0 ) beamA           = atoi(str[0].c_str());
      if( i ==  1 ) beamZ           = atoi(str[0].c_str());
      if( i ==  2 ) targetA         = atoi(str[0].c_str());
      if( i ==  3 ) targetZ         = atoi(str[0].c_str());
      if( i ==  4 ) recoilLightA    = atoi(str[0].c_str());
      if( i ==  5 ) recoilLightZ    = atoi(str[0].c_str());
      if( i ==  6 ) beamEnergy      = atof(str[0].c_str());
      if( i ==  7 ) beamEnergySigma = atof(str[0].c_str());
      if( i ==  8 ) beamAngle       = atof(str[0].c_str());
      if( i ==  9 ) beamAngleSigma  = atof(str[0].c_str());
      if( i == 10 ) beamX           = atof(str[0].c_str());
      if( i == 11 ) beamY           = atof(str[0].c_str());
      if( i == 12 ) numEvents       = atoi(str[0].c_str());
      if( i == 13 ) {
          if( str[0].compare("false") == 0 ) isTargetScattering = false;
          if( str[0].compare("true")  == 0 ) isTargetScattering = true;
      }
      if( i == 14 ) targetDensity   = atof(str[0].c_str());
      if( i == 15 ) targetThickness = atof(str[0].c_str());
      if( i == 16 ) beamStoppingPowerFile = str[0];
      if( i == 17 ) recoilLightStoppingPowerFile = str[0];
      if( i == 18 ) recoilHeavyStoppingPowerFile = str[0];
      if( i == 19 ) {
          if( str[0].compare("false") == 0 ) isDecay = false;
          if( str[0].compare("true")  == 0 ) isDecay = true;
      }
      if( i == 20 ) heavyDecayA = atoi(str[0].c_str());
      if( i == 21 ) heavyDecayZ = atoi(str[0].c_str());
      if( i == 22 ) {
          if( str[0].compare("false") == 0 ) isRedo = false;
          if( str[0].compare("true" ) == 0 ) isRedo = true;
      }
      
      if( i >= 23) {
          beamEx.push_back( atof(str[0].c_str()) );
      }
    }

      recoilHeavyA = beamA + targetA - recoilLightA;
      recoilHeavyZ = beamZ + targetZ - recoilLightZ;

  }

  void PrintReactionConfig(){

  printf("=====================================================\n");
  printf("   beam : A = %3d, Z = %2d \n", beamA, beamZ);
  printf(" target : A = %3d, Z = %2d \n", targetA, targetZ);
  printf("  light : A = %3d, Z = %2d \n", recoilLightA, recoilLightZ);

  printf(" beam Energy : %.2f +- %.2f MeV/u, dE/E = %5.2f %%\n", beamEnergy, beamEnergySigma, beamEnergySigma/beamEnergy);
  printf("       Angle : %.2f +- %.2f mrad\n", beamAngle, beamAngleSigma);
  printf("      offset : (x,y) = (%.2f, %.2f) mmm \n", beamX, beamY);

  printf("##### number of Simulation Events : %d \n", numEvents);
  
  printf("    is target scattering : %s \n", isTargetScattering ? "Yes" : "No");

  if(isTargetScattering){
    printf(" target density : %.f g/cm3\n", targetDensity);
    printf("      thickness : %.f cm\n", targetThickness);
    printf("         beam stopping file : %s \n", beamStoppingPowerFile.c_str());
    printf(" recoil light stopping file : %s \n", recoilLightStoppingPowerFile.c_str());
    printf(" recoil heavy stopping file : %s \n", recoilHeavyStoppingPowerFile.c_str());
  }
  
  printf("       is simulate decay : %s \n", isDecay ? "Yes" : "No");
  if( isDecay ){
    printf(" heavy decay : A = %d, Z = %d \n", heavyDecayA, heavyDecayZ);
  }
  printf(" is Redo until hit array : %s \n", isRedo ? "Yes" : "No");

  printf(" beam Ex : %.2f MeV \n", beamEx[0]);
  for( int i = 1; i < (int) beamEx.size(); i++){
    printf("         %.2f MeV \n", beamEx[i]);
  }
  
  printf("=====================================================\n");

  }

};

//=======================================================
//#######################################################
// Class for Transfer Reaction
// reaction notation A(a,b)B
// A = incident particle
// a = target
// b = light scattered particle
// B = heavy scattered particle
//======================================================= 
class TransferReaction {
public:
  TransferReaction();
  ~TransferReaction();

  void SetA(int A, int Z, double Ex);
  void Seta(int A, int Z);
  void Setb(int A, int Z);
  void SetB(int A, int Z);
  void SetIncidentEnergyAngle(double KEA, double theta, double phi);
  void SetExA(double Ex);
  void SetExB(double Ex);
  void SetReactionFromFile(string settingFile);
  
  TString GetReactionName();
  TString GetReactionName_Latex();
  
  ReactionConfig GetRectionConfig() { return reaction;}
  
  double GetMass_A(){return mA + ExA;}
  double GetMass_a(){return ma;}
  double GetMass_b(){return mb;}
  double GetMass_B(){return mB + ExB;}
  
  double GetCMTotalKE() {return Etot - mA - ma;}
  double GetQValue()    {return mA + ExA + ma - mb - mB - ExB;}
  double GetMaxExB()    {return Etot - mb - mB;}
  
  TLorentzVector GetPA(){return PA;}
  TLorentzVector GetPa(){return Pa;}
  TLorentzVector GetPb(){return Pb;}
  TLorentzVector GetPB(){return PB;}
  
  void CalReactionConstant();
  
  TLorentzVector * Event(double thetaCM, double phiCM);
  
  double GetEx(){return Ex;}
  double GetThetaCM(){return thetaCM;}
  
  double GetMomentumbCM()   {return p;}
  double GetReactionBeta()  {return beta;}
  double GetReactionGamma() {return gamma;}
  double GetCMTotalEnergy() {return Etot;}
  
private:

  ReactionConfig reaction;

  string nameA, namea, nameb, nameB;
  double thetaIN, phiIN;
  double mA, ma, mb, mB;

  double TA, T; // TA = KE of A pre u, T = total energy
  double ExA, ExB;
  double Ex, thetaCM; //calculated Ex using inverse mapping from e and z to thetaCM
  
  bool isReady;
  bool isBSet;
  
  double k; // CM Boost momentum
  double beta, gamma; //CM boost beta
  double Etot;
  double p; // CM frame momentum of b, B
  
  TLorentzVector PA, Pa, Pb, PB;

  TString format(TString name);
   
};

TransferReaction::TransferReaction(){
  
   thetaIN = 0.;
   phiIN = 0.;
   SetA(24, 12, 0);
   Seta(4,2);
   Setb(1,1);
   SetB(27,13);
   TA = 2.5;
   T = TA * reaction.beamA;
   
   ExA = 0;
   ExB = 0;
   
   Ex = TMath::QuietNaN();
   thetaCM = TMath::QuietNaN();
   
   CalReactionConstant();
   
   TLorentzVector temp (0,0,0,0);
   PA = temp;
   Pa = temp;
   Pb = temp;
   PB = temp;
   
}

TransferReaction::~TransferReaction(){

}

void TransferReaction::SetA(int A, int Z, double Ex = 0){
  Isotope temp (A, Z);
  mA = temp.Mass;
  reaction.beamA = A;
  reaction.beamZ = Z;
  ExA = Ex;
  nameA = temp.Name;
  isReady = false;
  isBSet = true;
  
}

void TransferReaction::Seta(int A, int Z){
  Isotope temp (A, Z);
  ma = temp.Mass;
  reaction.targetA = A;
  reaction.targetZ = Z;
  namea = temp.Name;
  isReady = false;
  isBSet = false;
}

void TransferReaction::Setb(int A, int Z){
  Isotope temp (A, Z);
  mb = temp.Mass;
  reaction.recoilLightA = A;
  reaction.recoilLightZ = Z;
  nameb = temp.Name;
  isReady = false;
  isBSet = false;
}
void TransferReaction::SetB(int A, int Z){
  Isotope temp (A, Z);
  mB = temp.Mass;
  reaction.recoilHeavyA = A;
  reaction.recoilHeavyZ = Z;
  nameB = temp.Name;
  isReady = false;
  isBSet = true;
}

void TransferReaction::SetIncidentEnergyAngle(double KEA, double theta, double phi){
  this->TA = KEA;
  this->T = TA * reaction.beamA;
  this->thetaIN = theta;
  this->phiIN = phi;
  isReady = false;
}

void TransferReaction::SetExA(double Ex){
  this->ExA = Ex;
  isReady = false;
}

void TransferReaction::SetExB(double Ex){
  this->ExB = Ex;
  isReady = false;
}

void TransferReaction::SetReactionFromFile(string settingFile){

   TMacro * haha = new TMacro();
   if( haha->ReadFile(settingFile.c_str()) > 0 ) {
      reaction.LoadReactionConfig(haha);

     SetA(reaction.beamA, reaction.beamZ);
     Seta(reaction.targetA, reaction.targetZ);
     Setb(reaction.recoilLightA, reaction.recoilLightZ);
     SetB(reaction.recoilHeavyA, reaction.recoilHeavyZ);

     SetIncidentEnergyAngle(reaction.beamEnergy, 0, 0);
     CalReactionConstant();
   }else{
     
     printf("cannot read file %s.\n", settingFile.c_str());
     isReady = false;
   }
   
}

TString TransferReaction::GetReactionName(){
  TString rName;
  rName.Form("%s(%s,%s)%s", nameA.c_str(), namea.c_str(), nameb.c_str(), nameB.c_str()); 
  return rName;
}

TString TransferReaction::format(TString name){
  if( name.IsAlpha() ) return name;
  int len = name.Length();
  TString temp = name;
  TString temp2 = name;
  if( temp.Remove(0, len-2).IsAlpha()){
     temp2.Remove(len-2);
  }else{
     temp = name;
     temp.Remove(0, len-1);
     temp2.Remove(len-1);
  }
  return "^{"+temp2+"}"+temp;
}

TString TransferReaction::GetReactionName_Latex(){
  TString rName;
  rName.Form("%s(%s,%s)%s", format(nameA).Data(), format(namea).Data(), format(nameb).Data(), format(nameB).Data()); 
  return rName;
}

void TransferReaction::CalReactionConstant(){
   if( !isBSet){
      reaction.recoilHeavyA = reaction.beamA + reaction.targetA - reaction.recoilLightA;
      reaction.recoilHeavyZ = reaction.beamZ + reaction.targetZ - reaction.recoilLightZ;
      Isotope temp (reaction.recoilHeavyA, reaction.recoilHeavyZ);
      mB = temp.Mass;
      isBSet = true;
   }
   
   k = TMath::Sqrt(TMath::Power(mA + ExA + T, 2) - (mA + ExA) * (mA + ExA)); 
   beta = k / (mA + ExA + ma + T);
   gamma = 1 / TMath::Sqrt(1- beta * beta);   
   Etot = TMath::Sqrt(TMath::Power(mA + ExA + ma + T,2) - k * k);
   p = TMath::Sqrt( (Etot*Etot - TMath::Power(mb + mB + ExB,2)) * (Etot*Etot - TMath::Power(mb - mB - ExB,2)) ) / 2 / Etot;
   
   PA.SetXYZM(0, 0, k, mA + ExA);
   PA.RotateY(thetaIN);
   PA.RotateZ(phiIN);
   
   Pa.SetXYZM(0,0,0,ma);
   
   isReady = true;
}

TLorentzVector * TransferReaction::Event(double thetaCM, double phiCM)
{
   if( isReady == false ){
      CalReactionConstant();
   }

   //TLorentzVector Pa(0, 0, 0, ma);
   
   //---- to CM frame
   TLorentzVector Pc = PA + Pa;
   TVector3 b = Pc.BoostVector();
   
   TVector3 vb(0,0,0);
   
   if( b.Mag() > 0 ){
     TVector3 v0 (0,0,0);
     TVector3 nb = v0 - b;
     
     TLorentzVector PAc = PA; 
     PAc.Boost(nb);
     TVector3 vA = PAc.Vect();
     
     TLorentzVector Pac = Pa;
     Pac.Boost(nb);
     TVector3 va = Pac.Vect();
     
     //--- construct vb
     vb = va;
     vb.SetMag(p);

     TVector3 ub = vb.Orthogonal();
     vb.Rotate(thetaCM, ub);
     vb.Rotate(phiCM + TMath::PiOver2(), va); // somehow, the calculation turn the vector 90 degree.
     //vb.Rotate(phiCM , va); // somehow, the calculation turn the vector 90 degree.
   }
   
   //--- from Pb
   TLorentzVector Pbc;
   Pbc.SetVectM(vb, mb);
   
   //--- from PB
   TLorentzVector PBc;
   //PBc.SetVectM(vB, mB + ExB);
   PBc.SetVectM(-vb, mB + ExB);
   
   //---- to Lab Frame
   TLorentzVector Pb = Pbc;
   Pb.Boost(b);
   TLorentzVector PB = PBc;
   PB.Boost(b);
   
   TLorentzVector * output = new TLorentzVector[4];
   output[0] = PA;
   output[1] = Pa;
   output[2] = Pb;
   output[3] = PB;
   
   this->Pb = Pb;
   this->PB = PB;
   
   
   return output;   
}

#endif 
