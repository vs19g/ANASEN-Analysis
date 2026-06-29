#define TrackRecon_cxx

#define RAW_HISTOS
// #define VTX_GATES
// #define AL_BEAM
// #define F_BEAM
// #define nA_Analysis

Int_t colors[40] = {
    kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange,
    kSpring, kTeal, kAzure, kViolet, kPink, kGray, kWhite,
    kRed + 2, kGreen + 2, kBlue + 2, kYellow + 2, kMagenta + 2, kCyan + 2, kOrange + 2,
    kSpring + 2, kTeal + 2, kAzure + 2, kViolet + 2, kPink + 2,
    kRed - 7, kGreen - 7, kBlue - 7, kYellow - 7, kMagenta - 7, kCyan - 7, kOrange - 7,
    kSpring - 7, kTeal - 7, kAzure - 7, kViolet - 7, kPink - 7, kGray + 2};

#include "TrackRecon.h"
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
#include <TVector3.h>
#include <TVector2.h>
#include <TRandom3.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <algorithm>

// --- Analysis Control Flags ---
bool process_alpha_proton_scattering = false,
     doMiscHistograms = true,
     doPCSX3ClusterAnalysis = true,
     doPCQQQClusterAnalysis = true,
     doOldAnalysis = false,
     do27AlapAnalysis = false,
     BenchMark = true,
     reactiondata = false;

// --- Geometry, Calibration, & Model Variables ---
double source_vertex = 53.0,
       z_entrance = -174.3 - 9.7 - 100.0,
       dither_sigma = 8.0,
       dither_sigma_c0 = 16.0,
       Gain = 1.0,
       cathode_gain = 1.0,
       a1c1_cfrac_split = 0.0,
       a1c1_missing_fmax = 2.0,
       a1c1_lowband_rfactor = 0.0,
       a1c1_z_scale_qqq = 0.0111081,
       a1c1_z_off_qqq = 34.501,
       a1c1_z_scale_sx3 = 0.0,
       a1c1_z_off_sx3 = 2.52614;

// --- Immutable Constants ---
const double qqq_z = 105.0,
             anode_gain = 1.5146e-5,
             sx3_phi_pitch = 6.5 * (M_PI / 180.0),
             qqq_wedge_pitch = (87.0 / 16.0) * (M_PI / 180.0),
             qqq_ring_pitch = 48.0 / 16.0;
std::string dataset;
int CO2percent;

TF1 pcfix_func("func", model_invert, -200, 200);

// --- A1C1 linear centre-fold model: cfrac = cfmin[cell] + k[cell]*fold ----------
// fold = |z - z_cellcentre| / halfcell in [0,1].
const double a1c1_zg[8] = {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};

static const double a1c1_cfmin_17F[7] = {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20};
static const double a1c1_k_17F[7] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
static const double a1c1_cfmin_27Al[7] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
static const double a1c1_k_27Al[7] = {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20};

// active per-cell set, populated by dataset in Begin()
double a1c1_cfmin_cell[7] = {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20};
double a1c1_k_cell[7] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

// --- Dead / missing PC wires -------------------------------------------------
// Some channels are unresponsive in a run (the white vertical gaps in
// PC_Index_Vs_Energy). A genuine two-wire track that straddles a dead wire
// collapses to a single fired wire: a "pseudo-1-wire" event. We still treat it
// as A1C1, but flag it so the excitation function can be split three ways:
//   _all       : every A1C1 event (unchanged from before)
//   _true1w    : neither neighbour of the fired anode/cathode is dead (genuine single)
//   _missingw  : a neighbouring wire is dead (suspected masked two-wire event)
// Index 0-23 within EACH plane; 1 = dead. Read the dead channels off the gaps
// in PC_Index_Vs_Energy (anode = index 0-23, cathode = index 24-47 -> w = idx-24)
// and fill the per-dataset arrays. Default all-alive => everything is _true1w
// and _missingw stays empty (no behaviour change until channels are entered).

static std::vector<int> a1c1_dead_anode_17F = {9, 12}; // 1 can be recovered
static std::vector<int> a1c1_dead_cathode_17F = {};    // 0,13,15 can be recovered
static std::vector<int> a1c1_dead_anode_27Al = {0, 12, 19};
static std::vector<int> a1c1_dead_cathode_27Al = {13};
std::vector<int> *a1c1_dead_anode = &a1c1_dead_anode_17F; // active set, chosen in Begin()
std::vector<int> *a1c1_dead_cathode = &a1c1_dead_cathode_17F;

// True if a neighbouring wire (index +/-1) of the fired anode OR cathode is
// dead -- i.e. this single-wire event may actually be a masked two-wire one.
inline bool a1c1_missing_neighbor(int awire, int cwire)
{
  auto deadAdj = [](const std::vector<int> &dead, int w)
  {
    if (w < 0 || w >= 24)
      return false;
    auto isDead = [&](int x)
    { return std::find(dead.begin(), dead.end(), x) != dead.end(); };
    return isDead(w - 1) || isDead(w + 1);
  };
  return deadAdj(*a1c1_dead_anode, awire) || deadAdj(*a1c1_dead_cathode, cwire);
}

// --- A1C1 LOW BAND (incomplete charge integration) -----------------------------
// 17F shows a second, parallel cfrac band (~0.10 -> 0.15 over the fold) where the
// cathode charge is only partially integrated. It still tracks z
// so it gets its own per-cell cfmin/k. Events with cfrac below a1c1_cfrac_split
// are reconstructed with this set instead of being rejected. a1c1_cfrac_split<=0
// disables the low band (e.g. 27Al, which shows no second band).
static const double a1c1_cfmin2_17F[7] = {0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10};
static const double a1c1_k2_17F[7] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
static const double a1c1_cfmin2_27Al[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // no low band
static const double a1c1_k2_27Al[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

double a1c1_cfmin2_cell[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double a1c1_k2_cell[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

inline double a1c1_zcorr(double z_a1c0, bool isQQQ)
{
  double scale = isQQQ ? a1c1_z_scale_qqq : a1c1_z_scale_sx3;
  double off = isQQQ ? a1c1_z_off_qqq : a1c1_z_off_sx3;
  return z_a1c0 * (1.0 - scale) - off;
}

// Sub-cell A1C1 z from cfrac (linear centre-fold). zf = crossover z (fired
// cathode, on the grid). Anchors on the fired wire, folds about the adjacent cell
// centre, inverts cfrac = cfmin+k*fold per cell. Picks the LOW band (band=1) when
// cfrac < a1c1_cfrac_split, else the MAIN band (band=0). The side (which adjacent
// cell) is the caller's job via a1c1_pick_side; both candidates are returned. Callers
// apply their own acceptance via the returned per-candidate flags:
//   inband  : f in [0,1]            (cfrac within the calibrated band)
//   pitchok : |pcz - zf| <= pitch   (consistent with the fired wire)
// A single candidate cell's complete charge-sharing solution.
struct A1C1CellSol
{
  int cell = -1;
  double pcz = -99999;
  double f = 0.0;
  double pitch = 0.0;
  bool inband = false;
  bool pitchok = false;
};

struct A1C1Sol
{
  int band;
  double cfrac_used; // cfrac actually inverted (after the low-band r-space fold)
  double pcz_lo;     // candidate z for the cell BELOW the fired wire (lower z) == lo.pcz
  double pcz_hi;     // candidate z for the cell ABOVE the fired wire (higher z) == hi.pcz
  A1C1CellSol hi;    // full solution for the cell ABOVE the fired wire
  A1C1CellSol lo;    // full solution for the cell BELOW the fired wire
};

// Charge-sharing inversion for ONE candidate cell. This contains ALL of the inversion
// mathematics: the sub-cell coordinate f, the position pcz, and the in-band / pitch
// acceptance. The SIDE (which of the two cells adjacent to the fired wire) is left to
// the caller -- a1c1_pick_side resolves it from the Si hit + beam axis. The sub-cell
// MAGNITUDE is side-independent; only the cell choice is ambiguous from the charge.
inline A1C1CellSol solve_cell(int cell, int wf, double zf, double cfrac,
                              const double *cfmin, const double *kk, int awire, int cwire)
{
  A1C1CellSol s;
  s.cell = cell;
  s.pcz = zf;

  if (cell < 0 || cell > 6)
    return s;

  double zc = 0.5 * (a1c1_zg[cell] + a1c1_zg[cell + 1]);   // cell centre
  double half = 0.5 * (a1c1_zg[cell] - a1c1_zg[cell + 1]); // half-cell width
  double pitch = a1c1_zg[cell] - a1c1_zg[cell + 1];        // full wire spacing

  if (half <= 0.0 || kk[cell] <= 0.0)
    return s;

  s.pitch = pitch;

  // f = 0 -> cell centre, f = 1 -> fired wire. Outside [0,1] = outside the band.
  s.f = (cfrac - cfmin[cell]) / kk[cell];

  // sign maps increasing f toward the fired cathode wire.
  double sgn = (a1c1_zg[wf] >= zc) ? +1.0 : -1.0;
  s.pcz = zc + sgn * s.f * half;

  // A dead neighbouring wire lets the fired wire collect the shared charge, so cfrac
  // (hence f) runs past the usual cfmin+k ceiling -- lift the in-band ceiling so these
  // events are accepted instead of rejected at f=1. pcz is never clamped.
  double fmax = a1c1_missing_neighbor(awire, cwire) ? a1c1_missing_fmax : 1.0;
  s.inband = (s.f >= 0.0 && s.f <= fmax);

  // Reconstructed position should remain within one cell pitch of the fired wire.
  s.pitchok = (TMath::Abs(s.pcz - zf) <= pitch);

  return s;
}

inline A1C1Sol a1c1_solve(double cfrac, double zf, int cwire = -1, double anodeE = -1, int awire = -1)
{
  A1C1Sol s{0, cfrac, zf, zf, {}, {}};
  const double *cfmin = a1c1_cfmin_cell;
  const double *kk = a1c1_k_cell;
  if (a1c1_cfrac_split > 0.0 && cfrac >= 0.0 && cfrac < a1c1_cfrac_split)
  {
    s.band = 1;
    if (a1c1_lowband_rfactor > 0.0 && cfrac > 0.0 && cfrac < 1.0)
    {
      double r = cfrac / (1.0 - cfrac);
      r *= a1c1_lowband_rfactor;
      cfrac = r / (1.0 + r);
    }
    else
    {
      cfmin = a1c1_cfmin2_cell;
      kk = a1c1_k2_cell;
    }
  }
  s.cfrac_used = cfrac;

  // Identify the fired cathode wire from the measured cathode position zf
  // which is assumed to lie on one side of the calibrated cathode wire positions a1c1_zg[].
  int wf = 0;
  for (int i = 1; i < 8; ++i)
    if (TMath::Abs(a1c1_zg[i] - zf) < TMath::Abs(a1c1_zg[wf] - zf))
      wf = i;

  // Full inversion for BOTH cells adjacent to the fired wire. The side is ambiguous
  // from the charge alone; a1c1_pick_side resolves it from the Si hit + beam axis.
  s.hi = solve_cell(wf - 1, wf, zf, cfrac, cfmin, kk, awire, cwire); // cell above (higher z)
  s.lo = solve_cell(wf, wf, zf, cfrac, cfmin, kk, awire, cwire);     // cell below (lower z)
  s.pcz_hi = s.hi.pcz;
  s.pcz_lo = s.lo.pcz;
  // Side selection (which candidate) is the caller's job: a1c1_pick_side uses the Si hit
  // + beam axis, then the caller reads the winning candidate (s.hi / s.lo).
  return s;
}

// Beam-axis 2-hypothesis side test. Given the two candidate PC z's (the cells either
// side of the fired wire), reconstruct the vertex (closest approach to the beam axis)
// for each and keep the beam-axis-consistent one. Acceptance is PERP-ONLY: a candidate
// is physical if its vertex sits within a1c1_side_perp_max of the beam axis.
// The old |vtx.Z - source_vertex| window was removed: vertex-Z consistency only has
// meaning for fixed-vertex SOURCE runs. For reaction/proton data the interaction
// vertex is distributed along the beam, so a window anchored on source_vertex rejected
// essentially every genuine event. status: 0 one side physical, 1 both physical,
// 2 neither (reject). When both candidates are physical the side is genuinely
// ambiguous without a vertex prior, so fall back to the smaller-Perp candidate.
double a1c1_side_perp_max = 20.0; // beam-axis Perp gate (mm)

// Which of the two candidate cells the beam-axis test selects.
enum class SideChoice
{
  High, // the cell ABOVE the fired wire (pcz_hi)
  Low   // the cell BELOW the fired wire (pcz_lo)
};

// Returns WHICH candidate wins (not the position) so the caller can copy the full
// winning solution. status: 0 one side physical, 1 both physical, 2 neither (reject).
inline SideChoice a1c1_pick_side(const TVector3 &si, double cx, double cy,
                                 double pcz_lo, double pcz_hi, int &status)
{
  auto vtxZP = [&](double pcz, double &z, double &perp)
  {
    TVector3 pc(cx, cy, pcz);
    TVector3 v = pc - si;
    double d = v.X() * v.X() + v.Y() * v.Y();
    double tm = (d > 0.0) ? -(si.X() * v.X() + si.Y() * v.Y()) / d : 0.0;
    TVector3 vtx = si + tm * v;
    z = vtx.Z();
    perp = vtx.Perp();
  };
  double zl, pl, zh, ph;
  vtxZP(pcz_lo, zl, pl);
  vtxZP(pcz_hi, zh, ph);
  bool okl = (pl <= a1c1_side_perp_max);
  bool okh = (ph <= a1c1_side_perp_max);
  status = (okl || okh) ? ((okl && okh) ? 1 : 0) : 2;
  if (okl && !okh)
    return SideChoice::Low;
  if (okh && !okl)
    return SideChoice::High;
  return (pl <= ph) ? SideChoice::Low : SideChoice::High; // both physical: smaller-Perp side
}

TGraph *MeV_to_cm = NULL, *cm_to_MeV = NULL;
TGraph *MeV_to_cm_p = NULL, *cm_to_MeVp = NULL;
TGraph *MeV_to_cm_27Al = NULL, *cm_to_MeV_27Al = NULL;
TGraph *MeV_to_cm_17F = NULL, *cm_to_MeV_17F = NULL;

// declaring masses for kinematics calculations
double mass_27Al = 26.981538;
double mass_4He = 4.002603254;
double mass_1H = 1.007825032;
double mass_30Si = 29.973770;
double mass_17F = 17.002095;
double mass_20Ne = 19.992440;

// new Parabola for 4wire shift
double z_to_crossover_rho(double z)
{
  return 1.65896E-4 * z * z + 4.61626E-8 * z + 32.067;
}

// Global instances
PW pwinstance;
TVector3 hitPos;
double qqqenergy, qqqtimestamp;
class Event
{
public:
  Event(TVector3 p, double e1, double e2, double t1, double t2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2) {}
  Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2) {}
  // Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2, int m1, int m2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2), multi1(m1), multi2(m2) {}
  Event(TVector3 p, double e1, double e2, double t1, double t2, int a, int c, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), Anodech(a), Cathodech(c), ch1(c1), ch2(c2) {}

  TVector3 pos;
  int ch1 = -1;        // int(ch1/16) gives qqq id, ch1%16 gives ring#
  int ch2 = -1;        // int(ch2/16) gives qqq id, ch2%16 gives wedge#
  double Energy1 = -1; // Front for QQQ, Anode for PC
  double Energy2 = -1; // Back for QQQ, Cathode for PC
  double Time1 = -1;
  double Time2 = -1;
  int Anodech = -1;
  int Cathodech = -1;

  // misc elements;
  int multi1 = -1, multi2 = -1;
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

bool PCQQQTimeCut = false;
bool PCSX3TimeCut = false, PCASX3TimeCut = false, PCCSX3TimeCut = false;
double anodeT = -99999, cathodeT = 99999;
int anodeIndex = -1, cathodeIndex = -1;

void protonAlphaHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events);
void miscHistograms_oneWire(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters);
void protonMiscHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events);
void protonMiscHistograms_sx3(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events);
void PCSX3ClusterAnalysis(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, const std::vector<std::vector<std::tuple<int, double, double>>> &cClusters);
void PCQQQClusterAnalysis(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, const std::vector<std::vector<std::tuple<int, double, double>>> &cClusters);

void TrackRecon::Begin(TTree * /*tree*/)
{
  pcfix_func.SetNpx(100000);
  ///// ---------Set Environment Variables--------- /////
  TString option = GetOption();
  if (option != "")
    plotter = new HistPlotter(option.Data(), "TFILE");
  else
    plotter = new HistPlotter("Analyzer_SX3.root", "TFILE");

  if (getenv("reactiondata"))
  {
    reactiondata = std::atoi(getenv("reactiondata"));
    std::cout << "Analyzing dataset as reactiondata" << std::endl;
  }

  if (getenv("DATASET"))
    dataset = std::string(getenv("DATASET"));
  if (getenv("source_vertex"))
    source_vertex = (double)std::atof(std::string(getenv("source_vertex")).c_str());

  if (getenv("CO2percent"))
    CO2percent = std::atoi(getenv("CO2percent"));
  std::cout << "CO2 percent set to  " << CO2percent << std::endl;

  if (getenv("source_vertex"))
    source_vertex = (double)std::atof(std::string(getenv("source_vertex")).c_str());

  if (getenv("DITHER_SIGMA"))
  {
    dither_sigma = std::atof(getenv("DITHER_SIGMA"));
    dither_sigma_c0 = dither_sigma;
    std::cout << "Dither Sigma set to " << dither_sigma << " mm" << std::endl;
  }

  if (getenv("Gain"))
    Gain = std::atof(getenv("Gain"));

  if (getenv("CATHODE_GAIN"))
    cathode_gain = std::atof(getenv("CATHODE_GAIN"));

  // --- A1C1 per-cell linear centre-fold constants: select the static, offline-
  // optimised set for this dataset (cfmin floats with the arbitrary gain, so the
  // values are dataset-specific). Re-fit them with fit_a1c1_cfrac.C on prebuilt
  // data and paste into the a1c1_{cfmin,k}_<dataset> arrays above. ---
  const double *cfmin_src = a1c1_cfmin_17F;
  const double *k_src = a1c1_k_17F;
  const double *cfmin2_src = a1c1_cfmin2_17F;
  const double *k2_src = a1c1_k2_17F;
  a1c1_cfrac_split = 0.15;    // 17F: split in the valley between low/main bands (cfrac<0.15 = low band)
  a1c1_lowband_rfactor = 7.0; // 17F: fold low band onto the main band (r-space).
                              // From the source-run cfrac_vs_sx3E (both bands at the
                              // alpha energy): g = r_main/r_low = (0.44/0.56)/(0.10/0.90)
                              // ~ 7.0, i.e. r_low*7 -> cfrac 0.10 maps to ~0.44.
  // a1c1_z_scale_qqq = 0.0;       // 17F: QQQ z scaling (REFIT from source runs)
  // a1c1_z_scale_sx3 = 0.0;       // 17F: SX3 z scaling (REFIT from source runs)
  // a1c1_z_off_qqq = 0.0;         // 17F: QQQ constant offset mm (REFIT)
  // a1c1_z_off_sx3 = 0.0;         // 17F: SX3 constant offset mm (REFIT)
  a1c1_dead_anode = &a1c1_dead_anode_17F;
  a1c1_dead_cathode = &a1c1_dead_cathode_17F;
  if (dataset == "27Al")
  {
    cfmin_src = a1c1_cfmin_17F;
    k_src = a1c1_k_17F;
    cfmin2_src = a1c1_cfmin2_27Al;
    k2_src = a1c1_k2_27Al;
    a1c1_cfrac_split = 0.0;     // 27Al: no second band, low band disabled
    a1c1_lowband_rfactor = 0.0; // 27Al: nothing to fold
    // a1c1_z_scale_qqq = 0.0;       // 27Al: QQQ z scaling (REFIT from source runs)
    // a1c1_z_scale_sx3 = 0.0;       // 27Al: SX3 z scaling (REFIT from source runs)
    // a1c1_z_off_qqq = 0.0;         // 27Al: QQQ constant offset mm (REFIT)
    // a1c1_z_off_sx3 = 0.0;         // 27Al: SX3 constant offset mm (REFIT)
    a1c1_dead_anode = &a1c1_dead_anode_27Al;
    a1c1_dead_cathode = &a1c1_dead_cathode_27Al;
  }
  if (getenv("A1C1_LOWBAND_RFACTOR"))
    a1c1_lowband_rfactor = std::atof(getenv("A1C1_LOWBAND_RFACTOR"));
  if (getenv("A1C1_Z_SCALE_QQQ"))
    a1c1_z_scale_qqq = std::atof(getenv("A1C1_Z_SCALE_QQQ"));
  if (getenv("A1C1_Z_SCALE_SX3"))
    a1c1_z_scale_sx3 = std::atof(getenv("A1C1_Z_SCALE_SX3"));
  if (getenv("A1C1_Z_OFF_QQQ"))
    a1c1_z_off_qqq = std::atof(getenv("A1C1_Z_OFF_QQQ"));
  if (getenv("A1C1_Z_OFF_SX3"))
    a1c1_z_off_sx3 = std::atof(getenv("A1C1_Z_OFF_SX3"));
  for (int i = 0; i < 7; ++i)
  {
    a1c1_cfmin_cell[i] = cfmin_src[i];
    a1c1_k_cell[i] = k_src[i];
    a1c1_cfmin2_cell[i] = cfmin2_src[i];
    a1c1_k2_cell[i] = k2_src[i];
  }
  std::cout << "A1C1 per-cell constants: using static " << (dataset.empty() ? "(default 17F)" : dataset)
            << " set; low-band split cfrac<" << a1c1_cfrac_split
            << "; low-band r-fold " << (a1c1_lowband_rfactor > 0.0 ? "ON x" : "OFF (")
            << a1c1_lowband_rfactor << (a1c1_lowband_rfactor > 0.0 ? "" : ")") << std::endl;

  pwinstance.ConstructGeo();

  for (int i = 0; i < 48; i++)
  {
    pcSlope[i] = 1.0;
    pcIntercept[i] = 0.0;
  }

  //  ------------Load PC Calibrations-------------- ///
  std::ifstream inputFile("slope_intercept_results_" + dataset + ".dat");
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

  // ------------Load QQQ Calibrations-------------- ///
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

  //  ------------Load SX3 Calibrations--------------- ///
  {
    std::ifstream infile("sx3cal/" + dataset + "/backgains.dat");
    std::string temp;
    int backpos, frontpos, clkpos;
    if (infile.is_open())
      while (infile >> clkpos >> temp >> frontpos >> temp >> backpos >> sx3BackGain[clkpos][frontpos][backpos])
        ; // std::cout << sx3BackGain[clkpos][frontpos][backpos] << std::endl;
    infile.close();

    infile.open("sx3cal/" + dataset + "/frontgains.dat");
    if (infile.is_open())
      while (infile >> clkpos >> temp >> temp >> frontpos >> sx3FrontOffset[clkpos][frontpos] >> sx3FrontGain[clkpos][frontpos])
        ; // std::cout << sx3FrontOffset[clkpos][frontpos] << " " << sx3FrontGain[clkpos][frontpos] << std::endl;
    infile.close();

    infile.open("sx3cal/" + dataset + "/rightgains.dat");
    if (infile.is_open())
      while (infile >> clkpos >> frontpos >> temp >> sx3RightGain[clkpos][frontpos])
      {
        sx3RightGain[clkpos][frontpos] = TMath::Abs(sx3RightGain[clkpos][frontpos]);
      }
    infile.close();
  }

  // ------------- ELOSS Correction read in from tables -------------

  if (dataset == "17F" && CO2percent == 3)
  {
    if (!MeV_to_cm)
      MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_20MeV_3pc.dat", "%lf %*lf %lf");
    if (!MeV_to_cm_p)
      MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV_3pc.dat", "%lf %*lf %lf");
    if (!MeV_to_cm_27Al)
      MeV_to_cm_27Al = new TGraph("eloss_calculations/aluminum_lookup_80MeV_3pc.dat", "%lf %*lf %lf");
    if (!MeV_to_cm_17F)
      MeV_to_cm_17F = new TGraph("eloss_calculations/fluorine_lookup_70MeV_3pc.dat", "%lf %*lf %lf");

    // MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_20MeV_3pc_350Torr.dat", "%lf %*lf %lf");
    // MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV_3pc_350Torr.dat", "%lf %*lf %lf");
    // MeV_to_cm_17F = new TGraph("eloss_calculations/fluorine_lookup_70MeV_3pc_350Torr.dat", "%lf %*lf %lf");
    // MeV_to_cm_27Al = new TGraph("eloss_calculations/aluminum_lookup_80MeV_3pc_350Torr.dat", "%lf %*lf %lf");
  }
  else if (dataset == "17F" && CO2percent == 4)
  {
    MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_20MeV_4pc.dat", "%lf %*lf %lf");
    MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV_4pc.dat", "%lf %*lf %lf");
    MeV_to_cm_27Al = new TGraph("eloss_calculations/aluminum_lookup_80MeV_4pc.dat", "%lf %*lf %lf");
    MeV_to_cm_17F = new TGraph("eloss_calculations/fluorine_lookup_70MeV_4pc.dat", "%lf %*lf %lf");

    // MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_20MeV_4pc_350Torr.dat", "%lf %*lf %lf");
    // MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV_4pc_350Torr.dat", "%lf %*lf %lf");
    // MeV_to_cm_17F = new TGraph("eloss_calculations/fluorine_lookup_70MeV_4pc_350Torr.dat", "%lf %*lf %lf");
    // MeV_to_cm_27Al = new TGraph("eloss_calculations/aluminum_lookup_80MeV_4pc_350Torr.dat", "%lf %*lf %lf");
  }
  else
  {
    MeV_to_cm = new TGraph("eloss_calculations/alpha_lookup_20MeV_3pc.dat", "%lf %*lf %lf");
    MeV_to_cm_p = new TGraph("eloss_calculations/proton_lookup_20MeV_3pc.dat", "%lf %*lf %lf");
    MeV_to_cm_27Al = new TGraph("eloss_calculations/aluminum_lookup_80MeV_3pc.dat", "%lf %*lf %lf");
    MeV_to_cm_17F = new TGraph("eloss_calculations/fluorine_lookup_70MeV_3pc.dat", "%lf %*lf %lf");
  }
  cm_to_MeV = new TGraph(MeV_to_cm->GetN(), MeV_to_cm->GetY(), MeV_to_cm->GetX());
  cm_to_MeVp = new TGraph(MeV_to_cm_p->GetN(), MeV_to_cm_p->GetY(), MeV_to_cm_p->GetX());
  cm_to_MeV_27Al = new TGraph(MeV_to_cm_27Al->GetN(), MeV_to_cm_27Al->GetY(), MeV_to_cm_27Al->GetX());
  cm_to_MeV_17F = new TGraph(MeV_to_cm_17F->GetN(), MeV_to_cm_17F->GetY(), MeV_to_cm_17F->GetX());
}

Bool_t TrackRecon::Process(Long64_t entry)
{
  hitPos.Clear();
  qqqenergy = -1;
  qqqtimestamp = -1;
  HitNonZero = false;
  PCQQQTimeCut = false;
  PCSX3TimeCut = false;
  PCASX3TimeCut = false;
  PCCSX3TimeCut = false;
  anodeT = -99999;
  cathodeT = 99999;
  anodeIndex = -1;
  cathodeIndex = -1;
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
  if (dataset == "17F" && reactiondata)
  {
    b_miscMulti->GetEntry(entry);
    b_miscID->GetEntry(entry);
    b_miscCh->GetEntry(entry);
    b_miscE->GetEntry(entry);
    b_miscT->GetEntry(entry);
    b_miscTf->GetEntry(entry);
  }
  // env vars are fixed for a run: read once, not per event
  static const double timecut_low = getenv("timecut_low") ? std::atof(getenv("timecut_low")) : 0;
  static const double timecut_high = getenv("timecut_high") ? std::atof(getenv("timecut_high")) : 1e15;

  if (pc.multi > 0)
  {
    for (int i = 0; i < pc.multi; i++)
    {
      if (pc.t[i] * 1e-9 < timecut_high && pc.t[i] * 1e-9 >= timecut_low)
      {
        // good, keep it moving
      }
      else
      {
        return kTRUE;
      }
    }
  }

  sx3.CalIndex();
  qqq.CalIndex();
  pc.CalIndex();

  std::vector<Event> SX3_Events;
  if (sx3.multi > 1)
  {
    std::array<sx3det, 24> Fsx3;
    // std::cout << "-----" << std::endl;
    bool found_upstream_sx3 = 0;
    for (int i = 0; i < sx3.multi; i++)
    {
      int id = sx3.id[i];
      if (id >= 12)
        continue;
      if (sx3.ch[i] >= 8)
      {
        int sx3ch = sx3.ch[i] - 8;
        sx3ch = (sx3ch + 3) % 4;
        if (id >= 12)
        {
          found_upstream_sx3 = 1;
          // std::cout << Form("f%d(",id) << sx3ch << "," << sx3.e[i] << ") " << std::flush;
        }
        // if(sx3ch==0 || sx3ch==3) continue;
        double value = sx3.e[i];
        int gch = sx3.id[i] * 4 + (sx3.ch[i] - 8);
        if (id < 12)
          Fsx3.at(id).fillevent("BACK", sx3ch, value);
        Fsx3.at(id).ts = static_cast<double>(sx3.t[i]);
#ifdef RAW_HISTOS
        plotter->Fill2D("sx3backs_all_raw", 100, 0, 100, 800, 0, 4096, gch, sx3.e[i]);
#endif
      }
      else
      {
        int sx3ch = sx3.ch[i] / 2;
        double value = sx3.e[i];
        if (id >= 12)
        {
          found_upstream_sx3 = 1;
          // std::cout << Form("b%d(",id) << sx3ch << "," << value << ") " << std::flush;
        }
        if (sx3.ch[i] % 2 == 0)
        {
          Fsx3.at(id).fillevent("FRONT_L", sx3ch, value * sx3RightGain[id][sx3ch]);
        }
        else
        {
          Fsx3.at(id).fillevent("FRONT_R", sx3ch, value);
        }
      }
    } // end for (i in sx3.multi)
    // if(found_upstream_sx3) std::cout << std::endl;

    for (int id = 0; id < 24; id++)
    {
      // std::cout << id << " " << Fsx3.at(id).valid_front_chans.size() << " " << Fsx3.at(id).valid_back_chans.size() << std::endl;;
      try
      {
        Fsx3.at(id).validate();
      }
      catch (std::exception exc)
      {
        std::cout << "oops! anyway " << std::endl;
        continue;
      }
      auto det = Fsx3.at(id);
      if (det.valid)
      {
        // std::cout << det.frontEL << " " << det.frontEL*sx3RightGain[id][det.stripF] << std::endl;
        // plotter->Fill2D("be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_b"+std::to_string(det.stripB),200,-1,1,800,0,8192,det.frontX,det.backE,"evsx");
        // plotter->Fill2D("unmatched_be_vs_x_sx3_id_" + std::to_string(id), 200, -1, 1, 800, 0, 4096, det.frontX, det.backE, "evsx");
        // plotter->Fill2D("unmatched_be_vs_x_sx3", 200, -1, 1, 800, 0, 4096, det.frontX, det.backE, "evsx");
        // plotter->Fill2D("matched_be_vs_x_sx3", 200, -60, 60, 800, 0, 8192, det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx");
        // plotter->Fill2D("matched_be_vs_x_sx3_id_" + std::to_string(id), 200, -60, 60, 800, 0, 8192, det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx");

        // plotter->Fill2D("matched_be_vs_x_sx3_id_" + std::to_string(id) + "_f" + std::to_string(det.stripF), 200, -60, 60, 800, 0, 8192,
        //                 det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx_matched");
        // plotter->Fill2D("fe_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_"+std::to_string(det.stripB),200,-1,1,800,0,4096,det.frontX,det.backE,"evsx");
        // plotter->Fill2D("l_vs_r_sx3_id_" + std::to_string(id) + "_f" + std::to_string(det.stripF), 800, 0, 4096, 800, 0, 4096, det.frontEL, det.frontER, "l_vs_r");
      }
      if (det.valid && (id == 9 || id == 7 || id == 1 || id == 3) && det.stripF != DEFAULT_NULL && det.stripB != DEFAULT_NULL)
      {
        double z = det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF];
        z = z + (75.0 / 2.0) - 3.0; // convert local sx3z to detector global coordinate system as indicated by measurements.
        // Note that this will be different for the upstream barrel, when it gets implemented
        double backE = det.backE * sx3BackGain[id][det.stripF][det.stripB];
        // if(backE<2000) continue;
        det.stripF = 3 - det.stripF;

        double strip_pitch = 40.30 / 4.0;                  // ~10.075 mm
        double y_local = (det.stripF - 1.5) * strip_pitch; // gives -15.11, -5.03, +5.03, +15.11
        double x_local = 88.0 * TMath::Cos(15.0 * M_PI / 180.0);
        double alpha_n = TMath::ATan2(y_local, x_local) * 180. / M_PI;
        double beta_n = 15.0 + alpha_n; // how much to add per strip to the starting position? this is the angle w.r.t an edge of the sx3, the above values run as (-10.08deg, -3.39deg, 3.39deg, 10.08deg)
        double phi_n = ((-id + 0.5) * 30 + beta_n);
        phi_n += 45;
        double rho_at_strip = 88.0 / TMath::Cos(alpha_n * M_PI / 180.0); //*TMath::Cos(15.0*M_PI/180.0) if the edge-length is 88mm

        // if(getenv("flip180")) {
        //	if(std::string(getenv("flip180"))=="1") {
        // if(dataset=="17F")
        //	phi_n+=180;//run 37 in 17F-->
        // }
        //}
        phi_n *= M_PI / 180.; // starting-position phi + strip contribution
        // Event sx3ev(TVector3(88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n),z),backE*0.001,-1,det.ts,-1,det.stripB+4*id,det.stripF+4*id);
        Event sx3ev(TVector3(rho_at_strip * TMath::Cos(phi_n), rho_at_strip * TMath::Sin(phi_n), z), backE * 0.001, -1, det.ts, -1, det.stripB + 4 * id, det.stripF + 4 * id);
        SX3_Events.push_back(sx3ev);
        plotter->Fill2D("sx3backs_gm", 100, 0, 100, 800, 0, 8192, det.stripB + 4 * id, backE);
        plotter->Fill1D("sx3backs_calib", 800, 0, 8192, backE);

        // plotter->Fill2D("SX3CartesianPlot", 200, -100, 100, 200, -100, 100, 88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n), "hCalSX3");
        plotter->Fill2D("SX3CartesianPlot" + std::to_string(id), 200, -100, 100, 200, -100, 100, 88.0 * TMath::Cos(phi_n), 88.0 * TMath::Sin(phi_n), "hCalSX3");
      }
#ifdef RAW_HISTOS
      if (det.valid && det.stripF != DEFAULT_NULL && det.stripB != DEFAULT_NULL)
      {
        plotter->Fill2D("sx3backs_raw", 100, 0, 100, 800, 0, 8192, det.stripB + 4 * id, det.backE);
      }
#endif
    }
  }
  // return kTRUE;
  //  QQQ Processing

  int qqqCount = 0;
  // REMOVE WHEN RERUNNING USING THE NEW CALIBRATION FILE
  std::vector<Event> QQQ_Events, PC_Events;
  // std::vector<Event> QQQ_Events_Raw, PC_Events_Raw;
  // std::vector<Event> QQQ_Events2; // clustering done

  // Check for muliplt hits in the same QQQ channel

  // -------Check Start ------

  std::unordered_map<int, std::tuple<int, int, double, double>> qvecr[4], qvecw[4];
  if (qqq.multi > 1)
  {
    // if(qqq.multi>=3) std::cout << "-----" << std::endl;
    for (int i = 0; i < qqq.multi; i++)
    {
      if (qqq.ch[i] / 16)
      {
        if (qvecr[qqq.id[i]].find(qqq.ch[i]) != qvecr[qqq.id[i]].end())
          std::cout << "mayday!" << std::endl;
        qvecr[qqq.id[i]][qqq.ch[i]] = std::tuple(qqq.id[i], qqq.ch[i], qqq.e[i], qqq.t[i]);
      }
      else
      {
        if (qvecw[qqq.id[i]].find(qqq.ch[i]) != qvecw[qqq.id[i]].end())
          std::cout << "mayday!" << std::endl;
        qvecw[qqq.id[i]][qqq.ch[i]] = std::tuple(qqq.id[i], qqq.ch[i], qqq.e[i], qqq.t[i]);
      }
    }
  }

  // -------Check End ------

  bool PCAQQQTimeCut = false;
  bool PCCQQQTimeCut = false;
  for (int i = 0; i < qqq.multi; i++)
  {
#ifdef RAW_HISTOS
    plotter->Fill2D("QQQ_Index_Vs_Energy", 16 * 8, 0, 16 * 8, 2000, 0, 8000, qqq.index[i], qqq.e[i], "hRawQQQ");

    for (int j = 0; j < qqq.multi; j++)
    {
      if (j == i)
        continue;
      plotter->Fill2D("QQQ_Coincidence_Matrix", 16 * 8, 0, 16 * 8, 16 * 8, 0, 16 * 8, qqq.index[i], qqq.index[j], "hRawQQQ");
    }

    for (int k = 0; k < pc.multi; k++)
    {
      if (pc.index[k] < 24 && pc.e[k] > 10)
      {
        plotter->Fill2D("QQQ_Vs_Anode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
        plotter->Fill2D("QQQ_Vs_PC_Index", 16 * 8, 0, 16 * 8, 24, 0, 24, qqq.index[i], pc.index[k], "hRawQQQ");
      }
      else if (pc.index[k] >= 24 && pc.e[k] > 10)
      {
        plotter->Fill2D("QQQ_Vs_Cathode_Energy", 400, 0, 4000, 1000, 0, 16000, qqq.e[i], pc.e[k], "hRawQQQ");
      }
    }
#endif
    for (int j = i + 1; j < qqq.multi; j++)
    {
      if (qqq.id[i] == qqq.id[j])
      {
        qqqCount++;

        int chWedge = -1;
        int chRing = -1;
        double eWedge = 0.0;
        double eWedgeMeV = 0.0;
        double eRing = 0.0;
        double eRingMeV = 0.0;
        double tRing = 0.0;
        double tWedge = 0.0;

        if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
        {
          chWedge = qqq.ch[i];
          eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
          chRing = qqq.ch[j] - 16;
          eRing = qqq.e[j];
          tRing = static_cast<double>(qqq.t[j]);
          tWedge = static_cast<double>(qqq.t[i]);
        }
        else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
        {
          chWedge = qqq.ch[j];
          eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
          chRing = qqq.ch[i] - 16;
          eRing = qqq.e[i];
          tRing = static_cast<double>(qqq.t[i]);
          tWedge = static_cast<double>(qqq.t[j]);
        }
        else
          continue;

        plotter->Fill1D("Wedgetime_Vs_Ringtime", 100, -1000, 1000, tWedge - tRing, "hTiming");
#ifdef RAW_HISTOS
        plotter->Fill2D("RingE_vs_Index", 16 * 4, 0, 16 * 4, 1000, 0, 16000, chRing + qqq.id[i] * 16, eRing, "hRawQQQ");
        plotter->Fill2D("WedgeE_vs_Index", 16 * 4, 0, 16 * 4, 1000, 0, 16000, chWedge + qqq.id[i] * 16, eWedge, "hRawQQQ");
#endif

        if (qqqCalibValid[qqq.id[i]][chWedge][chRing])
        {
          eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;
          eRingMeV = eRing * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;

          if (eRingMeV / eWedgeMeV > 3.0 || eRingMeV / eWedgeMeV < 1.0 / 3.0)
            continue;
          // if(eRingMeV<1.2 || eWedgeMeV<1.2) continue;

          // double theta = 2 * TMath::Pi() * (-qqq.id[i] * 16 + (15 - chWedge) + 0.5) / (16 * 4);

          double phi_qqq = (M_PI / 180.) * (-90 * qqq.id[i] + (87. / 16.) * ((15 - chWedge) + 0.5) + 3.0);
          double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
          // z used to be 75+30+23=128
          // we found a 12mm shift towards the vertex later --> 116
          Event qqqevent(TVector3(rho * TMath::Cos(phi_qqq), rho * TMath::Sin(phi_qqq), qqq_z), eRingMeV, eWedgeMeV, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);
          // Event qqqeventr(TVector3(rho * TMath::Cos(theta), rho * TMath::Sin(theta), qqq_z), eRing, eWedge, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);

          QQQ_Events.push_back(qqqevent);
          // QQQ_Events_Raw.push_back(qqqeventr);
          plotter->Fill2D("WedgeE_Vs_RingECal_selected", 1000, 0, 10, 1000, 0, 10, eWedgeMeV, eRingMeV, "hCalQQQ");

          plotter->Fill1D("QQQECal", 2048, 0, 10, eRingMeV);
          plotter->Fill1D("QQQECal", 2048, 0, 10, eWedgeMeV);

          const int channelsPerDetector = MAX_RING + MAX_WEDGE;
          int globalRingChannel = chRing + (qqq.id[i] * channelsPerDetector);
          int globalWedgeChannel = chWedge + (qqq.id[i] * channelsPerDetector) + MAX_RING;

          // Fill the histograms
          plotter->Fill2D("QQQ_CalibratedE_vs_Ch", 128, 0, 128, 1000, 0, 20, globalRingChannel, eRingMeV, "hCalQQQ");
          plotter->Fill2D("QQQ_CalibratedE_vs_Ch", 128, 0, 128, 1000, 0, 20, globalWedgeChannel, eWedgeMeV, "hCalQQQ");

          plotter->Fill2D("QQQCartesianPlot", 200, -100, 100, 200, -100, 100, rho * TMath::Cos(phi_qqq), rho * TMath::Sin(phi_qqq), "hCalQQQ");
          plotter->Fill2D("QQQCartesianPlot" + std::to_string(qqq.id[i]), 200, -100, 100, 200, -100, 100, rho * TMath::Cos(phi_qqq), rho * TMath::Sin(phi_qqq), "hCalQQQ");
          plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, rho * TMath::Cos(phi_qqq), rho * TMath::Sin(phi_qqq), "hPCQQQ");
        }
        else
          continue;

        for (int k = 0; k < pc.multi; k++)
        {
#ifdef RAW_HISTOS
          plotter->Fill2D("RingCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Anode_Index", 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Anode_Index" + std::to_string(qqq.id[i]), 16 * 4, 0, 16 * 4, 24, 0, 24, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("RingCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chRing + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
          plotter->Fill2D("WedgeCh_vs_Cathode_Index", 16 * 4, 0, 16 * 4, 24, 24, 48, chWedge + qqq.id[i] * 16, pc.index[k], "hRawQQQ");
#endif
          if (pc.index[k] < 24 && pc.e[k] > 10)
          {
            plotter->Fill2D("Timing_Difference_QQQ_PC", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
            plotter->Fill2D("DelT_Vs_QQQRingECal", 500, -2000, 2000, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hTiming");
            // if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
            if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
            {
              PCAQQQTimeCut = true;
              plotter->Fill2D("CalibratedQQQEvsPCE_R", 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
              plotter->Fill2D("CalibratedQQQEvsPCE_W", 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
            }
          }

          if (pc.index[k] >= 24 && pc.e[k] > 10)
          {
            if (tRing - static_cast<double>(pc.t[k]) < -200)
              PCCQQQTimeCut = true;
            // if (tRing - static_cast<double>(pc.t[k]) > 200) PCCQQQTimeCut = true;
            plotter->Fill2D("Timing_Difference_QQQ_PC_Cathode", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
          }
        } // end of pc k loop

        if (!HitNonZero)
        {
          // double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
          // double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
          double phi_qqq = (2 * M_PI) * (-90 * qqq.id[i] + (87. / 16.) * ((15 - chWedge) + 0.5) + 3.0);
          double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
          double x = rho * TMath::Cos(phi_qqq);
          double y = rho * TMath::Sin(phi_qqq);
          hitPos.SetXYZ(x, y, qqq_z);
          qqqenergy = eRingMeV;
          qqqtimestamp = tRing;
          HitNonZero = true;
        }
      } // if j==i
    } // j loop end
  } // i loop end

  PCQQQTimeCut = PCAQQQTimeCut && PCCQQQTimeCut;
#ifdef RAW_HISTOS
  plotter->Fill1D("QQQ_Multiplicity", 10, 0, 10, qqqCount, "hRawQQQ");
#endif
  aWireEvents.clear();
  aWireEvents.reserve(24);
  cWireEvents.clear();
  cWireEvents.reserve(24);
  // PC Gain Matching and Filling
  for (int i = 0; i < pc.multi; i++)
  {
    // std::cout << pc.index[i] << " " << pc.e[i] << " " << std::endl;
#ifdef RAW_HISTOS
    if (pc.e[i] > 50)
    {
      plotter->Fill2D("PC_Index_Vs_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], static_cast<double>(pc.e[i]), "hRawPC");
    }
#endif

    pc.e[i] = pcSlope[pc.index[i]] * pc.e[i] + pcIntercept[pc.index[i]];
    if (pc.e[i] > 50)
    {
      if (pc.index[i] >= 24)
        plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], pc.e[i] * cathode_gain, "hGMPC");
      else
        plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], pc.e[i], "hGMPC");
    }
    if (pc.e[i] > 50)
    {
      if (pc.index[i] < 24)
      {
        anodeT = static_cast<double>(pc.t[i]);
        anodeIndex = pc.index[i];
        aWireEvents[pc.index[i]] = std::tuple(pc.index[i], pc.e[i], static_cast<double>(pc.t[i]));
      }
      else
      {
        cathodeT = static_cast<double>(pc.t[i]);
        cathodeIndex = pc.index[i] - 24;
        // cWireEvents[pc.index[i] - 24] = std::tuple(pc.index[i] - 24, pc.e[i], static_cast<double>(pc.t[i]));
        cWireEvents[pc.index[i] - 24] = std::tuple(pc.index[i] - 24, pc.e[i] * cathode_gain, static_cast<double>(pc.t[i]));
      }
    }

    if (anodeT != -99999 && cathodeT != 99999)
    {
      for (int j = 0; j < qqq.multi; j++)
      {
        plotter->Fill1D("PC_Time_qqq", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
        plotter->Fill2D("PC_Time_Vs_QQQ_ch", 200, -2000, 2000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, qqq.ch[j], "hTiming");
        plotter->Fill2D("PC_Time_vs_AIndex_qqq", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
        plotter->Fill2D("PC_Time_vs_CIndex_qqq", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
        // plotter->Fill1D("PC_Time_A" + std::to_string(anodeIndex) + "_C" + std::to_string(cathodeIndex), 200, -1000, 1000, anodeT - cathodeT, "TimingPC");
      }

      for (int j = 0; j < sx3.multi; j++)
      {
        plotter->Fill1D("PC_Time_sx3", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
        // plotter->Fill2D("PC_Time_Vs_SX3_ch", 200, -2000, 2000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, sx3.ch[j], "hTiming");
        plotter->Fill2D("PC_Time_vs_AIndex_sx3", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
        plotter->Fill2D("PC_Time_vs_CIndex_sx3", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
      }
      for (const auto &sx3event : SX3_Events)
      {
        bool TCC = sx3event.Time1 - cathodeT < 0;
        bool TCA = sx3event.Time1 - anodeT < 0;
        // plotter->Fill2D("sx3_z_phi_awire"+std::to_string(anodeIndex)+"_TC"+std::to_string(TCA), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
        // plotter->Fill2D("sx3_z_phi_cwire"+std::to_string(cathodeIndex)+"_TC"+std::to_string(TCC), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
      }

      plotter->Fill1D("PC_Time", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
    }

    for (int j = i + 1; j < pc.multi; j++)
    {
#ifdef RAW_HISTOS
      plotter->Fill2D("PC_Coincidence_Matrix", 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
      plotter->Fill2D("PC_Coincidence_Matrix_anodeMinusCathode_lt_-200_" + std::to_string(anodeT - cathodeT < -200), 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
#endif
      plotter->Fill2D("Anode_V_Anode", 24, 0, 24, 24, 0, 24, pc.index[i], pc.index[j], "hGMPC");
    }
  }
  anodeHits.clear();
  cathodeHits.clear();
  corrcatMax.clear();

  for (int i = 0; i < pc.multi; i++)
  {
    // if (pc.e[i] > 100)
    {
      if (pc.index[i] < 24)
      {
        anodeHits.push_back(std::pair<int, double>(pc.index[i], pc.e[i]));
      }
      else if (pc.index[i] >= 24)
      {
        cathodeHits.push_back(std::pair<int, double>(pc.index[i] - 24, pc.e[i] * cathode_gain));
      }
    }
  }

  std::sort(anodeHits.begin(), anodeHits.end(), [](std::pair<int, double> a, std::pair<int, double> b)
            { return a.first < b.first; });

  std::sort(cathodeHits.begin(), cathodeHits.end(), [](std::pair<int, double> a, std::pair<int, double> b)
            { return a.first < b.first; });

  // clusters = collection of (collection of wires) where each wire is (index, energy, timestamp)
  std::vector<std::vector<std::tuple<int, double, double>>> aClusters = pwinstance.Make_Clusters(aWireEvents);
  std::vector<std::vector<std::tuple<int, double, double>>> cClusters = pwinstance.Make_Clusters(cWireEvents);

  for (const auto &aCluster : aClusters)
  {
    for (const auto &cCluster : cClusters)
    {
      if (aCluster.size() == 0)
        continue;
      if (cCluster.size() == 0)
        continue;
      // both have at least 1, here. Keep the a1, c1 events
      auto [crossover, alpha, apSumE, cpSumE, apMaxE, cpMaxE, apTSMaxE, cpTSMaxE] = pwinstance.FindCrossoverProperties(aCluster, cCluster);
      if (alpha != 9999999 && apSumE != -1)
      {
        // Event PCEvent(crossover,apMaxE,cpMaxE,apTSMaxE,cpTSMaxE);
        // Event PCEvent(crossover,apSumE,cpSumE,apTSMaxE,cpTSMaxE);
        Event PCEvent(crossover, apSumE, cpMaxE, apTSMaxE, cpTSMaxE); // run12 shows cathode-max and anode-sum provide best dE signals.
        // std::cout << apSumE << " " << crossover.Perp() << " " << apMaxE << " " << apTSMaxE << std::endl;
        PCEvent.multi1 = aCluster.size();
        PCEvent.multi2 = cCluster.size();
        PCEvent.Anodech = std::get<0>(aCluster[0]);
        PCEvent.Cathodech = std::get<0>(cCluster[0]);
        PC_Events.push_back(PCEvent);
      }
      else
      {
        ; // std::cout << "AAAA " << std::endl;
      }
    }
  }

  //////Timing stuff for F data

  static TRandom3 rnd(0); // seeded once (random seed via TUUID), not per event
  if (dataset == "17F" && reactiondata)
  {
    int ctr = 0;
    for (const auto &qqqevent : QQQ_Events)
    {
      double ts_rf = -987654321;
      double ts_needle = -987654321;
      double ts_mcp = -987654321;
      double ts_qqq = static_cast<double>(qqqevent.Time1) + (rnd.Uniform(16.0) - 8.0);
      bool found_rf = false;
      bool found_mcp = false;
      bool found_needle = false;
      bool qqq_inner_ring = (qqqevent.ch1 % 16) < 8;
      for (int j = 0; j < misc.multi; j++)
      {
        plotter->Fill1D("channels_misc_qqq", 20, 0, 20, misc.ch[j], "misc");
        if (misc.ch[j] == 2)
        { // Needle
          plotter->Fill2D("needle_vs_qqqE", 800, 0, 16384, 800, 0, 10, misc.e[j], qqqevent.Energy1, "misc");
          ts_needle = static_cast<double>(misc.t[j]) + static_cast<double>(misc.tf[j]);
          found_needle = 1;
          plotter->Fill1D("dt_qqq_needle", 800, -2000, 2000, ts_qqq - ts_needle, "misc");
        }
        if (misc.ch[j] == 3)
        { // RF
          ts_rf = static_cast<double>(misc.t[j]) + static_cast<double>(misc.tf[j]);
          found_rf = 1;
          plotter->Fill1D("dt_qqq_rf_innerring" + std::to_string(qqq_inner_ring), 800, -2000, 2000, ts_qqq - ts_rf, "misc");
        }
        if (misc.ch[j] == 4)
        { // mcp
          ts_mcp = static_cast<double>(misc.t[j]) + static_cast<double>(misc.tf[j]);
          found_mcp = 1;
          plotter->Fill1D("dt_qqq_mcp_innerring" + std::to_string(qqq_inner_ring), 800, -2000, 2000, ts_qqq - ts_mcp, "misc");
        }
      }
      if (found_rf && found_mcp)
      {
        if (ctr == 0)
          plotter->Fill1D("dt_rf_mcp_qqq_innerring" + std::to_string(qqq_inner_ring), 500, -1000, 1000, ts_rf - ts_mcp, "misc");
        double dt_rf_mcp = ts_rf - ts_mcp;
        double dt_qqq_rf = ts_qqq - ts_rf;
        double dt_qqq_mcp = ts_qqq - ts_mcp;
        plotter->Fill2D("dt(qqq,rf)_vs_(rf,mcp)_innerring" + std::to_string(qqq_inner_ring), 640, -2000, 2000, 640, -2000, 2000, dt_qqq_rf, dt_rf_mcp, "misc");
        plotter->Fill2D("dt_(qqq,mcp)_vs_(qqq,rf)_innerring" + std::to_string(qqq_inner_ring), 640, -1400, 2000, 640, -2000, 2000, dt_qqq_mcp, dt_qqq_rf, "misc");
        plotter->Fill2D("dt_(qqq,mcp)_vs_(rf,mcp)_innerring" + std::to_string(qqq_inner_ring), 640, -1400, -600, 640, -2000, 2000, dt_qqq_mcp, dt_rf_mcp, "misc");
      }
      ctr += 1;
    }

    for (const auto &sx3event : SX3_Events)
    {
      double ts_rf = -987654321;
      double ts_needle = -987654321;
      double ts_mcp = -987654321;
      double ts_sx3 = static_cast<double>(sx3event.Time1) + (rnd.Uniform(16.0) - 8.0);
      bool found_rf = false;
      bool found_mcp = false;
      bool found_needle = false;
      for (int j = 0; j < misc.multi; j++)
      {
        plotter->Fill1D("channels_misc_sx3", 20, 0, 20, misc.ch[j], "misc");
        if (misc.ch[j] == 2)
        { // Needle
          plotter->Fill2D("needle_vs_sx3E", 800, 0, 16384, 800, 0, 10, misc.e[j], sx3event.Energy1, "misc");
          ts_needle = static_cast<double>(misc.t[j]) + static_cast<double>(misc.tf[j]);
          found_needle = 1;
          plotter->Fill1D("dt_sx3_needle", 800, -2000, 2000, ts_sx3 - ts_needle, "misc");
        }
        if (misc.ch[j] == 3)
        { // RF
          ts_rf = static_cast<double>(misc.t[j]) + static_cast<double>(misc.tf[j]);
          found_rf = 1;
          plotter->Fill1D("dt_sx3_rf", 800, -2000, 2000, ts_sx3 - ts_rf, "misc");
        }
        if (misc.ch[j] == 4)
        { // mcp
          ts_mcp = static_cast<double>(misc.t[j]) + static_cast<double>(misc.tf[j]);
          found_mcp = 1;
          plotter->Fill1D("dt_sx3_mcp", 800, -2000, 2000, ts_sx3 - ts_mcp, "misc");
        }
      }
      if (found_rf && found_mcp)
      {
        if (ctr == 0)
          plotter->Fill1D("dt_rf_mcp_sx3", 500, -1000, 1000, ts_rf - ts_mcp, "misc");
        double dt_rf_mcp = ts_rf - ts_mcp;
        double dt_sx3_rf = ts_sx3 - ts_rf;
        double dt_sx3_mcp = ts_sx3 - ts_mcp;
        plotter->Fill2D("dt(sx3,rf)_vs_(rf,mcp)", 640, -2000, 2000, 640, -2000, 2000, dt_sx3_rf, dt_rf_mcp, "misc");
        plotter->Fill2D("dt_(sx3,mcp)_vs_(sx3,rf)", 640, -1400, 2000, 640, -2000, 2000, dt_sx3_mcp, dt_sx3_rf, "misc");
        plotter->Fill2D("dt_(sx3,mcp)_vs_(rf,mcp)", 640, -1400, -600, 640, -2000, 2000, dt_sx3_mcp, dt_rf_mcp, "misc");
      }
      ctr += 1;
    }
  }

  if (process_alpha_proton_scattering)
  {
    protonAlphaHistograms(plotter, QQQ_Events, SX3_Events, PC_Events);
    // return kTRUE;
  } // end if(process_alpha_proton_scattering)

  if (doMiscHistograms)
  {
    miscHistograms_oneWire(plotter, QQQ_Events, aClusters);
    protonMiscHistograms_sx3(plotter, QQQ_Events, SX3_Events, PC_Events);
    protonMiscHistograms(plotter, QQQ_Events, SX3_Events, PC_Events);
    // return kTRUE;
  }

#ifdef RAW_HISTOS
  if (QQQ_Events.size() && PC_Events.size())
    plotter->Fill2D("PCEv_vs_QQQEv", 20, 0, 20, 20, 0, 20, QQQ_Events.size(), PC_Events.size());

  plotter->Fill2D("ac_vs_cc", 20, 0, 20, 20, 0, 20, aClusters.size(), cClusters.size(), "wiremult");
  for (const auto &cluster : aClusters)
  {
    plotter->Fill1D("aClusters" + std::to_string(aClusters.size()), 20, -5, 15, cluster.size(), "wiremult");
  }
  for (const auto &cluster : cClusters)
  {
    plotter->Fill1D("cClusters" + std::to_string(cClusters.size()), 20, -5, 15, cluster.size(), "wiremult");
  }

  if (cClusters.size() && aClusters.size())
  {
    plotter->Fill2D("ac_vs_cc_ign0", 20, 0, 20, 20, 0, 20, aClusters.size(), cClusters.size(), "wiremult");
  }
#endif

  for (const auto &sx3event : SX3_Events)
  {
    for (int i = 0; i < 24; i++)
    {
      if (aWireEvents.find(i) != aWireEvents.end())
      {
        auto awire = aWireEvents[i];
        if (sx3event.Time1 - (double)std::get<2>(awire) < 150)
        {
          // plotter->Fill2D("sx3_z_phi2_awire"+std::to_string(std::get<0>(awire)), 400,-100,100, 100, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
          // plotter->Fill2D("sx3_z_strip#_awire"+std::to_string(std::get<0>(awire)), 400,-100,100, 100, -50,50,sx3event.pos.Z(), sx3event.ch2);
          plotter->Fill2D("onewire_dEa_Esx3_TC1_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, sx3event.Energy1, std::get<1>(awire), "1wire");
          plotter->Fill2D("onewire_aNum_sx3Phi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -200,200, i, sx3event.pos.Phi() * 180. / M_PI, "1wire");
        }
      }

      if (cWireEvents.find(i) != cWireEvents.end())
      {
        auto cwire = cWireEvents[i];
        if (sx3event.Time1 - (double)std::get<2>(cwire) < 150)
        {
          // plotter->Fill2D("sx3_z_phi2_cwire"+std::to_string(std::get<0>(cwire)),400,-100,100, 100, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
          // plotter->Fill2D("sx3_z_strip#_cwire"+std::to_string(std::get<0>(cwire)),400,-100,100, 100, -50,50,sx3event.pos.Z(), sx3event.ch2 );
          plotter->Fill2D("onewire_dEc_Esx3_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, sx3event.Energy1, std::get<1>(cwire), "1wire");
          plotter->Fill2D("onewire_cNum_sx3Phi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -200,200, i, sx3event.pos.Phi() * 180. / M_PI, "1wire");
        }
      }
    } // for 'i' loop

    for (const auto &acluster : aClusters)
    {
      auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster, "ANODE");
      int a_number = acluster.size();
      TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire, sx3event.pos.Phi());
      if (sx3event.Time1 - apTSMaxE < 150)
      {
        plotter->Fill2D("dEa_interp_Esx3_TC1_ignC" + std::to_string(acluster.size()), 400, 0, 10, 800, 0, 40000, sx3event.Energy1, apSumE, "ainterp_noc");
        plotter->Fill2D("aPhi_interp_sx3Phi_TC1_ignC" + std::to_string(acluster.size()), 120, -200,200, 120, -200,200, pc_closest.Phi() * 180. / M_PI, sx3event.pos.Phi() * 180. / M_PI, "ainterp_noc");
        plotter->Fill2D("aZ_interp_sx3Z_TC1_ignC" + std::to_string(acluster.size()), 400, -200, 200, 300, -100, 200, pc_closest.Z(), sx3event.pos.Z(), "ainterp_noc");
        TVector3 x2(pc_closest), x1(sx3event.pos);
        TVector3 v = x2 - x1;
        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
        TVector3 vector_closest_to_axis = x1 + t_minimum * v;
        plotter->Fill2D("vertexZ_interp_sx3Z_TC1_ignC" + std::to_string(acluster.size()), 400, -200, 200, 300, -100, 200, vector_closest_to_axis.Z(), sx3event.pos.Z(), "ainterp_noc");
        plotter->Fill2D("vertexXY_interp_TC1_ignC" + std::to_string(acluster.size()), 200, -100, 100, 200, -100, 100, vector_closest_to_axis.X(), vector_closest_to_axis.Y(), "ainterp_noc");
      }
    }
  }

  for (const auto &qqqevent : QQQ_Events)
  {
    for (int i = 0; i < 24; i++)
    {
      if (aWireEvents.find(i) != aWireEvents.end())
      {
        const auto &awire = aWireEvents[i];
        if (qqqevent.Time1 - (double)std::get<2>(awire) < 150)
        {
          plotter->Fill2D("onewire_dEa_Eqqq_TC1_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, std::get<1>(awire), "1wire");
          plotter->Fill2D("onewire_aNum_QQQPhi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -200,200, i, qqqevent.pos.Phi() * 180. / M_PI, "1wire");
        }
      }

      if (cWireEvents.find(i) != cWireEvents.end())
      {
        auto cwire = cWireEvents[i];
        if (qqqevent.Time1 - (double)std::get<2>(cwire) < 150)
        {
          plotter->Fill2D("onewire_dEc_Eqqq_TC1_fullev" + std::to_string(PC_Events.size() > 0), 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, std::get<1>(cwire), "1wire");
          plotter->Fill2D("onewire_cNum_QQQPhi_TC1_fullev" + std::to_string(PC_Events.size() > 0), 24, 0, 24, 120, -200,200, i, qqqevent.pos.Phi() * 180. / M_PI, "1wire");
        }
      }
    } // for 'i' loop
  }

  if (doPCSX3ClusterAnalysis)
  {
    PCSX3ClusterAnalysis(plotter, QQQ_Events, SX3_Events, PC_Events, aClusters, cClusters);
  }
  if (doPCQQQClusterAnalysis)
  {
    PCQQQClusterAnalysis(plotter, QQQ_Events, SX3_Events, PC_Events, aClusters, cClusters);
  }

///////////////////nA analysis using pseudo-wire (GetPseudoWire + getClosestWirePosAtWirePhi)///////////////////
#ifdef nA_Analysis
  if (aClusters.size() > 0)
  {
    // ---------------------------------------------------------
    // PROTON LOOP (SX3 BARREL)
    // ---------------------------------------------------------
    for (const auto &sx3event : SX3_Events)
    {
      // Pick the anode cluster closest in phi to this SX3 hit
      const std::vector<std::tuple<int, double, double>> *bestCluster = &aClusters[0];
      double bestDphi = 9999.0;

      for (const auto &acluster : aClusters)
      {
        auto [pw, sumE, maxE, tsMax] = pwinstance.GetPseudoWire(acluster, "ANODE");
        TVector3 pos = pwinstance.getClosestWirePosAtWirePhi(pw, sx3event.pos.Phi());
        double dphi = TMath::Abs(TVector2::Phi_mpi_pi(sx3event.pos.Phi() - pos.Phi()));
        if (dphi < bestDphi)
        {
          bestDphi = dphi;
          bestCluster = &acluster;
        }
      }

      // Extract the virtual wire specifically for the best cluster
      auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(*bestCluster, "ANODE");
      std::string nA_label = std::to_string(bestCluster->size()) + "A";

      TVector3 pcz_intersect = pwinstance.getClosestWirePosAtWirePhi(apwire, sx3event.pos.Phi());

      double deltaRho = sx3event.pos.Perp() - pcz_intersect.Perp();
      double deltaZ = sx3event.pos.Z() - pcz_intersect.Z();
      double vertex_recon = sx3event.pos.Z() - sx3event.pos.Perp() * (deltaZ / deltaRho);

      std::string vtx_gate = "";

#ifdef VTX_GATE
      if (vertex_recon >= -176.0 && vertex_recon < -100.0)
        vtx_gate = "_Z[-176_to_-100]";
      else if (vertex_recon >= -100.0 && vertex_recon < -50.0)
        vtx_gate = "_Z[-100_to_-50]";
      else if (vertex_recon >= -50.0 && vertex_recon < 0.0)
        vtx_gate = "_Z[-50_to_0]";
      else if (vertex_recon >= 0.0 && vertex_recon < 50.0)
        vtx_gate = "_Z[0_to_50]";
      else if (vertex_recon >= 50.0 && vertex_recon < 100.0)
        vtx_gate = "_Z[50_to_100]";
      else if (vertex_recon >= 100.0 && vertex_recon < 176.0)
        vtx_gate = "_Z[100_to_176]";
#endif

      // *0.1 converts mm to cm for the Eloss
      double beam_path_length = TMath::Abs(vertex_recon - z_entrance) * 0.1;

#ifdef AL_BEAM
      double beam_energy_at_vertex = cm_to_MeV_27Al->Eval(MeV_to_cm_27Al->Eval(65.195) - beam_path_length); // 72MeV is the energy of the 27Al beam right ebfore teh chamber
      Kinematics apkin_p(mass_27Al, mass_4He, mass_1H, mass_30Si, beam_energy_at_vertex / mass_27Al);
      Kinematics apkin_a(mass_27Al, mass_4He, mass_4He, mass_30Si, beam_energy_at_vertex / mass_27Al);
#endif

#ifdef F_BEAM
      double beam_energy_at_vertex = cm_to_MeV_17F->Eval(MeV_to_cm_17F->Eval(60.741) - beam_path_length); // 67.8 MeV iws the energy of the 17F beam right before the chamber
      Kinematics apkin_p(mass_17F, mass_4He, mass_1H, mass_20Ne, beam_energy_at_vertex / mass_17F);
      Kinematics apkin_a(mass_17F, mass_4He, mass_4He, mass_20Ne, beam_energy_at_vertex / mass_17F);
#endif

      double path_length = (sx3event.pos - TVector3(0, 0, vertex_recon)).Mag() * 0.1;
      double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1) - path_length);
      double sx3Efixalpha = cm_to_MeV->Eval(MeV_to_cm->Eval(sx3event.Energy1) - path_length);

      double theta_recon = (sx3event.pos - TVector3(0, 0, vertex_recon)).Theta();
      double sinTheta = TMath::Sin(theta_recon);

      // Fill standard un-gated plots
      plotter->Fill1D(nA_label + "_vertex_recon_SX3", 400, -200, 200, vertex_recon, nA_label);
      plotter->Fill1D(nA_label + "_vertex_recon", 400, -200, 200, vertex_recon, nA_label);

      plotter->Fill2D(nA_label + "_dE_Ecorr_Anode_SX3", 800, 0, 30, 800, 0, 30000, sx3Efix, apSumE * sinTheta, nA_label);
      plotter->Fill2D(nA_label + "_dE_Ecorr_Anode_SX3_alpha", 800, 0, 30, 800, 0, 30000, sx3Efixalpha, apSumE * sinTheta, nA_label);
      plotter->Fill2D(nA_label + "_sx3_E_vs_theta_raw_SX3", 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, sx3event.Energy1, nA_label);
      plotter->Fill2D(nA_label + "_sx3_E_vs_theta_corr_SX3", 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, sx3Efix, nA_label);

#if defined(F_BEAM) || defined(AL_BEAM)
      double Ex_from_proton = apkin_p.getExc(sx3Efix, theta_recon * 180. / M_PI);
      double Ex_from_alpha = apkin_a.getExc(sx3Efixalpha, theta_recon * 180. / M_PI);
      plotter->Fill1D(nA_label + "_Ex_from_protons_SX3", 1200, -30, 30, Ex_from_proton, nA_label);
      plotter->Fill1D(nA_label + "_Ex_from_alphas_SX3", 1200, -30, 30, Ex_from_alpha, nA_label);
#endif

      // Fill Gated Plots
      // if (vtx_gate != "")
      // {
      //   plotter->Fill2D("dE_Ecorr_Anode_SX3" + vtx_gate, 400, 0, 30, 800, 0, 40000, sx3Efix, apSumE * sinTheta, nA_label);
      //   plotter->Fill2D("dE_Ecorr_Anode_SX3_alpha" + vtx_gate, 400, 0, 30, 800, 0, 40000, sx3Efixalpha, apSumE * sinTheta, nA_label);
      //   plotter->Fill1D("twisted_pcz_recon_SX3" + vtx_gate, 600, -300, 300, pcz_intersect.Z(), nA_label);
      //   plotter->Fill1D("twisted_vertex_recon_SX3" + vtx_gate, 600, -300, 300, vertex_recon, nA_label);
      //   plotter->Fill1D("Ex_from_protons_SX3" + vtx_gate, 1200, -30, 30, Ex_from_proton, nA_label);
      //   plotter->Fill1D("Ex_from_alphas_SX3" + vtx_gate, 1200, -30, 30, Ex_from_alpha, nA_label);
      // }
    }

    // ---------------------------------------------------------
    // PROTON LOOP (QQQ ENDCAP)
    // ---------------------------------------------------------
    for (const auto &qqqevent : QQQ_Events)
    {
      const std::vector<std::tuple<int, double, double>> *bestCluster = nullptr;
      double bestDphi = 9999.0;

      for (const auto &acluster : aClusters)
      {
        auto [apw, sumE, maxE, tsMax] = pwinstance.GetPseudoWire(acluster, "ANODE");
        TVector3 pcPos = pwinstance.getClosestWirePosAtWirePhi(apw, qqqevent.pos.Phi());
        double dphi = TMath::Abs(TVector2::Phi_mpi_pi(qqqevent.pos.Phi() - pcPos.Phi()));
        if (dphi < bestDphi)
        {
          bestDphi = dphi;
          bestCluster = &acluster;
        }
      }
      if (!bestCluster)
        continue;

      // Extract the virtual wire specifically for the best cluster
      auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(*bestCluster, "ANODE");
      std::string nA_label = std::to_string(bestCluster->size()) + "A";

      TVector3 pcz_intersect = pwinstance.getClosestWirePosAtWirePhi(apwire, qqqevent.pos.Phi());

      double deltaRho = qqqevent.pos.Perp() - pcz_intersect.Perp();
      double deltaZ = qqqevent.pos.Z() - pcz_intersect.Z();
      double vertex_recon = qqqevent.pos.Z() - qqqevent.pos.Perp() * (deltaZ / deltaRho);

      std::string vtx_gate = "";

#ifdef VTX_GATE
      if (vertex_recon >= -176.0 && vertex_recon < -100.0)
        vtx_gate = "_Z[-176_to_-100]";
      else if (vertex_recon >= -100.0 && vertex_recon < -50.0)
        vtx_gate = "_Z[-100_to_-50]";
      else if (vertex_recon >= -50.0 && vertex_recon < 0.0)
        vtx_gate = "_Z[-50_to_0]";
      else if (vertex_recon >= 0.0 && vertex_recon < 50.0)
        vtx_gate = "_Z[0_to_50]";
      else if (vertex_recon >= 50.0 && vertex_recon < 100.0)
        vtx_gate = "_Z[50_to_100]";
      else if (vertex_recon >= 100.0 && vertex_recon < 176.0)
        vtx_gate = "_Z[100_to_176]";
#endif

      double beam_path_length = TMath::Abs(vertex_recon - z_entrance) * 0.1;

#ifdef AL_BEAM
      double beam_energy_at_vertex = cm_to_MeV_27Al->Eval(MeV_to_cm_27Al->Eval(72.0) - beam_path_length);
      Kinematics apkin_p(mass_27Al, mass_4He, mass_1H, mass_30Si, beam_energy_at_vertex / mass_27Al);
      Kinematics apkin_a(mass_27Al, mass_4He, mass_4He, mass_30Si, beam_energy_at_vertex / mass_27Al);
#endif

#ifdef F_BEAM
      double beam_energy_at_vertex = cm_to_MeV_17F->Eval(MeV_to_cm_17F->Eval(65.0) - beam_path_length);
      Kinematics apkin_p(mass_17F, mass_4He, mass_1H, mass_20Ne, beam_energy_at_vertex / mass_17F);
      Kinematics apkin_a(mass_17F, mass_4He, mass_4He, mass_20Ne, beam_energy_at_vertex / mass_17F);
#endif

      // Energy Loss Correction
      double path_length = (qqqevent.pos - TVector3(0, 0, vertex_recon)).Mag() * 0.1;
      double qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1) - path_length);
      double qqqEfixalpha = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length);

      double theta_recon = (qqqevent.pos - TVector3(0, 0, vertex_recon)).Theta();
      double sinTheta = TMath::Sin(theta_recon);

      // Fill standard un-gated plots
      plotter->Fill1D(nA_label + "_vertex_recon_QQQ", 400, -200, 200, vertex_recon, nA_label);
      plotter->Fill1D(nA_label + "_vertex_recon", 400, -200, 200, vertex_recon, nA_label);
      plotter->Fill2D(nA_label + "_dE_Ecorr_Anode_QQQ", 800, 0, 30, 800, 0, 30000, qqqEfix, apSumE * sinTheta, nA_label);
      plotter->Fill2D(nA_label + "_dE_Ecorr_Anode_QQQ_alpha", 800, 0, 30, 800, 0, 30000, qqqEfixalpha, apSumE * sinTheta, nA_label);
      plotter->Fill2D(nA_label + "_qqq_E_vs_theta_raw_QQQ", 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, qqqevent.Energy1, nA_label);
      plotter->Fill2D(nA_label + "_qqq_E_vs_theta_corr_QQQ", 180, 0, 180, 400, 0, 30, theta_recon * 180. / M_PI, qqqEfix, nA_label);

#if defined(F_BEAM) || defined(AL_BEAM)
      double Ex_from_alpha = apkin_a.getExc(qqqEfixalpha, theta_recon * 180. / M_PI);
      plotter->Fill1D(nA_label + "_Ex_from_alphas_QQQ", 1200, -30, 30, Ex_from_alpha, nA_label);
      plotter->Fill1D(nA_label + "_Ex_from_protons_QQQ", 1200, -30, 30, Ex_from_proton, nA_label);
#endif

      // Fill Gated Plots
      // if (vtx_gate != "")
      // {
      //   plotter->Fill1D(nA_label + "_twisted_pcz_recon_QQQ" + vtx_gate, 600, -300, 300, pcz_intersect.Z(), nA_label);
      //   plotter->Fill1D(nA_label + "_twisted_vertex_recon_QQQ" + vtx_gate, 600, -300, 300, vertex_recon, nA_label);
      //   plotter->Fill2D(nA_label + "_dE_Ecorr_Anode_QQQ" + vtx_gate, 400, 0, 30, 800, 0, 40000, qqqEfix, apSumE * sinTheta, nA_label);
      //   plotter->Fill2D(nA_label + "_dE_Ecorr_Anode_QQQ_alpha" + vtx_gate, 400, 0, 30, 800, 0, 40000, qqqEfixalpha, apSumE * sinTheta, nA_label);
      //   plotter->Fill1D(nA_label + "_Ex_from_alphas_QQQ" + vtx_gate, 1200, -30, 30, Ex_from_alpha, nA_label);
      //   plotter->Fill1D(nA_label + "_Ex_from_protons_QQQ" + vtx_gate, 1200, -30, 30, Ex_from_proton, nA_label);
      // }
    }
  }
#endif
  if (doOldAnalysis)
    OldAnalysis();
  return kTRUE;
}

void TrackRecon::Terminate()
{
  plotter->FlushToDisk(10);
}

void protonAlphaHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events)
{

  // Sidetrack for a(p,p)
  std::string aplabel = "a(p,p)";
  double initial_energy = 6.78;
  // if (dataset == "27Al")
  //   initial_energy = 6.78;
  // if (dataset == "17F")
  //   initial_energy = 6.78;

  Kinematics apkin_p(mass_1H, mass_4He, mass_1H, mass_4He, initial_energy); // m3 is proton
  Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, initial_energy); // m3 is alpha

  for (const auto &qqqevent : QQQ_Events)
  {
    for (const auto &sx3event : SX3_Events)
    {
      plotter->Fill1D("ap_qqq_sx3_dt", 800, -2000, 2000, qqqevent.Time1 - sx3event.Time1, aplabel);
      if (TMath::Abs(qqqevent.Time1 - sx3event.Time1) > 300)
        continue;
      // sx3event.pos.SetZ(sx3event.pos.Z()+5.0);
      plotter->Fill1D("ap_qqq_sx3_dt_timecut", 800, -2000, 2000, qqqevent.Time1 - sx3event.Time1, aplabel);
      plotter->Fill1D("ap_qqq_sx3_dphi", 180, -200,200, qqqevent.pos.Phi() * 180 / M_PI - sx3event.pos.Phi() * 180 / M_PI, aplabel);
      plotter->Fill2D("ap_qqq_sx3_dphi_vs_qqqphi", 180, -200,200, 180, -200,200, qqqevent.pos.Phi() * 180 / M_PI - sx3event.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, aplabel);
      plotter->Fill2D("ap_qqq_sx3_matrix", 400, 0, 10, 400, 0, 10, qqqevent.Energy1, sx3event.Energy1, aplabel);

      for (const auto &pcevent : PC_Events)
      {

        double pcz_fix = pcfix_func.Eval(pcevent.pos.Z()) - 5.0;
        TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
        TVector3 x1(qqqevent.pos);
        TVector3 v = x2f - x1;
        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
        TVector3 r_rhoMin_fix = x1 + t_minimum * v;
        double vertex_z = r_rhoMin_fix.Z();
        double theta_q = (qqqevent.pos - TVector3(0, 0, vertex_z)).Theta();
        // double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
        double sinTheta_customV = TMath::Sin(theta_q);
        double theta_s = (sx3event.pos - TVector3(0, 0, vertex_z)).Theta();
        // double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
        double sinTheta_s = TMath::Sin(theta_s);
        // if(vertex_z<0 || vertex_z>100) continue;

        // double sinTheta = TMath::Sin((qqqevent.pos - pcevent.pos).Theta());
        // plotter->Fill2D("sinTheta2_vs_sinTheta",80,-2,2,80,-2,2,sinTheta,sinTheta_customV,aplabel);

        plotter->Fill2D("ap_dE_E_Anodesx3B", 400, 0, 10, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, aplabel);
        plotter->Fill2D("ap_dE_E_Cathodesx3B", 400, 0, 10, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, aplabel);
        plotter->Fill2D("ap_dE_E_AnodeQQQ", 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, aplabel);
        plotter->Fill2D("ap_dE_E_CathodeQQQ", 400, 0, 10, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, aplabel);
        plotter->Fill2D("ap_dE3_E_AnodeQQQ", 400, 0, 10, 400, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV, aplabel);
        plotter->Fill2D("ap_dE3_E_CathodeQQQ", 400, 0, 10, 400, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV, aplabel);

        plotter->Fill2D("ap_dPhi_QQQ_PC", 180, -200,200, 180, -200,200, pcevent.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, aplabel);
        plotter->Fill2D("ap_dPhi_SX3_PC", 180, -200,200, 180, -200,200, pcevent.pos.Phi() * 180 / M_PI, sx3event.pos.Phi() * 180 / M_PI, aplabel);
        plotter->Fill1D("ap_dt_Anode_QQQ", 600, -2000, 2000, pcevent.Time1 - qqqevent.Time1, aplabel);
        plotter->Fill1D("ap_dt_Cathode_QQQ", 600, -2000, 2000, pcevent.Time2 - qqqevent.Time1, aplabel);
        plotter->Fill1D("ap_dt_Anode_SX3", 600, -2000, 2000, pcevent.Time1 - sx3event.Time1, aplabel);
        plotter->Fill1D("ap_dt_Cathode_SX3", 600, -2000, 2000, pcevent.Time2 - sx3event.Time1, aplabel);
        plotter->Fill1D("ap_pczfix", 600, -300, 300, pcz_fix, aplabel);
        plotter->Fill1D("ap_pcz", 600, -300, 300, pcevent.pos.Z(), aplabel);

        double path_length_q = (qqqevent.pos - TVector3(0, 0, vertex_z)).Mag() * 0.1;
        double path_length_s = (sx3event.pos - TVector3(0, 0, vertex_z)).Mag() * 0.1;
        // double path_length_q = (qqqevent.pos-r_rhoMin_fix).Mag()*0.1;
        // double path_length_s = (sx3event.pos-r_rhoMin_fix).Mag()*0.1;

        // We know that alphas predominantly are detected in QQQs, and protons in SX3s, and that protons don't leave much of a trace in dE layer.
        // Using the estimated path lengths, we correct alpha eloss in qqq, and protons in sx3. The result should (hopefully be) vertex independent.

        double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length_q);
        double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1) - path_length_s);
        // plotter->Fill2D("qqqEf_sx3E_matrix_all",400,0,10,400,0,10,qqqEfix,sx3event.Energy1,aplabel);
        plotter->Fill2D("ap_qqqEf_sx3Ef_matrix", 400, 0, 10, 400, 0, 10, qqqEfix, sx3Efix, aplabel);

        plotter->Fill2D("ap_Ef_vs_theta_qqq", 100, 0, 180, 400, 0, 10, theta_q * 180 / M_PI, qqqEfix, aplabel);
        plotter->Fill2D("ap_Ef_vs_theta_sx3", 100, 0, 180, 400, 0, 10, theta_s * 180 / M_PI, sx3Efix, aplabel);
        plotter->Fill2D("ap_theta_vs_theta_qqq_sx3", 100, 0, 180, 100, 0, 180, theta_q * 180 / M_PI, theta_s * 180 / M_PI, aplabel);
        plotter->Fill1D("ap_VertexReconZ", 400, -200, 200, vertex_z, aplabel);
        plotter->Fill2D("ap_VertexReconXY", 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), aplabel);
        plotter->Fill1D("ap_Ex_from_protons", 200, -10, 10, apkin_p.getExc(sx3Efix, theta_s * 180 / M_PI), aplabel);
        plotter->Fill1D("ap_Ex_from_alpha", 200, -10, 10, apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI), aplabel);

        if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
        { // one-anode, two-cathode events, as originally intended
          // std::cout << "Test" << std::endl;
          plotter->Fill1D("ap_VertexReconZ_a1c2", 400, -200, 200, vertex_z, aplabel);
          plotter->Fill2D("ap_VertexReconXY_a1c2", 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), aplabel);
          plotter->Fill2D("ap_theta_vs_theta_qqq_sx3_a1c2", 100, 0, 180, 100, 0, 180, theta_q * 180 / M_PI, theta_s * 180 / M_PI, aplabel);
          plotter->Fill2D("ap_Ef_vs_theta_qqq_a1c2", 100, 0, 180, 400, 0, 10, theta_q * 180 / M_PI, qqqEfix, aplabel);
          plotter->Fill1D("ap_Ex_from_protons_a1c2", 200, -10, 10, apkin_p.getExc(sx3Efix, theta_s * 180 / M_PI), aplabel);
          plotter->Fill1D("ap_Ex_from_alpha_a1c2", 200, -10, 10, apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI), aplabel);

          // std::cout << apkin_p.getExc(sx3Efix,theta_s*180/M_PI) << " " << apkin_a.getExc(qqqEfix,theta_q*180/M_PI)<<  std::endl;
          plotter->Fill2D("ap_Ef_vs_theta_sx3_a1c2", 100, 0, 180, 400, 0, 10, theta_s * 180 / M_PI, sx3Efix, aplabel);

          // plotter->Fill2D("qqqEf_sx3E_matrix",400,0,10,400,0,10,qqqEfix,sx3event.Energy1,aplabel);
          plotter->Fill2D("ap_qqq_sx3_matrix_a1c2", 400, 0, 10, 400, 0, 10, qqqevent.Energy1, sx3event.Energy1, aplabel);
          plotter->Fill2D("ap_qqqEf_sx3Ef_matrix_a1c2", 400, 0, 10, 400, 0, 10, qqqEfix, sx3Efix, aplabel);
          // std::cout << sx3event.Energy1 	<< " " <<  path_length_s << " " << sx3Efix << std::endl;

          // plotter->Fill2D("dE3_Ef_AnodeQQQ_a1c2",400,0,10,400,0,40000,qqqEfix,pcevent.Energy1*sinTheta_customV,aplabel);
          // plotter->Fill2D("dE3_Ef_CathodeQQQ_a1c2",400,0,10,400,0,10000,qqqEfix,pcevent.Energy2*sinTheta_customV,aplabel);

        } // end if(a1c2) loop
      } // end PC_Events for loop

    } // end SX3_Events for loop
  } // end QQQ_Events for loop, end sidetrack a(p,p)

  return;
}

void PCSX3ClusterAnalysis(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, const std::vector<std::vector<std::tuple<int, double, double>>> &cClusters)
{

  static TRandom3 rand(0);

  // --- GENUINE A1C0 events: one anode cluster and NO cathode cluster. These never
  // form an anode-cathode crossover, so the builder produces no PC_Events for them
  // and the pcevent loop below never sees them. Reconstruct them anode-only here, the
  // same way the A1C1 benchmark does: pair each Si hit with the anode pseudo-wire and
  // gate on anode-timestamp time coincidence + Si-anode phi + the beam-axis vertex.
  // Reference is the Si-only geometric guess pczguess (no cathode crossover exists). ---
  if (BenchMark && aClusters.size() == 1 && cClusters.size() == 0)
  {
    const auto &aCl = aClusters.front();
    auto vertexFrom = [](const TVector3 &si, const TVector3 &pcpoint)
    {
      TVector3 vf = pcpoint - si;
      double tm = -1.0 * (si.X() * vf.X() + si.Y() * vf.Y()) / (vf.X() * vf.X() + vf.Y() * vf.Y());
      return TVector3(si + tm * vf);
    };
    auto aPw = pwinstance.GetPseudoWire(aCl, "ANODE");
    auto apwire_bm = std::get<0>(aPw);
    double anodeTS = std::get<3>(aPw);
    for (const auto &sx3event : SX3_Events)
    {
      bool PCSX3TimeCut = (sx3event.Time1 - anodeTS < 150);
      TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, sx3event.pos.Phi());
      bool phicut = TMath::Abs(sx3event.pos.DeltaPhi(pc)) <= TMath::Pi() / 4.0;
      if (!(phicut && PCSX3TimeCut))
        continue;
      double smeared_phi = sx3event.pos.Phi() + rand.Uniform(-sx3_phi_pitch / 2.0, sx3_phi_pitch / 2.0);
      TVector3 smeared_sx3(sx3event.pos.Perp() * TMath::Cos(smeared_phi), sx3event.pos.Perp() * TMath::Sin(smeared_phi), sx3event.pos.Z());
      pc.SetZ(a1c1_zcorr(pc.Z(), true));
      double pc_dither = rand.Gaus(pc.Z(), dither_sigma);
      TVector3 vtx0 = vertexFrom(sx3event.pos, pc);
      TVector3 vtx1 = vertexFrom(smeared_sx3, TVector3(pc.X(), pc.Y(), pc_dither));

      if (!(vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6))
        continue;
      double sx3theta = TMath::ATan2(88.0, sx3event.pos.Z() - source_vertex);
      double pczguess = 37.0 / TMath::Tan(sx3theta) + source_vertex;
      plotter->Fill1D("Benchmark_SX3_VertexZ_trueA1C0", 800, -400, 400, vtx0.Z(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_VertexZ_trueA1C0_Hybrid", 800, -400, 400, vtx1.Z(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_VertexZ_trueA1C0_Hybrid_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 800, -400, 400, vtx1.Z(), "A1C0True_SX3");
      plotter->Fill2D("Benchmark_SX3_VertexXY_trueA1C0_Hybrid", 200, -100, 100, 200, -100, 100, vtx1.X(), vtx1.Y(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C0_Hybrid", 600, -200, 200, pc_dither, "A1C0True_SX3");
      plotter->Fill2D("Benchmark_SX3_PCZ_trueA1C0_Hybrid_vs_sx3pczguess", 400, -200, 200, 400, -200, 200, pczguess, pc_dither, "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C0_Hybrid_minus_sx3pczguess", 400, -100, 100, pc_dither - pczguess, "A1C0True_SX3");
    }
  }

  for (const auto &pcevent : PC_Events)
  {
    bool PCSX3TimeCut = false;
    bool PCASX3TimeCut = false;
    bool PCCSX3TimeCut = false;
    for (const auto &sx3event : SX3_Events)
    {
      plotter->Fill1D("dt_pcA_sx3B" + std::to_string(sx3event.ch2), 640, -2000, 2000, sx3event.Time1 - pcevent.Time1, "Timing");
      plotter->Fill1D("dt_pcC_sx3B" + std::to_string(sx3event.ch2), 640, -2000, 2000, sx3event.Time1 - pcevent.Time2, "Timing");
      if (sx3event.Time1 - pcevent.Time1 < 0) //-150 for alphas
        PCASX3TimeCut = 1;
      if (sx3event.Time1 - pcevent.Time2 < 0) //-200 for alphas
        PCCSX3TimeCut = 1;
      PCSX3TimeCut = PCASX3TimeCut && PCCSX3TimeCut;

      bool phicut = TMath::Abs(sx3event.pos.DeltaPhi(pcevent.pos)) <= TMath::Pi() / 4.0;

      plotter->Fill1D("dt_pcA_sx3B", 640, -2000, 2000, sx3event.Time1 - pcevent.Time1, "Timing");
      plotter->Fill1D("dt_pcC_sx3B", 640, -2000, 2000, sx3event.Time1 - pcevent.Time2, "Timing");
      plotter->Fill2D("dt_pcA_vs_sx3RE", 640, -2000, 2000, 400, 0, 30, sx3event.Time1 - pcevent.Time1, sx3event.Energy1, "Timing");
      plotter->Fill2D("dE_E_Anodesx3B", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, "PID_dE_E");
      plotter->Fill2D("dE_E_Cathodesx3B", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, "PID_dE_E");
      if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
        plotter->Fill2D("dE_E_Anodesx3B_a1c2", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, "PID_dE_E");
      if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
        plotter->Fill2D("dE_E_Cathodesx3B_a1c2", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, "PID_dE_E");
      if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
        plotter->Fill2D("dE_E_Anodesx3B_a2c1", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, "PID_dE_E");
      if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
        plotter->Fill2D("dE_E_Cathodesx3B_a2c1", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, "PID_dE_E");
      if (pcevent.multi1 == 1 && pcevent.multi2 == 1)
        plotter->Fill2D("dE_E_Anodesx3B_a1c1", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, "PID_dE_E");
      if (pcevent.multi1 == 1 && pcevent.multi2 == 1)
        plotter->Fill2D("dE_E_Cathodesx3B_a1c1", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, "PID_dE_E");
      if (pcevent.multi1 == 1 && pcevent.multi2 == 0)
        plotter->Fill2D("dE_E_Anodesx3B_a1c0", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, "PID_dE_E");
      if (pcevent.multi1 == 1 && pcevent.multi2 == 0)
        plotter->Fill2D("dE_E_Cathodesx3B_a1c0", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, "PID_dE_E");

      plotter->Fill2D("sx3phi_vs_pcphi" + std::to_string(sx3event.Time1 - pcevent.Time1 < -150), 100, -200,200, 100, -200,200, sx3event.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");
      plotter->Fill1D("sx3phi_minus_pcphi" + std::to_string(sx3event.Time1 - pcevent.Time1 < -150), 100, -180, 180, (sx3event.pos.DeltaPhi(pcevent.pos)) * 180 / M_PI, "Kinematics_Angles");

      if (PCSX3TimeCut)
      {
        plotter->Fill1D("dt_pcA_sx3B_timecut", 640, -2000, 2000, sx3event.Time1 - pcevent.Time1, "Timing");
        plotter->Fill1D("dt_pcC_sx3B_timecut", 640, -2000, 2000, sx3event.Time1 - pcevent.Time2, "Timing");
        plotter->Fill2D("xyplot_sx3" + std::to_string(sx3event.ch2 / 4), 100, -100, 100, 100, -100, 100, sx3event.pos.X(), sx3event.pos.Y(), "Vertex_Reconstruction");
        plotter->Fill2D("xyplot_sx3" + std::to_string(sx3event.ch2 / 4), 100, -100, 100, 100, -100, 100, pcevent.pos.X(), pcevent.pos.Y(), "Vertex_Reconstruction");
        plotter->Fill2D("pcz_vs_pcphi_TimeCut", 600, -200, 200, 120, -200,200, pcevent.pos.Z(), pcevent.pos.Phi() * 180 / M_PI, "PCZ_Recon");
      }

      double sx3rho = 88.0;
      double sx3z = sx3event.pos.Z();
      double pcz = pcevent.pos.Z();
      double calcsx3theta = TMath::ATan2(sx3rho - z_to_crossover_rho(pcz), sx3z - pcz);
      plotter->Fill2D("dE2_E_Anodesx3B", 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * TMath::Sin(calcsx3theta), "PID_dE_E");
      plotter->Fill2D("dE2_E_Cathodesx3B", 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * TMath::Sin(calcsx3theta), "PID_dE_E");

      double sx3theta = TMath::ATan2(sx3rho, sx3z - source_vertex);
      double pczguess = 37.0 / TMath::Tan(sx3theta) + source_vertex;
      double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z()) / TMath::Tan(sx3theta) + source_vertex;
      double sinTheta = TMath::Sin(sx3theta);

      TVector3 x2(pcevent.pos), x1(sx3event.pos);
      TVector3 v = x2 - x1;
      double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
      TVector3 vector_closest_to_z_sx3 = x1 + t_minimum * v;
      plotter->Fill1D("VertexReconZ_SX3" + std::to_string(PCSX3TimeCut), 600, -1300, 1300, vector_closest_to_z_sx3.Z(), "Vertex_Reconstruction");
      plotter->Fill1D("VertexReconZ_SX3", 600, -1300, 1300, vector_closest_to_z_sx3.Z(), "Vertex_Reconstruction");
      plotter->Fill2D("VertexReconXY_SX3" + std::to_string(PCSX3TimeCut), 100, -100, 100, 100, -100, 100, vector_closest_to_z_sx3.X(), vector_closest_to_z_sx3.Y(), "Vertex_Reconstruction");

      plotter->Fill2D("pcz_vs_time", 2000, 0, 2000, 600, -200, 200, pcevent.Time1 * 1e-9, pcevent.pos.Z(), "Timing");
      plotter->Fill2D("pcphi_vs_time", 2000, 0, 2000, 180, -200,200, pcevent.Time1 * 1e-9, pcevent.pos.Phi() * 180. / M_PI, "Timing");
      plotter->Fill2D("sx3phi_vs_time", 2000, 0, 2000, 180, -200,200, pcevent.Time1 * 1e-9, sx3event.pos.Phi() * 180. / M_PI, "Timing");

      plotter->Fill2D("pcz_vs_sx3pczguess", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");

      if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
      {
        plotter->Fill2D("pcz_vs_sx3pczguess_A1C2", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");
        double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());

        TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
        TVector3 v = x2f - x1;
        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
        TVector3 r_rhoMin_fix = x1 + t_minimum * v;
        plotter->Fill1D("VertexRecon_pczfix_sx3", 800, -300, 300, r_rhoMin_fix.Z(), "Vertex_Reconstruction");
        plotter->Fill1D("VertexRecon_pczfix", 800, -300, 300, r_rhoMin_fix.Z(), "Vertex_Reconstruction");
        plotter->Fill1D("pczfix_A1C2_1d_sx3", 600, -200, 200, pcz_fix, "PCZ_Recon");
        plotter->Fill2D("pczfix_vs_sx3pczguess_A1C2", 600, -200, 200, 600, -200, 200, pczguess, pcz_fix, "PCZ_Recon");
        plotter->Fill2D("pczfix_vs_sx3pczguess_int_A1C2", 600, -200, 200, 600, -200, 200, pcz_guess_int, pcz_fix, "PCZ_Recon");
        plotter->Fill2D("pczguess_vs_int", 600, -200, 200, 600, -200, 200, pcz_guess_int, pczguess, "PCZ_Recon");
        plotter->Fill1D("pczguess_vs_int_residualsx3", 200, -50, 50, pcz_guess_int - pczguess, "Residuals");
        plotter->Fill2D("pczfix_residual_vs_pczguess_A1C2", 600, -200, 200, 200, -100, 100, pczguess, pcz_fix - pczguess, "Residuals");
        plotter->Fill2D("pczfix_residual_vs_phi_A1C2", 200, 0, 6.28, 200, -100, 100, r_rhoMin_fix.Phi(), pcz_fix - pczguess, "Residuals");
        plotter->Fill2D("pczguess_vs_int_residual_vs_phi_A1C2", 200, 0, 6.28, 200, -100, 100, r_rhoMin_fix.Phi(), pcz_guess_int - pczguess, "Residuals");
        plotter->Fill1D("pczfix-sx3pczguess_A1C2", 200, -100, 100, pcz_fix - pczguess, "Residuals");
        plotter->Fill2D("pczfix_vs_sx3pczguess_A1C2_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");

        double sinTheta_customV = TMath::Sin((sx3event.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta());
        plotter->Fill2D("dE3_E_CathodeSX3_A1C2_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
        plotter->Fill2D("dE3_E_AnodeSX3_A1C2_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV, "PID_dE_E");

        if (TMath::Abs(r_rhoMin_fix.Z()) < 200.0)
        {
          plotter->Fill2D("dE3_E_AnodeSX3B_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV, "PID_dE_E");
          plotter->Fill2D("dE3_E_CathodeSX3B_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
        }

        if (pcevent.multi1 == 1 && pcevent.multi2 == 3)
        {
          plotter->Fill2D("pcz_vs_sx3pczguess_A1C3", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");
          plotter->Fill2D("pcz_vs_sx3pczguess_A1C3_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");
        }

        plotter->Fill2D("pcz_vs_sx3pczguess_int", 600, -200, 200, 600, -200, 200, pcz_guess_int, pcevent.pos.Z(), "PCZ_Recon");
        plotter->Fill2D("pcz_vs_sx3pczguess_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");

        bool sx3PhiCut = (TMath::Abs(sx3event.pos.Phi() - pcevent.pos.Phi()) < 45.0 * M_PI / 180.);

        plotter->Fill1D("pcz_sx3Coinc_phiCut" + std::to_string(sx3PhiCut) + "_TC" + std::to_string(PCSX3TimeCut), 300, 0, 200, sx3z, "PCZ_Recon");
        plotter->Fill2D("pcz_vs_sx3z_phiCut" + std::to_string(sx3PhiCut) + "_TC" + std::to_string(PCSX3TimeCut), 300, 0, 200, 600, -400, 400, sx3z, pcevent.pos.Z(), "PCZ_Recon");

        plotter->Fill2D("sx3E_vs_sx3z", 400, 0, 30, 300, 0, 200, sx3event.Energy1, sx3z, "Kinematics_Angles");

        plotter->Fill2D("pcdEA_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy1, sx3z, "Kinematics_Angles");
        plotter->Fill2D("pcdEA_vs_sx3pczguess", 800, 0, 20000, 600, -200, 200, pcevent.Energy1, pczguess, "Kinematics_Angles");
        plotter->Fill2D("pcdEA_vs_pczfix", 800, 0, 20000, 600, -200, 200, pcevent.Energy1, pcz_fix, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy2, sx3z, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_sx3pczguess", 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pczguess, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_pczfix", 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pcz_fix, "Kinematics_Angles");

        plotter->Fill2D("pcdEA_vs_sx3z" + std::to_string(sx3event.ch2), 800, 0, 20000, 300, 0, 200, pcevent.Energy1, sx3z, "Kinematics_Angles");
        plotter->Fill2D("pcdEA_vs_sx3pczguess" + std::to_string(sx3event.ch2), 800, 0, 20000, 600, -200, 200, pcevent.Energy1, pczguess, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_sx3z" + std::to_string(sx3event.ch2), 800, 0, 20000, 300, 0, 200, pcevent.Energy2, sx3z, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_sx3pczguess" + std::to_string(sx3event.ch2), 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pczguess, "Kinematics_Angles");

        plotter->Fill2D("pcdE2A_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy1 * sinTheta, sx3z, "Kinematics_Angles");
        plotter->Fill2D("pcdE2C_vs_sx3z", 800, 0, 20000, 300, 0, 200, pcevent.Energy2 * sinTheta, sx3z, "Kinematics_Angles");
        plotter->Fill2D("phi_vs_stripnum", 180, -180, 180, 48, 0, 48, pcevent.pos.Phi() * 180. / M_PI, sx3event.ch2, "Kinematics_Angles");
        plotter->Fill2D("E_theta_AnodeSX3", 400, -20, 180, 300, 0, 15, sx3theta * 180 / M_PI, sx3event.Energy1, "Kinematics_Angles");
      }
      if (PCSX3TimeCut)
      {
        plotter->Fill1D("PCZ_sx3", 800, -200, 200, pcevent.pos.Z(), "PCZ_Recon");
      }

      //-----------------------Benchmarking Method for Source Runs (SX3)------------------------//
      if (BenchMark && aClusters.size() == 1 && cClusters.size() == 1)
      {
        const auto &aCl = aClusters.front();
        const auto &cCl = cClusters.front();
        const std::string benchBranch = "Benchmark_SX3";
        auto vertexFrom = [](const TVector3 &si, const TVector3 &pcpoint)
        {
          TVector3 vf = pcpoint - si;
          double tm = -1.0 * (si.X() * vf.X() + si.Y() * vf.Y()) / (vf.X() * vf.X() + vf.Y() * vf.Y());
          return TVector3(si + tm * vf);
        };

        auto fillSuite = [&](const std::string &tag, double pcz_method, const TVector3 &vtx, const std::string &branch)
        {
          plotter->Fill1D("Benchmark_SX3_VertexZ_" + tag, 800, -400, 400, vtx.Z(), branch);
          plotter->Fill1D("Benchmark_SX3_VertexZ_" + tag + "_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 800, -400, 400, vtx.Z(), branch);
          plotter->Fill2D("Benchmark_SX3_VertexXY_" + tag, 200, -100, 100, 200, -100, 100, vtx.X(), vtx.Y(), branch);
          plotter->Fill1D("Benchmark_SX3_PCZ_" + tag, 600, -200, 200, pcz_method, branch);
        };

        auto fillVsRef = [&](const std::string &tag, double pcz_method, const TVector3 &vtx, double pcz_ref, const TVector3 &vtx_ref)
        {
          plotter->Fill2D("Benchmark_SX3_PCZ_" + tag + "_vs_ref", 400, -200, 200, 400, -200, 200, pcz_ref, pcz_method, "Benchmark_SX3_ref");
          plotter->Fill1D("Benchmark_SX3_PCZ_" + tag + "_minus_ref", 400, -100, 100, pcz_method - pcz_ref, "Benchmark_SX3_ref");
          plotter->Fill2D("Benchmark_SX3_PCZ_" + tag + "_vs_sx3pczguess", 400, -200, 200, 400, -200, 200, pczguess, pcz_method, "Benchmark_SX3_ref");
          plotter->Fill1D("Benchmark_SX3_PCZ_" + tag + "_minus_sx3pczguess", 400, -100, 100, pcz_method - pczguess, "Benchmark_SX3_ref");
        };

        double pcz_ref = pcfix_func.Eval(pcevent.pos.Z());
        TVector3 vtx_ref = vertexFrom(sx3event.pos, TVector3(pcevent.pos.X(), pcevent.pos.Y(), pcz_ref));

        auto pw_tuple = pwinstance.GetPseudoWire(aCl, "ANODE");
        std::pair<TVector3, TVector3> apwire_bm = std::get<0>(pw_tuple);

        auto cMaxWire = *std::max_element(cCl.begin(), cCl.end(), [](const auto &a, const auto &b)
                                          { return std::get<1>(a) < std::get<1>(b); });
        auto aMaxWire = *std::max_element(aCl.begin(), aCl.end(), [](const auto &a, const auto &b)
                                          { return std::get<1>(a) < std::get<1>(b); });
        std::vector<std::tuple<int, double, double>> cOne = {cMaxWire};

        auto xo_tuple = pwinstance.FindCrossoverProperties(aCl, cOne);
        TVector3 xo_a1c1 = std::get<0>(xo_tuple);
        double alpha_a1c1 = std::get<1>(xo_tuple);
        bool a1c1Good = (alpha_a1c1 != 9999999 && std::get<2>(xo_tuple) != -1);

        // --- A1C1 charge fraction (single max-E cathode vs anode, pseudo-wire sums) ---
        double aSumE_bm = std::get<1>(pw_tuple);
        double cSumE_bm = std::get<1>(cMaxWire);
        double ac_sum = aSumE_bm + cSumE_bm;
        double cfrac = (ac_sum > 0.0) ? cSumE_bm / ac_sum : -1.0;

        // Cmax/anode vs Si phi for ALL events reaching this block (not just the
        // A1C2 subset). Cmax = single max-E cathode wire (cSumE_bm), anode =
        // anode pseudo-wire sum. Ratio is unbounded above (arbitrary cathode gains).
        if (aSumE_bm > 0.0)
          plotter->Fill2D("Benchmark_SX3_CmaxOverAnode_vs_phi", 180, -180, 180, 250, 0, 5,
                          sx3event.pos.Phi() * 180. / M_PI, cSumE_bm / aSumE_bm, "Benchmark_SX3_ref");

        // Smearing setup
        double sx3_phi_pitch = 6.5 * (M_PI / 180.0);
        double smeared_phi = sx3event.pos.Phi() + rand.Uniform(-sx3_phi_pitch / 2.0, sx3_phi_pitch / 2.0);
        TVector3 smeared_sx3_pos(sx3event.pos.Perp() * TMath::Cos(smeared_phi), sx3event.pos.Perp() * TMath::Sin(smeared_phi), sx3event.pos.Z());

        // Plain (dither=false) shows the discrete wire spikes; _Hyb (Si phi smeared
        // + PC z dithered) is the smoothed cross-reference.
        auto doA1C1 = [&](const std::string &tag, const TVector3 &si_point, bool dither = true)
        {
          if (!a1c1Good)
            return;
          double pcz = dither ? rand.Gaus(xo_a1c1.Z(), dither_sigma) : xo_a1c1.Z();
          TVector3 vtx = vertexFrom(si_point, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz));
          fillSuite(tag, pcz, vtx, benchBranch);
          fillVsRef(tag, pcz, vtx, pcz_ref, vtx_ref);
        };

        auto doAnodeOnly = [&](const std::string &tag, double phi_use, const TVector3 &si_point, bool dither = true)
        {
          TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, phi_use);
          pc.SetZ(a1c1_zcorr(pc.Z(), false));
          TVector3 vtx0 = vertexFrom(si_point, pc);
          if (!(vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6))
            return;
          double pcz = dither ? rand.Gaus(pc.Z(), dither_sigma) : pc.Z();
          TVector3 vtx = vertexFrom(si_point, TVector3(pc.X(), pc.Y(), pcz));
          fillSuite(tag, pcz, vtx, benchBranch);
          fillVsRef(tag, pcz, vtx, pcz_ref, vtx_ref);
        };

        // --- A1C1 with the cfrac sub-cell model (linear centre-fold) ---
        // Anchored on the fired cathode (xo_a1c1.Z()); cell from the anode z;
        // |offset| from cfrac. REJECT (no fill) when out of band or inconsistent.
        auto doA1C1Model = [&](const std::string &tag, const TVector3 &si_point)
        {
          if (!a1c1Good || cfrac < 0.0)
            return;
          A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm);
          // side from the beam-axis 2-hypothesis test (the Si hit picks the cell)
          int side_status;
          SideChoice side = a1c1_pick_side(si_point, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
          const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
          double pcz_pick = best.pcz;
          if (!(best.inband && side_status != 2))
            return;
          TVector3 vtx = vertexFrom(si_point, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_pick));
          fillSuite(tag, pcz_pick, vtx, benchBranch);
          fillVsRef(tag, pcz_pick, vtx, pcz_ref, vtx_ref);
        };

        if (phicut && PCSX3TimeCut)
        {
          if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
          {
            fillSuite("A1C2", pcz_ref, vtx_ref, benchBranch);
            // XY-offset diagnostic: if the source/beam is off-axis, the vertex-z
            // residual modulates as cos(phi - phi0). A flat line = no XY offset.
            {
              double phi_deg = sx3event.pos.Phi() * 180.0 / M_PI;
              double vz_resid = vtx_ref.Z() - source_vertex;
              plotter->Fill2D("Diag_SX3_A1C2_vtxZ_resid_vs_phi", 180, -180, 180, 400, -100, 100,
                              phi_deg, vz_resid, "Diag_XYoffset");
              plotter->Fill2D("Diag_Combined_A1C2_vtxZ_resid_vs_phi", 180, -180, 180, 400, -100, 100,
                              phi_deg, vz_resid, "Diag_XYoffset");
              plotter->Fill2D("Diag_SX3_A1C2_vtxXY", 200, -15, 15, 200, -15, 15,
                              vtx_ref.X(), vtx_ref.Y(), "Diag_XYoffset");
              plotter->Fill2D("Diag_Combined_A1C2_time_vs_phi", 2000, 0, 2000, 180, -180, 180,
                              pcevent.Time1 * 1e-9, phi_deg, "Diag_XYoffset");
              plotter->Fill2D("Diag_SX3_A1C2_T_vs_vtxX", 2000, 0, 2000, 200, -15, 15,
                              pcevent.Time1 * 1e-9, vtx_ref.X(), "Diag_XYoffset");
              plotter->Fill2D("Diag_SX3_A1C2_T_vs_vtxY", 2000, 0, 2000, 200, -15, 15,
                              pcevent.Time1 * 1e-9, vtx_ref.Y(), "Diag_XYoffset");
            }

            doA1C1("A1C1", sx3event.pos, false);
            doAnodeOnly("A1C0", sx3event.pos.Phi(), sx3event.pos, false);
            doA1C1("A1C1_Hyb", smeared_sx3_pos);
            doAnodeOnly("A1C0_Hyb", smeared_phi, smeared_sx3_pos);

            // cfrac-model-corrected A1C1, directly comparable to A1C1 / A1C2
            doA1C1Model("A1C1_Cfrac", sx3event.pos);

            // A1C0 anode-only PCz bias vs track theta: test whether the z_a1c0 offset
            // is theta-driven (forward tracks cross the PC at a glancing angle). If
            // QQQ+SX3 collapse onto ONE curve (~0 at 90 deg, rising forward ~cot theta),
            // the offset is a function of theta, not per-detector. theta from vtx_ref.
            {
              double pcz_a1c0 = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, sx3event.pos.Phi()).Z();
              double theta_ref = (sx3event.pos - TVector3(0, 0, vtx_ref.Z())).Theta() * 180. / M_PI;
              plotter->Fill2D("Benchmark_SX3_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200,
                              theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_SX3_ref");
              plotter->Fill2D("Benchmark_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200,
                              theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_AnodeOnly");

              // A1C0 z-residual vs phi: if XY offset is the root cause, this should
              // show the same cos(phi-phi0) as the A1C2 vertex diagnostic above.
              double phi_deg_a = sx3event.pos.Phi() * 180.0 / M_PI;
              plotter->Fill2D("Diag_SX3_A1C0_zresid_vs_phi", 180, -180, 180, 200, -100, 100,
                              phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
              plotter->Fill2D("Diag_Combined_A1C0_zresid_vs_phi", 180, -180, 180, 200, -100, 100,
                              phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
            }

            // --- A1C1 charge-fraction diagnostics ---
            if (a1c1Good && cfrac >= 0.0)
            {
              plotter->Fill1D("Benchmark_SX3_A1C1_cfrac", 220, -0.05, 1.05, cfrac, "Benchmark_SX3_ref");
              plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_ref", 400, -200, 200, 220, -0.05, 1.05,
                              pcz_ref, cfrac, "Benchmark_SX3_ref");
              plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_sx3pczguess", 400, -200, 200, 220, -0.05, 1.05,
                              pczguess, cfrac, "Benchmark_SX3_ref");

              static const double zg[8] = {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};
              double zp = xo_a1c1.Z();
              auto fillCfracS = [&](const char *name, double truth)
              {
                double sgn = (truth >= zp) ? 1.0 : -1.0;
                double znb = (sgn > 0) ? 1.0e30 : -1.0e30;
                for (int i = 0; i < 8; ++i)
                {
                  if (sgn > 0 && zg[i] > zp + 1e-6 && zg[i] < znb)
                    znb = zg[i];
                  if (sgn < 0 && zg[i] < zp - 1e-6 && zg[i] > znb)
                    znb = zg[i];
                }
                if (TMath::Abs(znb) < 1e8 && TMath::Abs(znb - zp) > 0.0)
                  plotter->Fill2D(name, 240, -1.2, 1.2, 220, -0.05, 1.05,
                                  (truth - zp) / TMath::Abs(znb - zp), cfrac, "Benchmark_SX3_ref");
              };
              fillCfracS("Benchmark_SX3_A1C1_cfrac_vs_s", pcz_ref);
              fillCfracS("Benchmark_SX3_A1C1_cfrac_vs_s_sx3pczguess", pczguess);

              // folded about the cell CENTRE: f = |ref - z_center|/halfcell in [0,1]
              for (int i = 0; i < 7; ++i)
              {
                if (pcz_ref <= zg[i] && pcz_ref > zg[i + 1])
                {
                  double zc = 0.5 * (zg[i] + zg[i + 1]);
                  double half = 0.5 * (zg[i] - zg[i + 1]);
                  if (half > 0.0)
                    plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05,
                                    TMath::Abs(pcz_ref - zc) / half, cfrac, "Benchmark_SX3_ref");
                  break;
                }
              }

              // --- item 2: cfrac vs anode E (is cfrac a pure position variable?) ---
              // Ideal cfrac depends only on z, not on the anode signal. A flat band
              // here => independent of dE/anode gain; a tilt => dE/gain leakage into
              // the reconstructed z. The low band sitting at LOW anode E would also
              // confirm it as an incomplete-integration effect (feeds item 5).
              plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05,
                              aSumE_bm, cfrac, "Benchmark_SX3_ref");
              // r-space linearity test: if C = c*A + C_off (additive cathode offset),
              // then r = C/A = c + C_off*(1/A) -> a STRAIGHT line vs 1/anodeE (slope
              // C_off = removable contamination, intercept c = z-bearing ratio).
              if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                plotter->Fill2D("Benchmark_SX3_A1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0,
                                1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "Benchmark_SX3_ref");

              // --- item 3: cell/side misclassification rate ---
              // Truth cell from the precise A1C2 crossover (pcz_ref) vs the cell the
              // beam-axis side test (a1c1_pick_side) selects. A wrong cell = wrong side
              // of the fired wire = a full-cell (~40-80 mm) error.
              // Profile *_vs_fold: failures should pile up near the wires/edges.
              {
                A1C1Sol sm = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
                int sm_side_status;
                SideChoice sm_side = a1c1_pick_side(sx3event.pos, xo_a1c1.X(), xo_a1c1.Y(), sm.pcz_lo, sm.pcz_hi, sm_side_status);
                int sm_cell = (sm_side == SideChoice::High) ? sm.hi.cell : sm.lo.cell;
                int cell_truth = -1;
                for (int i = 0; i < 7; ++i)
                  if (pcz_ref <= a1c1_zg[i] && pcz_ref > a1c1_zg[i + 1])
                  {
                    cell_truth = i;
                    break;
                  }
                if (cell_truth >= 0)
                {
                  bool wrong = (sm_cell != cell_truth);
                  plotter->Fill2D("Benchmark_SX3_A1C1_cellsel_confusion", 7, 0, 7, 7, 0, 7,
                                  cell_truth + 0.5, sm_cell + 0.5, "Benchmark_SX3_ref");
                  plotter->Fill1D("Benchmark_SX3_A1C1_cellsel_misclass", 2, 0, 2, wrong ? 1.0 : 0.0, "Benchmark_SX3_ref");
                  plotter->Fill2D("Benchmark_SX3_A1C1_cellsel_misclass_vs_cell", 7, 0, 7, 2, 0, 2,
                                  cell_truth + 0.5, wrong ? 1.0 : 0.0, "Benchmark_SX3_ref");
                  double zc = 0.5 * (a1c1_zg[cell_truth] + a1c1_zg[cell_truth + 1]);
                  double half = 0.5 * (a1c1_zg[cell_truth] - a1c1_zg[cell_truth + 1]);
                  if (half > 0.0)
                  {
                    plotter->Fill2D("Benchmark_SX3_A1C1_cellsel_misclass_vs_fold", 120, 0, 1.2, 2, 0, 2,
                                    TMath::Abs(pcz_ref - zc) / half, wrong ? 1.0 : 0.0, "Benchmark_SX3_ref");
                    // item 5 validation: post-fold cfrac vs truth fold. With the fold
                    // ON (rfactor>0) the low band collapses onto the main band -> one
                    // line here instead of two.
                    plotter->Fill2D("Benchmark_SX3_A1C1_cfracUsed_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05,
                                    TMath::Abs(pcz_ref - zc) / half, sm.cfrac_used, "Benchmark_SX3_ref");
                    // anode-E correction validation: post-correction cfrac should be
                    // FLAT vs anodeE (the slope/curve of cfrac_vs_anodeE is removed).
                    if (aSumE_bm > 0.0)
                      plotter->Fill2D("Benchmark_SX3_A1C1_cfracUsed_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05,
                                      aSumE_bm, sm.cfrac_used, "Benchmark_SX3_ref");
                  }
                }
              }
            }
          }

          // ---- TRUE A1C1: genuine one-anode/one-cathode events. There is no
          // A1C2 crossover reference for these, so compare only to the Si
          // geometric guess (pczguess). xo_a1c1 is the real measured crossover. ----
          else if (pcevent.multi1 == 1 && pcevent.multi2 == 1 && a1c1Good)
          {
            double pcz_raw = xo_a1c1.Z();
            TVector3 vtx_raw = vertexFrom(sx3event.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_raw));
            fillSuite("trueA1C1", pcz_raw, vtx_raw, "A1C1True_SX3");
            plotter->Fill2D("Benchmark_SX3_PCZ_trueA1C1_vs_sx3pczguess", 400, -200, 200, 400, -200, 200, pczguess, pcz_raw, "A1C1True_SX3");
            plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C1_minus_sx3pczguess", 400, -100, 100, pcz_raw - pczguess, "A1C1True_SX3");

            // cfrac sub-cell reconstruction, anchored on the FIRED CATHODE.
            // Vertex-independent (no pczguess): the cell is the one adjacent to
            // the cathode that fired (zf), the side is taken from the anode z,
            // and the hit must land within one cell pitch of zf to be valid.
            if (cfrac >= 0.0)
            {
              A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
              // side from the beam-axis 2-hypothesis test (smaller-Perp vertex wins)
              int side_status;
              SideChoice side = a1c1_pick_side(sx3event.pos, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
              const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
              int cell = best.cell;
              double f = best.f;
              double pcz_cf = best.pcz;
              bool valid = (side_status != 2); // accept if at least one side is beam-axis consistent
              plotter->Fill1D("Benchmark_SX3_trueA1C1_sideStatus", 4, -1, 3, side_status + 0.5, "A1C1True_SX3");

              TVector3 vtx_cf = vertexFrom(sx3event.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_cf));
              fillSuite(valid ? "trueA1C1_Cfrac" : "trueA1C1_Cfrac_invalid", pcz_cf, vtx_cf, "A1C1True_SX3");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_cfrac", 220, -0.05, 1.05, cfrac, "A1C1True_SX3");
              // item 2: cfrac vs anode E for genuine A1C1 (no A1C2 ref here)
              plotter->Fill2D("Benchmark_SX3_trueA1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05, aSumE_bm, cfrac, "A1C1True_SX3");
              // r-space linearity test (r = c + C_off/anodeE should be a straight line)
              if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                plotter->Fill2D("Benchmark_SX3_trueA1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0,
                                1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "A1C1True_SX3");
              // reference-free per-cell cfrac (cell from geometry, no A1C2 ref): the
              // offline fitter reads per-cell edges/percentiles to gain-match cfmin/k.
              plotter->Fill2D("Benchmark_SX3_trueA1C1_cfrac_vs_cell", 7, 0, 7, 220, -0.05, 1.05, cell + 0.5, cfrac, "A1C1True_SX3");
              plotter->Fill2D("Benchmark_SX3_trueA1C1_f_vs_cell", 7, 0, 7, 260, -1.5, 2.5, cell + 0.5, f, "A1C1True_SX3");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_f", 260, -1.5, 2.5, f, "A1C1True_SX3");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_valid", 2, 0, 2, valid ? 1.0 : 0.0, "Benchmark_SX3_trueA1C1");
              // failure-reason breakdown (why the cfrac estimate is / isn't usable):
              //   0 valid & f in [0,1]   ideal, inside the calibrated band
              //   1 valid & f in [-1,0)  below band but within the pitch tolerance
              //   2 valid & f in (1,3]   above band but within the pitch tolerance
              //   3 invalid: f < -1      cfrac far below band (cfmin too high / cathode low?)
              //   4 invalid: f > 3       cfrac far above band (cfmin too low / cathode high?)
              //   5 cell k<=0            cell not autocalibrated
              int reason;
              if (a1c1_k_cell[cell] <= 0.0)
                reason = 5;
              else if (!valid)
                reason = (f < 0.0) ? 3 : 4;
              else if (f < 0.0)
                reason = 1;
              else if (f > 1.0)
                reason = 2;
              else
                reason = 0;
              plotter->Fill1D("Benchmark_SX3_trueA1C1_failreason", 6, 0, 6, reason + 0.5, "Benchmark_SX3_trueA1C1");
              // valid-only breakdown (reason is 0/1/2 here): how clean the accepted
              // sample is -- 0 ideal (f in [0,1]), 1 below-band, 2 above-band marginal.
              if (valid)
                plotter->Fill1D("Benchmark_SX3_trueA1C1_validreason", 3, 0, 3, reason + 0.5, "Benchmark_SX3_trueA1C1");
              // which cfrac band reconstructed this event: 0 = main, 1 = low (incomplete integration)
              plotter->Fill1D("Benchmark_SX3_trueA1C1_band", 2, 0, 2, s.band + 0.5, "Benchmark_SX3_trueA1C1");
              if (valid)
                plotter->Fill1D("Benchmark_SX3_trueA1C1_band_valid", 2, 0, 2, s.band + 0.5, "Benchmark_SX3_trueA1C1");
              // source-run-only diagnostics: pczguess is meaningful only when the vertex is fixed
              if (valid)
              {
                plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C1_Cfrac_minus_sx3pczguess_DIAG", 400, -100, 100, pcz_cf - pczguess, "Benchmark_SX3_trueA1C1");
                plotter->Fill2D("Benchmark_SX3_PCZ_trueA1C1_Cfrac_vs_sx3pczguess_DIAG", 400, -200, 200, 400, -200, 200, pczguess, pcz_cf, "Benchmark_SX3_trueA1C1");
              }
            }

            // anode-only (A1C0) estimate for the same events
            {
              TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, sx3event.pos.Phi());
              TVector3 vtx0 = vertexFrom(sx3event.pos, pc);
              if (vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6)
              {
                fillSuite("A1C1asA1C0", pc.Z(), vtx0, "A1C1True_SX3");
                plotter->Fill2D("Benchmark_SX3_PCZ_A1C1asA1C0_vs_sx3pczguess", 400, -200, 200, 400, -200, 200, pczguess, pc.Z(), "A1C1True_SX3");
              }
            }
          }
        }
      }
    }
  }
}

void PCQQQClusterAnalysis(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, const std::vector<std::vector<std::tuple<int, double, double>>> &cClusters)
{
  static TRandom3 rand(0);

  // --- GENUINE A1C0 events (QQQ twin): one anode cluster and NO cathode cluster.
  // See the SX3 block for rationale. Anode-timestamp time coincidence + Si-anode phi
  // + beam-axis vertex gate. Reference is the Si-only guess pcz_guess_37. ---
  if (BenchMark && aClusters.size() == 1 && cClusters.size() == 0)
  {
    const auto &aCl = aClusters.front();
    auto vertexFrom = [](const TVector3 &si, const TVector3 &pcpoint)
    {
      TVector3 vf = pcpoint - si;
      double tm = -1.0 * (si.X() * vf.X() + si.Y() * vf.Y()) / (vf.X() * vf.X() + vf.Y() * vf.Y());
      return TVector3(si + tm * vf);
    };
    auto aPw = pwinstance.GetPseudoWire(aCl, "ANODE");
    auto apwire_bm = std::get<0>(aPw);
    double anodeTS = std::get<3>(aPw);
    for (const auto &qqqevent : QQQ_Events)
    {
      bool timecut = (qqqevent.Time1 - anodeTS < 150);
      double smeared_phi = qqqevent.pos.Phi() + rand.Uniform(-qqq_wedge_pitch / 2.0, qqq_wedge_pitch / 2.0);
      TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, smeared_phi);

      bool phicut = TMath::Abs(qqqevent.pos.DeltaPhi(pc)) <= TMath::Pi() / 4.0;

      if (!(phicut && timecut))
        continue;
      double smeared_rho = qqqevent.pos.Perp() + rand.Uniform(-qqq_ring_pitch / 2.0, qqq_ring_pitch / 2.0);
      TVector3 smeared_qqq(smeared_rho * TMath::Cos(smeared_phi), smeared_rho * TMath::Sin(smeared_phi), qqqevent.pos.Z());
      pc.SetZ(a1c1_zcorr(pc.Z(), true));
      double pc_dither = rand.Gaus(pc.Z(), dither_sigma_c0 / 2.0);
      TVector3 vtx0 = vertexFrom(qqqevent.pos, pc);
      TVector3 vtx1 = vertexFrom(smeared_qqq, TVector3(pc.X(), pc.Y(), pc_dither));

      if (!(vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6))
        continue;
      double qqqTheta = (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta();
      double pcz_guess_37 = 37. / TMath::Tan(qqqTheta) + source_vertex;
      plotter->Fill1D("Benchmark_QQQ_VertexZ_trueA1C0", 800, -400, 400, vtx0.Z(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_VertexZ_trueA1C0_Hybrid", 800, -400, 400, vtx1.Z(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_VertexZ_trueA1C0_Hybrid_TC" + std::to_string(timecut) + "_PC" + std::to_string(phicut), 800, -400, 400, vtx1.Z(), "A1C0True_QQQ");
      plotter->Fill2D("Benchmark_QQQ_VertexXY_trueA1C0_Hybrid", 200, -100, 100, 200, -100, 100, vtx1.X(), vtx1.Y(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C0_Hybrid", 600, -200, 200, pc_dither, "A1C0True_QQQ");
      plotter->Fill2D("Benchmark_QQQ_PCZ_trueA1C0_Hybrid_vs_qqqpczguess", 400, -200, 200, 400, -200, 200, pcz_guess_37, pc_dither, "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C0_Hybrid_minus_qqqpczguess", 400, -100, 100, pc_dither - pcz_guess_37, "A1C0True_QQQ");
    }
  }

  for (const auto &pcevent : PC_Events)
  {
    for (const auto &qqqevent : QQQ_Events)
    {
      plotter->Fill1D("dt_pcA_qqqR", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1, "Timing");
      plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE", 640, -2000, 2000, 400, 0, 30, qqqevent.Time1 - pcevent.Time1, qqqevent.Energy1, "Timing");
      plotter->Fill1D("dt_pcC_qqqW", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2, "Timing");
      plotter->Fill2D("phiPC_vs_phiQQQ", 180, -200,200, 180, -200,200, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");

      double qqqTheta = (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta();
      double sinTheta = TMath::Sin(qqqTheta);

      TVector3 x2(pcevent.pos);
      TVector3 x1(qqqevent.pos);
      TVector3 v = x2 - x1;
      double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
      TVector3 r_rhoMin = x1 + t_minimum * v;

      bool timecut = (qqqevent.Time1 - pcevent.Time1 < 150);
      bool lowercut_cath = pcevent.Energy2 * sinTheta < 1 && (qqqevent.Energy2 < 5.0 || qqqevent.Energy1 < 5.0);

      bool phicut = TMath::Abs(qqqevent.pos.DeltaPhi(pcevent.pos)) <= TMath::Pi() / 4.0;

      if (lowercut_cath && phicut)
      {
        plotter->Fill1D("dt_pcA_qqqR_pidlow_PC1", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1, "Timing");
        plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE_pidlow_PC1", 640, -2000, 2000, 400, 0, 30, qqqevent.Time1 - pcevent.Time1, qqqevent.Energy1, "Timing");
        plotter->Fill1D("dt_pcC_qqqW_pidlow_PC1", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2, "Timing");
      }
      if (timecut)
      {
        plotter->Fill2D("dE_E_AnodeQQQR", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
        plotter->Fill2D("dE_E_CathodeQQQR", 400, 0, 30, 800, 0, 10000, qqqevent.Energy2, pcevent.Energy2, "PID_dE_E");
        if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
          plotter->Fill2D("dE_E_AnodeQQQR_a1c2", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
        if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
          plotter->Fill2D("dE_E_CathodeQQQR_a1c2", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
          plotter->Fill2D("dE_E_AnodeQQQR_a2c1", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
        if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
          plotter->Fill2D("dE_E_CathodeQQQR_a2c1", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        if (pcevent.multi1 == 1 && pcevent.multi2 == 1)
          plotter->Fill2D("dE_E_Anodesx3B_a1c1", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
        if (pcevent.multi1 == 1 && pcevent.multi2 == 1)
          plotter->Fill2D("dE_E_Cathodesx3B_a1c1", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        if (pcevent.multi1 == 1 && pcevent.multi2 == 0)
          plotter->Fill2D("dE_E_Anodesx3B_a1c0", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
        if (pcevent.multi1 == 1 && pcevent.multi2 == 0)
          plotter->Fill2D("dE_E_Cathodesx3B_a1c0", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");

        if (phicut)
        {
          plotter->Fill2D("dE2_E_AnodeQQQR_TC1PC1_pidlow" + std::to_string(lowercut_cath), 400, 0, 30, 800, 0, 4000, qqqevent.Energy1, pcevent.Energy1 * sinTheta, "PID_dE_E");
          plotter->Fill2D("dE2_E_CathodeQQQW_TC1PC1_pidlow" + std::to_string(lowercut_cath), 400, 0, 30, 800, 0, 1000, qqqevent.Energy2, pcevent.Energy2 * sinTheta, "PID_dE_E");
          plotter->Fill2D("E_theta_zoomin_AnodeQQQR_TC1PC1_pidlow" + std::to_string(lowercut_cath), 60, 0, 30, 300, 0, 15, qqqTheta * 180 / M_PI, qqqevent.Energy1, "Kinematics_Angles");
        }

        plotter->Fill2D("dE2_E_AnodeQQQR_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 4000, qqqevent.Energy1, pcevent.Energy1 * sinTheta, "PID_dE_E");
        plotter->Fill2D("dE2_E_CathodeQQQR_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 1000, qqqevent.Energy2, pcevent.Energy2 * sinTheta, "PID_dE_E");
        plotter->Fill2D("dEC_vs_dEA_TC1_PC" + std::to_string(phicut), 800, 0, 40000, 800, 0, 10000, pcevent.Energy1, pcevent.Energy2, "PID_dE_E");
        plotter->Fill2D("qqqphi_vs_time", 2000, 0, 2000, 180, -200,200, pcevent.Time1 * 1e-9, qqqevent.pos.Phi() * 180. / M_PI, "Timing");

        plotter->Fill1D("dt_pcA_qqqR_timecut", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1, "Timing");
        plotter->Fill1D("dt_pcC_qqqW_timecut", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2, "Timing");
        plotter->Fill2D("dE_theta_AnodeQQQR", 90, 0, 90, 400, 0, 20000, qqqTheta * 180 / M_PI, pcevent.Energy1, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_AnodeQQQR_zoomin", 60, 0, 30, 400, 0, 5000, qqqTheta * 180 / M_PI, pcevent.Energy1 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_AnodeQQQR", 90, 0, 90, 400, 0, 20000, qqqTheta * 180 / M_PI, pcevent.Energy1 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("phiPC_vs_phiQQQ_TimeCut", 180, -200,200, 180, -200,200, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");
        double pcz_guess_37 = 37. / TMath::Tan(qqqTheta) + source_vertex;

        if ((qqqevent.ch1) % 16 == 7)
          plotter->Fill2D("phiPC_vs_phiQQQ_TimeCut_allring8", 180, -200,200, 180, -200,200, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");

        plotter->Fill1D("phiQQQ_minus_phiPC_TimeCut_QQQ" + std::to_string(qqqevent.ch1 / 16), 180, -180, 180, qqqevent.pos.DeltaPhi(pcevent.pos) * 180 / M_PI, "Kinematics_Angles");

        plotter->Fill2D("Etot2_theta_AnodeQQQR", 75, 0, 90, 300, 0, 15, qqqTheta * 180 / M_PI, qqqevent.Energy1 + pcevent.Energy1 * anode_gain * sinTheta, "Kinematics_Angles");

        plotter->Fill2D("dE_theta_CathodeQQQR", 75, 0, 90, 800, 0, 10000, qqqTheta * 180 / M_PI, pcevent.Energy2, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_CathodeQQQR", 75, 0, 90, 800, 0, 10000, qqqTheta * 180 / M_PI, pcevent.Energy2 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_CathodeQQQR_zoomin", 60, 0, 30, 800, 0, 3000, qqqTheta * 180 / M_PI, pcevent.Energy2 * sinTheta, "Kinematics_Angles");

        plotter->Fill2D("dE_phi_AnodeQQQR", 100, -180, 180, 800, 0, 40000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Phi() * 180 / M_PI, pcevent.Energy1, "Kinematics_Angles");
        plotter->Fill2D("dE_phi_CathodeQQQR", 100, -180, 180, 800, 0, 40000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Phi() * 180 / M_PI, pcevent.Energy2, "Kinematics_Angles");
        plotter->Fill2D("pcdEA_vs_qqqpczguess", 800, 0, 20000, 600, -200, 200, pcevent.Energy1, pcz_guess_37, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_qqqpczguess", 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pcz_guess_37, "Kinematics_Angles");

        plotter->Fill1D("PCZ", 800, -200, 200, pcevent.pos.Z(), "PCZ_Recon");

        plotter->Fill2D("pczguess_vs_pc_37", 180, 0, 200, 150, 0, 200, pcz_guess_37, pcevent.pos.Z(), "PCZ_Recon");

        double pcz_guess_42 = 42. / TMath::Tan(qqqTheta) + source_vertex;
        plotter->Fill2D("pczguess_vs_pc_42", 180, 0, 200, 150, 0, 200, pcz_guess_42, pcevent.pos.Z(), "PCZ_Recon");

        double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z()) / TMath::Tan(qqqTheta) + source_vertex;
        plotter->Fill2D("pczguess_vs_pc_int", 400, -200, 200, 600, -400, 400, pcz_guess_int, pcevent.pos.Z(), "PCZ_Recon");

        if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
        {
          double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
          TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
          TVector3 v = x2f - x1;
          double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
          TVector3 r_rhoMin_fix = x1 + t_minimum * v;

          double sinTheta_customV = TMath::Sin((qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta());
          plotter->Fill2D("dE3_E_CathodeQQQW_A1C2_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, qqqevent.Energy2, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
          plotter->Fill2D("dE3_E_AnodeQQQR_A1C2_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV, "PID_dE_E");

          plotter->Fill1D("VertexRecon_pczfix_qqq", 800, -300, 300, r_rhoMin_fix.Z(), "Vertex_Reconstruction");
          plotter->Fill1D("VertexRecon_pczfix", 800, -300, 300, r_rhoMin_fix.Z(), "Vertex_Reconstruction");
          plotter->Fill1D("VertexRecon_pczfix_qqq_PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 800, -400, 400, r_rhoMin_fix.Z(), "Vertex_Reconstruction");

          if (TMath::Abs(r_rhoMin_fix.Z()) < 200.0)
          {
            plotter->Fill2D("dE3_E_AnodeQQQR_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV, "PID_dE_E");
            plotter->Fill2D("dE3_E_CathodeQQQR_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
          }

          plotter->Fill1D("pczfix_A1C2_1d_qqq", 600, -200, 200, pcz_fix, "PCZ_Recon");
          plotter->Fill2D("pczfix_vs_qqqpczguess_A1C2", 600, -200, 200, 600, -200, 200, pcz_guess_int, pcz_fix, "PCZ_Recon");
          plotter->Fill2D("pczguess_vs_pc_int_A1C2", 400, -200, 200, 600, -400, 400, pcz_guess_int, pcevent.pos.Z(), "PCZ_Recon");
          plotter->Fill2D("pczfix_residual_vs_pczguess_A1C2", 600, -200, 200, 200, -100, 100, pcz_guess_37, pcz_fix - pcz_guess_37, "Residuals");
          plotter->Fill2D("pczfix_residual_vs_phi_A1C2", 200, 0, 6.28, 200, -100, 100, r_rhoMin_fix.Phi(), pcz_fix - pcz_guess_37, "Residuals");
          plotter->Fill1D("pczfix-qqqpczguess_A1C2", 200, -100, 100, pcz_fix - pcz_guess_37, "Residuals");
          plotter->Fill1D("pczfix-qqqpczint_A1C2", 200, -100, 100, pcz_fix - pcz_guess_int, "Residuals");
          plotter->Fill1D("pczguess_vs_int_residualsqqq", 200, -50, 50, pcz_guess_int - pcz_guess_37, "Residuals");
          plotter->Fill2D("pcdEA_vs_pczfix", 800, 0, 20000, 600, -200, 200, pcevent.Energy1, pcz_fix, "Kinematics_Angles");
          plotter->Fill2D("pcdEC_vs_pczfix", 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pcz_fix, "Kinematics_Angles");

          double path_length = (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Mag() * 0.1;
          double qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length);
          double qqqEfix_p = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1) - path_length);

          plotter->Fill2D("E_thetaf_AnodeQQQR_TC1_PC" + std::to_string(phicut), 180, 0, 180, 600, 0, 15, (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta() * 180 / M_PI, qqqevent.Energy1, "Kinematics_Angles");
          if (lowercut_cath)
            plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 180, 0, 180, 600, 0, 15, (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta() * 180 / M_PI, qqqEfix_p, "Kinematics_Angles");
          else
          {
            std::string zcut = "_" + std::to_string((TMath::Abs(r_rhoMin_fix.Z()) < 180));
            plotter->Fill2D("Ef_thetaf_AnodeQQQR_TC1_PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath) + zcut, 180, 0, 180, 600, 0, 15, (qqqevent.pos - TVector3(0, 0, r_rhoMin_fix.Z())).Theta() * 180 / M_PI, qqqEfix, "Kinematics_Angles");
          }

          plotter->Fill2D("dE3_Ef_AnodeQQQR_TC1" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 600, 0, 15, 800, 0, 40000, qqqEfix, pcevent.Energy1 * sinTheta_customV, "PID_dE_E");
          plotter->Fill2D("dE3_Ef_CathodeQQQR_TC1PC" + std::to_string(phicut) + "_pidlow" + std::to_string(lowercut_cath), 600, 0, 15, 800, 0, 10000, qqqEfix, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
        }

        //-----------------------Benchmarking Method for Source Runs (QQQ)------------------------//
        if (BenchMark && aClusters.size() == 1 && cClusters.size() == 1)
        {
          const auto &aCl = aClusters.front();
          const auto &cCl = cClusters.front();
          const std::string benchBranch = "Benchmark_QQQ";

          auto vertexFrom = [](const TVector3 &si, const TVector3 &pcpoint)
          {
            TVector3 vf = pcpoint - si;
            double tm = -1.0 * (si.X() * vf.X() + si.Y() * vf.Y()) / (vf.X() * vf.X() + vf.Y() * vf.Y());
            return TVector3(si + tm * vf);
          };

          auto fillSuite = [&](const std::string &tag, double pcz_method, const TVector3 &vtx, const std::string &branch)
          {
            plotter->Fill1D("Benchmark_QQQ_VertexZ_" + tag, 800, -400, 400, vtx.Z(), branch);
            plotter->Fill1D("Benchmark_QQQ_VertexZ_" + tag + "_TC" + std::to_string(timecut) + "_PC" + std::to_string(phicut), 800, -400, 400, vtx.Z(), branch);
            plotter->Fill2D("Benchmark_QQQ_VertexXY_" + tag, 200, -100, 100, 200, -100, 100, vtx.X(), vtx.Y(), branch);
            plotter->Fill1D("Benchmark_QQQ_PCZ_" + tag, 600, -200, 200, pcz_method, branch);
          };

          auto fillVsRef = [&](const std::string &tag, double pcz_method, const TVector3 &vtx, double pcz_ref, const TVector3 &vtx_ref)
          {
            plotter->Fill2D("Benchmark_QQQ_PCZ_" + tag + "_vs_ref", 400, -200, 200, 400, -200, 200, pcz_ref, pcz_method, "Benchmark_QQQ_ref");
            plotter->Fill1D("Benchmark_QQQ_PCZ_" + tag + "_minus_ref", 400, -100, 100, pcz_method - pcz_ref, "Benchmark_QQQ_ref");
          };

          double pcz_ref = pcfix_func.Eval(pcevent.pos.Z());
          TVector3 vtx_ref = vertexFrom(qqqevent.pos, TVector3(pcevent.pos.X(), pcevent.pos.Y(), pcz_ref));

          auto pw_tuple = pwinstance.GetPseudoWire(aCl, "ANODE");
          std::pair<TVector3, TVector3> apwire_bm = std::get<0>(pw_tuple);

          auto cMaxWire = *std::max_element(cCl.begin(), cCl.end(), [](const auto &a, const auto &b)
                                            { return std::get<1>(a) < std::get<1>(b); });
          auto aMaxWire = *std::max_element(aCl.begin(), aCl.end(), [](const auto &a, const auto &b)
                                            { return std::get<1>(a) < std::get<1>(b); });
          std::vector<std::tuple<int, double, double>> cOne = {cMaxWire};

          auto xo_tuple = pwinstance.FindCrossoverProperties(aCl, cOne);
          TVector3 xo_a1c1 = std::get<0>(xo_tuple);
          double alpha_a1c1 = std::get<1>(xo_tuple);
          bool a1c1Good = (alpha_a1c1 != 9999999 && std::get<2>(xo_tuple) != -1);

          // --- A1C1 charge-fraction diagnostic ---
          // One cathode in A1C1, so the z-sensitive variable is cathode-vs-anode
          // charge (the V across each cell), not cathode/cathode. Use pseudo-wire
          // sums for gain-consistent energies.
          double aSumE_bm = std::get<1>(pw_tuple); // anode sum (reuse pw_tuple)
          double cSumE_bm = std::get<1>(cMaxWire); // Extract the energy directly!
          double ac_sum = aSumE_bm + cSumE_bm;
          double cfrac = (ac_sum > 0.0) ? cSumE_bm / ac_sum : -1.0; // bounded [0,1]

          // Cmax/anode vs Si phi for ALL events reaching this block (not just the
          // A1C2 subset). Cmax = single max-E cathode wire (cSumE_bm), anode =
          // anode pseudo-wire sum. Ratio is unbounded above (arbitrary cathode gains).
          if (aSumE_bm > 0.0)
            plotter->Fill2D("Benchmark_QQQ_CmaxOverAnode_vs_phi", 180, -180, 180, 250, 0, 5,
                            qqqevent.pos.Phi() * 180. / M_PI, cSumE_bm / aSumE_bm, "Benchmark_QQQ_ref");

          // Smearing Setup
          double qqq_wedge_pitch = (87.0 / 16.0) * (M_PI / 180.0);
          double qqq_ring_pitch = 48.0 / 16.0;
          double smeared_phi = qqqevent.pos.Phi() + rand.Uniform(-qqq_wedge_pitch / 2.0, qqq_wedge_pitch / 2.0);
          double smeared_rho = qqqevent.pos.Perp() + rand.Uniform(-qqq_ring_pitch / 2.0, qqq_ring_pitch / 2.0);

          TVector3 smeared_qqq_pos(smeared_rho * TMath::Cos(smeared_phi), smeared_rho * TMath::Sin(smeared_phi), qqqevent.pos.Z());

          // Plain (dither=false) shows the discrete wire spikes; _Hyb (Si phi smeared
          // + PC z dithered) is the smoothed cross-reference.
          auto doA1C1 = [&](const std::string &tag, const TVector3 &si_point, bool dither = true)
          {
            if (!a1c1Good)
              return;
            double pcz = dither ? rand.Gaus(xo_a1c1.Z(), dither_sigma) : xo_a1c1.Z();
            TVector3 vtx = vertexFrom(si_point, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz));
            fillSuite(tag, pcz, vtx, benchBranch);
            fillVsRef(tag, pcz, vtx, pcz_ref, vtx_ref);
          };

          auto doAnodeOnly = [&](const std::string &tag, double phi_use, const TVector3 &si_point, bool dither = true)
          {
            TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, phi_use);
            pc.SetZ(a1c1_zcorr(pc.Z(), true));
            TVector3 vtx0 = vertexFrom(si_point, pc);
            if (!(vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6))
              return;
            double pcz = dither ? rand.Gaus(pc.Z(), dither_sigma_c0 / 2.0) : pc.Z();
            TVector3 vtx = vertexFrom(si_point, TVector3(pc.X(), pc.Y(), pcz));
            fillSuite(tag, pcz, vtx, benchBranch);
            fillVsRef(tag, pcz, vtx, pcz_ref, vtx_ref);
          };

          // --- A1C1 with the cfrac sub-cell model (QQQ, linear centre-fold) ---
          // Anchored on the fired cathode (xo_a1c1.Z()); cell from the anode-only
          // z (z_a1c0, independent of the Si guess); |offset| from cfrac. REJECT
          // when out of band or inconsistent with the fired wire.
          auto doA1C1Model = [&](const std::string &tag, const TVector3 &si_point)
          {
            if (!a1c1Good || cfrac < 0.0)
              return;
            A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
            // side from the beam-axis 2-hypothesis test (the Si hit picks the cell)
            int side_status;
            SideChoice side = a1c1_pick_side(si_point, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
            const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
            double pcz_pick = best.pcz;
            if (!(best.inband && side_status != 2))
              return;
            TVector3 vtx = vertexFrom(si_point, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_pick));
            fillSuite(tag, pcz_pick, vtx, benchBranch);
            fillVsRef(tag, pcz_pick, vtx, pcz_ref, vtx_ref);
          };

          if (phicut && timecut)
          {
            if (pcevent.multi1 == 1 && pcevent.multi2 == 2)
            {
              fillSuite("A1C2", pcz_ref, vtx_ref, benchBranch);
              // XY-offset diagnostic (QQQ twin of SX3 block above)
              {
                double phi_deg = qqqevent.pos.Phi() * 180.0 / M_PI;
                double vz_resid = vtx_ref.Z() - source_vertex;
                plotter->Fill2D("Diag_QQQ_A1C2_vtxZ_resid_vs_phi", 180, -180, 180, 400, -100, 100,
                                phi_deg, vz_resid, "Diag_XYoffset");
                plotter->Fill2D("Diag_Combined_A1C2_vtxZ_resid_vs_phi", 180, -180, 180, 400, -100, 100,
                                phi_deg, vz_resid, "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_vtxXY", 200, -15, 15, 200, -15, 15,
                                vtx_ref.X(), vtx_ref.Y(), "Diag_XYoffset");
                plotter->Fill2D("Diag_Combined_A1C2_time_vs_phi", 2000, 0, 2000, 180, -180, 180,
                                pcevent.Time1 * 1e-9, phi_deg, "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_time_vs_phi", 2000, 0, 2000, 180, -180, 180,
                                pcevent.Time1 * 1e-9, phi_deg, "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_T_vs_vtxX", 2000, 0, 2000, 200, -15, 15,
                                pcevent.Time1 * 1e-9, vtx_ref.X(), "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_T_vs_vtxY", 2000, 0, 2000, 200, -15, 15,
                                pcevent.Time1 * 1e-9, vtx_ref.Y(), "Diag_XYoffset");
              }

              doA1C1("A1C1", qqqevent.pos, false);
              doAnodeOnly("A1C0", qqqevent.pos.Phi(), qqqevent.pos, false);
              doA1C1("A1C1_Hyb", smeared_qqq_pos);
              doAnodeOnly("A1C0_Hyb", smeared_phi, smeared_qqq_pos);

              // --- Execute the cfrac-model for QQQ ---
              doA1C1Model("A1C1_Cfrac", qqqevent.pos);

              // A1C0 anode-only PCz bias vs track theta (see SX3 twin): if QQQ and SX3
              // fall on the same curve, the z_a1c0 offset is theta-driven. theta from vtx_ref.
              {
                double pcz_a1c0 = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, qqqevent.pos.Phi()).Z();
                double theta_ref = (qqqevent.pos - TVector3(0, 0, vtx_ref.Z())).Theta() * 180. / M_PI;
                plotter->Fill2D("Benchmark_QQQ_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200,
                                theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_QQQ_ref");
                plotter->Fill2D("Benchmark_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200,
                                theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_AnodeOnly");

                double phi_deg_a = qqqevent.pos.Phi() * 180.0 / M_PI;
                plotter->Fill2D("Diag_QQQ_A1C0_zresid_vs_phi", 180, -180, 180, 200, -100, 100,
                                phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
                plotter->Fill2D("Diag_Combined_A1C0_zresid_vs_phi", 180, -180, 180, 200, -100, 100,
                                phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
              }

              // --- A1C1 charge-fraction diagnostics ---
              if (a1c1Good && cfrac >= 0.0)
              {
                plotter->Fill1D("Benchmark_QQQ_A1C1_cfrac", 220, -0.05, 1.05, cfrac, "Benchmark_QQQ_ref");
                plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_ref", 400, -200, 200, 220, -0.05, 1.05,
                                pcz_ref, cfrac, "Benchmark_QQQ_ref");
                plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_qqqpczguess", 400, -200, 200, 220, -0.05, 1.05,
                                pcz_guess_37, cfrac, "Benchmark_QQQ_ref");

                // --- cfrac vs normalized cell position (collapses all cells) ---
                static const double zg[8] = {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};
                double zp = xo_a1c1.Z();
                auto fillCfracS = [&](const char *name, double truth)
                {
                  double sgn = (truth >= zp) ? 1.0 : -1.0;
                  double znb = (sgn > 0) ? 1.0e30 : -1.0e30;
                  for (int i = 0; i < 8; ++i)
                  {
                    if (sgn > 0 && zg[i] > zp + 1e-6 && zg[i] < znb)
                      znb = zg[i];
                    if (sgn < 0 && zg[i] < zp - 1e-6 && zg[i] > znb)
                      znb = zg[i];
                  }
                  if (TMath::Abs(znb) < 1e8 && TMath::Abs(znb - zp) > 0.0)
                    plotter->Fill2D(name, 240, -1.2, 1.2, 220, -0.05, 1.05,
                                    (truth - zp) / TMath::Abs(znb - zp), cfrac, "Benchmark_QQQ_ref");
                };
                fillCfracS("Benchmark_QQQ_A1C1_cfrac_vs_s", pcz_ref);
                fillCfracS("Benchmark_QQQ_A1C1_cfrac_vs_s_qqqpczguess", pcz_guess_37);

                // --- cfrac vs cell-centre fold: f = |ref - z_center|/halfcell in [0,1] ---
                for (int i = 0; i < 7; ++i)
                {
                  if (pcz_ref <= zg[i] && pcz_ref > zg[i + 1])
                  {
                    double zc = 0.5 * (zg[i] + zg[i + 1]);
                    double half = 0.5 * (zg[i] - zg[i + 1]);
                    if (half > 0.0)
                      plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05,
                                      TMath::Abs(pcz_ref - zc) / half, cfrac, "Benchmark_QQQ_ref");
                    break;
                  }
                }

                // --- item 2: cfrac vs anode E (pure position variable check) ---
                plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05,
                                aSumE_bm, cfrac, "Benchmark_QQQ_ref");
                // r-space linearity test (r = c + C_off/anodeE should be a straight line)
                if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                  plotter->Fill2D("Benchmark_QQQ_A1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0,
                                  1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "Benchmark_QQQ_ref");

                // --- item 3: cell/side misclassification rate (truth = pcz_ref) ---
                {
                  A1C1Sol sm = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
                  int sm_side_status;
                  SideChoice sm_side = a1c1_pick_side(qqqevent.pos, xo_a1c1.X(), xo_a1c1.Y(), sm.pcz_lo, sm.pcz_hi, sm_side_status);
                  int sm_cell = (sm_side == SideChoice::High) ? sm.hi.cell : sm.lo.cell;
                  int cell_truth = -1;
                  for (int i = 0; i < 7; ++i)
                    if (pcz_ref <= a1c1_zg[i] && pcz_ref > a1c1_zg[i + 1])
                    {
                      cell_truth = i;
                      break;
                    }
                  if (cell_truth >= 0)
                  {
                    bool wrong = (sm_cell != cell_truth);
                    plotter->Fill2D("Benchmark_QQQ_A1C1_cellsel_confusion", 7, 0, 7, 7, 0, 7,
                                    cell_truth + 0.5, sm_cell + 0.5, "Benchmark_QQQ_ref");
                    plotter->Fill1D("Benchmark_QQQ_A1C1_cellsel_misclass", 2, 0, 2, wrong ? 1.0 : 0.0, "Benchmark_QQQ_ref");
                    plotter->Fill2D("Benchmark_QQQ_A1C1_cellsel_misclass_vs_cell", 7, 0, 7, 2, 0, 2,
                                    cell_truth + 0.5, wrong ? 1.0 : 0.0, "Benchmark_QQQ_ref");
                    double zc = 0.5 * (a1c1_zg[cell_truth] + a1c1_zg[cell_truth + 1]);
                    double half = 0.5 * (a1c1_zg[cell_truth] - a1c1_zg[cell_truth + 1]);
                    if (half > 0.0)
                    {
                      plotter->Fill2D("Benchmark_QQQ_A1C1_cellsel_misclass_vs_fold", 120, 0, 1.2, 2, 0, 2,
                                      TMath::Abs(pcz_ref - zc) / half, wrong ? 1.0 : 0.0, "Benchmark_QQQ_ref");
                      // item 5 validation: post-fold cfrac vs truth fold (bands merge)
                      plotter->Fill2D("Benchmark_QQQ_A1C1_cfracUsed_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05,
                                      TMath::Abs(pcz_ref - zc) / half, sm.cfrac_used, "Benchmark_QQQ_ref");
                      // anode-E correction validation: post-correction cfrac vs anodeE
                      // should be FLAT.
                      if (aSumE_bm > 0.0)
                        plotter->Fill2D("Benchmark_QQQ_A1C1_cfracUsed_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05,
                                        aSumE_bm, sm.cfrac_used, "Benchmark_QQQ_ref");
                    }
                  }
                }
              }
            }
          }

          // ---- TRUE A1C1 (QQQ): genuine one-anode/one-cathode events. No A1C2
          // crossover reference exists, so compare only to the geometric guess
          // (pcz_guess_int). xo_a1c1 is the real measured crossover. ----
          else if (pcevent.multi1 >= 1 && pcevent.multi2 == 1 && a1c1Good)
          {
            double pcz_raw = xo_a1c1.Z();
            TVector3 vtx_raw = vertexFrom(qqqevent.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_raw));
            fillSuite("trueA1C1", pcz_raw, vtx_raw, "A1C1True_QQQ");
            plotter->Fill2D("Benchmark_QQQ_PCZ_trueA1C1_vs_qqqpczguess", 400, -200, 200, 400, -200, 200, pcz_guess_int, pcz_raw, "Benchmark_QQQ_trueA1C1");
            plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C1_minus_qqqpczguess", 400, -100, 100, pcz_raw - pcz_guess_int, "Benchmark_QQQ_trueA1C1");

            // cfrac sub-cell reconstruction, anchored on the FIRED CATHODE.
            // Vertex-independent (no pcz_guess_int): the cell is the one adjacent
            // to the cathode that fired (zf), the side is taken from the anode z,
            // and the hit must land within one cell pitch of zf to be valid.
            if (cfrac >= 0.0)
            {
              A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
              // side from the beam-axis 2-hypothesis test (smaller-Perp vertex wins)
              int side_status;
              SideChoice side = a1c1_pick_side(qqqevent.pos, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
              const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
              int cell = best.cell;
              double f = best.f;
              double pcz_cf = best.pcz;
              bool valid = (side_status != 2); // accept if at least one side is beam-axis consistent
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_sideStatus", 4, -1, 3, side_status + 0.5, "Benchmark_QQQ_trueA1C1");

              TVector3 vtx_cf = vertexFrom(qqqevent.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_cf));
              fillSuite(valid ? "trueA1C1_Cfrac" : "trueA1C1_Cfrac_invalid", pcz_cf, vtx_cf, "A1C1True_QQQ");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_cfrac", 220, -0.05, 1.05, cfrac, "Benchmark_QQQ_trueA1C1");
              // item 2: cfrac vs anode E for genuine A1C1 (no A1C2 ref here)
              plotter->Fill2D("Benchmark_QQQ_trueA1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05, aSumE_bm, cfrac, "Benchmark_QQQ_trueA1C1");
              // r-space linearity test (r = c + C_off/anodeE should be a straight line)
              if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                plotter->Fill2D("Benchmark_QQQ_trueA1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0,
                                1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "Benchmark_QQQ_trueA1C1");
              // plotter->Fill1D("Benchmark_QQQ_trueA1C1_cfrac_cathode" + std::to_string(pcevent.Cathodech), 220, -0.05, 1.05, cfrac, "Benchmark_trueA1C1_cathode");
              // reference-free per-cell cfrac (cell from geometry, no A1C2 ref): the
              // offline fitter reads per-cell edges/percentiles to gain-match cfmin/k.
              plotter->Fill2D("Benchmark_QQQ_trueA1C1_cfrac_vs_cell", 7, 0, 7, 220, -0.05, 1.05, cell + 0.5, cfrac, "Benchmark_QQQ_trueA1C1");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_f", 260, -1.5, 2.5, f, "Benchmark_QQQ_trueA1C1");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_valid", 2, 0, 2, valid ? 1.0 : 0.0, "Benchmark_QQQ_trueA1C1");
              // failure-reason breakdown (why the cfrac estimate is / isn't usable):
              //   0 valid & f in [0,1]   ideal, inside the calibrated band
              //   1 valid & f in [-1,0)  below band but within the pitch tolerance
              //   2 valid & f in (1,3]   above band but within the pitch tolerance
              //   3 invalid: f < -1      cfrac far below band (cfmin too high / cathode low?)
              //   4 invalid: f > 3       cfrac far above band (cfmin too low / cathode high?)
              //   5 cell k<=0            cell not autocalibrated
              int reason;
              if (a1c1_k_cell[cell] <= 0.0)
                reason = 5;
              else if (!valid)
                reason = (f < 0.0) ? 3 : 4;
              else if (f < 0.0)
                reason = 1;
              else if (f > 1.0)
                reason = 2;
              else
                reason = 0;
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_failreason", 6, 0, 6, reason + 0.5, "Benchmark_QQQ_trueA1C1");
              // valid-only breakdown (reason is 0/1/2 here): how clean the accepted
              // sample is -- 0 ideal (f in [0,1]), 1 below-band, 2 above-band marginal.
              if (valid)
                plotter->Fill1D("Benchmark_QQQ_trueA1C1_validreason", 3, 0, 3, reason + 0.5, "Benchmark_QQQ_trueA1C1");
              // which cfrac band reconstructed this event: 0 = main, 1 = low (incomplete integration)
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_band", 2, 0, 2, s.band + 0.5, "Benchmark_QQQ_trueA1C1");
              if (valid)
                plotter->Fill1D("Benchmark_QQQ_trueA1C1_band_valid", 2, 0, 2, s.band + 0.5, "Benchmark_QQQ_trueA1C1");
              // source-run-only diagnostics: pcz_guess_int is meaningful only when the vertex is fixed
              if (valid)
              {
                plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C1_Cfrac_minus_qqqpczguess_DIAG", 400, -100, 100, pcz_cf - pcz_guess_int, "Benchmark_QQQ_trueA1C1");
                plotter->Fill2D("Benchmark_QQQ_PCZ_trueA1C1_Cfrac_vs_qqqpczguess_DIAG", 400, -200, 200, 400, -200, 200, pcz_guess_int, pcz_cf, "Benchmark_QQQ_trueA1C1");
              }
            }

            // anode-only (A1C0) estimate for the same events
            {
              TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, qqqevent.pos.Phi());
              TVector3 vtx0 = vertexFrom(qqqevent.pos, pc);
              if (vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6)
              {
                fillSuite("A1C1asA1C0", pc.Z(), vtx0, "A1C1True_QQQ");
                plotter->Fill2D("Benchmark_QQQ_PCZ_A1C1asA1C0_vs_qqqpczguess", 400, -200, 200, 400, -200, 200, pcz_guess_int, pc.Z(), "A1C1True_QQQ");
              }
            }
          }
        }
        double qqqrho = qqqevent.pos.Perp();
        double qqqz = (qqqevent.pos - TVector3(0, 0, source_vertex)).Z();
        double tan_theta = qqqrho / qqqz;
        double pcz_guess_int2 = z_to_crossover_rho(pcevent.pos.Z()) / tan_theta + source_vertex;
        plotter->Fill2D("pczguess_vs_pc_int2", 180, 0, 200, 150, 0, 200, pcz_guess_int2, pcevent.pos.Z(), "PCZ_Recon");

        double qqqz2 = (qqqevent.pos - r_rhoMin).Z();
        double tan_theta2 = qqqrho / qqqz2;
        double pcz_guess_int3 = z_to_crossover_rho(pcevent.pos.Z()) / tan_theta2 + r_rhoMin.Z();
        plotter->Fill2D("pczguess_vs_pc_int3", 180, 0, 200, 150, 0, 200, pcz_guess_int3, pcevent.pos.Z(), "PCZ_Recon");

        double pcz_guess = pcz_guess_int;
        plotter->Fill2D("pctheta_vs_qqqtheta_sv", 180, -200,200, 180, -200,200, qqqTheta * 180 / M_PI, (pcevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, "Kinematics_Angles");
        plotter->Fill2D("pctheta_vs_qqqtheta_rmz", 180, -200,200, 180, -200,200, (qqqevent.pos - TVector3(0, 0, r_rhoMin.Z())).Theta() * 180 / M_PI, (pcevent.pos - TVector3(0, 0, r_rhoMin.Z())).Theta() * 180 / M_PI, "Kinematics_Angles");
        plotter->Fill2D("pctheta_vs_qqqtheta_rm", 180, -200,200, 180, -200,200, (qqqevent.pos - r_rhoMin).Theta() * 180 / M_PI, (pcevent.pos - r_rhoMin).Theta() * 180 / M_PI, "Kinematics_Angles");
        plotter->Fill2D("pczguess_vs_pc_phi=" + std::to_string(qqqevent.pos.Phi() * 180. / M_PI), 300, 0, 200, 150, 0, 200, pcz_guess, pcevent.pos.Z(), "Z_Reconstruction");
      }
    }
  }
}

void TrackRecon::OldAnalysis()
{
  int aID = 0, cID = 0;
  double aE = 0, cE = 0;
  double aESum = 0, cESum = 0;
  double aEMax = 0, cEMax = 0;
  int aIDMax = 0, cIDMax = 0;

  if (anodeHits.size() >= 1 && cathodeHits.size() >= 1)
  {
    // 2. CRITICAL FIX: Define reference vector 'a'

    {
      for (const auto &anode : anodeHits)
      {
        aID = anode.first;
        aE = anode.second;
        aESum += aE;
        if (aE > aEMax)
        {
          aEMax = aE;
          aIDMax = aID;
        }
      }

      for (const auto &cathode : cathodeHits)
      {
        cID = cathode.first;
        cE = cathode.second;
        plotter->Fill2D("AnodeMax_Vs_Cathode_Coincidence_Matrix", 24, 0, 24, 24, 0, 24, aIDMax, cID, "hRawPC");
        plotter->Fill2D("Anode_Vs_Cathode_Coincidence_Matrix", 24, 0, 24, 24, 0, 24, aID, cID, "hRawPC");
        plotter->Fill2D("Anode_Vs_Cathode_Coincidence_Matrix_qqq" + std::to_string(HitNonZero), 24, 0, 24, 24, 0, 24, aID, cID, "hRawPC");
        plotter->Fill2D("Anode_vs_CathodeE", 2000, 0, 30000, 2000, 0, 30000, aE, cE, "hGMPC");
        plotter->Fill2D("CathodeMult_V_CathodeE", 6, 0, 6, 2000, 0, 30000, cathodeHits.size(), cE, "hGMPC");
        if (((aIDMax + cID) % 24) >= 20 || ((aIDMax + cID) % 24) <= 3)
        {
          corrcatMax.push_back(std::pair<int, double>(cID, cE));
          cESum += cE;
          if (cE > cEMax)
          {
            cEMax = cE;
            cIDMax = cID;
          }
        }
      }
    }
  }

  TVector3 anodeIntersection, vector_closest_to_z;
  anodeIntersection.Clear();
  vector_closest_to_z.Clear();
  if (corrcatMax.size() > 0)
  {
    double x = 0, y = 0, z = 0;
    for (const auto &corr : corrcatMax)
    {
      if (pwinstance.Crossover[aIDMax][corr.first][0].z > 9000000)
        continue;
      if (cESum > 0)
      {
        x += (corr.second) / cESum * pwinstance.Crossover[aIDMax][corr.first][0].x;
        y += (corr.second) / cESum * pwinstance.Crossover[aIDMax][corr.first][0].y;
        z += (corr.second) / cESum * pwinstance.Crossover[aIDMax][corr.first][0].z;
      }
    }
    if (x == 0 && y == 0 && z == 0)
      ;
    // to ignore events with no valid crossover points
    else
    {
      anodeIntersection = TVector3(x, y, z);
      // std::cout << "Anode Intersection: " << anodeIntersection.X() << ", " << anodeIntersection.Y() << ", " << anodeIntersection.Z() << " " << aIDMax <<  std::endl;
    }
  }
  bool PCQQQPhiCut = false;
  // flip the algorithm for cathode 1 multi anode events
  if ((hitPos.Phi() > (anodeIntersection.Phi() - TMath::PiOver4())) && (hitPos.Phi() < (anodeIntersection.Phi() + TMath::PiOver4())))
  {
    PCQQQPhiCut = true;
  }

  if (anodeIntersection.Z() != 0 && anodeIntersection.Perp() > 0 && HitNonZero)
  {
    plotter->Fill1D("PC_Z_Projection", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("Z_Proj_VsDelTime", 600, -300, 300, 200, -2000, 2000, anodeIntersection.Z(), anodeT - cathodeT, "hPCzQQQ");
    plotter->Fill2D("IntPhi_vs_QQQphi", 100, -200, 200, 80, -200, 200, anodeIntersection.Phi() * 180. / TMath::Pi(), hitPos.Phi() * 180. / TMath::Pi(), "hPCQQQ");
    // plotter->Fill2D("Inttheta_vs_QQQtheta", 90, 0, 180, 20, 0, 45, anodeIntersection.Theta() * 180. / TMath::Pi(), hitPos.Theta() * 180. / TMath::Pi(), "hPCQQQ");
    // plotter->Fill2D("Inttheta_vs_QQQtheta_TC" + std::to_string(PCQQQTimeCut)+ "_PC"+std::to_string(PCQQQPhiCut), 90, 0, 180, 20, 0, 45, anodeIntersection.Theta() * 180. / TMath::Pi(), hitPos.Theta() * 180. / TMath::Pi(), "hPCQQQ");
    plotter->Fill2D("IntPhi_vs_QQQphi_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 100, -200, 200, 80, -200, 200, anodeIntersection.Phi() * 180. / TMath::Pi(), hitPos.Phi() * 180. / TMath::Pi(), "hPCQQQ");
  }

  if (anodeIntersection.Z() != 0 && anodeIntersection.Perp() > 0 && PCSX3TimeCut)
  {
    plotter->Fill1D("PC_Z_Projection_sx3", 600, -200, 200, anodeIntersection.Z(), "hPCZSX3");
  }
  if (anodeIntersection.Z() != 0 && cathodeHits.size() >= 2)
    plotter->Fill1D("PC_Z_Projection_TC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");

  if (anodeIntersection.Z() != 0 && cathodeHits.size() == 1)
  {
    plotter->Fill1D("PC_Z_proj_1C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("IntersectionPhi_vs_AnodeZ_1C", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hPCzQQQ");
  }

  if (anodeIntersection.Z() != 0 && cathodeHits.size() == 2)
  {
    plotter->Fill1D("PC_Z_proj_2C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("IntersectionPhi_vs_AnodeZ_2C", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hGMPC");
  }
  if (anodeIntersection.Z() != 0 && cathodeHits.size() > 2)
  {
    plotter->Fill1D("PC_Z_proj_nC", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
    plotter->Fill2D("IntersectionPhi_vs_AnodeZ_nC", 400, -200, 200, 600, -300, 300, anodeIntersection.Phi() * 180. / TMath::Pi(), anodeIntersection.Z(), "hGMPC");
  }
  if (anodeHits.size() > 0 && cathodeHits.size() > 0)
    plotter->Fill2D("AHits_vs_CHits", 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");

  // make another plot with nearest neighbour constraint
  bool hasNeighbourAnodes = false;
  bool hasNeighbourCathodes = false;

  // 1. Check Anodes for neighbours (including wrap-around 0-23)
  for (size_t i = 0; i < anodeHits.size(); i++)
  {
    for (size_t j = i + 1; j < anodeHits.size(); j++)
    {
      int diff = std::abs(anodeHits[i].first - anodeHits[j].first);
      if (diff == 1 || diff == 23)
      { // 23 handles the cylindrical wrap
        hasNeighbourAnodes = true;
        break;
      }
    }
    if (hasNeighbourAnodes)
      break;
  }

  // 2. Check Cathodes for neighbours (including wrap-around 0-23)
  for (size_t i = 0; i < cathodeHits.size(); i++)
  {
    for (size_t j = i + 1; j < cathodeHits.size(); j++)
    {
      int diff = std::abs(cathodeHits[i].first - cathodeHits[j].first);
      if (diff == 1 || diff == 23)
      {
        hasNeighbourCathodes = true;
        break;
      }
    }
    if (hasNeighbourCathodes)
      break;
  }

  // ---------------------------------------------------------
  // FILL PLOTS
  // ---------------------------------------------------------
  if (anodeHits.size() > 0 && cathodeHits.size() > 0)
  {
#ifdef RAW_HISTOS
    plotter->Fill2D("AHits_vs_CHits_NA" + std::to_string(hasNeighbourAnodes), 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");
    plotter->Fill2D("AHits_vs_CHits_NC" + std::to_string(hasNeighbourCathodes), 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");

    // Constraint Plot: Only fill if BOTH planes have adjacent hits
    // This effectively removes events with only isolated single-wire hits (noise)
    if (hasNeighbourAnodes && hasNeighbourCathodes)
    {
      plotter->Fill2D("AHits_vs_CHits_NN", 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");
    }
#endif
  }

  if (HitNonZero && anodeIntersection.Z() != 0)
  {
    pwinstance.CalTrack2(hitPos, anodeIntersection);
    plotter->Fill1D("VertexRecon", 600, -1300, 1300, pwinstance.GetZ0());
    plotter->Fill1D("VertexRecon_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, pwinstance.GetZ0());

    if (cathodeHits.size() == 2)
      plotter->Fill1D("VertexRecon_2c_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, pwinstance.GetZ0());

    TVector3 x2(anodeIntersection), x1(hitPos);

    TVector3 v = x2 - x1;
    double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
    vector_closest_to_z = x1 + t_minimum * v;

    plotter->Fill1D("VertexRecon_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");

    if (qqqenergy < 4.0)
      plotter->Fill1D("VertexRecon_Z(qqqE<4.0MeV)_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");

    if (vector_closest_to_z.Perp() < 20)
    {
      plotter->Fill1D("VertexRecon_RadialCut_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");
    }

    plotter->Fill2D("VertexRecon_XY_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 100, -100, 100, 100, -100, 100, vector_closest_to_z.X(), vector_closest_to_z.Y(), "customVertex");
    if (cathodeHits.size() == 2)
    {
      plotter->Fill1D("VertexRecon2C_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");
      if (vector_closest_to_z.Perp() < 20)
      {
        plotter->Fill1D("VertexRecon2C_RadialCut_Z_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, vector_closest_to_z.Z(), "customVertex");
      }
      plotter->Fill2D("VertexRecon2C_XY_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 100, -100, 100, 100, -100, 100, vector_closest_to_z.X(), vector_closest_to_z.Y(), "customVertex");
      plotter->Fill2D("VertexRecon2C_RhoZ_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 100, -100, 100, 600, -1300, 1300, vector_closest_to_z.Perp(), vector_closest_to_z.Z(), "customVertex");
      plotter->Fill2D("VertexRecon2C_Z_vs_QQQE_TC" + std::to_string(PCQQQTimeCut) + "_PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, 800, 0, 20, vector_closest_to_z.Z(), qqqenergy, "customVertex");
    }
  }

  for (int i = 0; i < qqq.multi; i++)
  {
    if (anodeIntersection.Perp() > 0)
    { // suppress x,y=0,0 events
      if (PCQQQTimeCut)
      {
        plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");
        plotter->Fill2D("PC_XY_Projection_QQQ_TimeCut" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, hitPos.X(), hitPos.Y(), "hPCQQQ");
      }
      plotter->Fill2D("PC_XY_Projection_QQQ" + std::to_string(qqq.id[i]), 400, -100, 100, 400, -100, 100, anodeIntersection.X(), anodeIntersection.Y(), "hPCQQQ");
    }
    for (int j = i + 1; j < qqq.multi; j++)
    {
      if (qqq.id[i] == qqq.id[j])
      {
        int chWedge = -1;
        int chRing = -1;
        double eWedge = 0.0;
        double eWedgeMeV = 0.0;
        double eRing = 0.0;
        double eRingMeV = 0.0;
        double tRing = 0.0;
        int qqqID = -1;
        if (qqq.ch[i] < 16 && qqq.ch[j] >= 16 && qqqGainValid[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16])
        {
          chWedge = qqq.ch[i];
          eWedge = qqq.e[i] * qqqGain[qqq.id[i]][qqq.ch[i]][qqq.ch[j] - 16];
          chRing = qqq.ch[j] - 16;
          eRing = qqq.e[j];
          tRing = static_cast<double>(qqq.t[j]);
          qqqID = qqq.id[i];
        }
        else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
        {
          chWedge = qqq.ch[j];
          eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
          chRing = qqq.ch[i] - 16;
          tRing = static_cast<double>(qqq.t[i]);
          eRing = qqq.e[i];
          qqqID = qqq.id[i];
        }
        else
          continue;

        if (qqqCalibValid[qqq.id[i]][chWedge][chRing])
        {
          eWedgeMeV = eWedge * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;
          eRingMeV = eRing * qqqCalib[qqq.id[i]][chWedge][chRing] / 1000;
        }
        else
          continue;

        // if (anodeIntersection.Z() != 0)
        {
          plotter->Fill2D("PC_Z_vs_QQQRing", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCzQQQ");
          plotter->Fill2D("PC_Z_vs_QQQRho", 600, -300, 300, 40, 40, 110, anodeIntersection.Z(), hitPos.Perp(), "hPCzQQQ");
        }

        if (anodeIntersection.Z() != 0 && cathodeHits.size() == 2)
        {
          plotter->Fill2D("PC_Z_vs_QQQRing_2C", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCzQQQ");
          plotter->Fill2D("PC_Z_vs_QQQRing_2C" + std::to_string(qqq.id[i]), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCzQQQ");
          plotter->Fill2D("PC_Z_vs_QQQWedge_2C", 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chWedge, "hPCzQQQ");
        }
        plotter->Fill2D("VertexRecon_QQQRingTC" + std::to_string(PCQQQTimeCut) + "PhiC" + std::to_string(PCQQQPhiCut), 600, -1300, 1300, 16, 0, 16, vector_closest_to_z.Z(), chRing, "hPCQQQ");
        double phi = TMath::ATan2(anodeIntersection.Y(), anodeIntersection.X()) * 180. / TMath::Pi();
        plotter->Fill2D("PolarAngle_Vs_QQQWedge" + std::to_string(qqqID), 360, -200,200, 16, 0, 16, phi, chWedge, "hPCQQQ");
        // plotter->Fill2D("EdE_PC_vs_QQQ_timegate_ls1000"+std::to_string())

        plotter->Fill2D("PC_Z_vs_QQQRing_Det" + std::to_string(qqqID), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCQQQ");
        // double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
        // double rho = 50. + 40. / 16. * (chRing + 0.5);

        for (int k = 0; k < pc.multi; k++)
        {
          if (pc.index[k] >= 24)
            continue;

          //          double sinTheta = TMath::Sin((hitPos-vector_closest_to_z).Theta());
          double sinTheta = TMath::Sin((anodeIntersection - TVector3(0, 0, 90.0)).Theta());
          //           double sinTheta = TMath::Sin((anodeIntersection-vector_closest_to_z).Theta());
          //          double sinTheta = TMath::Sin((hitPos-TVector3(0,0,30.0)).Theta());
          //          double sinTheta = TMath::Sin(hitPos.Theta());

          if (cathodeHits.size() == 2 && PCQQQPhiCut)
          {
            plotter->Fill2D("CalibratedQQQE_RvsCPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eRingMeV, pc.e[k] * sinTheta, "hPCQQQ");
            plotter->Fill2D("CalibratedQQQE_WvsCPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eWedgeMeV, pc.e[k] * sinTheta, "hPCQQQ");
            plotter->Fill2D("CalibratedQQQE_RvsPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
            plotter->Fill2D("CalibratedQQQE_WvsPCE_TC" + std::to_string(PCQQQTimeCut), 400, 0, 10, 400, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
            plotter->Fill2D("PCQQQ_dTimevsdPhi", 200, -2000, 2000, 80, -200, 200, tRing - static_cast<double>(pc.t[k]), (hitPos.Phi() - anodeIntersection.Phi()) * 180. / TMath::Pi(), "hTiming");
          }
        }
      } /// qqq i==j case end
    } // j loop end
  } // qqq i loop end

  for (int i = 0; i < sx3.multi; i++)
  {
    // plotting sx3 strip hits vs anode phi
    if (sx3.ch[i] < 8 && anodeIntersection.Perp() > 0)
      plotter->Fill2D("PCPhi_vs_SX3Strip", 100, -200, 200, 8 * 24, 0, 8 * 24, anodeIntersection.Phi() * 180. / TMath::Pi(), sx3.id[i] * 8 + sx3.ch[i]);
  }

  if (anodeIntersection.Z() != 0 && cathodeHits.size() == 3)
  {
    plotter->Fill1D("PC_Z_proj_3C", 600, -300, 300, anodeIntersection.Z(), "hPCzQQQ");
  }

  if (anodeIntersection.Perp() != 0)
  {
    plotter->Fill2D("AnodeMaxE_Vs_Cathode_Sum_Energy", 2000, 0, 20000, 2000, 0, 10000, aEMax, cESum, "hGMPC");
    plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy", 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
    plotter->Fill2D("AnodeMaxE_Vs_Cathode_Max_Energy", 800, 0, 20000, 800, 0, 10000, aEMax, cEMax, "hGMPC");
    // double sinTheta = TMath::Sin((anodeIntersection - TVector3(0,0,source_vertex)).Theta());///TMath::Sin((TVector3(51.5,0,128.) - TVector3(0,0,85)).Theta());
    // plotter->Fill2D("AnodeMaxE_Vs_Cathode_Max_Energy_path_corrected", 800, 0, 20000, 800, 0, 10000, aEMax*sinTheta, cEMax*sinTheta, "hGMPC");
    plotter->Fill2D("AnodeSumE_Vs_Cathode_Sum_Energy", 800, 0, 20000, 800, 0, 10000, aESum, cESum, "hGMPC");
    plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_TC" + std::to_string(PCQQQTimeCut) + "_PC" + std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
    // plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_path_corrected"+std::to_string(PCQQQTimeCut)+"_PC"+std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cEMax*sinTheta, "hGMPC");
    // plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_path_corrected", 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cEMax*sinTheta, "hGMPC");
    if (aEMax > 0)
    {
      double ratio = cEMax / aEMax;
      std::string folder = "Diagnostics_CMax";

      // 1. Summary 2D Plots
      plotter->Fill2D("CMax_over_Anode_vs_Z", 600, -300, 300, 200, 0, 2.0, anodeIntersection.Z(), ratio, folder);
      plotter->Fill2D("CMax_over_Anode_vs_AnodeID", 24, 0, 24, 200, 0, 2.0, aIDMax, ratio, folder);
      plotter->Fill2D("CMax_over_Anode_vs_CathodeID", 24, 0, 24, 200, 0, 2.0, cIDMax, ratio, folder);

      // 2. Individual 1D Histogram for this SPECIFIC Anode-Cathode Pair
      std::string pairName = "Ratio_A" + std::to_string(aIDMax) + "_C" + std::to_string(cIDMax);
      plotter->Fill1D(pairName, 200, 0, 2.0, ratio, folder + "/Pairs");

      // (Optional) If you also still want the independent ones:
      plotter->Fill1D("Ratio_A" + std::to_string(aIDMax), 200, 0, 2.0, ratio, folder + "/PerAnode");
      plotter->Fill1D("Ratio_C" + std::to_string(cIDMax), 200, 0, 2.0, ratio, folder + "/PerCathode");
    }

    if (PCQQQTimeCut && PCQQQPhiCut)
    {
      plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_TC" + std::to_string(PCQQQTimeCut) + "_PC" + std::to_string(PCQQQPhiCut) + "_cMax" + std::to_string(cIDMax), 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
    }
    // plotter->Fill2D("AnodeSumE_Vs_CathodeSum_Energy_path_corrected", 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cESum*sinTheta, "hGMPC");
    // plotter->Fill2D("AnodeSumE_Vs_CathodeSum_Energy_path_corrected_TC"+std::to_string(PCQQQTimeCut)+"_PC"+std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum*sinTheta, cESum*sinTheta, "hGMPC");        */
  }
  plotter->Fill1D("Correlated_Cathode_MaxAnode", 6, 0, 5, corrcatMax.size(), "hGMPC");
  plotter->Fill2D("Correlated_Cathode_VS_MaxAnodeEnergy", 6, 0, 5, 2000, 0, 30000, corrcatMax.size(), aEMax, "hGMPC");
  plotter->Fill1D("AnodeHits", 12, 0, 11, anodeHits.size(), "hGMPC");
  plotter->Fill2D("AnodeMaxE_vs_AnodeHits", 12, 0, 11, 2000, 0, 30000, anodeHits.size(), aEMax, "hGMPC");

  if (anodeHits.size() < 1)
  {
    plotter->Fill1D("NoAnodeHits_CathodeHits", 6, 0, 5, cathodeHits.size(), "hGMPC");
  }

  for (const auto &cwevent : cWireEvents)
  {
    // plotter->Fill1D("cwdtqqq_vs_cw"+std::to_string(PCQQQTimeCut),800,-2000,2000,24,0,24,std::get<2>(cwevent)-qqqtimestamp,std::get<0>(cwevent));
    for (const auto &awevent : aWireEvents)
    {
      plotter->Fill2D("aw_vs_cw", 24, 0, 24, 24, 0, 24, std::get<0>(awevent), std::get<0>(cwevent));
      plotter->Fill2D("aw_vs_cw_dtq" + std::to_string(PCQQQTimeCut), 24, 0, 24, 24, 0, 24, std::get<0>(awevent), std::get<0>(cwevent));
    }
  }
  for (const auto &awevent : aWireEvents)
  {
    // plotter->Fill1D("awdtqqq_vs_aw"+std::to_string(PCQQQTimeCut),800,-2000,2000,24,0,24,std::get<2>(awevent)-qqqtimestamp,std::get<0>(awevent));
  }
}

void miscHistograms_oneWire(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters)
{
  // consider the 'proton-like' QQQ branch seen in a,p data
  static TRandom3 rand(0); // seeded once (random seed via TUUID), not per call
  double initial_energy = 6.78;
  // if (dataset == "27Al") /// m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.6, 1 gain)
  //   initial_energy = 6.79;
  // if (dataset == "17F")
  //   initial_energy = 6.78; // m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.6, 350 gain)
  // initial_energy = 6.32; // m3 is alpha, 6.411 MeV is 7.0 MeV proton energy after havar+mylar+kapton+100mm 4He gas (molar mass 5.3, 1 gain)

  Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, initial_energy);
  for (const auto &qqqevent : QQQ_Events)
  {
    if (qqqevent.Energy1 < 0.6)
      continue; // coarse gating
    // if(qqqevent.Energy1 > 5.0) continue; //coarse gating
    for (const auto &acluster : aClusters)
    {
      auto [apwire, apSumE, apMaxE, apTSMaxE] = pwinstance.GetPseudoWire(acluster, "ANODE");
      // if(apSumE<6000) continue;
      int a_number = acluster.size();
      TVector3 pc_closest = pwinstance.getClosestWirePosAtWirePhi(apwire, qqqevent.pos.Phi());
      plotter->Fill1D("dt_anode_interp_qqq", 800, -2000, 2000, qqqevent.Time1 - apTSMaxE, "ainterp_noc");
      if (qqqevent.Time1 - apTSMaxE < 150)
      {
        bool phicut = TMath::Abs(qqqevent.pos.DeltaPhi(pc_closest)) <= TMath::Pi() / 4.0;
        TVector3 x2(pc_closest), x1(qqqevent.pos);
        TVector3 v = x2 - x1;
        double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
        TVector3 r_rhoMin_fix = x1 + t_minimum * v;

        double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
        double sinTheta2 = TMath::Sin(theta_q);

        if (r_rhoMin_fix.Perp() > 6.0)
          continue;
        if (r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100)
          continue;
        if (!phicut)
          continue;
        plotter->Fill1D("dt_anode_ainterp_qqq_gated", 800, -2000, 2000, qqqevent.Time1 - apTSMaxE, "ainterp_noc");
        plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_qqqE", 800, -2000, 2000, 800, 0, 10, qqqevent.Time1 - apTSMaxE, qqqevent.Energy1, "ainterp_noc");
        plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a" + std::to_string(acluster.size()), 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, apSumE, "ainterp_noc");
        plotter->Fill2D("pcPhi_ainterp_qqqPhi_TC1_ignC_a" + std::to_string(acluster.size()), 120, -200,200, 120, -200,200, pc_closest.Phi() * 180. / M_PI, qqqevent.pos.Phi() * 180. / M_PI, "ainterp_noc");
        plotter->Fill2D("pcZ_ainterp_qqqZ_TC1_ignC_a" + std::to_string(acluster.size()) + "_PC" + std::to_string(phicut), 300, -100, 200, 400, -200, 200, qqqevent.pos.Z(), pc_closest.Z(), "ainterp_noc");

        // plotter->Fill2D("pcZ_ainterp_qqqpczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_closest.Z(),"ainterp_noc");
        plotter->Fill2D("dEa3_ainterp_Eqqq_TC1_ignC_a" + std::to_string(acluster.size()) + "_PC" + std::to_string(phicut), 1200, 0, 30, 800, 0, 30000, qqqevent.Energy1, apSumE * sinTheta2 * 3., "ainterp_noc");

        plotter->Fill2D("vertexZ_ainterp_qqqZ_TC1_ignC_a" + std::to_string(acluster.size()), 300, -100, 200, 800, -400, 400, qqqevent.pos.Z(), r_rhoMin_fix.Z(), "ainterp_noc");
        plotter->Fill1D("vertexZ1d_ainterp_qqqZ_TC1_ignC_a" + std::to_string(acluster.size()), 800, -400, 400, r_rhoMin_fix.Z(), "ainterp_noc");
        plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a" + std::to_string(acluster.size()), 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), "ainterp_noc");

        double path_length_q = (qqqevent.pos - r_rhoMin_fix).Mag() * 0.1;
        double qqqEfix;
        qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length_q);
        plotter->Fill1D("pmisc_ow_Ex_from_alpha", 200, -10, 10, apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI), "ainterp_noc");
        plotter->Fill2D("pmisc_ow_Ef_vs_theta_qqq", 100, 0, 180, 800, 0, 20, theta_q * 180 / M_PI, qqqEfix, "ainterp_noc");
        plotter->Fill2D("pmisc_ow_VertexReconZ_vs_Ef", 800, -400, 400, 800, 0, 20, r_rhoMin_fix.Z(), qqqEfix, "ainterp_noc");
      }
    }
  } // end QQQEvents loop
}

void protonMiscHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events)
{
  // consider the 'proton-like' QQQ branch seen in a,p data
  static TRandom3 rand(0); // seeded once (random seed via TUUID), not per call
  double initial_energy = 6.78;
  // if (dataset == "27Al")
  //   initial_energy = 6.79; // m3 is alpha, 6.79 MeV is 7.0 MeV proton energy after kapton+100mm 4He gas (molar mass 5.2, 1 gain)
  // if (dataset == "17F")
  //   initial_energy = 6.32; // m3 is alpha, 6.411 MeV is 7.0 MeV proton energy after Havar+kapton+100mm 4He gas (molar mass 5.2, 1 gain)

  Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, initial_energy); // m3 is alpha

  for (const auto &qqqevent : QQQ_Events)
  {
    if (qqqevent.Energy1 < 0.6)
      continue; // coarse gating
    // if(qqqevent.Energy1 > 5.0) continue; //coarse gating
    for (const auto &pcevent : PC_Events)
    {
      if (!(pcevent.multi1 == 1 && pcevent.multi2 <= 2))
        continue;
      // if(pcevent.Energy1 > 11000) continue; //coarse gating

      bool phicut = TMath::Abs(qqqevent.pos.DeltaPhi(pcevent.pos)) <= TMath::Pi() / 4.0;
      if (!phicut)
        continue;
      // if(pcevent.Time1-qqqevent.Time1<-150 || pcevent.Time1-qqqevent.Time1 >850) continue;

      double pcz_fix, pcz_dith = pcevent.pos.Z();
      if (pcevent.multi2 == 2)
        pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
      else
      {
        pcz_fix = rand.Gaus(pcevent.pos.Z(), 8.0); // dither for a1c1 events
        pcz_dith = pcz_fix;
      }

      // --- A1C1 dither vs cfrac comparison (keep BOTH methods) ---
      // For genuine A1C1 events, reconstruct the PC z two ways and fill a parallel
      // set of excitation observables (suffix _dither / _cfrac) so the methods can
      // be overlaid directly. This runs independently of the main-flow vertex cut
      // below; each method applies its own vertex gate. The dither value also
      // feeds the main (untagged) histograms exactly as before. cfrac =
      // cpMax/(apSum+cpMax) = Energy2/(Energy1+Energy2): cell anchored on the fired
      // cathode, side from the anode-only z, offset from the per-cell autocal.
      if (pcevent.multi2 == 1 && pcevent.Energy2 > 1400)
      {
        // wire-topology category: _true1w (no dead neighbour) vs _missingw
        // (a neighbouring wire is dead -> possible masked two-wire event).
        const std::string wcat = a1c1_missing_neighbor(pcevent.Anodech, pcevent.Cathodech) ? "_missingw" : "_true1w";
        auto fillCmp = [&](double pcz, const std::string &m)
        {
          TVector3 x2(pcevent.pos.X(), pcevent.pos.Y(), pcz);
          TVector3 vv = x2 - qqqevent.pos;
          double tm = -1.0 * (qqqevent.pos.X() * vv.X() + qqqevent.pos.Y() * vv.Y()) / (vv.X() * vv.X() + vv.Y() * vv.Y());
          TVector3 rv = qqqevent.pos + tm * vv;
          if (rv.Perp() > 6.0 || rv.Z() < -173.6 || rv.Z() > 100)
            return;
          double th = (qqqevent.pos - rv).Theta();
          double pl = (qqqevent.pos - rv).Mag() * 0.1;
          double Ef = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - pl);
          double Ex = apkin_a.getExc(Ef, th * 180 / M_PI);
          std::string lbl = "proton+misc_a1c1cmp";
          // fill "all" (existing names) plus the wire-topology split (_true1w/_missingw)
          for (const std::string &w : {std::string(""), wcat})
          {
            plotter->Fill1D("pmisc_a1c1cmp_pcz_" + m + w, 600, -300, 300, pcz, lbl);
            plotter->Fill1D("pmisc_a1c1cmp_Ex_" + m + w, 200, -10, 10, Ex, lbl);
            plotter->Fill1D("pmisc_a1c1cmp_VertexZ_" + m + w, 800, -400, 400, rv.Z(), lbl);
            plotter->Fill2D("pmisc_a1c1cmp_VertexZ_vs_Ef_" + m + w, 800, -400, 400, 800, 0, 20, rv.Z(), Ef, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_phi_vs_Ef_" + m + w, 180, -180, 180, 800, 0, 20, qqqevent.pos.Phi() * 180 / M_PI, Ef, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_phi_vs_Ex_" + m + w, 180, -180, 180, 800, -5, 5, qqqevent.pos.Phi() * 180 / M_PI, Ex, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_Ef_vs_theta_" + m + w, 100, 0, 180, 800, 0, 20, th * 180 / M_PI, Ef, lbl);
          }
        };

        fillCmp(pcz_dith, "dither"); // method 1: Gaussian dither (main-flow value)

        // method 2: cfrac sub-cell linear centre-fold (side ref = anode-only z
        // rebuilt from the fired anode wire)
        double cfrac = (pcevent.Energy1 > 0.0 && pcevent.Energy2 > 0.0) ? pcevent.Energy2 / (pcevent.Energy1 + pcevent.Energy2) : -1.0;

        if (cfrac >= 0.0)
        {
          std::vector<std::tuple<int, double, double>> aOne = {std::make_tuple(pcevent.Anodech, 1.0, 0.0)};
          auto apw = pwinstance.GetPseudoWire(aOne, "ANODE");
          A1C1Sol s = a1c1_solve(cfrac, pcevent.pos.Z(), pcevent.Cathodech, pcevent.Energy1, pcevent.Anodech);
          // beam-axis 2-hypothesis side test (crossover = PC point, Si = qqq hit).
          int side_status;
          SideChoice side = a1c1_pick_side(qqqevent.pos, pcevent.pos.X(), pcevent.pos.Y(), s.pcz_lo, s.pcz_hi, side_status);
          double pcz_pick = (side == SideChoice::High) ? s.pcz_hi : s.pcz_lo;
          // cfrac_all = beam-axis pick for ALL events; "cfrac" = beam-axis pick for
          // events where at least one side is on-axis (side_status != 2).
          fillCmp(pcz_pick, "cfrac_all");
          if (side_status != 2)
          {
            fillCmp(pcz_pick, "cfrac");
            plotter->Fill2D("pmisc_a1c1cmp_pcz_cfrac_vs_dither", 600, -300, 300, 600, -300, 300, pcz_dith, pcz_pick, "proton+misc_a1c1cmp");
          }
        }
      }

      TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
      TVector3 x1(qqqevent.pos);
      TVector3 v = x2f - x1;
      double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
      TVector3 r_rhoMin_fix = x1 + t_minimum * v;
      double vertex_z = r_rhoMin_fix.Z();
      // double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
      double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
      double sinTheta_customV = TMath::Sin(theta_q);
      // if(r_rhoMin_fix.Perp()>6) continue;
      bool cathode_alpha_select = (pcevent.Energy2 > 1400);
      if (vertex_z < -173.6 || vertex_z > 100)
        continue;

      // What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
      auto plot_with_tag = [&](std::string tag = "")
      {
        std::string pmlabel = "proton+misc" + tag;
        plotter->Fill2D("pmisc_dE_E_AnodeQQQ" + tag, 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, pmlabel);
        plotter->Fill2D("pmisc_dE_E_CathodeQQQ" + tag, 400, 0, 10, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, pmlabel);
        plotter->Fill2D("pmisc_dE3_E_AnodeQQQ" + tag, 400, 0, 10, 400, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV * 3., pmlabel);
        plotter->Fill2D("pmisc_dE3_E_CathodeQQQ" + tag, 400, 0, 10, 400, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV, pmlabel);
        plotter->Fill2D("pmisc_dPhi_QQQ_PC" + tag, 180, -200,200, 180, -200,200, pcevent.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, pmlabel);
        plotter->Fill1D("pmisc_dt_Anode_QQQ_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, pcevent.Time1 - qqqevent.Time1, pmlabel);
        plotter->Fill1D("pmisc_dt_Cathode_QQQ" + tag, 600, -2000, 2000, pcevent.Time2 - qqqevent.Time1, pmlabel);
        plotter->Fill2D("pmisc_dt_Anode_E_QQQ_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time1 - qqqevent.Time1, qqqevent.Energy1, pmlabel);
        plotter->Fill2D("pmisc_dt_AnodeQQQ_vsPCPhi" + tag, 600, -2000, 2000, 180, -200,200, pcevent.Time1 - qqqevent.Time1, pcevent.pos.Phi() * 180. / M_PI, pmlabel);
        plotter->Fill2D("pmisc_dt_Cathode_E_QQQ" + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time2 - qqqevent.Time1, qqqevent.Energy1, pmlabel);
        plotter->Fill2D("pmisc_dt_CathodeQQQ_vsPCPhi" + tag, 600, -2000, 2000, 180, -200,200, pcevent.Time2 - qqqevent.Time1, pcevent.pos.Phi() * 180. / M_PI, pmlabel);
        plotter->Fill1D("pmisc_pczfix" + tag, 600, -300, 300, pcz_fix, pmlabel);
        if (pcevent.multi2 == 2)
        {
          plotter->Fill1D("pmisc_pcz" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);
          plotter->Fill1D("pmisc_pcz2" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);
        }
        if (pcevent.multi2 == 1)
        {
          plotter->Fill1D("pmisc_pcz" + tag, 600, -300, 300, pcz_fix, pmlabel);
          plotter->Fill1D("pmisc_pcz1" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);
          plotter->Fill1D("pmisc_pcz_dith" + tag, 600, -300, 300, pcz_dith, pmlabel);
        }

        // double path_length_q = (qqqevent.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
        // double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
        double path_length_q = (qqqevent.pos - r_rhoMin_fix).Mag() * 0.1;
        double qqqEfix;
        if (tag == "_cathode_alphas")
        { // satisfied when find succeeds
          qqqEfix = cm_to_MeV->Eval(MeV_to_cm->Eval(qqqevent.Energy1) - path_length_q);
          plotter->Fill1D("pmisc_Ex_from_alpha", 200, -10, 10, apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI), pmlabel);
        }
        else
          qqqEfix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(qqqevent.Energy1) - path_length_q);
        // plotter->Fill2D("qqqEf_sx3E_matrix_all"+tag,400,0,10,400,0,10,qqqEfix,sx3event.Energy1,pmlabel);
        plotter->Fill2D("pmisc_dE3_Ef_AnodeQQQ" + tag, 400, 0, 10, 400, 0, 40000, qqqEfix, pcevent.Energy1 * sinTheta_customV * 3, pmlabel);
        plotter->Fill2D("pmisc_dE3_Ef_CathodeQQQ" + tag, 400, 0, 10, 400, 0, 10000, qqqEfix, pcevent.Energy2 * sinTheta_customV, pmlabel);

        plotter->Fill1D("pmisc_VertexReconZ" + tag, 800, -400, 400, vertex_z, pmlabel);
        plotter->Fill2D("pmisc_VertexReconXY" + tag, 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), pmlabel);
        plotter->Fill2D("pmisc_VertexReconZ_vs_Ef" + tag, 800, -400, 400, 800, 0, 20, vertex_z, qqqEfix, pmlabel);
        plotter->Fill2D("pmisc_VertexReconZ_vs_Ef" + tag + "_a" + std::to_string(pcevent.multi1), 800, -400, 400, 800, 0, 20, vertex_z, qqqEfix, pmlabel);

        plotter->Fill2D("pmisc_Ef_vs_theta_qqq" + tag, 100, 0, 180, 800, 0, 20, theta_q * 180 / M_PI, qqqEfix, pmlabel);
        if (pcevent.multi2 == 1)
        {
          plotter->Fill2D("pmisc_Ef_vs_theta_qqq_a1c1" + tag, 100, 0, 180, 800, 0, 20, theta_q * 180 / M_PI, qqqEfix, pmlabel);
          plotter->Fill2D("pmisc_VertexReconZ_vs_Ef_a1c1" + tag, 800, -400, 400, 800, 0, 20, vertex_z, qqqEfix, pmlabel);
        }
      };

      if (cathode_alpha_select)
        plot_with_tag("_cathode_alphas");
      // else
      //	plot_with_tag("_cathode_protons");
      // plot_with_tag();

      // plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(qqqEfix,theta_s*180/M_PI),pmlabel);

    } // end PCEvents loop
  } // end QQQEvents loop
}

void protonMiscHistograms_sx3(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events)
{
  // consider the 'proton-like' QQQ branch seen in a,p data
  static TRandom3 rand(0); // seeded once (random seed via TUUID), not per call
  for (const auto &sx3event : SX3_Events)
  {
    if (sx3event.Energy1 < 1.2)
      continue; // coarse gating
    // if(sx3event.Energy1 > 5.0) continue; //coarse gating
    for (const auto &pcevent : PC_Events)
    {
      if (!(pcevent.multi1 == 1 && pcevent.multi2 == 2))
        continue;
      // if(pcevent.Energy1 > 11000) continue; //coarse gating

      bool phicut = TMath::Abs(sx3event.pos.DeltaPhi(pcevent.pos)) <= TMath::Pi() / 3.0;

      if (!phicut)
        continue;
      // if(pcevent.Time1-sx3event.Time1<-150 || pcevent.Time1-sx3event.Time1 >850) continue;

      double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
      TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
      TVector3 x1(sx3event.pos);
      TVector3 v = x2f - x1;
      double t_minimum = -1.0 * (x1.X() * v.X() + x1.Y() * v.Y()) / (v.X() * v.X() + v.Y() * v.Y());
      TVector3 r_rhoMin_fix = x1 + t_minimum * v;
      double vertex_z = r_rhoMin_fix.Z();
      // double theta_q = (sx3event.pos - TVector3(0,0,vertex_z)).Theta();

      if (r_rhoMin_fix.Perp() > 10.0)
        continue;
      double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
      double sinTheta_customV = TMath::Sin(theta_s);
      bool cathode_alpha_select = (pcevent.Energy2 > 1400);

      // What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
      auto plot_with_tag = [&](std::string tag = "")
      {
        std::string pmlabel = "proton+miscsx3" + tag;
        plotter->Fill2D("pmiscs_dE_E_Anodesx3" + tag, 400, 0, 10, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, pmlabel);
        plotter->Fill2D("pmiscs_dE_E_Cathodesx3" + tag, 400, 0, 10, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, pmlabel);
        plotter->Fill2D("pmiscs_dE3_E_Anodesx3" + tag, 400, 0, 10, 400, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV * 3., pmlabel);
        plotter->Fill2D("pmiscs_dE3_E_Cathodesx3" + tag, 400, 0, 10, 400, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV, pmlabel);
        plotter->Fill2D("pmiscs_dPhi_sx3_PC" + tag, 180, -200,200, 180, -200,200, pcevent.pos.Phi() * 180 / M_PI, sx3event.pos.Phi() * 180 / M_PI, pmlabel);
        plotter->Fill1D("pmiscs_dt_Anode_sx3_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, pcevent.Time1 - sx3event.Time1, pmlabel);
        plotter->Fill1D("pmiscs_dt_Cathode_sx3" + tag, 600, -2000, 2000, pcevent.Time2 - sx3event.Time1, pmlabel);
        plotter->Fill2D("pmiscs_dt_Anode_E_sx3_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time1 - sx3event.Time1, sx3event.Energy1, pmlabel);
        plotter->Fill2D("pmiscs_dt_Cathode_E_sx3" + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time2 - sx3event.Time1, sx3event.Energy1, pmlabel);
        plotter->Fill2D("pmiscs_dt_Cathodesx3_vsPCPhi" + tag, 600, -2000, 2000, 180, -200,200, pcevent.Time2 - sx3event.Time1, pcevent.pos.Phi() * 180. / M_PI, pmlabel);
        plotter->Fill1D("pmiscs_pczfix" + tag, 600, -300, 300, pcz_fix, pmlabel);
        plotter->Fill1D("pmiscs_pcz" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);

        // double path_length_q = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
        // double path_length_s = (sx3event.pos-TVector3(0,0,vertex_z)).Mag()*0.1;
        double path_length_s = (sx3event.pos - r_rhoMin_fix).Mag() * 0.1;
        double sx3Efix = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1) - path_length_s);

        // plotter->Fill2D("sx3Ef_sx3E_matrix_all"+tag,400,0,10,400,0,10,sx3Efix,sx3event.Energy1,pmlabel);
        plotter->Fill2D("pmiscs_dE3_Ef_Anodesx3" + tag, 400, 0, 10, 400, 0, 40000, sx3Efix, pcevent.Energy1 * sinTheta_customV * 3, pmlabel);
        plotter->Fill2D("pmiscs_dE3_Ef_Cathodesx3" + tag, 400, 0, 10, 400, 0, 10000, sx3Efix, pcevent.Energy2 * sinTheta_customV, pmlabel);

        plotter->Fill2D("pmiscs_Ef_vs_theta_sx3" + tag, 100, 0, 180, 800, 0, 20, theta_s * 180 / M_PI, sx3Efix, pmlabel);
        plotter->Fill1D("pmiscs_VertexReconZ" + tag, 800, -400, 400, vertex_z, pmlabel);
        plotter->Fill2D("pmiscs_VertexReconXY" + tag, 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), pmlabel);
        plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef" + tag, 800, -400, 400, 800, 0, 20, vertex_z, sx3Efix, pmlabel);
        plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef" + tag + "_a" + std::to_string(pcevent.multi1), 800, -400, 400, 800, 0, 20, vertex_z, sx3Efix, pmlabel);
      };

      plot_with_tag();
      if (cathode_alpha_select)
        plot_with_tag("_cathode_alphas");
      else
        plot_with_tag("_cathode_protons");

      // plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),pmlabel);

    } // end PCEvents loop (A1C2 main flow)

    // --- A1C1 dither vs cfrac comparison (keep BOTH methods) ---
    // The A1C2 loop above does not touch A1C1 events. Here we reconstruct the PC z
    // for genuine A1C1 events two ways and fill a parallel set of excitation
    // observables (suffix _dither / _cfrac) so the methods can be overlaid. The
    // dither is a Gaussian-smeared crossover z baseline; cfrac is the linear
    // centre-fold sub-cell model (cell anchored on the fired cathode, side from
    // the anode-only z at the SX3 phi, offset from the per-cell autocal).
    for (const auto &pcevent : PC_Events)
    {
      if (!(pcevent.multi1 == 1 && pcevent.multi2 == 1))
        continue;
      bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi() + TMath::Pi() / 3. && sx3event.pos.Phi() >= pcevent.pos.Phi() - TMath::Pi() / 3.;
      if (!phicut)
        continue;
      if (!(pcevent.Energy2 > 1400))
        continue;

      // wire-topology category: _true1w (no dead neighbour) vs _missingw
      // (a neighbouring wire is dead -> possible masked two-wire event).
      const std::string wcat = a1c1_missing_neighbor(pcevent.Anodech, pcevent.Cathodech) ? "_missingw" : "_true1w";
      auto fillCmp = [&](double pcz, const std::string &m)
      {
        TVector3 x2(pcevent.pos.X(), pcevent.pos.Y(), pcz);
        TVector3 vv = x2 - sx3event.pos;
        double tm = -1.0 * (sx3event.pos.X() * vv.X() + sx3event.pos.Y() * vv.Y()) / (vv.X() * vv.X() + vv.Y() * vv.Y());
        TVector3 rv = sx3event.pos + tm * vv;
        if (rv.Perp() > 10.0)
          return;
        double th = (sx3event.pos - rv).Theta();
        double pl = (sx3event.pos - rv).Mag() * 0.1;
        double Ef = cm_to_MeVp->Eval(MeV_to_cm_p->Eval(sx3event.Energy1) - pl);
        std::string lbl = "proton+miscsx3_a1c1cmp";
        // fill "all" (existing names) plus the wire-topology split (_true1w/_missingw)
        for (const std::string &w : {std::string(""), wcat})
        {
          plotter->Fill1D("pmiscs_a1c1cmp_pcz_" + m + w, 600, -300, 300, pcz, lbl);
          plotter->Fill1D("pmiscs_a1c1cmp_VertexZ_" + m + w, 800, -400, 400, rv.Z(), lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_VertexZ_vs_Ef_" + m + w, 800, -400, 400, 800, 0, 20, rv.Z(), Ef, lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_Ef_vs_theta_" + m + w, 100, 0, 180, 800, 0, 20, th * 180 / M_PI, Ef, lbl);
          plotter->Fill2D("pmisc_a1c1cmp_phi_vs_Ef_" + m + w, 180, -180, 180, 800, 0, 20, sx3event.pos.Phi() * 180 / M_PI, Ef, lbl);
        }
      };

      fillCmp(rand.Gaus(pcevent.pos.Z(), 8.0), "dither"); // method 1: Gaussian dither baseline

      // method 2: cfrac sub-cell linear centre-fold (side ref = anode-only z
      // rebuilt from the fired anode wire)
      double ac = pcevent.Energy1 + pcevent.Energy2;
      double cfrac = (ac > 0.0) ? pcevent.Energy2 / ac : -1.0;
      if (cfrac >= 0.0)
      {
        std::vector<std::tuple<int, double, double>> aOne = {std::make_tuple(pcevent.Anodech, 1.0, 0.0)};
        auto apw = pwinstance.GetPseudoWire(aOne, "ANODE");
        A1C1Sol s = a1c1_solve(cfrac, pcevent.pos.Z(), pcevent.Cathodech, pcevent.Energy1, pcevent.Anodech);
        // beam-axis 2-hypothesis side test (crossover = PC point, Si = sx3 hit).
        int side_status;
        SideChoice side = a1c1_pick_side(sx3event.pos, pcevent.pos.X(), pcevent.pos.Y(), s.pcz_lo, s.pcz_hi, side_status);
        double pcz_pick = (side == SideChoice::High) ? s.pcz_hi : s.pcz_lo;
        // cfrac_all = beam-axis pick for ALL events; "cfrac" = accepted (side_status != 2).
        fillCmp(pcz_pick, "cfrac_all");
        if (side_status != 2)
          fillCmp(pcz_pick, "cfrac");
      }
    } // end A1C1 comparison loop
  } // end sx3Events loop
}