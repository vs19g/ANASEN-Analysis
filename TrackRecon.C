#define TrackRecon_cxx

#define RAW_HISTOS

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
#include <TSpline.h>
#include <TSystem.h> // gSystem->mkdir for the pc_calib_raw/ output directory

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <unistd.h> // getpid(), for a unique per-process pc_calib_raw/ filename
#include <map>
#include <utility>
#include <algorithm>

// --- Analysis Control Flags ---
bool process_alpha_proton_scattering = false,
     doMiscHistograms = true,
     doPCSX3ClusterAnalysis = true,
     doPCQQQClusterAnalysis = true,
     doOldAnalysis = false,
     BenchMark = false,
     onwire_analysis = true,
     diagnostic_eplots = false,
     diagnostic_tplots = true,
     reactiondata = false,
     doPCEnergyCalibration = true,
     ta_foil_run = false,
     source_run = false;

// --- Geometry, Calibration, & Model Variables ---
double source_vertex = 53.0,
       //  z_entrance = -174.3 - 9.7 - 100.0,
    z_entrance = -174.3 - 9.7 - 270.0, // new measurement of the chamber length puts the chamber at
                                       // 1175mm instead of 1105, plus some part of the window actually lies outside the chamber
    dither_sigma = 8.0,
       dither_sigma_c0 = 16.0,
       cathode_gain = 1.0,
       a1c1_cfrac_split = 0.0,
       a1c1_missing_fmax = 2.0,
       a1c1_lowband_rfactor = 0.0,
       a1c1_z_scale_qqq = 0.0111081,
       a1c1_z_off_qqq = 34.501,
       a1c1_z_scale_sx3 = 0.0,
       a1c1_z_off_sx3 = 2.52614,
       beam_axis_x = 0.0,
       beam_axis_y = 0.0,
       ta_foil_z_mm = 0.0,
       pc_calib_alpha_source_mev = 5.486;

// --- Immutable Constants ---
const double qqq_z = 105.0,
             anode_gain = 1.5146e-5,
             sx3_phi_pitch = 6.5 * (M_PI / 180.0),
             qqq_wedge_pitch = (87.0 / 16.0) * (M_PI / 180.0),
             qqq_ring_pitch = 48.0 / 16.0;
std::string dataset;
int co2pc = 3;      // default to 3% CO2; also selects the Eloss table pc suffix.
int pressure = 250; // gas pressure (torr) for the Eloss-table filenames;
                    // overridable via the pressure_in_torr env var.

TF1 pcfix_func("func", model_invert, -200, 200);

inline TVector3 beamVertex(const TVector3 &si, const TVector3 &dir)
{
  double d = dir.X() * dir.X() + dir.Y() * dir.Y();
  double t = (d > 0.0) ? -((si.X() - beam_axis_x) * dir.X() + (si.Y() - beam_axis_y) * dir.Y()) / d : 0.0;
  return si + t * dir;
}
inline double beamPerp(const TVector3 &p)
{
  return TMath::Sqrt((p.X() - beam_axis_x) * (p.X() - beam_axis_x) + (p.Y() - beam_axis_y) * (p.Y() - beam_axis_y));
}

struct PCPath
{
  bool ok;
  double gap_cm;     // anode-cathode gap traversed in the PC gas
  double anode_cm;   // Si -> anode surface
  double cathode_cm; // Si -> cathode surface
};
inline PCPath pcPath(const TVector3 &vtx, const TVector3 &si)
{
  auto [cint, aint, dl] = find_PC_PathLength(vtx, si);
  if (dl >= 54321.0)
    return {false, (si - vtx).Mag() * 0.1, 0.0, 0.0};
  double a = (si - aint).Mag() * 0.1;
  return {true, dl, a, a - dl};
}

// Charge-collection region, guard wires -> cathode. The track crosses these outward-bound
// (vtx -> guard -> anode -> cathode -> si), so measured back from the Si end the guard is
// the farther surface: guard_cm > cathode_cm and thick_cm = guard_cm - cathode_cm.
struct PCCollect
{
  bool ok;
  double thick_cm;   // guard -> cathode, cm
  double guard_cm;   // Si -> guard surface
  double cathode_cm; // Si -> cathode surface
};
inline PCCollect pcCollectionPath(const TVector3 &vtx, const TVector3 &si)
{
  auto [gint, cint, dl] = find_PC_CollectionPath(vtx, si);
  if (dl >= 54321.0)
    return {false, 0.0, 0.0, 0.0};
  double g = (si - gint).Mag() * 0.1;
  return {true, dl, g, g - dl};
}

double ejectile_z_deut = -220.0;  // vertex_z below this -> deuteron
double ejectile_e_alpha = 2000.0; // PC anode energy above this -> alpha
// Per-ejectile (ejectile mass, recoil mass) for a reaction's (a,a)/(a,d)/(a,p)
// channels, passed into reaction_aa_core so the right Kinematics is built.
struct AAEjectileMasses
{
  double m_a, m_ra; // alpha ejectile, recoil
  double m_d, m_rd; // deuteron ejectile, recoil
  double m_p, m_rp; // proton ejectile, recoil
};

const double a1c1_zg[8] = {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};

static const double a1c1_cfmin_17F[7] = {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20};
static const double a1c1_k_17F[7] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
static const double a1c1_cfmin_27Al[7] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
static const double a1c1_k_27Al[7] = {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20};

// active per-cell set, populated by dataset in Begin()
double a1c1_cfmin_cell[7] = {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20};
double a1c1_k_cell[7] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};

static std::vector<int> a1c1_dead_anode_17F = {9, 12}; // 1 can be recovered
static std::vector<int> a1c1_dead_cathode_17F = {};    // 0,13,15 can be recovered
static std::vector<int> a1c1_dead_anode_27Al = {0, 12, 19};
static std::vector<int> a1c1_dead_cathode_27Al = {13};

std::vector<std::pair<double, double>> pcCalibData[48];

std::vector<int> *a1c1_dead_anode = &a1c1_dead_anode_17F; // active set, chosen in Begin()
std::vector<int> *a1c1_dead_cathode = &a1c1_dead_cathode_17F;

bool a1c1_is_anode_dead[24] = {false};
bool a1c1_is_cathode_dead[24] = {false};

inline void a1c1_rebuild_dead_masks()
{
  std::fill(std::begin(a1c1_is_anode_dead), std::end(a1c1_is_anode_dead), false);
  std::fill(std::begin(a1c1_is_cathode_dead), std::end(a1c1_is_cathode_dead), false);
  for (int w : *a1c1_dead_anode)
    if (w >= 0 && w < 24)
      a1c1_is_anode_dead[w] = true;
  for (int w : *a1c1_dead_cathode)
    if (w >= 0 && w < 24)
      a1c1_is_cathode_dead[w] = true;
}

// True if a neighbouring wire (index +/-1) of the fired anode OR cathode is
// dead -- i.e. this single-wire event may actually be a masked two-wire one.
inline bool a1c1_missing_neighbor(int awire, int cwire)
{
  auto deadAdj = [](const bool *deadArr, int w)
  {
    if (w < 0 || w >= 24)
      return false;
    return (w > 0 && deadArr[w - 1]) || (w < 23 && deadArr[w + 1]);
  };
  return deadAdj(a1c1_is_anode_dead, awire) || deadAdj(a1c1_is_cathode_dead, cwire);
}

inline double pathLengthCm(const TVector3 &a, const TVector3 &b)
{
  double dx = a.X() - b.X(), dy = a.Y() - b.Y(), dz = a.Z() - b.Z();
  return std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1;
}

constexpr double kTaFoilElossMeV = 0.04;

struct TaFoilRun
{
  int run;
  double z_mm;
};
static const TaFoilRun kTaFoilRuns[] = {
    // 27Al proton-scattering campaign (run_tr.sh block 3, runs 15, 17-22)
    {15, -57.28},
    {17, -135.68},
    {18, -27.88},
    {19, -8.28},
    {20, 11.32},
    {21, 30.92},
    {22, 70.12},
    // 17F proton-scattering campaign (run_tr.sh block 6, runs 38-48)
    {38, 11.32},
    {39, 30.92},
    {40, 50.52},
    {41, -184.68},
    {42, 70.12},
    {43, 109.32},
    {44, 50.52},
    {45, 30.92},
    {46, -8.28},
    {47, -8.28},
    {48, -57.28},
};

inline double applyTaFoilEloss(double beam_energy_at_vertex, double vertex_z)
{
  if (!ta_foil_run || vertex_z <= ta_foil_z_mm)
    return beam_energy_at_vertex;
  return beam_energy_at_vertex - kTaFoilElossMeV;
}

inline double evalEloss(TSpline3 *fwd, TSpline3 *inv, double E, double pathlen)
{
  if (!fwd || !inv || !std::isfinite(E) || !std::isfinite(pathlen))
    return 0.0;
  double residual = fwd->Eval(E) - pathlen;
  if (!std::isfinite(residual))
    return 0.0;
  if (residual <= 0.0)
    return 0.0;
  double e = inv->Eval(residual);
  return std::isfinite(e) ? e : 0.0;
}

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
  double cfrac_used;
  double pcz_lo;
  double pcz_hi;
  A1C1CellSol hi;
  A1C1CellSol lo;
};

inline A1C1CellSol solve_cell(int cell, int wf, double zf, double cfrac,
                              const double *cfmin, const double *kk, bool dead_neighbor)
{
  A1C1CellSol s;
  s.cell = cell;
  s.pcz = zf; // safe sentinel: fired-wire position so edge-wire defaults don't read as z=0

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

  double fmax = dead_neighbor ? a1c1_missing_fmax : 1.0;
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

  int wf = 0;
  for (int i = 1; i < 8; ++i)
    if (TMath::Abs(a1c1_zg[i] - zf) < TMath::Abs(a1c1_zg[wf] - zf))
      wf = i;

  bool dead_neighbor = a1c1_missing_neighbor(awire, cwire);           // same for both cells; hoist to avoid double scan
  s.hi = solve_cell(wf - 1, wf, zf, cfrac, cfmin, kk, dead_neighbor); // cell above (higher z)
  s.lo = solve_cell(wf, wf, zf, cfrac, cfmin, kk, dead_neighbor);     // cell below (lower z)
  s.pcz_hi = s.hi.pcz;
  s.pcz_lo = s.lo.pcz;
  return s;
}

double a1c1_side_perp_max = 20.0; // beam-axis Perp gate (mm)

// Which of the two candidate cells the beam-axis test selects.
enum class SideChoice
{
  High, // the cell ABOVE the fired wire (pcz_hi)
  Low   // the cell BELOW the fired wire (pcz_lo)
};

inline SideChoice a1c1_pick_side(const TVector3 &si, double cx, double cy, double pcz_lo, double pcz_hi, int &status)
{
  auto vtxZP = [&](double pcz, double &z, double &perp)
  {
    TVector3 pc(cx, cy, pcz);
    TVector3 vtx = beamVertex(si, pc - si);
    z = vtx.Z();
    perp = beamPerp(vtx);
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
TGraph *MeV_to_cm_d = NULL, *cm_to_MeVd = NULL;
TGraph *MeV_to_cm_27Al = NULL, *cm_to_MeV_27Al = NULL;
TGraph *MeV_to_cm_17F = NULL, *cm_to_MeV_17F = NULL;

TSpline3 *MeV_to_cm_spl = NULL, *cm_to_MeV_spl = NULL;
TSpline3 *MeV_to_cm_p_spl = NULL, *cm_to_MeVp_spl = NULL;
TSpline3 *MeV_to_cm_d_spl = NULL, *cm_to_MeVd_spl = NULL;
TSpline3 *MeV_to_cm_27Al_spl = NULL, *cm_to_MeV_27Al_spl = NULL;
TSpline3 *MeV_to_cm_17F_spl = NULL, *cm_to_MeV_17F_spl = NULL;

// declaring masses for kinematics calculations
double mass_27Al = 26.981538;
double mass_4He = 4.002603254;
double mass_1H = 1.007825032;
double mass_30Si = 29.973770;
double mass_17F = 17.002095;
double mass_20Ne = 19.992440;
double mass_2H = 2.014101778; // deuteron, for (a,d) ejectile kinematics
// Recoil masses for the (a,X) ejectile channels (from the MakeVertex branch).
double mass_19Ne_rec = 19.001880903; // 17F(a,d) recoil
double mass_29Si_rec = 28.976494664; // 27Al(a,d) recoil

// new Parabola for 4wire shift
double z_to_crossover_rho(double z)
{
  return 1.65896E-4 * z * z + 4.61626E-8 * z + 32.067;
}

// Global instances
PW pwinstance;

inline TVector3 a1c0_hybrid_pcz(const std::pair<TVector3, TVector3> &apwire, double phi,
                                bool isQQQ, double sigma, TRandom3 &rand)
{
  TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire, phi);
  pc.SetZ(a1c1_zcorr(pc.Z(), isQQQ));
  pc.SetZ(rand.Gaus(pc.Z(), sigma));
  return pc;
}

TVector3 hitPos;
double qqqenergy, qqqtimestamp;
class Event
{
public:
  Event(TVector3 p, double e1, double e2, double t1, double t2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2) {}
  Event(TVector3 p, double e1, double e2, double esum, double t1, double t2) : pos(p), Energy1(e1), Energy2(e2), EnergySum(esum), Time1(t1), Time2(t2) {}
  Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2) {}
  // Event(TVector3 p, double e1, double e2, double t1, double t2, int c1, int c2, int m1, int m2) : pos(p), Energy1(e1), Energy2(e2), Time1(t1), Time2(t2), ch1(c1), ch2(c2), multi1(m1), multi2(m2) {}
  Event(TVector3 p, double e1, double e2, double esum, double t1, double t2, int a, int c, int c1, int c2) : pos(p), Energy1(e1), Energy2(e2), EnergySum(esum), Time1(t1), Time2(t2), Anodech(a), Cathodech(c), ch1(c1), ch2(c2) {}

  TVector3 pos;
  int ch1 = -1;        // int(ch1/16) gives qqq id, ch1%16 gives ring#
  int ch2 = -1;        // int(ch2/16) gives qqq id, ch2%16 gives wedge#
  double Energy1 = -1; // Front for QQQ, Anode for PC
  double Energy2 = -1; // Back for QQQ, Cathode for PC
  double EnergySum = -1;
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
double pcEnergySlope[48];
double pcEnergyIntercept[48];
bool pcEnergyCalibLoaded = false;

HistPlotter *plotter;

bool HitNonZero;
bool sx3ecut;
bool qqqEcut;

bool PCQQQTimeCut = false;
bool PCSX3TimeCut = false, PCASX3TimeCut = false, PCCSX3TimeCut = false;
double anodeT = -99999, cathodeT = 99999;
int anodeIndex = -1, cathodeIndex = -1;

double a1c1_cfrac_pcz(const Event &pcevent, const TVector3 &si, bool &inband);
void protonAlphaHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events);
void pcCalibratedHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events_calibrated);
void miscHistograms_oneWire(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters);
void protonMiscHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events);
void protonMiscHistograms_sx3(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events);
void miscHistograms_17Fax(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, std::string globaltag = "");
void miscHistograms_27Alax(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                           const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, std::string globaltag = "");
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

  plotter->set_barrier_limit(getenv("FLUSH_BARRIER") ? std::atoll(getenv("FLUSH_BARRIER")) : 50000);

  if (getenv("reactiondata"))
  {
    reactiondata = std::atoi(getenv("reactiondata"));
    std::cout << "Analyzing dataset as reactiondata" << std::endl;
  }

  // Run classification is by OUT_DIR, not RUN_NUMBER: run numbers collide
  // across datasets/blocks (e.g. 17F's alpha+gas block and 27Al's proton
  // block both use runs 18-21), so a RUN_NUMBER-keyed lookup misclassifies
  // whichever dataset didn't originally populate kTaFoilRuns.
  std::string outdir = getenv("OUT_DIR") ? getenv("OUT_DIR") : "";
  ta_foil_run = (outdir == "Output_p");
  source_run = (outdir == "Output_a");
  if (ta_foil_run && getenv("RUN_NUMBER"))
  {
    int run_number = std::atoi(getenv("RUN_NUMBER"));
    for (const auto &r : kTaFoilRuns)
    {
      if (r.run == run_number)
      {
        ta_foil_z_mm = r.z_mm;
        break;
      }
    }
  }
  std::cout << "OUT_DIR=" << outdir << " -> ta_foil_run=" << ta_foil_run
            << " (z=" << ta_foil_z_mm << " mm), source_run=" << source_run << std::endl;

  // if (getenv("PC_ENERGY_CALIBRATION"))
  //   doPCEnergyCalibration = std::atoi(getenv("PC_ENERGY_CALIBRATION")) != 0;

  if (getenv("DATASET"))
    dataset = std::string(getenv("DATASET"));
  if (getenv("source_vertex"))
    source_vertex = (double)std::atof(std::string(getenv("source_vertex")).c_str());

  if (doPCEnergyCalibration)
    std::cout << "PC energy calibration ON: alpha source = " << pc_calib_alpha_source_mev
              << " MeV, source position = (" << beam_axis_x << ", " << beam_axis_y << ", " << source_vertex
              << ") mm -- appends raw calibration points to pc_calib_raw/ in Terminate()" << std::endl;

  if (getenv("CO2percent"))
    co2pc = std::atoi(getenv("CO2percent"));
  std::cout << "CO2 percent set to  " << co2pc << std::endl;

  if (getenv("source_vertex"))
    source_vertex = (double)std::atof(std::string(getenv("source_vertex")).c_str());

  if (getenv("DITHER_SIGMA"))
  {
    dither_sigma = std::atof(getenv("DITHER_SIGMA"));
    dither_sigma_c0 = dither_sigma;
    std::cout << "Dither Sigma set to " << dither_sigma << " mm" << std::endl;
  }

  if (getenv("CATHODE_GAIN"))
    cathode_gain = std::atof(getenv("CATHODE_GAIN"));

  const double *cfmin_src = a1c1_cfmin_17F;
  const double *k_src = a1c1_k_17F;
  const double *cfmin2_src = a1c1_cfmin2_17F;
  const double *k2_src = a1c1_k2_17F;
  a1c1_cfrac_split = 0.15;
  a1c1_lowband_rfactor = 7.0;
  a1c1_dead_anode = &a1c1_dead_anode_17F;
  a1c1_dead_cathode = &a1c1_dead_cathode_17F;
  if (dataset == "27Al")
  {
    cfmin_src = a1c1_cfmin_17F;
    k_src = a1c1_k_17F;
    cfmin2_src = a1c1_cfmin2_27Al;
    k2_src = a1c1_k2_27Al;
    a1c1_cfrac_split = 0.0;
    a1c1_lowband_rfactor = 0.0;
    a1c1_dead_anode = &a1c1_dead_anode_27Al;
    a1c1_dead_cathode = &a1c1_dead_cathode_27Al;
  }
  a1c1_rebuild_dead_masks();
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

  if (dataset == "27Al" && reactiondata)
  {
    beam_axis_x = -15.0;
    beam_axis_y = -5.0;
  }
  if (getenv("BEAM_AXIS_X"))
    beam_axis_x = std::atof(getenv("BEAM_AXIS_X"));
  if (getenv("BEAM_AXIS_Y"))
    beam_axis_y = std::atof(getenv("BEAM_AXIS_Y"));
  std::cout << "Beam-axis origin (x,y) = (" << beam_axis_x << ", " << beam_axis_y << ") mm" << std::endl;
  if (doPCEnergyCalibration)
    std::cout << "PC energy calibration ON: alpha source = " << pc_calib_alpha_source_mev
              << " MeV, source position = (" << beam_axis_x << ", " << beam_axis_y << ", " << source_vertex
              << ") mm -- appends raw calibration points to pc_calib_raw/ in Terminate()"
              << " (run pccal/fit_pc_energy_calibration.C afterward to produce pc_energy_calibration.dat)" << std::endl;
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

  //  ------------Load independent PC energy calibration (ADC -> dE_gas MeV)-------------- ///
  for (int i = 0; i < 48; i++)
  {
    pcEnergySlope[i] = 1.0;
    pcEnergyIntercept[i] = 0.0;
  }
  {
    std::ifstream pcEnergyFile("pc_energy_calibration.dat");
    if (pcEnergyFile.is_open())
    {
      std::string line;
      int index;
      double slope, intercept;
      while (std::getline(pcEnergyFile, line))
      {
        std::stringstream ss(line);
        ss >> index >> slope >> intercept;
        if (index >= 0 && index <= 47)
        {
          pcEnergySlope[index] = slope;
          pcEnergyIntercept[index] = intercept;
        }
      }
      pcEnergyFile.close();
      pcEnergyCalibLoaded = true;
      std::cout << "Loaded independent PC energy calibration from pc_energy_calibration.dat"
                << " -- populating PC_Events_calibrated" << std::endl;
    }
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

  if (getenv("pressure_in_torr"))
    pressure = std::atoi(getenv("pressure_in_torr"));
  std::cout << "Loading Eloss tables: alpha/proton/deutron/aluminum/fluorine at " << pressure
            << " torr, " << co2pc << "% CO2" << std::endl;
  MeV_to_cm = new TGraph(Form("eloss_calculations/alpha_lookup_50MeV_%dtorr_%dpc.dat", pressure, co2pc), "%lf %*lf %lf");
  MeV_to_cm_p = new TGraph(Form("eloss_calculations/proton_lookup_30MeV_%dtorr_%dpc.dat", pressure, co2pc), "%lf %*lf %lf");
  MeV_to_cm_d = new TGraph(Form("eloss_calculations/deutron_lookup_30MeV_%dtorr_%dpc.dat", pressure, co2pc), "%lf %*lf %lf");
  MeV_to_cm_27Al = new TGraph(Form("eloss_calculations/aluminum_lookup_80MeV_%dtorr_%dpc.dat", pressure, co2pc), "%lf %*lf %lf");
  MeV_to_cm_17F = new TGraph(Form("eloss_calculations/fluorine_lookup_70MeV_%dtorr_%dpc.dat", pressure, co2pc), "%lf %*lf %lf");

  auto invert = [](TGraph *g) -> TGraph *
  {
    return (g && g->GetN() > 0) ? new TGraph(g->GetN(), g->GetY(), g->GetX()) : new TGraph();
  };
  cm_to_MeV = invert(MeV_to_cm);
  cm_to_MeVp = invert(MeV_to_cm_p);
  cm_to_MeVd = invert(MeV_to_cm_d);
  cm_to_MeV_27Al = invert(MeV_to_cm_27Al);
  cm_to_MeV_17F = invert(MeV_to_cm_17F);

  auto buildSpline = [](const char *name, TGraph *g) -> TSpline3 *
  {
    if (g && g->GetN() >= 2)
    {
      TGraph sorted(*g);
      sorted.Sort();
      return new TSpline3(name, &sorted);
    }
    TGraph empty;
    empty.SetPoint(0, 0.0, 0.0);
    empty.SetPoint(1, 1.0, 0.0);
    return new TSpline3(name, &empty);
  };
  MeV_to_cm_spl = buildSpline("MeV_to_cm_spl", MeV_to_cm);
  cm_to_MeV_spl = buildSpline("cm_to_MeV_spl", cm_to_MeV);
  MeV_to_cm_p_spl = buildSpline("MeV_to_cm_p_spl", MeV_to_cm_p);
  cm_to_MeVp_spl = buildSpline("cm_to_MeVp_spl", cm_to_MeVp);
  MeV_to_cm_d_spl = buildSpline("MeV_to_cm_d_spl", MeV_to_cm_d);
  cm_to_MeVd_spl = buildSpline("cm_to_MeVd_spl", cm_to_MeVd);
  MeV_to_cm_27Al_spl = buildSpline("MeV_to_cm_27Al_spl", MeV_to_cm_27Al);
  cm_to_MeV_27Al_spl = buildSpline("cm_to_MeV_27Al_spl", cm_to_MeV_27Al);
  MeV_to_cm_17F_spl = buildSpline("MeV_to_cm_17F_spl", MeV_to_cm_17F);
  cm_to_MeV_17F_spl = buildSpline("cm_to_MeV_17F_spl", cm_to_MeV_17F);
}

inline double evalElossForward(TSpline3 *fwd, TSpline3 *inv, double E, double pathlen)
{
  if (!fwd || !inv || !std::isfinite(E) || !std::isfinite(pathlen))
    return 0.0;
  double depth0 = fwd->Eval(E);
  if (!std::isfinite(depth0))
    return 0.0;
  double depth = depth0 + pathlen;
  if (depth >= inv->GetXmax())
    return 0.0; // path length exceeds the tabulated range -> particle has fully stopped
  double e = inv->Eval(depth);
  if (!std::isfinite(e) || e < 0.0 || e > E)
    return 0.0; // extrapolated past the tabulated stopping point -> treat as fully stopped
  return e;
}

inline void pcEnergyCalibrationAccumulate(const std::vector<Event> &PC_Events,
                                          const std::vector<Event> &QQQ_Events,
                                          const std::vector<Event> &SX3_Events)
{
  const TVector3 source_pos(beam_axis_x, beam_axis_y, source_vertex);
  for (const auto &pcevent : PC_Events)
  {
    if (!(pcevent.multi1 >= 1 && pcevent.multi2 >= 1))
      continue;

    TVector3 interaction = pcevent.pos;
    if (pcevent.multi1 == 1 && pcevent.multi2 == 1)
    {
      bool inband = false;
      double pcz = a1c1_cfrac_pcz(pcevent, source_pos, inband);
      if (!inband)
        continue;
      interaction.SetZ(pcz);
    }

    // Extend the source->interaction track to find where it enters the collection region
    // (guard wires) and where it leaves it (cathode).
    TVector3 trackVec = interaction - source_pos;
    if (trackVec.Mag() < 0.01)
      continue;                                                // degenerate -- source and hit coincide
    TVector3 farPoint = source_pos + 2000.0 * trackVec.Unit(); // well beyond Si at ~88 mm
    auto [gint_s, cint_s, dl_s] = find_PC_CollectionPath(source_pos, farPoint);
    if (dl_s >= 54321.0)
      continue;                                               // geometry intersection failed
    double dist_to_entry = (gint_s - source_pos).Mag() * 0.1; // source -> guard wires, cm
    double dist_to_exit = (cint_s - source_pos).Mag() * 0.1;  // source -> cathode, cm
    if (!std::isfinite(dist_to_entry) || dist_to_entry <= 0.0 ||
        !std::isfinite(dist_to_exit) || dist_to_exit <= 0.0 ||
        dist_to_entry >= dist_to_exit)
      continue;
    double E_entry = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, pc_calib_alpha_source_mev, dist_to_entry);
    double E_exit = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, pc_calib_alpha_source_mev, dist_to_exit);
    if (!std::isfinite(E_entry) || E_entry <= 0.0 ||
        !std::isfinite(E_exit) || E_exit < 0.0 || E_entry <= E_exit)
      continue;
    if (pcevent.Anodech >= 0 && pcevent.Anodech < 24)
      pcCalibData[pcevent.Anodech].push_back({pcevent.Energy1, E_entry - E_exit});

    // Si-coincidence supplement: for each matching Si event, project the pcz using the
    // z-dependent anode radius (FIX 2: z_to_crossover_rho, not a flat 37 mm) as a sanity
    // gate, then take the collection-region dE via pcCollectionPath(source, si).
    auto considerSi = [&](const Event &sievent, double phi_win)
    {
      if (TMath::Abs(sievent.pos.DeltaPhi(pcevent.pos)) > phi_win)
        return;
      double theta = (sievent.pos - source_pos).Theta();
      if (theta <= 0.0 || !std::isfinite(theta))
        return;
      // Use z-dependent anode crossover radius for the projected pcz validity check.
      double z = z_to_crossover_rho(pcevent.pos.Z()) / TMath::Tan(theta) + source_vertex;
      if (!std::isfinite(z) || TMath::Abs(z) > 200)
        return;
      // pcCollectionPath: guard_cm = si->guard, cathode_cm = si->cathode (both from the si end).
      // Crossing order from the beam axis: source -> guard -> cathode -> si, so measured from
      // the source, dist_to_entry = total - guard_cm < dist_to_exit = total - cathode_cm.
      PCCollect pc = pcCollectionPath(source_pos, sievent.pos);
      if (!pc.ok)
        return;
      double tot = pathLengthCm(source_pos, sievent.pos);
      double d_en = tot - pc.guard_cm;
      double d_ex = tot - pc.cathode_cm;
      if (!std::isfinite(d_en) || d_en <= 0.0 || !std::isfinite(d_ex) || d_ex <= 0.0 ||
          d_en >= d_ex)
        return;
      double Ee = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, pc_calib_alpha_source_mev, d_en);
      double Ex = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, pc_calib_alpha_source_mev, d_ex);
      if (!std::isfinite(Ee) || Ee <= 0.0 || !std::isfinite(Ex) || Ex < 0.0 || Ee <= Ex)
        return;
      if (pcevent.Anodech >= 0 && pcevent.Anodech < 24)
        pcCalibData[pcevent.Anodech].push_back({pcevent.Energy1, Ee - Ex});
    };
    for (const auto &qqqevent : QQQ_Events)
      considerSi(qqqevent, TMath::Pi() / 4.0);
    for (const auto &sx3event : SX3_Events)
      considerSi(sx3event, TMath::Pi() / 3.0);
  }
}

inline double invertBeamEnergyMeV(double m1, double m2, double m3, double m4, double t3, double angle3_deg, double assumedEx = 0.0,
                                  double ebeamMeV_lo = 10.0, double ebeamMeV_hi = 100.0, int iters = 60)
{
  Kinematics kin;
  auto excAtBeamMeV = [&](double ebeamMeV)
  {
    kin.setValues(m1, m2, m3, m4, ebeamMeV / m1); // Kinematics wants E/u, not total E
    return kin.getExc(t3, angle3_deg) - assumedEx;
  };
  double f_lo = excAtBeamMeV(ebeamMeV_lo);
  double f_hi = excAtBeamMeV(ebeamMeV_hi);
  if (!std::isfinite(f_lo) || !std::isfinite(f_hi) || f_lo * f_hi > 0.0)
    return -1.0; // no root in range
  for (int i = 0; i < iters; ++i)
  {
    double mid = 0.5 * (ebeamMeV_lo + ebeamMeV_hi);
    double f_mid = excAtBeamMeV(mid);
    if (!std::isfinite(f_mid))
      return -1.0;
    if (f_mid * f_lo <= 0.0)
      ebeamMeV_hi = mid;
    else
    {
      ebeamMeV_lo = mid;
      f_lo = f_mid;
    }
  }
  return 0.5 * (ebeamMeV_lo + ebeamMeV_hi); // total MeV
}

inline double predictElasticEnergy(Kinematics &kin, double angle3_deg, double t3_lo = 0.001, double t3_hi = 60.0, int iters = 60)
{
  const int N = 200;
  double dt = (t3_hi - t3_lo) / N;
  int n_sign_changes = 0;
  double seg_lo = t3_lo, seg_hi = t3_hi;
  double prev = kin.getExc(t3_lo, angle3_deg);
  for (int k = 1; k <= N; ++k)
  {
    double t = t3_lo + k * dt;
    double cur = kin.getExc(t, angle3_deg);
    if (std::isfinite(prev) && std::isfinite(cur) && prev * cur < 0.0)
    {
      ++n_sign_changes;
      seg_lo = t - dt;
      seg_hi = t;
    }
    if (std::isfinite(cur))
      prev = cur;
  }
  if (n_sign_changes == 0)
    return -1.0; // no root in range (e.g. kinematically forbidden angle)
  if (n_sign_changes > 1)
    return -1.0; // ambiguous (multi-valued) locus -> reject
  // Single sign change: bisect within [seg_lo, seg_hi] only.
  double f_lo = kin.getExc(seg_lo, angle3_deg);
  for (int i = 0; i < iters; ++i)
  {
    double t3_mid = 0.5 * (seg_lo + seg_hi);
    double f_mid = kin.getExc(t3_mid, angle3_deg);
    if (!std::isfinite(f_mid))
      return -1.0;
    if (f_mid * f_lo <= 0.0)
      seg_hi = t3_mid;
    else
    {
      seg_lo = t3_mid;
      f_lo = f_mid;
    }
  }
  return 0.5 * (seg_lo + seg_hi);
}

inline void pcEnergyCalibrationAccumulateProton(const std::vector<Event> &PC_Events, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events)
{
  if (!ta_foil_run)
    return; // only meaningful for the proton-scattering campaign
  static const double initial_energy = 6.89;
  Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, initial_energy / mass_1H);

  auto tryEvent = [&](const Event &pcevent, const Event &sievent, double perp_max, double phi_win)
  {
    if (!(pcevent.multi1 >= 1 && pcevent.multi2 >= 1))
      return;
    if (!(pcevent.Energy2 > 1400)) // cathode-tagged alpha, same cut as protonMiscHistograms
      return;
    if (TMath::Abs(sievent.pos.DeltaPhi(pcevent.pos)) > phi_win)
      return;

    double pcz;
    if (pcevent.multi2 == 2)
      pcz = pcfix_func.Eval(pcevent.pos.Z());
    else
    {
      bool inband;
      pcz = a1c1_cfrac_pcz(pcevent, sievent.pos, inband);
      if (!inband)
        return; // only trust in-band A1C1 solutions for calibration
    }

    TVector3 x2(pcevent.pos.X(), pcevent.pos.Y(), pcz);
    TVector3 vertex = beamVertex(sievent.pos, x2 - sievent.pos);
    if (beamPerp(vertex) > perp_max || vertex.Z() < z_entrance || vertex.Z() > 100)
      return;

    double theta = (sievent.pos - vertex).Theta();
    double beam_path_length = TMath::Abs(vertex.Z() - z_entrance) * 0.1;
    double beam_energy_at_vertex = evalElossForward(MeV_to_cm_p_spl, cm_to_MeVp_spl, initial_energy, beam_path_length);
    beam_energy_at_vertex = applyTaFoilEloss(beam_energy_at_vertex, vertex.Z());
    if (beam_energy_at_vertex <= 0.0)
      return;
    apkin_a.setValues(mass_1H, mass_4He, mass_4He, mass_1H, beam_energy_at_vertex);

    double predicted_alpha_E = predictElasticEnergy(apkin_a, theta * 180.0 / M_PI);
    if (predicted_alpha_E <= 0.0)
      return;

    // pcCollectionPath: guard_cm = si->guard, cathode_cm = si->cathode (both from the si end).
    // Crossing order from the beam axis: vertex -> guard -> cathode -> si, so measured from
    // the vertex, dist_to_entry = total - guard_cm < dist_to_exit = total - cathode_cm.
    PCCollect pc = pcCollectionPath(vertex, sievent.pos);
    if (!pc.ok)
      return;
    double total_cm = pathLengthCm(vertex, sievent.pos);
    double dist_to_entry = total_cm - pc.guard_cm;  // vertex -> guard wires, cm
    double dist_to_exit = total_cm - pc.cathode_cm; // vertex -> cathode, cm
    if (!std::isfinite(dist_to_entry) || dist_to_entry <= 0.0 ||
        !std::isfinite(dist_to_exit) || dist_to_exit <= 0.0 ||
        dist_to_entry >= dist_to_exit)
      return;
    double E_entry = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, predicted_alpha_E, dist_to_entry);
    double E_exit = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, predicted_alpha_E, dist_to_exit);
    if (!std::isfinite(E_entry) || E_entry <= 0.0 ||
        !std::isfinite(E_exit) || E_exit < 0.0 || E_entry <= E_exit)
      return;
    if (pcevent.multi1 == 1 && pcevent.Anodech >= 0 && pcevent.Anodech < 24)
      pcCalibData[pcevent.Anodech].push_back({pcevent.Energy1, E_entry - E_exit});
  };

  for (const auto &pcevent : PC_Events)
  {
    for (const auto &qqqevent : QQQ_Events)
      tryEvent(pcevent, qqqevent, 6.0, TMath::Pi() / 4.0);
    for (const auto &sx3event : SX3_Events)
      tryEvent(pcevent, sx3event, 10.0, TMath::Pi() / 3.0);
  }
}

// Reads VmRSS (resident memory, MB) for this process from /proc/self/status.
// Returns -1.0 if unavailable (e.g. non-Linux) so callers can skip the check.
inline double currentRSS_MB()
{
  std::ifstream statusFile("/proc/self/status");
  std::string line;
  while (std::getline(statusFile, line))
  {
    if (line.compare(0, 6, "VmRSS:") == 0)
    {
      std::istringstream iss(line.substr(6));
      double kb = -1.0;
      iss >> kb;
      return kb > 0.0 ? kb / 1024.0 : -1.0;
    }
  }
  return -1.0;
}

Bool_t TrackRecon::Process(Long64_t entry)
{

  static const double maxRSS_MB = getenv("MAX_RSS_MB") ? std::atof(getenv("MAX_RSS_MB")) : 0.0;
  static const Long64_t checkStride = getenv("MEMCHECK_STRIDE") ? std::atoll(getenv("MEMCHECK_STRIDE")) : 5000;
  static Long64_t processedCount = 0;
  ++processedCount;

  if (maxRSS_MB > 0.0 && (processedCount % checkStride == 0))
  {
    double rss = currentRSS_MB();
    if (rss > 0.0 && rss > maxRSS_MB)
    {
      std::cout << "MAX_RSS_MB (" << maxRSS_MB << ") exceeded (RSS=" << rss
                << " MB) at entry " << entry << " -- forcing a cache flush and continuing." << std::endl;
      plotter->force_flush_caches();
    }
  }

  plotter->barrier_increment();

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

  static TRandom3 rnd_qqq(0);
  static TRandom3 rnd_sx3(0);

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
        Fsx3.at(id).ts = static_cast<double>(sx3.t[i]) + (rnd_sx3.Uniform(16.0) - 8.0);
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
      // if (det.valid)
      // {
      //   // std::cout << det.frontEL << " " << det.frontEL*sx3RightGain[id][det.stripF] << std::endl;
      //   // plotter->Fill2D("be_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_b"+std::to_string(det.stripB),200,-1,1,800,0,8192,det.frontX,det.backE,"evsx");
      //   // plotter->Fill2D("unmatched_be_vs_x_sx3_id_" + std::to_string(id), 200, -1, 1, 800, 0, 4096, det.frontX, det.backE, "evsx");
      //   // plotter->Fill2D("unmatched_be_vs_x_sx3", 200, -1, 1, 800, 0, 4096, det.frontX, det.backE, "evsx");
      //   // plotter->Fill2D("matched_be_vs_x_sx3", 200, -60, 60, 800, 0, 8192, det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx");
      //   // plotter->Fill2D("matched_be_vs_x_sx3_id_" + std::to_string(id), 200, -60, 60, 800, 0, 8192, det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx");

      //   // plotter->Fill2D("matched_be_vs_x_sx3_id_" + std::to_string(id) + "_f" + std::to_string(det.stripF), 200, -60, 60, 800, 0, 8192,
      //   //                 det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF], det.backE * sx3BackGain[id][det.stripF][det.stripB], "evsx_matched");
      //   // plotter->Fill2D("fe_vs_x_sx3_id_"+std::to_string(id)+"_f"+std::to_string(det.stripF)+"_"+std::to_string(det.stripB),200,-1,1,800,0,4096,det.frontX,det.backE,"evsx");
      //   // plotter->Fill2D("l_vs_r_sx3_id_" + std::to_string(id) + "_f" + std::to_string(det.stripF), 800, 0, 4096, 800, 0, 4096, det.frontEL, det.frontER, "l_vs_r");
      // }
      if (det.valid && (id == 9 || id == 7 || id == 1 || id == 3) && det.stripF != DEFAULT_NULL && det.stripB != DEFAULT_NULL)
      {
        double z = det.frontX * sx3FrontGain[id][det.stripF] + sx3FrontOffset[id][det.stripF];
        z = z + (75.0 / 2.0) - 3.0; // convert local sx3z to detector global coordinate system as indicated by measurements.
        // Note that this will be different for the upstream barrel, when it gets implemented
        double backE = det.backE * sx3BackGain[id][det.stripF][det.stripB];
        // det.stripF = 3 - det.stripF;
        if (id == 9 && backE < 2000)
          continue; // SX3 id 9 has an elevated low-energy noise floor (ported from MakeVertex.C)

        double alpha_n = TMath::ATan2((2 * (3 - det.stripF) - 3) * 40.30, 8.0 * 88.0 * TMath::Cos(15.0 * M_PI / 180.0)) * 180. / M_PI; // angle subtended w.r.t the radial perpendicular bisector of each sx3
        double beta_n = 15.0 + alpha_n;                                                                                                // how much to add per strip to the starting position? this is the angle w.r.t an edge of the sx3, the above values run as (-10.08deg, -3.39deg, 3.39deg, 10.08deg)
        double phi_n = ((-id + 0.5) * 30 + beta_n);
        phi_n += 45;
        double rho_at_strip = 88.0 / TMath::Cos(alpha_n * M_PI / 180.0); // TMath::Cos(15.0*M_PI/180.0) if the edge-length is 88mm
        phi_n *= M_PI / 180.;                                            // starting-position phi + strip contribution
        // Event sx3ev(TVector3(88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n),z),backE*0.001,-1,det.ts,-1,det.stripB+4*id,det.stripF+4*id);
        Event sx3ev(TVector3(rho_at_strip * TMath::Cos(phi_n), rho_at_strip * TMath::Sin(phi_n), z), backE * 0.001, -1, det.ts, -1, det.stripB + 4 * id, det.stripF + 4 * id);
        SX3_Events.push_back(sx3ev);
        if (diagnostic_eplots)
        {
          plotter->Fill2D("sx3backs_gm", 100, 0, 100, 800, 0, 8192, det.stripB + 4 * id, backE, "hCalSX3");
          plotter->Fill1D("sx3backs_calib", 800, 0, 8192, backE, "hCalSX3");

          // plotter->Fill2D("SX3CartesianPlot", 200, -100, 100, 200, -100, 100, 88.0*TMath::Cos(phi_n),88.0*TMath::Sin(phi_n), "hCalSX3");
          plotter->Fill2D("SX3CartesianPlot" + std::to_string(id), 200, -100, 100, 200, -100, 100, 88.0 * TMath::Cos(phi_n), 88.0 * TMath::Sin(phi_n), "hCalSX3");
        }
      }
    }
  }
  // return kTRUE;
  //  QQQ Processing

  int qqqCount = 0;
  std::vector<Event> QQQ_Events, PC_Events;
  std::vector<Event> PC_Events_calibrated; // independent of PC_Events; ADC->MeV via pcEnergySlope/Intercept
  // std::vector<Event> QQQ_Events_Raw, PC_Events_Raw;
  // std::vector<Event> QQQ_Events2; // clustering done

  bool PCAQQQTimeCut = false;
  bool PCCQQQTimeCut = false;
  for (int i = 0; i < qqq.multi; i++)
  {
    if (qqq.index[i] == 112)
      continue; // known-bad QQQ channel (ported from MakeVertex.C)
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
          tRing = static_cast<double>(qqq.t[j]) + (rnd_qqq.Uniform(16.0) - 8.0);
          tWedge = static_cast<double>(qqq.t[i]) + (rnd_qqq.Uniform(16.0) - 8.0);
        }
        else if (qqq.ch[j] < 16 && qqq.ch[i] >= 16 && qqqGainValid[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16])
        {
          chWedge = qqq.ch[j];
          eWedge = qqq.e[j] * qqqGain[qqq.id[j]][qqq.ch[j]][qqq.ch[i] - 16];
          chRing = qqq.ch[i] - 16;
          eRing = qqq.e[i];
          tRing = static_cast<double>(qqq.t[i]) + (rnd_qqq.Uniform(16.0) - 8.0);
          tWedge = static_cast<double>(qqq.t[j]) + (rnd_qqq.Uniform(16.0) - 8.0);
        }
        else
          continue;

        // known-bad QQQ wedge/ring channels (ported from MakeVertex.C)
        if (chWedge + qqq.id[i] * 16 == 49 || chWedge + qqq.id[i] * 16 == 48)
          continue;
        if (chRing + qqq.id[i] * 16 == 63)
          continue;

        if (diagnostic_tplots)
        {
          plotter->Fill1D("Wedgetime_Vs_Ringtime", 100, -1000, 1000, tWedge - tRing, "hTiming");
        }
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

          // double theta = 2 * TMath::Pi() * (-qqq.id[i] * 16 + (15 - chWedge) + 0.5)/(16*4);

          double phi_qqq = (M_PI / 180.) * (-90 * qqq.id[i] + (87. / 16.) * ((15 - chWedge) + 0.5) + 3.0);
          double rho = 50. + (50. / 16.) * (chRing + 0.5); //"?"
          // z used to be 75+30+23=128
          // we found a 12mm shift towards the vertex later --> 116
          Event qqqevent(TVector3(rho * TMath::Cos(phi_qqq), rho * TMath::Sin(phi_qqq), qqq_z), eRingMeV, eWedgeMeV, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);
          // Event qqqeventr(TVector3(rho * TMath::Cos(theta), rho * TMath::Sin(theta), qqq_z), eRing, eWedge, tRing, tWedge, chRing + qqq.id[i] * 16, chWedge + qqq.id[i] * 16);

          QQQ_Events.push_back(qqqevent);
          // QQQ_Events_Raw.push_back(qqqeventr);
          if (diagnostic_eplots)
          {
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
            if (diagnostic_tplots)
            {
              plotter->Fill2D("Timing_Difference_QQQ_PC", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
              plotter->Fill2D("DelT_Vs_QQQRingECal", 500, -2000, 2000, 1000, 0, 10, tRing - static_cast<double>(pc.t[k]), eRingMeV, "hTiming");
            }
            if (diagnostic_eplots)
            {
              // if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
              if (tRing - static_cast<double>(pc.t[k]) < -150) // proton tests, 27Al
              {
                PCAQQQTimeCut = true;
                plotter->Fill2D("CalibratedQQQEvsPCE_R", 1000, 0, 10, 2000, 0, 30000, eRingMeV, pc.e[k], "hPCQQQ");
                plotter->Fill2D("CalibratedQQQEvsPCE_W", 1000, 0, 10, 2000, 0, 30000, eWedgeMeV, pc.e[k], "hPCQQQ");
              }
            }
          }

          if (pc.index[k] >= 24 && pc.e[k] > 10)
          {
            if (tRing - static_cast<double>(pc.t[k]) < -200)
              PCCQQQTimeCut = true;
            if (diagnostic_tplots)
            {
              // if (tRing - static_cast<double>(pc.t[k]) > 200) PCCQQQTimeCut = true;
              plotter->Fill2D("Timing_Difference_QQQ_PC_Cathode", 500, -2000, 2000, 16, 0, 16, tRing - static_cast<double>(pc.t[k]), chRing, "hTiming");
            }
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
    if (diagnostic_eplots)
    {
      if (pc.e[i] > 50)
      {
        if (pc.index[i] >= 24)
          plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], pc.e[i] * cathode_gain, "hGMPC");
        else
          plotter->Fill2D("PC_Index_VS_GainMatched_Energy", 48, 0, 48, 2000, 0, 30000, pc.index[i], pc.e[i], "hGMPC");
      }
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
        if (pc.index[i] - 24 >= 15 && pc.index[i] - 24 <= 23)
          pc.t[i] += 300;
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
        if (diagnostic_tplots)
        {
          plotter->Fill1D("PC_Time_qqq", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
          plotter->Fill2D("PC_Time_Vs_QQQ_ch", 200, -2000, 2000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, qqq.ch[j], "hTiming");
          plotter->Fill2D("PC_Time_vs_AIndex_qqq", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
          plotter->Fill2D("PC_Time_vs_CIndex_qqq", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
          // plotter->Fill1D("PC_Time_A" + std::to_string(anodeIndex) + "_C" + std::to_string(cathodeIndex), 200, -1000, 1000, anodeT - cathodeT, "TimingPC");
        }
      }

      for (int j = 0; j < sx3.multi; j++)
      {
        if (diagnostic_tplots)
        {
          plotter->Fill1D("PC_Time_sx3", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
          // plotter->Fill2D("PC_Time_Vs_SX3_ch", 200, -2000, 2000, 16 * 8, 0, 16 * 8, anodeT - cathodeT, sx3.ch[j], "hTiming");
          plotter->Fill2D("PC_Time_vs_AIndex_sx3", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, anodeIndex, "hTiming");
          plotter->Fill2D("PC_Time_vs_CIndex_sx3", 200, -2000, 2000, 24, 0, 24, anodeT - cathodeT, cathodeIndex, "hTiming");
        }
      }
      for (const auto &sx3event : SX3_Events)
      {
        bool TCC = sx3event.Time1 - cathodeT < 0;
        bool TCA = sx3event.Time1 - anodeT < 0;
        // plotter->Fill2D("sx3_z_phi_awire"+std::to_string(anodeIndex)+"_TC"+std::to_string(TCA), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
        // plotter->Fill2D("sx3_z_phi_cwire"+std::to_string(cathodeIndex)+"_TC"+std::to_string(TCC), 400,-100,100, 200, -200,200,sx3event.pos.Z(), sx3event.pos.Phi()*180/M_PI );
      }

      if (diagnostic_tplots)
      {
        plotter->Fill1D("PC_Time", 200, -2000, 2000, anodeT - cathodeT, "hTiming");
      }
    }

    for (int j = i + 1; j < pc.multi; j++)
    {
#ifdef RAW_HISTOS
      plotter->Fill2D("PC_Coincidence_Matrix", 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
      plotter->Fill2D("PC_Coincidence_Matrix_anodeMinusCathode_lt_-200_" + std::to_string(anodeT - cathodeT < -200), 48, 0, 48, 48, 0, 48, pc.index[i], pc.index[j], "hRawPC");
#endif

      if (diagnostic_eplots)
      {
        plotter->Fill2D("Anode_V_Anode", 24, 0, 24, 24, 0, 24, pc.index[i], pc.index[j], "hGMPC");
      }
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
    if (aCluster.size() == 2)
    {
      double ae0 = std::get<1>(aCluster[0]);
      double ae1 = std::get<1>(aCluster[1]);
      double alo = std::min(ae0, ae1);
      double ahi = std::max(ae0, ae1);
      if (ahi > 0.0)
      {
        double aratio = alo / ahi;
        plotter->Fill1D("A2_anode_ratio", 120, 0, 1.2, aratio, "hGMPC");
        plotter->Fill2D("A2_anode_ratio_vs_sum", 800, 0, 40000, 120, 0, 1.2, ae0 + ae1, aratio, "hGMPC");
        plotter->Fill2D("A2_anode_ratio_vs_lowerIndex", 24, 0, 24, 120, 0, 1.2,
                        std::min(std::get<0>(aCluster[0]), std::get<0>(aCluster[1])), aratio, "hGMPC");
      }

      plotter->Fill1D("Raw_A2_AnodeSum", 800, 0, 40000, ae0 + ae1, "hGMPC");
    }
    else if (aCluster.size() == 1)
    {
      plotter->Fill1D("Raw_A1_AnodeSum", 800, 0, 40000, std::get<1>(aCluster[0]), "hGMPC");
    }
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
        Event PCEvent(crossover, apSumE, cpMaxE, cpSumE, apTSMaxE, cpTSMaxE); // run12 shows cathode-max and anode-sum provide best dE signals.
        // std::cout << apSumE << " " << crossover.Perp() << " " << apMaxE << " " << apTSMaxE << std::endl;
        PCEvent.multi1 = aCluster.size();
        PCEvent.multi2 = cCluster.size();
        PCEvent.Anodech = std::get<0>(aCluster[0]);
        PCEvent.Cathodech = std::get<0>(cCluster[0]);
        PC_Events.push_back(PCEvent);

        if (pcEnergyCalibLoaded)
        {
          Event PCEventCalibrated = PCEvent;
          double anodeCalibSum = 0.0;
          for (const auto &w : aCluster)
          {
            int wi = std::get<0>(w);
            if (wi >= 0 && wi < 24)
              anodeCalibSum += pcEnergySlope[wi] * std::get<1>(w) + pcEnergyIntercept[wi];
          }
          PCEventCalibrated.Energy1 = anodeCalibSum;
          // Cathode uses the single max wire (cpMaxE) -- indexed by z, so it's
          // already phi-consistent; leave it as-is.
          PCEventCalibrated.Energy2 = pcEnergySlope[24 + PCEvent.Cathodech] * cpMaxE + pcEnergyIntercept[24 + PCEvent.Cathodech];
          PC_Events_calibrated.push_back(PCEventCalibrated);
        }
      }
      else
      {
        ; // std::cout << "AAAA " << std::endl;
      }
    }
  }

  if ((pcEnergyCalibLoaded || doPCEnergyCalibration) && cClusters.empty())
  {
    const TVector3 source_pos_a1c0(beam_axis_x, beam_axis_y, source_vertex);
    for (const auto &aCl : aClusters)
    {
      if (aCl.empty())
        continue;
      auto aPw = pwinstance.GetPseudoWire(aCl, "ANODE");
      auto apwire = std::get<0>(aPw);
      double apSumE = std::get<1>(aPw);
      double apTSMaxE = std::get<3>(aPw);
      int anodeIdx = std::get<0>(aCl[0]);
      if (anodeIdx < 0 || anodeIdx >= 24)
        continue;

      const Event *bestSi = nullptr;
      bool bestIsQQQ = true;
      double bestDphi = 1e9;
      auto consider = [&](const std::vector<Event> &sis, bool isQQQ)
      {
        for (const auto &si : sis)
        {
          if (!(si.Time1 - apTSMaxE < 150)) // loose time coincidence (benchmark sign)
            continue;
          TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire, si.pos.Phi());
          double dphi = TMath::Abs(si.pos.DeltaPhi(pc));
          if (dphi <= TMath::Pi() / 4.0 && dphi < bestDphi)
          {
            bestDphi = dphi;
            bestSi = &si;
            bestIsQQQ = isQQQ;
          }
        }
      };
      consider(QQQ_Events, true);
      consider(SX3_Events, false);
      if (!bestSi)
        continue;

      TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire, bestSi->pos.Phi());
      pc.SetZ(a1c1_zcorr(pc.Z(), bestIsQQQ)); // same A1C0 z reference as the benchmark

      if (pcEnergyCalibLoaded)
      {
        double anodeCalibSum = 0.0;
        for (const auto &w : aCl)
        {
          int wi = std::get<0>(w);
          if (wi >= 0 && wi < 24)
            anodeCalibSum += pcEnergySlope[wi] * std::get<1>(w) + pcEnergyIntercept[wi];
        }
        Event ev(pc, anodeCalibSum, -1.0, apTSMaxE, -1.0);
        ev.multi1 = static_cast<int>(aCl.size());
        ev.multi2 = 0; // no cathode -> a{n}c0 topology in pcCalibratedHistograms
        ev.Anodech = anodeIdx;
        ev.Cathodech = -1;
        PC_Events_calibrated.push_back(ev);
      }

      // Anode-wire calibration point -- source runs only. The fixed alpha-source
      // energy is only valid there; proton-run A1C0 has no elastic tag to predict
      // its energy, so it contributes to the display but not the fit.
      if (doPCEnergyCalibration && !ta_foil_run)
      {
        TVector3 ray_dir = (pc - source_pos_a1c0).Unit();
        TVector3 virt_out = source_pos_a1c0 + ray_dir * 120.0;
        PCCollect pcc = pcCollectionPath(source_pos_a1c0, virt_out);
        if (pcc.ok)
        {
          double tc = pathLengthCm(source_pos_a1c0, virt_out);
          double de = tc - pcc.guard_cm;   // source -> guard wires, cm
          double dx = tc - pcc.cathode_cm; // source -> cathode, cm
          if (std::isfinite(de) && de > 0.0 && std::isfinite(dx) && dx > de)
          {
            double Ee = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, pc_calib_alpha_source_mev, de);
            double Ex = evalElossForward(MeV_to_cm_spl, cm_to_MeV_spl, pc_calib_alpha_source_mev, dx);
            if (std::isfinite(Ee) && std::isfinite(Ex) && Ee > Ex && Ex >= 0.0)
              pcCalibData[anodeIdx].push_back({apSumE, Ee - Ex});
          }
        }
      }
    }
  }

  if (doPCEnergyCalibration)
  {
    pcEnergyCalibrationAccumulate(PC_Events, QQQ_Events, SX3_Events);
    pcEnergyCalibrationAccumulateProton(PC_Events, QQQ_Events, SX3_Events);
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

  if (pcEnergyCalibLoaded)
    pcCalibratedHistograms(plotter, QQQ_Events, SX3_Events, PC_Events_calibrated);

  auto hasPCCoincidence = [&](const TVector3 &pos)
  {
    for (const auto &pcevent : PC_Events)
    {
      if (pcevent.multi1 < 1)
        continue;
      if (TMath::Abs(pos.DeltaPhi(pcevent.pos)) <= TMath::Pi() / 4.0)
        return true;
    }
    return false;
  };
  for (const auto &qqqevent : QQQ_Events)
  {
    plotter->Fill1D("siE_qqq_calibrated_all", 800, 0, 15, qqqevent.Energy1, "siE");
    bool coinc = hasPCCoincidence(qqqevent.pos);
    plotter->Fill1D(coinc ? "siE_qqq_calibrated_withPC" : "siE_qqq_calibrated_noPC", 800, 0, 15, qqqevent.Energy1, "siE");
  }
  for (const auto &sx3event : SX3_Events)
  {
    plotter->Fill1D("siE_sx3_calibrated_all", 800, 0, 15, sx3event.Energy1, "siE");
    bool coinc = hasPCCoincidence(sx3event.pos);
    plotter->Fill1D(coinc ? "siE_sx3_calibrated_withPC" : "siE_sx3_calibrated_noPC", 800, 0, 15, sx3event.Energy1, "siE");
  }

  if (doMiscHistograms && ta_foil_run)
  {
    if (onwire_analysis)
      miscHistograms_oneWire(plotter, QQQ_Events, aClusters);
    protonMiscHistograms_sx3(plotter, QQQ_Events, SX3_Events, PC_Events);
    protonMiscHistograms(plotter, QQQ_Events, SX3_Events, PC_Events);
  }

  if (reactiondata)
  {
    if (dataset == "17F")
      miscHistograms_17Fax(plotter, QQQ_Events, SX3_Events, PC_Events, aClusters);
    if (dataset == "27Al")
      miscHistograms_27Alax(plotter, QQQ_Events, SX3_Events, PC_Events, aClusters);
  }
  // return kTRUE;

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
  if (doPCSX3ClusterAnalysis)
  {
    PCSX3ClusterAnalysis(plotter, QQQ_Events, SX3_Events, PC_Events, aClusters, cClusters);
  }
  if (doPCQQQClusterAnalysis)
  {
    PCQQQClusterAnalysis(plotter, QQQ_Events, SX3_Events, PC_Events, aClusters, cClusters);
  }

  if (doOldAnalysis)
    OldAnalysis();
  return kTRUE;
}

void TrackRecon::Terminate()
{
  plotter->FlushToDisk(10);

  if (doPCEnergyCalibration)
  {
    gSystem->mkdir("pc_calib_raw", kTRUE);
    std::string tag = getenv("RUN_NUMBER") ? std::string("run") + getenv("RUN_NUMBER")
                                           : dataset + "_pid" + std::to_string(getpid());
    std::string outname = "pc_calib_raw/points_" + tag + ".dat";
    std::ofstream outfile(outname);
    outfile << std::scientific << std::setprecision(6);
    long long nPoints = 0;
    for (int wire = 0; wire < 48; ++wire)
    {
      for (const auto &p : pcCalibData[wire])
      {
        outfile << wire << " " << p.first << " " << p.second << "\n";
        ++nPoints;
      }
    }
    outfile.close();
    std::cout << "PC energy calibration: appended " << nPoints << " raw points to " << outname
              << " -- run pccal/fit_pc_energy_calibration.C once all calibration runs are done"
              << " to (re)produce pc_energy_calibration.dat" << std::endl;
  }
}

void protonAlphaHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events)
{

  std::string aplabel = "a(p,p)";
  double initial_energy = 6.89;

  Kinematics apkin_p(mass_1H, mass_4He, mass_1H, mass_4He, initial_energy / mass_1H); // m3 is proton
  Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, initial_energy / mass_1H); // m3 is alpha

  for (const auto &qqqevent : QQQ_Events)
  {
    for (const auto &sx3event : SX3_Events)
    {
      plotter->Fill1D("ap_qqq_sx3_dt", 800, -2000, 2000, qqqevent.Time1 - sx3event.Time1, aplabel);
      if (TMath::Abs(qqqevent.Time1 - sx3event.Time1) > 300)
        continue;
      // sx3event.pos.SetZ(sx3event.pos.Z()+5.0);
      plotter->Fill1D("ap_qqq_sx3_dt_timecut", 800, -2000, 2000, qqqevent.Time1 - sx3event.Time1, aplabel);
      plotter->Fill1D("ap_qqq_sx3_dphi", 100, -200, 200, qqqevent.pos.Phi() * 180 / M_PI - sx3event.pos.Phi() * 180 / M_PI, aplabel);
      plotter->Fill2D("ap_qqq_sx3_dphi_vs_qqqphi", 100, -200, 200, 100, -200, 200, qqqevent.pos.Phi() * 180 / M_PI - sx3event.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, aplabel);
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

        plotter->Fill2D("ap_dPhi_QQQ_PC", 100, -200, 200, 100, -200, 200, pcevent.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, aplabel);
        plotter->Fill2D("ap_dPhi_SX3_PC", 100, -200, 200, 100, -200, 200, pcevent.pos.Phi() * 180 / M_PI, sx3event.pos.Phi() * 180 / M_PI, aplabel);
        plotter->Fill1D("ap_dt_Anode_QQQ", 600, -2000, 2000, pcevent.Time1 - qqqevent.Time1, aplabel);
        plotter->Fill1D("ap_dt_Cathode_QQQ", 600, -2000, 2000, pcevent.Time2 - qqqevent.Time1, aplabel);
        plotter->Fill1D("ap_dt_Anode_SX3", 600, -2000, 2000, pcevent.Time1 - sx3event.Time1, aplabel);
        plotter->Fill1D("ap_dt_Cathode_SX3", 600, -2000, 2000, pcevent.Time2 - sx3event.Time1, aplabel);
        plotter->Fill1D("ap_pczfix", 600, -300, 300, pcz_fix, aplabel);
        plotter->Fill1D("ap_pcz", 600, -300, 300, pcevent.pos.Z(), aplabel);

        double dzq = qqqevent.pos.Z() - vertex_z;
        double dzs = sx3event.pos.Z() - vertex_z;
        double path_length_q = std::sqrt(qqqevent.pos.Perp2() + dzq * dzq) * 0.1;
        double path_length_s = std::sqrt(sx3event.pos.Perp2() + dzs * dzs) * 0.1;

        double qqqEfix = evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, qqqevent.Energy1, path_length_q);
        double sx3Efix = evalEloss(MeV_to_cm_p_spl, cm_to_MeVp_spl, sx3event.Energy1, path_length_s);
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

void pcCalibratedHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events_calibrated)
{
  for (const auto &pcevent : PC_Events_calibrated)
  {
    const std::string topo = "_a" + std::to_string(pcevent.multi1) + "c" + std::to_string(pcevent.multi2);
    const bool hasCathode = (pcevent.Cathodech >= 0);
        if (hasCathode)
      plotter->Fill2D("Calib_AnodeE_vs_CathodeE_a1c1andup", 800, 0, 3, 800, 0, 3, pcevent.Energy1, pcevent.Energy2, "hCalibPC");
    for (const std::string &t : {std::string(""), topo})
    {
      plotter->Fill2D("Calib_AnodeE_vs_AnodeIndex" + t, 24, 0, 24, 800, 0, 3, pcevent.Anodech, pcevent.Energy1, "hCalibPC");
      plotter->Fill1D("Calib_AnodeE" + t, 800, 0, 3, pcevent.Energy1, "hCalibPC");
      if (hasCathode)
      {
        plotter->Fill2D("Calib_CathodeE_vs_CathodeIndex" + t, 24, 0, 24, 800, 0, 3, pcevent.Cathodech, pcevent.Energy2, "hCalibPC");
        plotter->Fill1D("Calib_CathodeE" + t, 800, 0, 3, pcevent.Energy2, "hCalibPC");
        plotter->Fill2D("Calib_AnodeE_vs_CathodeE" + t, 800, 0, 3, 800, 0, 3, pcevent.Energy1, pcevent.Energy2, "hCalibPC");
      }

      for (const auto &qqqevent : QQQ_Events)
      {
        plotter->Fill2D("Calib_dE_AnodeE_vs_QQQE" + t, 400, 0, 10, 800, 0, 3, qqqevent.Energy1, pcevent.Energy1, "hCalibPC");
       if (hasCathode)
          plotter->Fill2D("Calib_dE_CathodeE_vs_QQQE" + t, 400, 0, 10, 800, 0, 3, qqqevent.Energy1, pcevent.Energy2, "hCalibPC");
      }
      for (const auto &sx3event : SX3_Events)
      {
        plotter->Fill2D("Calib_dE_AnodeE_vs_SX3E" + t, 400, 0, 10, 800, 0, 3, sx3event.Energy1, pcevent.Energy1, "hCalibPC");
        if (hasCathode)
          plotter->Fill2D("Calib_dE_CathodeE_vs_SX3E" + t, 400, 0, 10, 800, 0, 3, sx3event.Energy1, pcevent.Energy2, "hCalibPC");
      }
    }
  }
}

void PCSX3ClusterAnalysis(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, const std::vector<std::vector<std::tuple<int, double, double>>> &cClusters)
{

  static TRandom3 rand(0);

  // --- GENUINE A1C0 events:
  if (BenchMark && aClusters.size() == 1 && cClusters.size() == 0)
  {
    const auto &aCl = aClusters.front();
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
      // A1C0 hybrid z (shared with the QQQ twin block + miscHistograms_oneWire).
      TVector3 pc_hybrid = a1c0_hybrid_pcz(apwire_bm, sx3event.pos.Phi(), true, dither_sigma, rand);
      TVector3 vtx0 = beamVertex(sx3event.pos, pc - sx3event.pos);
      TVector3 vtx1 = beamVertex(smeared_sx3, pc_hybrid - smeared_sx3);

      if (!(vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6))
        continue;
      double sx3theta = TMath::ATan2(88.0, sx3event.pos.Z() - source_vertex);
      double pczguess = 37.0 / TMath::Tan(sx3theta) + source_vertex;
      plotter->Fill1D("Benchmark_SX3_VertexZ_trueA1C0", 800, -400, 400, vtx0.Z(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_VertexZ_trueA1C0_Hybrid", 800, -400, 400, vtx1.Z(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_VertexZ_trueA1C0_Hybrid_TC" + std::to_string(PCSX3TimeCut) + "_PC" + std::to_string(phicut), 800, -400, 400, vtx1.Z(), "A1C0True_SX3");
      plotter->Fill2D("Benchmark_SX3_VertexXY_trueA1C0_Hybrid", 200, -100, 100, 200, -100, 100, vtx1.X(), vtx1.Y(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C0_Hybrid", 600, -200, 200, pc_hybrid.Z(), "A1C0True_SX3");
      plotter->Fill2D("Benchmark_SX3_PCZ_trueA1C0_Hybrid_vs_sx3pczguess", 400, -200, 200, 400, -200, 200, pczguess, pc_hybrid.Z(), "A1C0True_SX3");
      plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C0_Hybrid_minus_sx3pczguess", 400, -100, 100, pc_hybrid.Z() - pczguess, "A1C0True_SX3");
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

      plotter->Fill2D("sx3phi_vs_pcphi" + std::to_string(sx3event.Time1 - pcevent.Time1 < -150), 100, -200, 200, 100, -200, 200, sx3event.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");
      plotter->Fill1D("sx3phi_minus_pcphi" + std::to_string(sx3event.Time1 - pcevent.Time1 < -150), 100, -180, 180, (sx3event.pos.DeltaPhi(pcevent.pos)) * 180 / M_PI, "Kinematics_Angles");

      if (PCSX3TimeCut)
      {
        plotter->Fill1D("dt_pcA_sx3B_timecut", 640, -2000, 2000, sx3event.Time1 - pcevent.Time1, "Timing");
        plotter->Fill1D("dt_pcC_sx3B_timecut", 640, -2000, 2000, sx3event.Time1 - pcevent.Time2, "Timing");
        plotter->Fill2D("xyplot_sx3" + std::to_string(sx3event.ch2 / 4), 100, -100, 100, 100, -100, 100, sx3event.pos.X(), sx3event.pos.Y(), "Vertex_Reconstruction");
        plotter->Fill2D("xyplot_sx3" + std::to_string(sx3event.ch2 / 4), 100, -100, 100, 100, -100, 100, pcevent.pos.X(), pcevent.pos.Y(), "Vertex_Reconstruction");
        plotter->Fill2D("pcz_vs_pcphi_TimeCut", 600, -200, 200, 120, -200, 200, pcevent.pos.Z(), pcevent.pos.Phi() * 180 / M_PI, "PCZ_Recon");
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
      plotter->Fill2D("pcphi_vs_time", 2000, 0, 2000, 100, -200, 200, pcevent.Time1 * 1e-9, pcevent.pos.Phi() * 180. / M_PI, "Timing");
      plotter->Fill2D("sx3phi_vs_time", 2000, 0, 2000, 100, -200, 200, pcevent.Time1 * 1e-9, sx3event.pos.Phi() * 180. / M_PI, "Timing");

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

        // if (TMath::Abs(r_rhoMin_fix.Z()) < 200.0)
        // {
        //   plotter->Fill2D("dE3_E_AnodeSX3B_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV, "PID_dE_E");
        //   plotter->Fill2D("dE3_E_CathodeSX3B_A1C2_(vertex_fix_z/100)=" + std::to_string(floor(r_rhoMin_fix.Z() / 100.0)), 400, 0, 30, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
        // }

        if (pcevent.multi1 == 1 && pcevent.multi2 == 3)
        {
          plotter->Fill2D("pcz_vs_sx3pczguess_A1C3", 600, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");
        }

        plotter->Fill2D("pcz_vs_sx3pczguess_int", 600, -200, 200, 600, -200, 200, pcz_guess_int, pcevent.pos.Z(), "PCZ_Recon");
        // plotter->Fill2D("pcz_vs_sx3pczguess_strip" + std::to_string(sx3event.ch2), 300, -200, 200, 600, -200, 200, pczguess, pcevent.pos.Z(), "PCZ_Recon");

        bool sx3PhiCut = (TMath::Abs(sx3event.pos.Phi() - pcevent.pos.Phi()) < 45.0 * M_PI / 180.);

        plotter->Fill1D("pcz_sx3Coinc_phiCut" + std::to_string(sx3PhiCut) + "_TC" + std::to_string(PCSX3TimeCut), 300, 0, 200, sx3z, "PCZ_Recon");
        plotter->Fill2D("pcz_vs_sx3z_phiCut" + std::to_string(sx3PhiCut) + "_TC" + std::to_string(PCSX3TimeCut), 300, 0, 200, 600, -400, 400, sx3z, pcevent.pos.Z(), "PCZ_Recon");

        plotter->Fill2D("sx3E_vs_sx3z", 400, 0, 30, 300, 0, 200, sx3event.Energy1, sx3z, "Kinematics_Angles");

        // plotter->Fill2D("pcdEA_vs_sx3z", 300, 0, 200, 800, 0, 20000, sx3z, pcevent.Energy1, "Kinematics_Angles");
        // plotter->Fill2D("pcdEA_vs_sx3pczguess", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy1, "Kinematics_Angles");
        plotter->Fill2D("pcdEA_vs_pczfix", 600, -200, 200, 800, 0, 20000, pcz_fix, pcevent.Energy1, "PCdE_vs_Z");
        // plotter->Fill2D("pcdEC_vs_sx3z", 300, 0, 200, 800, 0, 20000, sx3z, pcevent.Energy2, "Kinematics_Angles");
        // plotter->Fill2D("pcdEC_vs_sx3pczguess", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "Kinematics_Angles");
        plotter->Fill2D("pcdEC_vs_pczfix", 600, -200, 200, 800, 0, 20000, pcz_fix, pcevent.Energy2, "PCdE_vs_Z");
        plotter->Fill2D("pcdEACSum_vs_pczfix", 600, -200, 200, 800, 0, 20000, pcz_fix, (pcevent.Energy1 + pcevent.Energy2) / 2, "PCdE_vs_Z");

        // plotter->Fill2D("pcdEA_vs_sx3z" + std::to_string(sx3event.ch2), 300, 0, 200, 800, 0, 20000, sx3z, pcevent.Energy1, "Kinematics_Angles");
        // plotter->Fill2D("pcdEA_vs_sx3pczguess" + std::to_string(sx3event.ch2), 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy1, "Kinematics_Angles");
        // plotter->Fill2D("pcdEC_vs_sx3z" + std::to_string(sx3event.ch2), 300, 0, 200, 800, 0, 20000, sx3z, pcevent.Energy2, "Kinematics_Angles");
        // plotter->Fill2D("pcdEC_vs_sx3pczguess" + std::to_string(sx3event.ch2), 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "Kinematics_Angles");

        plotter->Fill2D("pcdE2A_vs_sx3z", 300, 0, 200, 800, 0, 20000, sx3z, pcevent.Energy1 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("pcdE2C_vs_sx3z", 300, 0, 200, 800, 0, 20000, sx3z, pcevent.Energy2 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("phi_vs_stripnum", 180, -180, 180, 48, 0, 48, pcevent.pos.Phi() * 180. / M_PI, sx3event.ch2, "Kinematics_Angles");
        plotter->Fill2D("E_theta_AnodeSX3", 300, 0, 15, 400, -20, 180, sx3event.Energy1, sx3theta * 180 / M_PI, "Kinematics_Angles");
      }

      // plotter->Fill2D("pcdEA_vs_sx3pczguess_A" + std::to_string(pcevent.multi1) + "C" + std::to_string(pcevent.multi2), 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy1, "PCdE_vs_Z");
      plotter->Fill2D("pcdEA_vs_sx3pczguess", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy1, "PCdE_vs_Z");
      plotter->Fill2D("pcdEA_vs_pczguess", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy1, "PCdE_vs_Z");
      // plotter->Fill2D("pcdEA_vs_pczfix" + std::to_string(pcevent.multi1) + "A" + std::to_string(pcevent.multi1) + "C", 600, -200, 200, 800, 0, 20000, pcz_fix, pcevent.Energy1, "PCdE_vs_Z");
      // plotter->Fill2D("pcdEC_vs_sx3pczguess_A" + std::to_string(pcevent.multi1) + "C" + std::to_string(pcevent.multi2), 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "PCdE_vs_Z");
      if (pcevent.multi1 == 1)
      {
        plotter->Fill2D("pcdEA_vs_sx3pczguess_A1", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy1, "PCdE_vs_Z");
      }
      if (pcevent.multi2 == 1)
      {
        plotter->Fill2D("pcdEC_vs_sx3pczguess_C1", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "PCdE_vs_Z");
      }
      if (pcevent.multi2 == 2)
      {
        plotter->Fill2D("pcdEC_vs_sx3pczguess_C2", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "PCdE_vs_Z");
      }
      plotter->Fill2D("pcdEC_vs_sx3pczguess", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "PCdE_vs_Z");
      plotter->Fill2D("pcdEC_vs_pczguess", 600, -200, 200, 800, 0, 20000, pczguess, pcevent.Energy2, "PCdE_vs_Z");
      plotter->Fill2D("pcdEACSum_vs_pczguess", 600, -200, 200, 800, 0, 20000, pczguess, (pcevent.Energy2 + pcevent.Energy1) / 2, "PCdE_vs_Z");

      // plotter->Fill2D("pcdEC_vs_pczfix" + std::to_string(pcevent.multi1) + "A" + std::to_string(pcevent.multi1) + "C", 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pcz_fix, "PCdE_vs_Z");

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

        if (aSumE_bm > 0.0)
          plotter->Fill2D("Benchmark_SX3_CmaxOverAnode_vs_phi", 90, -180, 180, 250, 0, 5,
                          sx3event.pos.Phi() * 180. / M_PI, cSumE_bm / aSumE_bm, "Benchmark_SX3_ref");

        double sx3_phi_pitch = 6.5 * (M_PI / 180.0);
        double smeared_phi = sx3event.pos.Phi() + rand.Uniform(-sx3_phi_pitch / 2.0, sx3_phi_pitch / 2.0);
        TVector3 smeared_sx3_pos(sx3event.pos.Perp() * TMath::Cos(smeared_phi), sx3event.pos.Perp() * TMath::Sin(smeared_phi), sx3event.pos.Z());

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

        auto doA1C1Model = [&](const std::string &tag, const TVector3 &si_point)
        {
          if (!a1c1Good || cfrac < 0.0)
            return;
          A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
          int side_status;
          SideChoice side = a1c1_pick_side(si_point, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
          const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
          double pcz_pick = best.pcz;
          if (!(best.inband && best.pitchok && side_status != 2))
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
            {
              double phi_deg = sx3event.pos.Phi() * 180.0 / M_PI;
              double vz_resid = vtx_ref.Z() - source_vertex;
              plotter->Fill2D("Diag_SX3_A1C2_vtxZ_resid_vs_phi", 90, -180, 180, 400, -100, 100, phi_deg, vz_resid, "Diag_XYoffset");
              plotter->Fill2D("Diag_Combined_A1C2_vtxZ_resid_vs_phi", 90, -180, 180, 400, -100, 100, phi_deg, vz_resid, "Diag_XYoffset");
              plotter->Fill2D("Diag_SX3_A1C2_vtxXY", 200, -15, 15, 200, -15, 15, vtx_ref.X(), vtx_ref.Y(), "Diag_XYoffset");
              plotter->Fill2D("Diag_Combined_A1C2_time_vs_phi", 2000, 0, 2000, 90, -180, 180, pcevent.Time1 * 1e-9, phi_deg, "Diag_XYoffset");
              plotter->Fill2D("Diag_SX3_A1C2_T_vs_vtxX", 2000, 0, 2000, 200, -15, 15, pcevent.Time1 * 1e-9, vtx_ref.X(), "Diag_XYoffset");
              plotter->Fill2D("Diag_SX3_A1C2_T_vs_vtxY", 2000, 0, 2000, 200, -15, 15, pcevent.Time1 * 1e-9, vtx_ref.Y(), "Diag_XYoffset");
            }

            doA1C1("A1C1", sx3event.pos, false);
            doAnodeOnly("A1C0", sx3event.pos.Phi(), sx3event.pos, false);
            doA1C1("A1C1_Hyb", smeared_sx3_pos);
            doAnodeOnly("A1C0_Hyb", smeared_phi, smeared_sx3_pos);

            doA1C1Model("A1C1_Cfrac", sx3event.pos);

            {
              double pcz_a1c0 = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, sx3event.pos.Phi()).Z();
              double theta_ref = (sx3event.pos - TVector3(0, 0, vtx_ref.Z())).Theta() * 180. / M_PI;
              plotter->Fill2D("Benchmark_SX3_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200, theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_SX3_ref");
              plotter->Fill2D("Benchmark_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200, theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_AnodeOnly");

              double phi_deg_a = sx3event.pos.Phi() * 180.0 / M_PI;
              plotter->Fill2D("Diag_SX3_A1C0_zresid_vs_phi", 90, -180, 180, 200, -100, 100, phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
              plotter->Fill2D("Diag_Combined_A1C0_zresid_vs_phi", 90, -180, 180, 200, -100, 100, phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
            }

            if (a1c1Good && cfrac >= 0.0)
            {
              plotter->Fill1D("Benchmark_SX3_A1C1_cfrac", 220, -0.05, 1.05, cfrac, "Benchmark_SX3_ref");
              plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_ref", 400, -200, 200, 220, -0.05, 1.05, pcz_ref, cfrac, "Benchmark_SX3_ref");
              plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_sx3pczguess", 400, -200, 200, 220, -0.05, 1.05, pczguess, cfrac, "Benchmark_SX3_ref");

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

              for (int i = 0; i < 7; ++i)
              {
                if (pcz_ref <= zg[i] && pcz_ref > zg[i + 1])
                {
                  double zc = 0.5 * (zg[i] + zg[i + 1]);
                  double half = 0.5 * (zg[i] - zg[i + 1]);
                  if (half > 0.0)
                    plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05, TMath::Abs(pcz_ref - zc) / half, cfrac, "Benchmark_SX3_ref");
                  break;
                }
              }

              plotter->Fill2D("Benchmark_SX3_A1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05, aSumE_bm, cfrac, "Benchmark_SX3_ref");
              if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                plotter->Fill2D("Benchmark_SX3_A1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0, 1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "Benchmark_SX3_ref");

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
                  plotter->Fill2D("Benchmark_SX3_A1C1_cellsel_confusion", 7, 0, 7, 7, 0, 7, cell_truth + 0.5, sm_cell + 0.5, "Benchmark_SX3_ref");
                  plotter->Fill1D("Benchmark_SX3_A1C1_cellsel_misclass", 2, 0, 2, wrong ? 1.0 : 0.0, "Benchmark_SX3_ref");
                  plotter->Fill2D("Benchmark_SX3_A1C1_cellsel_misclass_vs_cell", 7, 0, 7, 2, 0, 2, cell_truth + 0.5, wrong ? 1.0 : 0.0, "Benchmark_SX3_ref");

                  double zc = 0.5 * (a1c1_zg[cell_truth] + a1c1_zg[cell_truth + 1]);
                  double half = 0.5 * (a1c1_zg[cell_truth] - a1c1_zg[cell_truth + 1]);

                  plotter->Fill2D("AnodeEnergy_vs_CellSX3", 120, 0, 1.2, 800, 0, 40000, 1 - TMath::Abs(pcz_ref - zc) / half, pcevent.Energy1);
                  plotter->Fill2D("CathodeEnergy_vs_CellSX3", 120, 0, 1.2, 800, 0, 40000, TMath::Abs(pcz_ref - zc) / half, pcevent.Energy2);
                  plotter->Fill2D("FracEnergy_vs_CellSX3", 120, 0, 1.2, 800, 0, 10, TMath::Abs(pcz_ref - zc) / half, pcevent.Energy2 / pcevent.Energy1);
                  plotter->Fill2D("SumEnergy_vs_CellSX3", 120, 0, 1.2, 800, 0, 10, TMath::Abs(pcz_ref - zc) / half, (pcevent.Energy1 + pcevent.Energy2) / 2);

                  if (half > 0.0)
                  {
                    plotter->Fill2D("Benchmark_SX3_A1C1_cellsel_misclass_vs_fold", 120, 0, 1.2, 2, 0, 2, TMath::Abs(pcz_ref - zc) / half, wrong ? 1.0 : 0.0, "Benchmark_SX3_ref");
                    plotter->Fill2D("Benchmark_SX3_A1C1_cfracUsed_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05, TMath::Abs(pcz_ref - zc) / half, sm.cfrac_used, "Benchmark_SX3_ref");
                    if (aSumE_bm > 0.0)
                    {
                      plotter->Fill2D("Benchmark_SX3_A1C1_cfracUsed_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05,
                                      aSumE_bm, sm.cfrac_used, "Benchmark_SX3_ref");
                    }
                  }
                }
              }
            }
          }

          else if (pcevent.multi1 == 1 && pcevent.multi2 == 1 && a1c1Good)
          {
            double pcz_raw = xo_a1c1.Z();
            TVector3 vtx_raw = vertexFrom(sx3event.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_raw));
            fillSuite("trueA1C1", pcz_raw, vtx_raw, "A1C1True_SX3");
            plotter->Fill2D("Benchmark_SX3_PCZ_trueA1C1_vs_sx3pczguess", 400, -200, 200, 400, -200, 200, pczguess, pcz_raw, "A1C1True_SX3");
            plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C1_minus_sx3pczguess", 400, -100, 100, pcz_raw - pczguess, "A1C1True_SX3");

            if (cfrac >= 0.0)
            {
              A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
              int side_status;
              SideChoice side = a1c1_pick_side(sx3event.pos, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
              const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
              int cell = best.cell;
              double f = best.f;
              double pcz_cf = best.pcz;
              bool valid = (side_status != 2);
              plotter->Fill1D("Benchmark_SX3_trueA1C1_sideStatus", 4, -1, 3, side_status + 0.5, "A1C1True_SX3");

              TVector3 vtx_cf = vertexFrom(sx3event.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_cf));
              fillSuite(valid ? "trueA1C1_Cfrac" : "trueA1C1_Cfrac_invalid", pcz_cf, vtx_cf, "A1C1True_SX3");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_cfrac", 220, -0.05, 1.05, cfrac, "A1C1True_SX3");
              plotter->Fill2D("Benchmark_SX3_trueA1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05, aSumE_bm, cfrac, "A1C1True_SX3");
              if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                plotter->Fill2D("Benchmark_SX3_trueA1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0, 1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "A1C1True_SX3");
              plotter->Fill2D("Benchmark_SX3_trueA1C1_cfrac_vs_cell", 7, 0, 7, 220, -0.05, 1.05, cell + 0.5, cfrac, "A1C1True_SX3");
              plotter->Fill2D("Benchmark_SX3_trueA1C1_f_vs_cell", 7, 0, 7, 260, -1.5, 2.5, cell + 0.5, f, "A1C1True_SX3");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_f", 260, -1.5, 2.5, f, "A1C1True_SX3");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_valid", 2, 0, 2, valid ? 1.0 : 0.0, "Benchmark_SX3_trueA1C1");
              int reason;
              if (cell < 0 || cell > 6 || a1c1_k_cell[cell] <= 0.0)
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
              if (valid)
                plotter->Fill1D("Benchmark_SX3_trueA1C1_validreason", 3, 0, 3, reason + 0.5, "Benchmark_SX3_trueA1C1");
              plotter->Fill1D("Benchmark_SX3_trueA1C1_band", 2, 0, 2, s.band + 0.5, "Benchmark_SX3_trueA1C1");
              if (valid)
                plotter->Fill1D("Benchmark_SX3_trueA1C1_band_valid", 2, 0, 2, s.band + 0.5, "Benchmark_SX3_trueA1C1");
              if (valid)
              {
                plotter->Fill1D("Benchmark_SX3_PCZ_trueA1C1_Cfrac_minus_sx3pczguess_DIAG", 400, -100, 100, pcz_cf - pczguess, "Benchmark_SX3_trueA1C1");
                plotter->Fill2D("Benchmark_SX3_PCZ_trueA1C1_Cfrac_vs_sx3pczguess_DIAG", 400, -200, 200, 400, -200, 200, pczguess, pcz_cf, "Benchmark_SX3_trueA1C1");
              }
            }

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

  // --- GENUINE A1C0 events (QQQ twin):
  if (BenchMark && aClusters.size() == 1 && cClusters.size() == 0)
  {
    const auto &aCl = aClusters.front();
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
      // A1C0 hybrid z (shared with the SX3 twin block + miscHistograms_oneWire).
      TVector3 pc_hybrid = a1c0_hybrid_pcz(apwire_bm, smeared_phi, true, dither_sigma_c0 / 2.0, rand);
      TVector3 vtx0 = beamVertex(qqqevent.pos, pc - qqqevent.pos);
      TVector3 vtx1 = beamVertex(smeared_qqq, pc_hybrid - smeared_qqq);

      if (!(vtx0.Perp() <= 6.0 && vtx0.Z() >= -173.6))
        continue;
      double qqqTheta = (qqqevent.pos - TVector3(0, 0, source_vertex)).Theta();
      double pcz_guess_37 = 37. / TMath::Tan(qqqTheta) + source_vertex;
      plotter->Fill1D("Benchmark_QQQ_VertexZ_trueA1C0", 800, -400, 400, vtx0.Z(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_VertexZ_trueA1C0_Hybrid", 800, -400, 400, vtx1.Z(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_VertexZ_trueA1C0_Hybrid_TC" + std::to_string(timecut) + "_PC" + std::to_string(phicut), 800, -400, 400, vtx1.Z(), "A1C0True_QQQ");
      plotter->Fill2D("Benchmark_QQQ_VertexXY_trueA1C0_Hybrid", 200, -100, 100, 200, -100, 100, vtx1.X(), vtx1.Y(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C0_Hybrid", 600, -200, 200, pc_hybrid.Z(), "A1C0True_QQQ");
      plotter->Fill2D("Benchmark_QQQ_PCZ_trueA1C0_Hybrid_vs_qqqpczguess", 400, -200, 200, 400, -200, 200, pcz_guess_37, pc_hybrid.Z(), "A1C0True_QQQ");
      plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C0_Hybrid_minus_qqqpczguess", 400, -100, 100, pc_hybrid.Z() - pcz_guess_37, "A1C0True_QQQ");
    }
  }

  for (const auto &pcevent : PC_Events)
  {

    plotter->Fill2D("pcdEACAvg_vs_anodechannel", 24, 0, 23, 800, 0, 20000, pcevent.Anodech, (pcevent.Energy1 + pcevent.Energy2) / 2);
    plotter->Fill2D("pcdEACAvg_vs_cathodechannel", 24, 0, 23, 800, 0, 20000, pcevent.Cathodech, (pcevent.Energy1 + pcevent.Energy2) / 2);

    for (const auto &qqqevent : QQQ_Events)
    {
      plotter->Fill1D("dt_pcA_qqqR", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1, "Timing");
      plotter->Fill2D("dt_pcA_qqqR_vs_qqqRE", 640, -2000, 2000, 400, 0, 30, qqqevent.Time1 - pcevent.Time1, qqqevent.Energy1, "Timing");
      plotter->Fill1D("dt_pcC_qqqW", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2, "Timing");
      plotter->Fill2D("phiPC_vs_phiQQQ", 100, -200, 200, 100, -200, 200, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");

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
        {
          plotter->Fill2D("dE_E_AnodeQQQR_a1c2", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
          plotter->Fill2D("dE_E_CathodeQQQR_a1c2", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        }
        if (pcevent.multi1 == 2 && pcevent.multi2 == 1)
        {
          plotter->Fill2D("dE_E_AnodeQQQR_a2c1", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
          plotter->Fill2D("dE_E_CathodeQQQR_a2c1", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        }
        if (pcevent.multi1 == 1 && pcevent.multi2 == 1)
        {
          plotter->Fill2D("dE_E_Anodesx3B_a1c1", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
          plotter->Fill2D("dE_E_Cathodesx3B_a1c1", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        }
        if (pcevent.multi1 == 1 && pcevent.multi2 == 0)
        {
          plotter->Fill2D("dE_E_Anodesx3B_a1c0", 400, 0, 30, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, "PID_dE_E");
          plotter->Fill2D("dE_E_Cathodesx3B_a1c0", 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, "PID_dE_E");
        }
        if (phicut)
        {
          plotter->Fill2D("dE2_E_AnodeQQQR_TC1PC1_pidlow" + std::to_string(lowercut_cath), 400, 0, 30, 800, 0, 4000, qqqevent.Energy1, pcevent.Energy1 * sinTheta, "PID_dE_E");
          plotter->Fill2D("dE2_E_CathodeQQQW_TC1PC1_pidlow" + std::to_string(lowercut_cath), 400, 0, 30, 800, 0, 1000, qqqevent.Energy2, pcevent.Energy2 * sinTheta, "PID_dE_E");
          plotter->Fill2D("E_theta_zoomin_AnodeQQQR_TC1PC1_pidlow" + std::to_string(lowercut_cath), 60, 0, 30, 300, 0, 15, qqqTheta * 180 / M_PI, qqqevent.Energy1, "Kinematics_Angles");
        }

        plotter->Fill2D("dE2_E_AnodeQQQR_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 4000, qqqevent.Energy1, pcevent.Energy1 * sinTheta, "PID_dE_E");
        plotter->Fill2D("dE2_E_CathodeQQQR_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 1000, qqqevent.Energy2, pcevent.Energy2 * sinTheta, "PID_dE_E");
        plotter->Fill2D("dEC_vs_dEA_TC1_PC" + std::to_string(phicut), 800, 0, 40000, 800, 0, 10000, pcevent.Energy1, pcevent.Energy2, "PID_dE_E");
        plotter->Fill2D("qqqphi_vs_time", 2000, 0, 2000, 100, -200, 200, pcevent.Time1 * 1e-9, qqqevent.pos.Phi() * 180. / M_PI, "Timing");

        plotter->Fill1D("dt_pcA_qqqR_timecut", 640, -2000, 2000, qqqevent.Time1 - pcevent.Time1, "Timing");
        plotter->Fill1D("dt_pcC_qqqW_timecut", 640, -2000, 2000, qqqevent.Time2 - pcevent.Time2, "Timing");
        plotter->Fill2D("dE_theta_AnodeQQQR", 90, 0, 90, 400, 0, 20000, qqqTheta * 180 / M_PI, pcevent.Energy1, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_AnodeQQQR_zoomin", 60, 0, 30, 400, 0, 5000, qqqTheta * 180 / M_PI, pcevent.Energy1 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_AnodeQQQR", 90, 0, 90, 400, 0, 20000, qqqTheta * 180 / M_PI, pcevent.Energy1 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("phiPC_vs_phiQQQ_TimeCut", 100, -200, 200, 100, -200, 200, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");
        double pcz_guess_37 = 37. / TMath::Tan(qqqTheta) + source_vertex;

        if ((qqqevent.ch1) % 16 == 7)
          plotter->Fill2D("phiPC_vs_phiQQQ_TimeCut_allring8", 100, -200, 200, 100, -200, 200, qqqevent.pos.Phi() * 180 / M_PI, pcevent.pos.Phi() * 180 / M_PI, "Kinematics_Angles");

        plotter->Fill1D("phiQQQ_minus_phiPC_TimeCut_QQQ" + std::to_string(qqqevent.ch1 / 16), 180, -180, 180, qqqevent.pos.DeltaPhi(pcevent.pos) * 180 / M_PI, "Kinematics_Angles");

        plotter->Fill2D("Etot2_theta_AnodeQQQR", 75, 0, 90, 300, 0, 15, qqqTheta * 180 / M_PI, qqqevent.Energy1 + pcevent.Energy1 * anode_gain * sinTheta, "Kinematics_Angles");

        plotter->Fill2D("dE_theta_CathodeQQQR", 75, 0, 90, 800, 0, 10000, qqqTheta * 180 / M_PI, pcevent.Energy2, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_CathodeQQQR", 75, 0, 90, 800, 0, 10000, qqqTheta * 180 / M_PI, pcevent.Energy2 * sinTheta, "Kinematics_Angles");
        plotter->Fill2D("dE2_theta_CathodeQQQR_zoomin", 60, 0, 30, 800, 0, 3000, qqqTheta * 180 / M_PI, pcevent.Energy2 * sinTheta, "Kinematics_Angles");

        plotter->Fill2D("dE_phi_AnodeQQQR", 90, -180, 180, 800, 0, 40000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Phi() * 180 / M_PI, pcevent.Energy1, "Kinematics_Angles");
        plotter->Fill2D("dE_phi_CathodeQQQR", 90, -180, 180, 800, 0, 40000, (qqqevent.pos - TVector3(0, 0, source_vertex)).Phi() * 180 / M_PI, pcevent.Energy2, "Kinematics_Angles");

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
          plotter->Fill2D("dE3_E_CathodeQQQW_A1C2_TC1_PC" + std::to_string(phicut), 400, 0, 30, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV, "PID_dE_E");
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
          plotter->Fill2D("pcdEA_vs_pczfix", 600, -200, 200, 800, 0, 20000, pcz_fix, pcevent.Energy1, "PCdE_vs_Z");
          plotter->Fill2D("pcdEC_vs_pczfix", 600, -200, 200, 800, 0, 20000, pcz_fix, pcevent.Energy2, "PCdE_vs_Z");

          plotter->Fill2D("pcdEACSum_vs_pczfix", 600, -200, 200, 800, 0, 20000, pcz_fix, (pcevent.Energy1 + pcevent.Energy2) / 2, "PCdE_vs_Z");

          double path_length = pathLengthCm(qqqevent.pos, r_rhoMin_fix);
          double qqqEfix = evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, qqqevent.Energy1, path_length);
          double qqqEfix_p = evalEloss(MeV_to_cm_p_spl, cm_to_MeVp_spl, qqqevent.Energy1, path_length);

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
          plotter->Fill2D("pcdEACAvg_A1C2_vs_qqqE", 400, 0, 10, 1600, 0, 40000, qqqevent.Energy1, (pcevent.Energy1 + pcevent.Energy2) / 2);
        }

        // plotter->Fill2D("pcdEA_vs_qqqpczguess_A" + std::to_string(pcevent.multi1) + "C" + std::to_string(pcevent.multi2), 600, -200, 200, 800, 0, 20000, pcz_guess_37, pcevent.Energy1, "PCdE_vs_Z");
        plotter->Fill2D("pcdEA_vs_qqqpczguess", 600, -200, 200, 800, 0, 20000, pcz_guess_37, pcevent.Energy1, "PCdE_vs_Z");
        plotter->Fill2D("pcdEA_vs_pczguess", 600, -200, 200, 800, 0, 20000, pcz_guess_37, pcevent.Energy1, "PCdE_vs_Z");
        // plotter->Fill2D("pcdEA_vs_pczfix" + std::to_string(pcevent.multi1) + "A" + std::to_string(pcevent.multi1) + "C", 600, -200, 200, 800, 0, 20000, pcz_fix, pcevent.Energy1, "PCdE_vs_Z");
        // plotter->Fill2D("pcdEC_vs_qqqpczguess_A" + std::to_string(pcevent.multi1) + "C" + std::to_string(pcevent.multi2), 600, -200, 200, 800, 0, 20000, pcz_guess_37, pcevent.Energy2, "PCdE_vs_Z");
        plotter->Fill2D("pcdEC_vs_qqqpczguess", 600, -200, 200, 800, 0, 20000, pcz_guess_37, pcevent.Energy2, "PCdE_vs_Z");
        plotter->Fill2D("pcdEC_vs_pczguess", 600, -200, 200, 800, 0, 20000, pcz_guess_37, pcevent.Energy2, "PCdE_vs_Z");

        plotter->Fill2D("pcdEACAvg_vs_qqqpczguess", 600, -200, 200, 800, 0, 20000, pcz_guess_37, (pcevent.Energy1 + pcevent.Energy2) / 2, "PCdE_vs_Z");
        plotter->Fill2D("pcdEACAvg_vs_pczguess", 600, -200, 200, 800, 0, 20000, pcz_guess_37, (pcevent.Energy1 + pcevent.Energy2) / 2, "PCdE_vs_Z");

        plotter->Fill2D("pcdEACAvg_vs_qqqE", 400, 0, 10, 1600, 0, 40000, qqqevent.Energy1, (pcevent.Energy1 + pcevent.Energy2) / 2);

        // plotter->Fill2D("pcdEC_vs_pczfix" + std::to_string(pcevent.multi1) + "A" + std::to_string(pcevent.multi1) + "C", 800, 0, 20000, 600, -200, 200, pcevent.Energy2, pcz_fix, "PCdE_vs_Z");

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

          double aSumE_bm = std::get<1>(pw_tuple);
          double cSumE_bm = std::get<1>(cMaxWire);
          double ac_sum = aSumE_bm + cSumE_bm;
          double cfrac = (ac_sum > 0.0) ? cSumE_bm / ac_sum : -1.0;

          if (aSumE_bm > 0.0)
            plotter->Fill2D("Benchmark_QQQ_CmaxOverAnode_vs_phi", 90, -180, 180, 250, 0, 5, qqqevent.pos.Phi() * 180. / M_PI, cSumE_bm / aSumE_bm, "Benchmark_QQQ_ref");

          double qqq_wedge_pitch = (87.0 / 16.0) * (M_PI / 180.0);
          double qqq_ring_pitch = 48.0 / 16.0;
          double smeared_phi = qqqevent.pos.Phi() + rand.Uniform(-qqq_wedge_pitch / 2.0, qqq_wedge_pitch / 2.0);
          double smeared_rho = qqqevent.pos.Perp() + rand.Uniform(-qqq_ring_pitch / 2.0, qqq_ring_pitch / 2.0);

          TVector3 smeared_qqq_pos(smeared_rho * TMath::Cos(smeared_phi), smeared_rho * TMath::Sin(smeared_phi), qqqevent.pos.Z());

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

          auto doA1C1Model = [&](const std::string &tag, const TVector3 &si_point)
          {
            if (!a1c1Good || cfrac < 0.0)
              return;
            A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
            int side_status;
            SideChoice side = a1c1_pick_side(si_point, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
            const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
            double pcz_pick = best.pcz;
            if (!(best.inband && best.pitchok && side_status != 2))
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
              {
                double phi_deg = qqqevent.pos.Phi() * 180.0 / M_PI;
                double vz_resid = vtx_ref.Z() - source_vertex;
                plotter->Fill2D("Diag_QQQ_A1C2_vtxZ_resid_vs_phi", 180, -180, 180, 400, -100, 100, phi_deg, vz_resid, "Diag_XYoffset");
                plotter->Fill2D("Diag_Combined_A1C2_vtxZ_resid_vs_phi", 90, -180, 180, 400, -100, 100, phi_deg, vz_resid, "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_vtxXY", 200, -15, 15, 200, -15, 15, vtx_ref.X(), vtx_ref.Y(), "Diag_XYoffset");
                plotter->Fill2D("Diag_Combined_A1C2_time_vs_phi", 2000, 0, 2000, 90, -180, 180, pcevent.Time1 * 1e-9, phi_deg, "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_time_vs_phi", 2000, 0, 2000, 90, -180, 180, pcevent.Time1 * 1e-9, phi_deg, "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_T_vs_vtxX", 2000, 0, 2000, 200, -15, 15, pcevent.Time1 * 1e-9, vtx_ref.X(), "Diag_XYoffset");
                plotter->Fill2D("Diag_QQQ_A1C2_T_vs_vtxY", 2000, 0, 2000, 200, -15, 15, pcevent.Time1 * 1e-9, vtx_ref.Y(), "Diag_XYoffset");
              }

              doA1C1("A1C1", qqqevent.pos, false);
              doAnodeOnly("A1C0", qqqevent.pos.Phi(), qqqevent.pos, false);
              doA1C1("A1C1_Hyb", smeared_qqq_pos);
              doAnodeOnly("A1C0_Hyb", smeared_phi, smeared_qqq_pos);

              doA1C1Model("A1C1_Cfrac", qqqevent.pos);

              {
                double pcz_a1c0 = pwinstance.getClosestWirePosAtWirePhi(apwire_bm, qqqevent.pos.Phi()).Z();
                double theta_ref = (qqqevent.pos - TVector3(0, 0, vtx_ref.Z())).Theta() * 180. / M_PI;
                plotter->Fill2D("Benchmark_QQQ_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200, theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_QQQ_ref");
                plotter->Fill2D("Benchmark_PCZ_A1C0_minus_ref_vs_theta", 180, 0, 180, 400, -200, 200, theta_ref, pcz_a1c0 - pcz_ref, "Benchmark_AnodeOnly");

                double phi_deg_a = qqqevent.pos.Phi() * 180.0 / M_PI;
                plotter->Fill2D("Diag_QQQ_A1C0_zresid_vs_phi", 90, -180, 180, 200, -100, 100, phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
                plotter->Fill2D("Diag_Combined_A1C0_zresid_vs_phi", 90, -180, 180, 200, -100, 100, phi_deg_a, pcz_a1c0 - pcz_ref, "Diag_XYoffset");
              }

              if (a1c1Good && cfrac >= 0.0)
              {
                plotter->Fill1D("Benchmark_QQQ_A1C1_cfrac", 220, -0.05, 1.05, cfrac, "Benchmark_QQQ_ref");
                plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_ref", 400, -200, 200, 220, -0.05, 1.05, pcz_ref, cfrac, "Benchmark_QQQ_ref");
                plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_qqqpczguess", 400, -200, 200, 220, -0.05, 1.05, pcz_guess_37, cfrac, "Benchmark_QQQ_ref");

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
                    plotter->Fill2D(name, 240, -1.2, 1.2, 220, -0.05, 1.05, (truth - zp) / TMath::Abs(znb - zp), cfrac, "Benchmark_QQQ_ref");
                };
                fillCfracS("Benchmark_QQQ_A1C1_cfrac_vs_s", pcz_ref);
                fillCfracS("Benchmark_QQQ_A1C1_cfrac_vs_s_qqqpczguess", pcz_guess_37);

                for (int i = 0; i < 7; ++i)
                {
                  if (pcz_ref <= zg[i] && pcz_ref > zg[i + 1])
                  {
                    double zc = 0.5 * (zg[i] + zg[i + 1]);
                    double half = 0.5 * (zg[i] - zg[i + 1]);
                    if (half > 0.0)
                      plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05, TMath::Abs(pcz_ref - zc) / half, cfrac, "Benchmark_QQQ_ref");
                    break;
                  }
                }

                plotter->Fill2D("Benchmark_QQQ_A1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05, aSumE_bm, cfrac, "Benchmark_QQQ_ref");
                if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                  plotter->Fill2D("Benchmark_QQQ_A1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0, 1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "Benchmark_QQQ_ref");

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
                    plotter->Fill2D("Benchmark_QQQ_A1C1_cellsel_confusion", 7, 0, 7, 7, 0, 7, cell_truth + 0.5, sm_cell + 0.5, "Benchmark_QQQ_ref");
                    plotter->Fill1D("Benchmark_QQQ_A1C1_cellsel_misclass", 2, 0, 2, wrong ? 1.0 : 0.0, "Benchmark_QQQ_ref");
                    plotter->Fill2D("Benchmark_QQQ_A1C1_cellsel_misclass_vs_cell", 7, 0, 7, 2, 0, 2, cell_truth + 0.5, wrong ? 1.0 : 0.0, "Benchmark_QQQ_ref");

                    double zc = 0.5 * (a1c1_zg[cell_truth] + a1c1_zg[cell_truth + 1]);
                    double half = 0.5 * (a1c1_zg[cell_truth] - a1c1_zg[cell_truth + 1]);

                    plotter->Fill2D("AnodeEnergy_vs_CellQQQ", 120, 0, 1.2, 800, 0, 40000, 1 - TMath::Abs(pcz_ref - zc) / half, pcevent.Energy1);
                    plotter->Fill2D("CathodeEnergy_vs_CellQQQ", 120, 0, 1.2, 800, 0, 40000, TMath::Abs(pcz_ref - zc) / half, pcevent.Energy2);
                    plotter->Fill2D("FracEnergy_vs_CellQQQ", 120, 0, 1.2, 1200, 0, 20, TMath::Abs(pcz_ref - zc) / half, pcevent.Energy2 / pcevent.Energy1);
                    plotter->Fill2D("SumEnergy_vs_CellQQQ", 120, 0, 1.2, 800, 0, 40000, TMath::Abs(pcz_ref - zc) / half, (pcevent.Energy2 + pcevent.Energy1) / 2);

                    if (half > 0.0)
                    {
                      plotter->Fill2D("Benchmark_QQQ_A1C1_cellsel_misclass_vs_fold", 120, 0, 1.2, 2, 0, 2, TMath::Abs(pcz_ref - zc) / half, wrong ? 1.0 : 0.0, "Benchmark_QQQ_ref");
                      plotter->Fill2D("Benchmark_QQQ_A1C1_cfracUsed_vs_fold", 120, 0, 1.2, 220, -0.05, 1.05, TMath::Abs(pcz_ref - zc) / half, sm.cfrac_used, "Benchmark_QQQ_ref");
                      if (aSumE_bm > 0.0)
                      {
                        plotter->Fill2D("Benchmark_QQQ_A1C1_cfracUsed_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05,
                                        aSumE_bm, sm.cfrac_used, "Benchmark_QQQ_ref");
                      }
                    }
                  }
                }
              }
            }
          }

          else if (pcevent.multi1 >= 1 && pcevent.multi2 == 1 && a1c1Good)
          {
            double pcz_raw = xo_a1c1.Z();
            TVector3 vtx_raw = vertexFrom(qqqevent.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_raw));
            fillSuite("trueA1C1", pcz_raw, vtx_raw, "A1C1True_QQQ");
            plotter->Fill2D("Benchmark_QQQ_PCZ_trueA1C1_vs_qqqpczguess", 400, -200, 200, 400, -200, 200, pcz_guess_int, pcz_raw, "Benchmark_QQQ_trueA1C1");
            plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C1_minus_qqqpczguess", 400, -100, 100, pcz_raw - pcz_guess_int, "Benchmark_QQQ_trueA1C1");

            if (cfrac >= 0.0)
            {
              A1C1Sol s = a1c1_solve(cfrac, xo_a1c1.Z(), std::get<0>(cMaxWire), aSumE_bm, std::get<0>(aMaxWire));
              int side_status;
              SideChoice side = a1c1_pick_side(qqqevent.pos, xo_a1c1.X(), xo_a1c1.Y(), s.pcz_lo, s.pcz_hi, side_status);
              const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
              int cell = best.cell;
              double f = best.f;
              double pcz_cf = best.pcz;
              bool valid = (side_status != 2);
              bool cfrac_valid = (valid && best.inband && best.pitchok);
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_sideStatus", 4, -1, 3, side_status + 0.5, "Benchmark_QQQ_trueA1C1");

              TVector3 vtx_cf = vertexFrom(qqqevent.pos, TVector3(xo_a1c1.X(), xo_a1c1.Y(), pcz_cf));
              fillSuite(valid ? "trueA1C1_Cfrac" : "trueA1C1_Cfrac_invalid", pcz_cf, vtx_cf, "A1C1True_QQQ");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_cfrac", 220, -0.05, 1.05, cfrac, "Benchmark_QQQ_trueA1C1");
              plotter->Fill2D("Benchmark_QQQ_trueA1C1_cfrac_vs_anodeE", 400, 0, 40000, 220, -0.05, 1.05, aSumE_bm, cfrac, "Benchmark_QQQ_trueA1C1");
              if (aSumE_bm > 0.0 && cfrac > 0.0 && cfrac < 1.0)
                plotter->Fill2D("Benchmark_QQQ_trueA1C1_r_vs_invAnodeE", 200, 0, 0.0004, 200, 0, 2.0,
                                1.0 / aSumE_bm, cfrac / (1.0 - cfrac), "Benchmark_QQQ_trueA1C1");
              plotter->Fill2D("Benchmark_QQQ_trueA1C1_cfrac_vs_cell", 7, 0, 7, 220, -0.05, 1.05, cell + 0.5, cfrac, "Benchmark_QQQ_trueA1C1");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_f", 260, -1.5, 2.5, f, "Benchmark_QQQ_trueA1C1");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_valid", 2, 0, 2, valid ? 1.0 : 0.0, "Benchmark_QQQ_trueA1C1");
              int reason;
              if (cell < 0 || cell > 6 || a1c1_k_cell[cell] <= 0.0)
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
              if (valid)
                plotter->Fill1D("Benchmark_QQQ_trueA1C1_validreason", 3, 0, 3, reason + 0.5, "Benchmark_QQQ_trueA1C1");
              plotter->Fill1D("Benchmark_QQQ_trueA1C1_band", 2, 0, 2, s.band + 0.5, "Benchmark_QQQ_trueA1C1");
              if (valid)
                plotter->Fill1D("Benchmark_QQQ_trueA1C1_band_valid", 2, 0, 2, s.band + 0.5, "Benchmark_QQQ_trueA1C1");
              if (valid)
              {
                plotter->Fill1D("Benchmark_QQQ_PCZ_trueA1C1_Cfrac_minus_qqqpczguess_DIAG", 400, -100, 100, pcz_cf - pcz_guess_int, "Benchmark_QQQ_trueA1C1");
                plotter->Fill2D("Benchmark_QQQ_PCZ_trueA1C1_Cfrac_vs_qqqpczguess_DIAG", 400, -200, 200, 400, -200, 200, pcz_guess_int, pcz_cf, "Benchmark_QQQ_trueA1C1");
              }
            }

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
        plotter->Fill2D("pctheta_vs_qqqtheta_sv", 180, -200, 200, 180, -200, 200, qqqTheta * 180 / M_PI, (pcevent.pos - TVector3(0, 0, source_vertex)).Theta() * 180 / M_PI, "Kinematics_Angles");
        plotter->Fill2D("pctheta_vs_qqqtheta_rmz", 180, -200, 200, 180, -200, 200, (qqqevent.pos - TVector3(0, 0, r_rhoMin.Z())).Theta() * 180 / M_PI, (pcevent.pos - TVector3(0, 0, r_rhoMin.Z())).Theta() * 180 / M_PI, "Kinematics_Angles");
        plotter->Fill2D("pctheta_vs_qqqtheta_rm", 180, -200, 200, 180, -200, 200, (qqqevent.pos - r_rhoMin).Theta() * 180 / M_PI, (pcevent.pos - r_rhoMin).Theta() * 180 / M_PI, "Kinematics_Angles");
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

  if (anodeHits.size() > 0 && cathodeHits.size() > 0)
  {
#ifdef RAW_HISTOS
    plotter->Fill2D("AHits_vs_CHits_NA" + std::to_string(hasNeighbourAnodes), 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");
    plotter->Fill2D("AHits_vs_CHits_NC" + std::to_string(hasNeighbourCathodes), 12, 0, 11, 6, 0, 5, anodeHits.size(), cathodeHits.size(), "hRawPC");

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
        plotter->Fill2D("PolarAngle_Vs_QQQWedge" + std::to_string(qqqID), 360, -200, 200, 16, 0, 16, phi, chWedge, "hPCQQQ");
        // plotter->Fill2D("EdE_PC_vs_QQQ_timegate_ls1000"+std::to_string())

        plotter->Fill2D("PC_Z_vs_QQQRing_Det" + std::to_string(qqqID), 600, -300, 300, 16, 0, 16, anodeIntersection.Z(), chRing, "hPCQQQ");
        // double theta = -TMath::Pi() / 2 + 2 * TMath::Pi() / 16 / 4. * (qqq.id[i] * 16 + chWedge + 0.5);
        // double rho = 50. + 40. / 16. * (chRing + 0.5);

        for (int k = 0; k < pc.multi; k++)
        {
          if (pc.index[k] >= 24)
            continue;

          double sinTheta = TMath::Sin((anodeIntersection - TVector3(0, 0, 90.0)).Theta());

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
    plotter->Fill2D("AnodeSumE_Vs_Cathode_Sum_Energy", 800, 0, 20000, 800, 0, 10000, aESum, cESum, "hGMPC");
    plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_TC" + std::to_string(PCQQQTimeCut) + "_PC" + std::to_string(PCQQQPhiCut), 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
    if (aEMax > 0)
    {
      double ratio = cEMax / aEMax;
      std::string folder = "Diagnostics_CMax";

      plotter->Fill2D("CMax_over_Anode_vs_Z", 600, -300, 300, 200, 0, 2.0, anodeIntersection.Z(), ratio, folder);
      plotter->Fill2D("CMax_over_Anode_vs_AnodeID", 24, 0, 24, 200, 0, 2.0, aIDMax, ratio, folder);
      plotter->Fill2D("CMax_over_Anode_vs_CathodeID", 24, 0, 24, 200, 0, 2.0, cIDMax, ratio, folder);
    }

    if (PCQQQTimeCut && PCQQQPhiCut)
    {
      plotter->Fill2D("AnodeSumE_Vs_Cathode_Max_Energy_TC" + std::to_string(PCQQQTimeCut) + "_PC" + std::to_string(PCQQQPhiCut) + "_cMax" + std::to_string(cIDMax), 800, 0, 20000, 800, 0, 10000, aESum, cEMax, "hGMPC");
    }
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
  double initial_energy = 6.89;

  Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, initial_energy / mass_1H);
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
        TVector3 pc_hybrid = a1c0_hybrid_pcz(apwire, qqqevent.pos.Phi(), true, dither_sigma_c0 / 2.0, rand);
        TVector3 r_rhoMin_fix = beamVertex(qqqevent.pos, pc_hybrid - qqqevent.pos);

        double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
        double sinTheta2 = TMath::Sin(theta_q);

        if (beamPerp(r_rhoMin_fix) > 6.0)
          continue;
        if (r_rhoMin_fix.Z() < -173.6 || r_rhoMin_fix.Z() > 100)
          continue;
        if (!phicut)
          continue;
        plotter->Fill1D("dt_anode_ainterp_qqq_gated", 800, -2000, 2000, qqqevent.Time1 - apTSMaxE, "ainterp_noc");
        plotter->Fill2D("dt_anode_ainterp_qqq_gated_vs_qqqE", 800, -2000, 2000, 800, 0, 10, qqqevent.Time1 - apTSMaxE, qqqevent.Energy1, "ainterp_noc");
        // plotter->Fill2D("dEa_ainterp_Eqqq_TC1_ignC_a" + std::to_string(acluster.size()), 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, apSumE, "ainterp_noc");
        // plotter->Fill2D("pcPhi_ainterp_qqqPhi_TC1_ignC_a" + std::to_string(acluster.size()), 120, -200, 200, 120, -200, 200, pc_closest.Phi() * 180. / M_PI, qqqevent.pos.Phi() * 180. / M_PI, "ainterp_noc");
        // plotter->Fill2D("pcZ_ainterp_qqqZ_TC1_ignC_a" + std::to_string(acluster.size()) + "_PC" + std::to_string(phicut), 300, -100, 200, 400, -200, 200, qqqevent.pos.Z(), pc_hybrid.Z(), "ainterp_noc");

        // plotter->Fill2D("pcZ_ainterp_qqqpczguess_TC1_ignC_a"+std::to_string(acluster.size()),300,-100,200,400,-200,200,pczguess,pc_hybrid.Z(),"ainterp_noc");
        // plotter->Fill2D("dEa3_ainterp_Eqqq_TC1_ignC_a" + std::to_string(acluster.size()) + "_PC" + std::to_string(phicut), 1200, 0, 30, 800, 0, 30000, qqqevent.Energy1, apSumE * sinTheta2 * 3., "ainterp_noc");

        // plotter->Fill2D("vertexZ_ainterp_qqqZ_TC1_ignC_a" + std::to_string(acluster.size()), 300, -100, 200, 800, -400, 400, qqqevent.pos.Z(), r_rhoMin_fix.Z(), "ainterp_noc");
        // plotter->Fill1D("vertexZ1d_ainterp_qqqZ_TC1_ignC_a" + std::to_string(acluster.size()), 800, -400, 400, r_rhoMin_fix.Z(), "ainterp_noc");
        // plotter->Fill2D("vertexXY_ainterp_TC1_ignC_a" + std::to_string(acluster.size()), 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), "ainterp_noc");

        double path_length_q = pathLengthCm(qqqevent.pos, r_rhoMin_fix);
        double qqqEfix = evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, qqqevent.Energy1, path_length_q);
        double qqqEx = apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI);
        plotter->Fill1D("pmisc_ow_Ex_from_alpha", 600, -10, 10, qqqEx, "ainterp_noc");
        plotter->Fill1D("pmisc_ow_Ef_from_alpha", 600, 0, 20, qqqEfix, "ainterp_noc");
        plotter->Fill2D("pmisc_ow_Ex_vs_theta_qqq", 100, 0, 180, 800, 0, 20, theta_q * 180 / M_PI, qqqEx, "ainterp_noc");
        plotter->Fill2D("pmisc_ow_Ef_vs_theta_qqq", 100, 0, 180, 800, 0, 20, theta_q * 180 / M_PI, qqqEfix, "ainterp_noc");
        plotter->Fill2D("pmisc_ow_VertexReconZ_vs_Ef", 800, -400, 400, 800, 0, 20, r_rhoMin_fix.Z(), qqqEfix, "ainterp_noc");
      }
    }
  } // end QQQEvents loop
}

void protonMiscHistograms(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events)
{
  // consider the 'proton-like' QQQ branch seen in a,p data
  static TRandom3 rand(0);
  double initial_energy = 6.89;

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

      if (pcevent.multi2 == 1 && pcevent.Energy2 > 1400)
      {
        const std::string wcat = a1c1_missing_neighbor(pcevent.Anodech, pcevent.Cathodech) ? "_missingw" : "_true1w";
        auto fillCmp = [&](double pcz, const std::string &m)
        {
          TVector3 x2(pcevent.pos.X(), pcevent.pos.Y(), pcz);
          TVector3 rv = beamVertex(qqqevent.pos, x2 - qqqevent.pos);
          if (beamPerp(rv) > 6.0 || rv.Z() < -173.6 || rv.Z() > 100)
            return;
          double th = (qqqevent.pos - rv).Theta();
          double pl = pathLengthCm(qqqevent.pos, rv);
          double Ef = evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, qqqevent.Energy1, pl);
          double beam_pl_cmp = TMath::Abs(rv.Z() - z_entrance) * 0.1;
          double beam_E_cmp = evalElossForward(MeV_to_cm_p_spl, cm_to_MeVp_spl, initial_energy, beam_pl_cmp);
          beam_E_cmp = applyTaFoilEloss(beam_E_cmp, rv.Z());
          if (beam_E_cmp <= 0.0)
            beam_E_cmp = 0.001;
          Kinematics apkin_a_cmp(mass_1H, mass_4He, mass_4He, mass_1H, beam_E_cmp / mass_1H);
          double Ex = apkin_a_cmp.getExc(Ef, th * 180 / M_PI);
          std::string lbl = "proton+misc_a1c1cmp";
          // fill "all" (existing names) plus the wire-topology split (_true1w/_missingw)
          for (const std::string &w : {std::string(""), wcat})
          {
            plotter->Fill1D("pmisc_a1c1cmp_pcz_" + m + w, 600, -300, 300, pcz, lbl);
            plotter->Fill1D("pmisc_a1c1cmp_Ex_" + m + w, 200, -10, 10, Ex, lbl);
            plotter->Fill1D("pmisc_a1c1cmp_VertexZ_" + m + w, 800, -400, 400, rv.Z(), lbl);
            plotter->Fill2D("pmisc_a1c1cmp_VertexZ_vs_Ef_" + m + w, 800, -400, 400, 800, 0, 20, rv.Z(), Ef, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_VertexZ_vs_Ex_" + m + w, 800, -400, 400, 400, -10, 10, rv.Z(), Ex, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_phi_vs_Ef_" + m + w, 90, -180, 180, 800, 0, 20, qqqevent.pos.Phi() * 180 / M_PI, Ef, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_phi_vs_Ex_" + m + w, 90, -180, 180, 800, -10, 10, qqqevent.pos.Phi() * 180 / M_PI, Ex, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_Ef_vs_theta_" + m + w, 100, 0, 180, 800, 0, 20, th * 180 / M_PI, Ef, lbl);
            plotter->Fill2D("pmisc_a1c1cmp_Ex_vs_theta_" + m + w, 100, 0, 180, 800, -10, 10, th * 180 / M_PI, Ex, lbl);
          }
        };

        fillCmp(pcz_dith, "dither"); // method 1: Gaussian dither (main-flow value)
        double ac = pcevent.Energy1 + pcevent.Energy2;
        double cfrac = (ac > 0.0) ? pcevent.Energy2 / ac : -1.0;
        if (cfrac >= 0.0)
        {
          std::vector<std::tuple<int, double, double>> aOne = {std::make_tuple(pcevent.Anodech, 1.0, 0.0)};
          auto apw = pwinstance.GetPseudoWire(aOne, "ANODE");
          A1C1Sol s = a1c1_solve(cfrac, pcevent.pos.Z(), pcevent.Cathodech, pcevent.Energy1, pcevent.Anodech);
          // beam-axis 2-hypothesis side test (crossover = PC point, Si = qqq hit).
          int side_status;
          SideChoice side = a1c1_pick_side(qqqevent.pos, pcevent.pos.X(), pcevent.pos.Y(), s.pcz_lo, s.pcz_hi, side_status);
          const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
          double pcz_pick = best.pcz;
          // cfrac_all = beam-axis pick for ALL events; "cfrac" = inband + on-axis.
          fillCmp(pcz_pick, "cfrac_all");
          if (best.inband && side_status != 2)
          {
            fillCmp(pcz_pick, "cfrac");
            plotter->Fill2D("pmisc_a1c1cmp_pcz_cfrac_vs_dither", 600, -300, 300, 600, -300, 300, pcz_dith, pcz_pick, "proton+misc_a1c1cmp");
          }
        }
      }

      TVector3 x2f(pcevent.pos.X(), pcevent.pos.Y(), pcz_fix);
      TVector3 x1(qqqevent.pos);
      TVector3 r_rhoMin_fix = beamVertex(x1, x2f - x1);
      double vertex_z = r_rhoMin_fix.Z();
      // double theta_q = (qqqevent.pos - TVector3(0,0,vertex_z)).Theta();
      double theta_q = (qqqevent.pos - r_rhoMin_fix).Theta();
      double sinTheta_customV = TMath::Sin(theta_q);
      // if(beamPerp(r_rhoMin_fix)>6) continue;
      bool cathode_alpha_select = (pcevent.Energy2 > 1400);
      if (vertex_z < -173.6 || vertex_z > 100)
        continue;

      double beam_path_length_q = TMath::Abs(vertex_z - z_entrance) * 0.1;
      double beam_energy_at_vertex_q = evalElossForward(MeV_to_cm_p_spl, cm_to_MeVp_spl, initial_energy, beam_path_length_q);
      beam_energy_at_vertex_q = applyTaFoilEloss(beam_energy_at_vertex_q, vertex_z);
      plotter->Fill2D("pmisc_BeamEnergy_vs_VertexZ", 800, -400, 400, 400, 0, initial_energy, vertex_z, beam_energy_at_vertex_q, "qqq");
      if (beam_energy_at_vertex_q <= 0.0)
        beam_energy_at_vertex_q = 0.001;
      Kinematics apkin_a(mass_1H, mass_4He, mass_4He, mass_1H, beam_energy_at_vertex_q / mass_1H);

      PCPath pa_pp = pcPath(r_rhoMin_fix, qqqevent.pos);
      bool pa_have_seg = pa_pp.ok;
      double pa_anode_cm = pa_pp.anode_cm, pa_cathode_cm = pa_pp.cathode_cm;
      double pa_dl_cm = pa_have_seg ? (pa_anode_cm - pa_cathode_cm) : 0.0;
      double pa_dist_mm = (qqqevent.pos - r_rhoMin_fix).Mag();
      double pa_pathfraction = (pa_dist_mm > 0.0) ? pa_dl_cm * 10.0 / pa_dist_mm : 0.0;
      double pcz_guess_int = z_to_crossover_rho(pcevent.pos.Z()) /
                                 TMath::Tan((qqqevent.pos - TVector3(0, 0, source_vertex)).Theta()) +
                             source_vertex;

      // What's below: radial cut, time coincident, phi-correlated events with possible energy selection applied to both E_si and dE_Anodes
      auto plot_with_tag = [&](std::string tag = "")
      {
        std::string pmlabel = "proton+misc" + tag;
        plotter->Fill2D("pmisc_dE_E_AnodeQQQ" + tag, 400, 0, 10, 800, 0, 40000, qqqevent.Energy1, pcevent.Energy1, pmlabel);
        plotter->Fill2D("pmisc_dE_E_CathodeQQQ" + tag, 400, 0, 10, 800, 0, 10000, qqqevent.Energy1, pcevent.Energy2, pmlabel);
        plotter->Fill2D("pmisc_dE3_E_AnodeQQQ" + tag, 400, 0, 10, 400, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * sinTheta_customV * 3., pmlabel);
        plotter->Fill2D("pmisc_dE3_E_CathodeQQQ" + tag, 400, 0, 10, 400, 0, 10000, qqqevent.Energy1, pcevent.Energy2 * sinTheta_customV, pmlabel);
        plotter->Fill2D("pmisc_dPhi_QQQ_PC" + tag, 100, -200, 200, 100, -200, 200, pcevent.pos.Phi() * 180 / M_PI, qqqevent.pos.Phi() * 180 / M_PI, pmlabel);
        plotter->Fill1D("pmisc_dt_Anode_QQQ_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, pcevent.Time1 - qqqevent.Time1, pmlabel);
        plotter->Fill1D("pmisc_dt_Cathode_QQQ" + tag, 600, -2000, 2000, pcevent.Time2 - qqqevent.Time1, pmlabel);
        plotter->Fill2D("pmisc_dt_Anode_E_QQQ_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time1 - qqqevent.Time1, qqqevent.Energy1, pmlabel);
        plotter->Fill2D("pmisc_dt_AnodeQQQ_vsPCPhi" + tag, 600, -2000, 2000, 100, -200, 200, pcevent.Time1 - qqqevent.Time1, pcevent.pos.Phi() * 180. / M_PI, pmlabel);
        plotter->Fill2D("pmisc_dt_Cathode_E_QQQ" + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time2 - qqqevent.Time1, qqqevent.Energy1, pmlabel);
        plotter->Fill2D("pmisc_dt_CathodeQQQ_vsPCPhi" + tag, 600, -2000, 2000, 100, -200, 200, pcevent.Time2 - qqqevent.Time1, pcevent.pos.Phi() * 180. / M_PI, pmlabel);
        plotter->Fill1D("pmisc_pczfix" + tag, 600, -300, 300, pcz_fix, pmlabel);

        double path_length_q = pathLengthCm(qqqevent.pos, r_rhoMin_fix);
        double qqqEfix = evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, qqqevent.Energy1, path_length_q);
        double qqqEx = apkin_a.getExc(qqqEfix, theta_q * 180 / M_PI);

        if (pcevent.multi2 == 2)
        {
          plotter->Fill1D("pmisc_pcz" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);
          plotter->Fill1D("pmisc_pcz2" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);
        }
        if (pcevent.multi2 == 1)
        {
          plotter->Fill1D("pmisc_pcz" + tag, 600, -300, 300, pcz_fix, pmlabel);
          plotter->Fill1D("pmisc_pcz1" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);
        }

        if (tag == "_cathode_alphas")
        {
          plotter->Fill1D("pmisc_Ex_from_alpha", 800, -10, 10, qqqEx, pmlabel);
          plotter->Fill2D("pmisc_Ex_vs_theta_qqq", 100, 0, 180, 800, -10, 10, theta_q * 180 / M_PI, qqqEx, pmlabel);
          plotter->Fill2D("pmisc_VertexReconZ_vs_Ex", 800, -400, 400, 800, -10, 10, vertex_z, qqqEx, pmlabel);
        }
        else
          qqqEfix = evalEloss(MeV_to_cm_p_spl, cm_to_MeVp_spl, qqqevent.Energy1, path_length_q);
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

        if (pa_have_seg)
        {
          // Per-electrode Eloss-corrected dE across the PC gas (proton table).
          double E_an = evalEloss(MeV_to_cm_p_spl, cm_to_MeVp_spl, qqqevent.Energy1, pa_anode_cm);
          double E_ca = evalEloss(MeV_to_cm_p_spl, cm_to_MeVp_spl, qqqevent.Energy1, pa_cathode_cm);
          plotter->Fill2D("pmisc_dEa_guess_vs_dEa" + tag, 400, 0, 5, 800, 0, 40000, E_an - E_ca, pcevent.Energy1, pmlabel);
        }
        if (pa_pathfraction > 0.0)
        {
          plotter->Fill2D("pmisc_dEapf_E_AnodeQQQ" + tag, 400, 0, 10, 400, 0, 40000, qqqevent.Energy1, pcevent.Energy1 / (20 * pa_pathfraction), pmlabel);
          plotter->Fill2D("pmisc_dEapf_Theta_TC1" + tag, 180, 0, 180, 800, 0, 40000, theta_q * 180 / M_PI, pcevent.Energy1 / (20 * pa_pathfraction), pmlabel);
        }
        if (pa_dl_cm > 0.0)
        {
          plotter->Fill2D("pmisc_dE4_E_AnodeQQQ" + tag, 400, 0, 10, 400, 0, 40000, qqqevent.Energy1, pcevent.Energy1 * 1.72 / pa_dl_cm, pmlabel);
          plotter->Fill2D("pmisc_dE4_Theta_TC1_" + tag, 180, 0, 180, 800, 0, 40000, theta_q * 180 / M_PI, pcevent.Energy1 * 1.72 / pa_dl_cm, pmlabel);
          plotter->Fill2D("pmisc_dE4_Rho_TC1_" + tag, 100, 0, 40, 400, 0, 40000, r_rhoMin_fix.Perp(), pcevent.Energy1 * 1.72 / pa_dl_cm, pmlabel);
        }
        plotter->Fill2D("pmisc_pcz_vs_pczguess" + tag, 600, -300, 300, 600, -300, 300, pcz_guess_int, pcevent.pos.Z(), pmlabel);
      };

      plot_with_tag();
      if (cathode_alpha_select)
        plot_with_tag("_cathode_alphas");
      else
        plot_with_tag("_cathode_protons");

      // plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(qqqEfix,theta_s*180/M_PI),pmlabel);

    } // end PCEvents loop
  } // end QQQEvents loop
}

void protonMiscHistograms_sx3(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events)
{
  // consider the 'proton-like' QQQ branch seen in a,p data
  static TRandom3 rand(0); // seeded once (random seed via TUUID), not per call
  double initial_energy = 6.89;

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
      TVector3 r_rhoMin_fix = beamVertex(x1, x2f - x1);
      double vertex_z = r_rhoMin_fix.Z();
      // double theta_q = (sx3event.pos - TVector3(0,0,vertex_z)).Theta();

      if (beamPerp(r_rhoMin_fix) > 10.0)
        continue;
      if (vertex_z < -173.6 || vertex_z > 100)
        continue; // same beam-region acceptance as the QQQ branch
      double theta_s = (sx3event.pos - r_rhoMin_fix).Theta();
      double sinTheta_customV = TMath::Sin(theta_s);
      bool cathode_alpha_select = (pcevent.Energy2 > 1400);
      double beam_path_length_s = TMath::Abs(vertex_z - z_entrance) * 0.1;
      double beam_energy_at_vertex_s = evalElossForward(MeV_to_cm_p_spl, cm_to_MeVp_spl, initial_energy, beam_path_length_s);
      beam_energy_at_vertex_s = applyTaFoilEloss(beam_energy_at_vertex_s, vertex_z);
      plotter->Fill2D("pmiscs_BeamEnergy_vs_VertexZ", 800, -400, 400, 400, 0, initial_energy, vertex_z, beam_energy_at_vertex_s, "sx3");
      if (beam_energy_at_vertex_s <= 0.0)
        beam_energy_at_vertex_s = 0.001;
      Kinematics apkin_a_s(mass_1H, mass_4He, mass_4He, mass_1H, beam_energy_at_vertex_s / mass_1H);

      auto plot_with_tag = [&](std::string tag = "")
      {
        std::string pmlabel = "proton+miscsx3" + tag;
        plotter->Fill2D("pmiscs_dE_E_Anodesx3" + tag, 400, 0, 10, 800, 0, 40000, sx3event.Energy1, pcevent.Energy1, pmlabel);
        plotter->Fill2D("pmiscs_dE_E_Cathodesx3" + tag, 400, 0, 10, 800, 0, 10000, sx3event.Energy1, pcevent.Energy2, pmlabel);
        plotter->Fill2D("pmiscs_dE3_E_Anodesx3" + tag, 400, 0, 10, 400, 0, 40000, sx3event.Energy1, pcevent.Energy1 * sinTheta_customV * 3., pmlabel);
        plotter->Fill2D("pmiscs_dE3_E_Cathodesx3" + tag, 400, 0, 10, 400, 0, 10000, sx3event.Energy1, pcevent.Energy2 * sinTheta_customV, pmlabel);
        plotter->Fill2D("pmiscs_dPhi_sx3_PC" + tag, 100, -200, 200, 100, -200, 200, pcevent.pos.Phi() * 180 / M_PI, sx3event.pos.Phi() * 180 / M_PI, pmlabel);
        plotter->Fill1D("pmiscs_dt_Anode_sx3_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, pcevent.Time1 - sx3event.Time1, pmlabel);
        plotter->Fill1D("pmiscs_dt_Cathode_sx3" + tag, 600, -2000, 2000, pcevent.Time2 - sx3event.Time1, pmlabel);
        plotter->Fill2D("pmiscs_dt_Anode_E_sx3_PC" + std::to_string(phicut) + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time1 - sx3event.Time1, sx3event.Energy1, pmlabel);
        plotter->Fill2D("pmiscs_dt_Cathode_E_sx3" + tag, 600, -2000, 2000, 400, 0, 10, pcevent.Time2 - sx3event.Time1, sx3event.Energy1, pmlabel);
        plotter->Fill2D("pmiscs_dt_Cathodesx3_vsPCPhi" + tag, 600, -2000, 2000, 100, -200, 200, pcevent.Time2 - sx3event.Time1, pcevent.pos.Phi() * 180. / M_PI, pmlabel);
        plotter->Fill1D("pmiscs_pczfix" + tag, 600, -300, 300, pcz_fix, pmlabel);
        plotter->Fill1D("pmiscs_pcz" + tag, 600, -300, 300, pcevent.pos.Z(), pmlabel);

        double path_length_s = pathLengthCm(sx3event.pos, r_rhoMin_fix);
        // alpha Eloss table for cathode-alpha events, proton otherwise (matches QQQ).
        double sx3Efix = cathode_alpha_select
                             ? evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, sx3event.Energy1, path_length_s)
                             : evalEloss(MeV_to_cm_p_spl, cm_to_MeVp_spl, sx3event.Energy1, path_length_s);

        // plotter->Fill2D("sx3Ef_sx3E_matrix_all"+tag,400,0,10,400,0,10,sx3Efix,sx3event.Energy1,pmlabel);
        plotter->Fill2D("pmiscs_dE3_Ef_Anodesx3" + tag, 400, 0, 10, 400, 0, 40000, sx3Efix, pcevent.Energy1 * sinTheta_customV * 3, pmlabel);
        plotter->Fill2D("pmiscs_dE3_Ef_Cathodesx3" + tag, 400, 0, 10, 400, 0, 10000, sx3Efix, pcevent.Energy2 * sinTheta_customV, pmlabel);

        plotter->Fill2D("pmiscs_Ef_vs_theta_sx3" + tag, 100, 0, 180, 800, 0, 20, theta_s * 180 / M_PI, sx3Efix, pmlabel);
        plotter->Fill1D("pmiscs_VertexReconZ" + tag, 800, -400, 400, vertex_z, pmlabel);
        plotter->Fill2D("pmiscs_VertexReconXY" + tag, 200, -100, 100, 200, -100, 100, r_rhoMin_fix.X(), r_rhoMin_fix.Y(), pmlabel);
        plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef" + tag, 800, -400, 400, 800, 0, 20, vertex_z, sx3Efix, pmlabel);
        plotter->Fill2D("pmiscs_VertexReconZ_vs_Ef" + tag + "_a" + std::to_string(pcevent.multi1), 800, -400, 400, 800, 0, 20, vertex_z, sx3Efix, pmlabel);
        if (tag == "_cathode_alphas")
          plotter->Fill1D("pmiscs_Ex_from_alpha", 200, -10, 10, apkin_a_s.getExc(sx3Efix, theta_s * 180 / M_PI), pmlabel);
      };

      plot_with_tag();
      if (cathode_alpha_select)
        plot_with_tag("_cathode_alphas");
      else
        plot_with_tag("_cathode_protons");

      // plotter->Fill1D("pmisc_Ex_from_protons",200,-10,10,apkin_p.getExc(sx3Efix,theta_s*180/M_PI),pmlabel);

    } // end PCEvents loop (A1C2 main flow)
    for (const auto &pcevent : PC_Events)
    {
      if (!(pcevent.multi1 == 1 && pcevent.multi2 == 1))
        continue;
      bool phicut = sx3event.pos.Phi() <= pcevent.pos.Phi() + TMath::Pi() / 3. && sx3event.pos.Phi() >= pcevent.pos.Phi() - TMath::Pi() / 3.;
      if (!phicut)
        continue;
      if (!(pcevent.Energy2 > 1400))
        continue;

      const std::string wcat = a1c1_missing_neighbor(pcevent.Anodech, pcevent.Cathodech) ? "_missingw" : "_true1w";
      auto fillCmp = [&](double pcz, const std::string &m)
      {
        TVector3 x2(pcevent.pos.X(), pcevent.pos.Y(), pcz);
        TVector3 rv = beamVertex(sx3event.pos, x2 - sx3event.pos);
        if (beamPerp(rv) > 10.0)
          return;
        double th = (sx3event.pos - rv).Theta();
        double pl = pathLengthCm(sx3event.pos, rv);
        double Ef = evalEloss(MeV_to_cm_spl, cm_to_MeV_spl, sx3event.Energy1, pl);
        double beam_pl_cmp = TMath::Abs(rv.Z() - z_entrance) * 0.1;
        double beam_E_cmp = evalElossForward(MeV_to_cm_p_spl, cm_to_MeVp_spl, initial_energy, beam_pl_cmp);
        beam_E_cmp = applyTaFoilEloss(beam_E_cmp, rv.Z());
        if (beam_E_cmp <= 0.0)
          beam_E_cmp = 0.001;
        Kinematics apkin_a_cmp(mass_1H, mass_4He, mass_4He, mass_1H, beam_E_cmp / mass_1H);
        double Ex = apkin_a_cmp.getExc(Ef, th * 180 / M_PI);
        std::string lbl = "proton+miscsx3_a1c1cmp";
        for (const std::string &w : {std::string(""), wcat})
        {
          plotter->Fill1D("pmiscs_a1c1cmp_pcz_" + m + w, 600, -300, 300, pcz, lbl);
          plotter->Fill1D("pmiscs_a1c1cmp_Ex_" + m + w, 200, -10, 10, Ex, lbl);
          plotter->Fill1D("pmiscs_a1c1cmp_VertexZ_" + m + w, 800, -400, 400, rv.Z(), lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_VertexZ_vs_Ef_" + m + w, 800, -400, 400, 800, 0, 20, rv.Z(), Ef, lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_VertexZ_vs_Ex_" + m + w, 800, -400, 400, 800, -10, 10, rv.Z(), Ex, lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_phi_vs_Ef_" + m + w, 90, -180, 180, 800, 0, 20, sx3event.pos.Phi() * 180 / M_PI, Ef, lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_phi_vs_Ex_" + m + w, 90, -180, 180, 800, -10, 10, sx3event.pos.Phi() * 180 / M_PI, Ex, lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_Ef_vs_theta_" + m + w, 100, 0, 180, 800, 0, 20, th * 180 / M_PI, Ef, lbl);
          plotter->Fill2D("pmiscs_a1c1cmp_Ex_vs_theta_" + m + w, 100, 0, 180, 800, -10, 10, th * 180 / M_PI, Ex, lbl);
        }
      };

      double pcz_dith_s = rand.Gaus(pcevent.pos.Z(), 8.0);
      fillCmp(pcz_dith_s, "dither");
      double ac = pcevent.Energy1 + pcevent.Energy2;
      double cfrac = (ac > 0.0) ? pcevent.Energy2 / ac : -1.0;
      if (cfrac >= 0.0)
      {
        std::vector<std::tuple<int, double, double>> aOne = {std::make_tuple(pcevent.Anodech, 1.0, 0.0)};
        auto apw = pwinstance.GetPseudoWire(aOne, "ANODE");
        A1C1Sol s = a1c1_solve(cfrac, pcevent.pos.Z(), pcevent.Cathodech, pcevent.Energy1, pcevent.Anodech);
        int side_status;
        SideChoice side = a1c1_pick_side(sx3event.pos, pcevent.pos.X(), pcevent.pos.Y(), s.pcz_lo, s.pcz_hi, side_status);
        const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
        double pcz_pick = best.pcz;
        fillCmp(pcz_pick, "cfrac_all");
        if (best.inband && side_status != 2)
        {
          fillCmp(pcz_pick, "cfrac");
          plotter->Fill2D("pmiscs_a1c1cmp_pcz_cfrac_vs_dither", 600, -300, 300, 600, -300, 300, pcz_dith_s, pcz_pick, "proton+miscsx3_a1c1cmp");
        }
      }
    } // end A1C1 comparison loop
  } // end sx3Events loop
}

inline double a1c1_cfrac_pcz(const Event &pcevent, const TVector3 &si, bool &inband)
{
  inband = false;
  double ac = pcevent.Energy1 + pcevent.Energy2;
  double cfrac = (ac > 0.0) ? pcevent.Energy2 / ac : -1.0;
  if (cfrac < 0.0)
    return pcevent.pos.Z();
  A1C1Sol s = a1c1_solve(cfrac, pcevent.pos.Z(), pcevent.Cathodech, pcevent.Energy1, pcevent.Anodech);
  int side_status;
  SideChoice side = a1c1_pick_side(si, pcevent.pos.X(), pcevent.pos.Y(), s.pcz_lo, s.pcz_hi, side_status);
  const A1C1CellSol &best = (side == SideChoice::High) ? s.hi : s.lo;
  inband = (best.inband && side_status != 2);
  return best.pcz;
}

static const std::vector<double> levels_30Si_MeV = {
    0.0, 2.235, 3.498, 6.550, 6.870, 7.043, 7.220, 8.660, 9.110,
    9.370, 9.840, 10.430, 10.720, 10.910, 11.220, 11.460, 12.340,
    14.240, 14.690, 15.950, 16.700};
// 27Al levels (27Al(a,a')27Al* inelastic recoil), from Adopted Levels.
static const std::vector<double> levels_27Al_MeV = {
    0.0, 6.1584, 6.4773, 6.6513, 7.2272, 7.4771, 7.935, 7.948};

inline double snapToNearestLevel(double ex, const std::vector<double> &levels, double &residual)
{
  double best = levels.front();
  residual = std::abs(ex - best);
  for (double lvl : levels)
  {
    double r = std::abs(ex - lvl);
    if (r < residual)
    {
      residual = r;
      best = lvl;
    }
  }
  return best;
}

static void reaction_ax_core(HistPlotter *plotter, const std::vector<Event> &Si_Events, const std::vector<Event> &PC_Events,
                             const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, bool isQQQ,
                             const std::string &rx, const std::string &det, double si_ecut, double perp_cut, double phi_win,
                             double dEa_max, double dEc_max, double ef_max,
                             double beamE0, TSpline3 *beam_MeV_to_cm, TSpline3 *beam_cm_to_MeV, double m_beam,
                             const AAEjectileMasses &ej_m, const std::string &globaltag)
{
  const std::string sfx = "_" + det + globaltag;

  for (const auto &sievent : Si_Events)
  {
    if (sievent.Energy1 < si_ecut)
      continue; // coarse Si energy cut

    auto reconstructAndFill = [&](double pcz_fix, const TVector3 &pcXY, double anodeE, double cathodeE, double anodeE_MeV, double cathodeE_MeV,
                                  const std::string &topo1, const std::string &topo2 = "")
    {
      TVector3 x2f(pcXY.X(), pcXY.Y(), pcz_fix);
      TVector3 r_rhoMin_fix = beamVertex(sievent.pos, x2f - sievent.pos);
      double vertex_z = r_rhoMin_fix.Z();
      if (beamPerp(r_rhoMin_fix) > perp_cut || vertex_z < z_entrance)
        return;

      double theta = (sievent.pos - r_rhoMin_fix).Theta();

      double beam_path_length = TMath::Abs(vertex_z - z_entrance) * 0.1; // mm -> cm
      double beam_energy_at_vertex = evalElossForward(beam_MeV_to_cm, beam_cm_to_MeV, beamE0, beam_path_length);
      if (beam_energy_at_vertex <= 0.0)
        return;

      plotter->Fill2D(rx + "_BeamEnergy_vs_VertexZ" + sfx, 800, -400, 400, 400, 0, beamE0, vertex_z, beam_energy_at_vertex, globaltag + "_" + rx + "+misc_" + det);

      bool trueProton = (beam_energy_at_vertex < 10.0);
      bool isAlpha = (!trueProton && anodeE >= 2200);

      double m3, m4;
      TSpline3 *ej_fwd, *ej_inv;
      std::string ejtag;

      if (isAlpha)
      {
        m3 = ej_m.m_a;
        m4 = ej_m.m_ra;
        ej_fwd = MeV_to_cm_spl;
        ej_inv = cm_to_MeV_spl;
        ejtag = "_a";
      }
      else
      {
        m3 = ej_m.m_p;
        m4 = ej_m.m_rp;
        ej_fwd = MeV_to_cm_p_spl;
        ej_inv = cm_to_MeVp_spl;
        ejtag = "_p";
      }

      Kinematics kin(m_beam, mass_4He, m3, m4, beam_energy_at_vertex / m_beam); // beamE given as E/u
      double path_length = pathLengthCm(sievent.pos, r_rhoMin_fix);
      double Efix = evalEloss(ej_fwd, ej_inv, sievent.Energy1, path_length);
      double Ex = kin.getExc(Efix, theta * 180 / M_PI);
      std::string pmlabel = globaltag + "_" + rx + "+misc_" + det + ejtag;

      const double ex_gate_MeV = 3;
      const std::vector<double> &levels = (ejtag == "_a") ? levels_27Al_MeV : levels_30Si_MeV;
      double level_residual = 0.0;
      // double snapped_level = snapToNearestLevel(Ex, levels, level_residual);
      double ebeam_kin_MeV = (Ex < ex_gate_MeV)
                                 ? invertBeamEnergyMeV(m_beam, mass_4He, m3, m4, Efix, theta * 180 / M_PI, 0.0)
                                 : -1.0;
      // double ebeam_kin_MeV = invertBeamEnergyMeV(m_beam, mass_4He, m3, m4, Efix, theta * 180 / M_PI, snapped_level);

      auto plot_with_tag = [&](const std::string &topo)
      {
        std::string t = topo.empty() ? "" : ("_" + topo);
        plotter->Fill1D(rx + "_Ex_from" + ejtag + t + sfx, 400, -20, 20, Ex, pmlabel);
        plotter->Fill2D(rx + "_VertexReconZ_vs_Ef" + ejtag + t + sfx, 800, -400, 400, 800, 0, ef_max, vertex_z, Efix, pmlabel);
        plotter->Fill2D(rx + "_VertexReconZ_vs_Ex" + ejtag + t + sfx, 800, -400, 400, 800, -20, 20, vertex_z, Ex, pmlabel);
        if (ebeam_kin_MeV > 0.0)
          plotter->Fill2D(rx + "_BeamEnergy_Eloss_vs_Kin" + ejtag + t + sfx, 400, 0, beamE0 * 1.5, 400, 0, beamE0 * 1.5,
                          beam_energy_at_vertex, ebeam_kin_MeV, pmlabel);
      };

      plotter->Fill2D(rx + "_dE_E_Anode" + sfx, 400, 0, dEa_max, 800, 0, 40000, sievent.Energy1, anodeE, pmlabel);
      plotter->Fill2D(rx + "_dE_E_Anode" + sfx + "_E<10MeV" + std::to_string(beam_energy_at_vertex < 10), 400, 0, dEa_max, 800, 0, 40000, sievent.Energy1, anodeE, pmlabel);
      if (cathodeE >= 0.0)
        plotter->Fill2D(rx + "_dE_E_Cathode" + sfx, 400, 0, dEa_max, 800, 0, dEc_max, sievent.Energy1, cathodeE, pmlabel);
      plotter->Fill1D(rx + "_pczfix" + sfx, 600, -300, 300, pcz_fix, pmlabel);
      plotter->Fill2D(rx + "_Ef_vs_theta" + ejtag + sfx, 100, 0, 180, 800, 0, ef_max, theta * 180 / M_PI, Efix, pmlabel);
      plotter->Fill2D(rx + "_Ex_vs_theta" + ejtag + sfx, 100, 0, 180, 800, -20, 20, theta * 180 / M_PI, Ex, pmlabel);
      plotter->Fill1D(rx + "_VertexReconZ" + sfx, 800, -400, 400, vertex_z, pmlabel);

      plot_with_tag(""); // unconditional master fill, all topologies
      plot_with_tag(topo1);
      if (!topo2.empty())
        plot_with_tag(topo2);
      if (trueProton)
        plot_with_tag("trueProton"); // clean, alpha-free proton sub-sample

      // Gas segmentation validation
      PCPath pp = pcPath(r_rhoMin_fix, sievent.pos);
      if (pp.ok)
      {
        double E_an = evalEloss(ej_fwd, ej_inv, sievent.Energy1, pp.anode_cm);
        double E_ca = evalEloss(ej_fwd, ej_inv, sievent.Energy1, pp.cathode_cm);
        plotter->Fill2D(rx + "_dEgas_vs_Ef" + ejtag + sfx, 400, 0, ef_max, 400, 0, 1, Efix, E_an - E_ca, pmlabel);
        if (anodeE_MeV >= 0.0 && cathodeE_MeV >= 0.0)
        {
          plotter->Fill2D(rx + "_dEgasCalib_vs_Ef" + ejtag + sfx, 400, 0, ef_max, 800, -2, 2, Efix, anodeE_MeV , pmlabel);
          plotter->Fill2D(rx + "_dEgasCalib_vs_E" + ejtag + sfx, 400, 0, ef_max, 800, -2, 2, sievent.Energy1, anodeE_MeV , pmlabel);
          plotter->Fill2D(rx + "_dEgasPred_vs_dEgasCalib" + ejtag + sfx, 800, -2, 2, 400, 0, 2, anodeE_MeV , E_an - E_ca, pmlabel);
        }
      }
    };

    for (const auto &pcevent : PC_Events)
    {
      if (!(pcevent.multi1 == 1 && (pcevent.multi2 == 1 || pcevent.multi2 == 2)))
        continue;
      if (TMath::Abs(sievent.pos.DeltaPhi(pcevent.pos)) > phi_win)
        continue;
      double anodeE_MeV = (pcevent.Anodech >= 0 && pcevent.Anodech < 24)
                              ? pcEnergySlope[pcevent.Anodech] * pcevent.Energy1 + pcEnergyIntercept[pcevent.Anodech]
                              : -1.0;
      double cathodeE_MeV = (pcevent.Cathodech >= 0 && pcevent.Cathodech < 24)
                                ? pcEnergySlope[24 + pcevent.Cathodech] * pcevent.Energy2 + pcEnergyIntercept[24 + pcevent.Cathodech]
                                : -1.0;

      if (pcevent.multi2 == 1) // A1C1
      {
        bool a1c1_inband = false;
        double pcz_fix = a1c1_cfrac_pcz(pcevent, sievent.pos, a1c1_inband);

        double ac = pcevent.Energy1 + pcevent.Energy2;
        double cfrac = (ac > 0.0) ? pcevent.Energy2 / ac : -1.0;
        if (cfrac >= 0.0)
        {
          std::string pmlabel = globaltag + "_" + rx + "+misc_" + det + "_a1c1cfrac";
          plotter->Fill1D(rx + "_a1c1_cfrac" + sfx, 220, -0.05, 1.05, cfrac, pmlabel);
          plotter->Fill2D(rx + "_a1c1_cfrac_vs_anodeE" + sfx, 400, 0, 40000, 220, -0.05, 1.05, pcevent.Energy1, cfrac, pmlabel);
          plotter->Fill1D(rx + "_a1c1_cfrac_inband" + sfx, 220, -0.05, 1.05, a1c1_inband ? cfrac : -1.0, pmlabel);
        }

        reconstructAndFill(pcz_fix, pcevent.pos, pcevent.Energy1, pcevent.Energy2, anodeE_MeV, cathodeE_MeV,
                           "a1c1", a1c1_inband ? "a1c1_inband" : "");
      }
      else // A1C2 (multi2 == 2)
      {
        double pcz_fix = pcfix_func.Eval(pcevent.pos.Z());
        reconstructAndFill(pcz_fix, pcevent.pos, pcevent.Energy1, pcevent.Energy2, anodeE_MeV, cathodeE_MeV, "a1c2fix");
      }
    }

    for (const auto &aCl : aClusters)
    {
      if (aCl.empty())
        continue;
      auto aPw = pwinstance.GetPseudoWire(aCl, "ANODE");
      auto apwire = std::get<0>(aPw);
      double apSumE = std::get<1>(aPw);

      TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire, sievent.pos.Phi());
      pc.SetZ(a1c1_zcorr(pc.Z(), isQQQ));

      if (TMath::Abs(sievent.pos.DeltaPhi(pc)) > phi_win)
        continue;

      std::string pmlabel = globaltag + "_" + rx + "+misc_" + det + "_a1c0";
      plotter->Fill2D(rx + "_dE_E_Anode_a1c0" + sfx, 400, 0, dEa_max, 800, 0, 40000, sievent.Energy1, apSumE, pmlabel);
      TVector3 r_rhoMin_a1c0 = beamVertex(sievent.pos, pc - sievent.pos);
      double beam_path_length_a1c0 = TMath::Abs(r_rhoMin_a1c0.Z() - z_entrance) * 0.1;
      double beam_energy_at_vertex_a1c0 = evalElossForward(beam_MeV_to_cm, beam_cm_to_MeV, beamE0, beam_path_length_a1c0);
      plotter->Fill2D(rx + "_dE_E_Anode_a1c0" + sfx + "_10MeV" + std::to_string(beam_energy_at_vertex_a1c0 < 10), 400, 0, dEa_max, 800, 0, 40000, sievent.Energy1, apSumE, pmlabel);
      plotter->Fill2D(rx + "_dPhi_a1c0" + sfx, 100, -200, 200, 100, -200, 200, pc.Phi() * 180 / M_PI, sievent.pos.Phi() * 180 / M_PI, pmlabel);
      plotter->Fill1D(rx + "_rawZ_a1c0" + sfx, 600, -300, 300, pc.Z(), pmlabel);

      reconstructAndFill(pc.Z(), pc, apSumE, -1.0, -1.0, -1.0, "a1c0");
    }
  }
}

void miscHistograms_17Fax(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                          const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, std::string globaltag)
{

  // 17F beam energy at the gas target, after the entrance-window foils:
  // 67.8 MeV -> Mylar (MCP, 4.2426 um-equiv) -> 64.0305 MeV -> Kapton (7.8 um) -> 56.7173 MeV.
  double ebeam_17F_MeV = 56.7173;
  // 17F(a,a)/(a,d)/(a,p): ejectile + recoil masses per channel.
  AAEjectileMasses ej17F{mass_4He, mass_17F, mass_2H, mass_19Ne_rec, mass_1H, mass_20Ne};
  reaction_ax_core(plotter, QQQ_Events, PC_Events, aClusters, true, "m17Fax", "qqq", 0.6, 6.0, TMath::Pi() / 4.0,
                   30.0, 40000.0, 30.0, ebeam_17F_MeV, MeV_to_cm_17F_spl, cm_to_MeV_17F_spl, mass_17F, ej17F, globaltag);
  reaction_ax_core(plotter, SX3_Events, PC_Events, aClusters, false, "m17Fax", "sx3", 1.2, 10.0, TMath::Pi() / 3.0,
                   30.0, 40000.0, 30.0, ebeam_17F_MeV, MeV_to_cm_17F_spl, cm_to_MeV_17F_spl, mass_17F, ej17F, globaltag);
}

// 27Al(a,a) excitation functions for BOTH silicon branches (QQQ + SX3), with the
// 27Al beam table. Same consistently-named histogram set as 17F.
void miscHistograms_27Alax(HistPlotter *plotter, const std::vector<Event> &QQQ_Events, const std::vector<Event> &SX3_Events, const std::vector<Event> &PC_Events,
                           const std::vector<std::vector<std::tuple<int, double, double>>> &aClusters, std::string globaltag)
{
  // 27Al(a,a)/(a,d)/(a,p): ejectile + recoil masses per channel.
  AAEjectileMasses ej27Al{mass_4He, mass_27Al, mass_2H, mass_29Si_rec, mass_1H, mass_30Si};
  reaction_ax_core(plotter, QQQ_Events, PC_Events, aClusters, true, "m27Alax", "qqq", 0.6, 6.0, TMath::Pi() / 4.0,
                   10.0, 10000.0, 20.0, 56.103, MeV_to_cm_27Al_spl, cm_to_MeV_27Al_spl, mass_27Al, ej27Al, globaltag);
  reaction_ax_core(plotter, SX3_Events, PC_Events, aClusters, false, "m27Alax", "sx3", 1.2, 10.0, TMath::Pi() / 3.0,
                   10.0, 10000.0, 20.0, 56.103, MeV_to_cm_27Al_spl, cm_to_MeV_27Al_spl, mass_27Al, ej27Al, globaltag);
}
