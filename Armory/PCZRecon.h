#ifndef PCZRecon_h
#define PCZRecon_h

// PC Z-position reconstruction, one section per anode/cathode topology:
// A1C0 (single anode wire only), A2C0 (two-wire anode cluster, no cathode --
// same math as A1C0, see that section), A1C1 (anode + single cathode, charge
// division), A1C2 (anode + two cathodes, "step ladder" correction). Each topology gets one
// well-defined entry point instead of the math being split across files by
// historical accident (A1C0/A1C1 used to live in TrackRecon.C itself; A1C2's
// underlying model lives in the separately-shared PC_StepLadder_Correction.h
// -- see the A1C2 section below for why that one isn't just moved in).
//
// This header holds the reconstruction MATH only. The per-dataset tuning
// constants it reads (cfrac fit parameters, dead-wire lists, Z calibration,
// beam-axis origin) are still owned and set by TrackRecon.C's Begin() --
// this header just declares them `extern` so the same single translation
// unit (TrackRecon.C is compiled as one .C file via ACLiC) can see them.
// Moving the constants themselves out is a separate, riskier change and is
// deliberately NOT done here.
//
// Relocated verbatim from TrackRecon.C (no logic changes): A1C1CellSol,
// A1C1Sol, solve_cell, a1c1_solve, SideChoice, a1c1_pick_side, a1c1_zcorr,
// a1c0_hybrid_pcz (split into a1c0_wirePos + a1c0_hybrid_pcz).
//
// New in this header (see call-site migration notes where each is used):
// a1c1_solve_pick, a1c1_cfrac_pcz, a1c2_zfix.

#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>
#include "ClassPW.h"

// --- Per-dataset tuning constants, defined and set in TrackRecon.C ---
extern PW pwinstance;
extern const double a1c1_zg[8];
extern double a1c1_cfmin_cell[7];
extern double a1c1_k_cell[7];
extern double a1c1_cfmin2_cell[7];
extern double a1c1_k2_cell[7];
extern double a1c1_cfrac_split;
extern double a1c1_lowband_rfactor;
extern double a1c1_missing_fmax;
extern double a1c1_z_scale_qqq;
extern double a1c1_z_off_qqq;
extern double a1c1_z_scale_sx3;
extern double a1c1_z_off_sx3;
bool a1c1_missing_neighbor(int awire, int cwire); // TrackRecon.C: dead-wire-adjacency check
TVector3 beamVertex(const TVector3 &si, const TVector3 &dir);
double beamPerp(const TVector3 &p);

// ---------------------------------------------------------------------
// A1C0: single-wire (anode only) position reconstruction
// ---------------------------------------------------------------------

// Anode-wire Z, corrected by the a1c1-derived scale+offset (a1c1_zcorr),
// same reference frame the a1c1 solve below reports in.
inline double a1c1_zcorr(double z_a1c0, bool isQQQ)
{
  double scale = isQQQ ? a1c1_z_scale_qqq : a1c1_z_scale_sx3;
  double off = isQQQ ? a1c1_z_off_qqq : a1c1_z_off_sx3;
  return z_a1c0 * (1.0 - scale) - off;
}

// The raw (undithered) A1C0 wire position: nearest-wire XY at the given
// track phi, Z corrected into the same frame a1c1 uses. This is what
// _rawZ_a1c0 / phi cuts / _dPhi_a1c0-style diagnostics should read --
// anywhere you want the true discretized wire position, not a smoothed one.
inline TVector3 a1c0_wirePos(const std::pair<TVector3, TVector3> &apwire, double phi, bool isQQQ)
{
  TVector3 pc = pwinstance.getClosestWirePosAtWirePhi(apwire, phi);
  pc.SetZ(a1c1_zcorr(pc.Z(), isQQQ));
  return pc;
}

// a1c0_wirePos, then Gaussian-dithered in Z to hide wire-pitch quantization.
// sigma is the caller's choice (dither_sigma or dither_sigma_c0/2.0 etc in
// TrackRecon.C) -- this function doesn't know which convention is "correct"
// for a given call site, only how to apply whichever sigma it's given.
inline TVector3 a1c0_hybrid_pcz(const std::pair<TVector3, TVector3> &apwire, double phi,
                                bool isQQQ, double sigma, TRandom3 &rand)
{
  TVector3 pc = a1c0_wirePos(apwire, phi, isQQQ);
  pc.SetZ(rand.Gaus(pc.Z(), sigma));
  return pc;
}

// ---------------------------------------------------------------------
// A2C0: two-wire anode-cluster (charge-shared), no-cathode position
// reconstruction
// ---------------------------------------------------------------------
//
// Same math as A1C0 above -- pseudowire + phi-minimization, no dither --
// just handed a genuine 2-wire energy-weighted pseudowire (GetPseudoWire
// over a 2-wire anode cluster) instead of a single real wire. a1c0_wirePos
// already treats its `apwire` argument as an opaque pseudowire pair
// regardless of how many physical wires went into it, so this is a
// documented, named entry point rather than new math: call sites can say
// what topology they mean instead of reusing a1c0_wirePos silently and
// trusting a comment to explain why.
//
// Deliberately no dithered twin (no a2c0_hybrid_pcz): A2C0 is meant to feed
// the reaction-analysis plots at its raw, undithered resolution, not stand
// in for a1c0_hybrid_pcz's benchmark-truth-comparison role. If a "genuine
// A2C0" BenchMark validation block is ever wanted (mirroring the existing
// aClusters.size()==1 && cClusters.size()==0 A1C0 block in TrackRecon.C),
// add one there rather than adding dithering here.
inline TVector3 a2c0_wirePos(const std::pair<TVector3, TVector3> &apwire, double phi, bool isQQQ)
{
  return a1c0_wirePos(apwire, phi, isQQQ);
}

// ---------------------------------------------------------------------
// A1C1: single-anode + single-cathode charge-division position reconstruction
// ---------------------------------------------------------------------

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

// Which of the two candidate cells the beam-axis test selects.
enum class SideChoice
{
  High, // the cell ABOVE the fired wire (pcz_hi)
  Low   // the cell BELOW the fired wire (pcz_lo)
};

// Picks between the two candidate cells (lo/hi) a1c1_solve returns.
//
// This used to compare beamVertex(si, ...)'s beamPerp for each candidate, but
// that check is a no-op: both candidates share the same (cx, cy), so
// beamVertex's t-parameter -- computed purely from dir.X()/dir.Y() -- comes out
// identical for both, and therefore so does the resulting vertex X/Y and
// beamPerp. The two disambiguating branches below were consequently dead code;
// the function always fell through to returning Low regardless of the actual
// physics.
//
// inband/pitchok DO differ between the two cells (each depends on cfrac's
// position within that cell's own band), so they're used here for status
// instead. The final tie-break prefers whichever candidate sits closer to the
// fired wire itself, rather than farther into the next cell over.
inline SideChoice a1c1_pick_side(const A1C1CellSol &lo, const A1C1CellSol &hi, double zf, int &status)
{
  bool okl = lo.inband && lo.pitchok;
  bool okh = hi.inband && hi.pitchok;
  status = (okl || okh) ? ((okl && okh) ? 1 : 0) : 2;
  if (okl && !okh)
    return SideChoice::Low;
  if (okh && !okl)
    return SideChoice::High;
  double dl = TMath::Abs(lo.pcz - zf);
  double dh = TMath::Abs(hi.pcz - zf);
  return (dl <= dh) ? SideChoice::Low : SideChoice::High;
}

// a1c1_solve() + a1c1_pick_side() together, with the picked cell already
// resolved. Every call site that needs more than just the final pcz (side
// status, cell index, f, inband, pitchok -- e.g. for benchmark/diagnostic
// histograms) was previously re-deriving this same 4-line pattern by hand;
// this is that pattern, named once.
struct A1C1PickedSol
{
  A1C1Sol sol;
  SideChoice side = SideChoice::Low;
  int side_status = -1;

  const A1C1CellSol &best() const { return (side == SideChoice::High) ? sol.hi : sol.lo; }
};

inline A1C1PickedSol a1c1_solve_pick(double cfrac, double zf, const TVector3 &si, double cx, double cy,
                                     int cwire = -1, double anodeE = -1, int awire = -1)
{
  // si/cx/cy are no longer used by the pick itself (see a1c1_pick_side above)
  // but are kept in the signature so the ~10 existing call sites throughout
  // TrackRecon.C don't need to change.
  A1C1PickedSol out;
  out.sol = a1c1_solve(cfrac, zf, cwire, anodeE, awire);
  out.side = a1c1_pick_side(out.sol.lo, out.sol.hi, zf, out.side_status);
  return out;
}

// Full pipeline: raw anode/cathode energies -> cfrac -> solve -> pick side ->
// picked Z, plus whether the pick landed in-band. This is what most call
// sites actually want (they don't need the intermediate A1C1Sol/side_status
// unless they're doing benchmark diagnostics -- for that, call
// a1c1_solve_pick directly instead).
//
// Takes primitives rather than TrackRecon.C's `Event` type so this header
// has no dependency on TrackRecon.C's class definitions; see the thin
// `Event`-taking overload kept in TrackRecon.C next to the `Event` class
// for the short call-site spelling existing code uses.
inline double a1c1_cfrac_pcz(double pcz_raw, double energyAnode, double energyCathode,
                             double cx, double cy, int cathodeCh, int anodeCh,
                             const TVector3 &si, bool &inband)
{
  inband = false;
  double ac = energyAnode + energyCathode;
  double cfrac = (ac > 0.0) ? energyCathode / ac : -1.0;
  if (cfrac < 0.0)
    return pcz_raw;
  A1C1PickedSol picked = a1c1_solve_pick(cfrac, pcz_raw, si, cx, cy, cathodeCh, energyAnode, anodeCh);
  const A1C1CellSol &best = picked.best();
  inband = (best.inband && picked.side_status != 2);
  return best.pcz;
}

// ---------------------------------------------------------------------
// A1C2: two-cathode "step ladder" Z reconstruction
// ---------------------------------------------------------------------
//
// Same step-ladder, pivot-about-cell-midpoint model as Armory/
// PC_StepLadder_Correction.h's model_invert. That file is shared with
// MakeVertex.C (a separate, parallel branch of development, out of scope
// here) and a few scratch macros, so rather than include it or edit it,
// this is a standalone copy -- keeping every TrackRecon.C reconstruction
// path (A1C0/A1C1/A1C2) self-contained in this one header rather than
// reaching into a file with unrelated consumers. If the underlying model
// ever changes, both copies need the same edit; there are exactly two.
//
// Rewritten as a plain scalar function rather than kept in model_invert's
// original `double f(double *y, double *p)` shape -- that signature only
// existed to match TF1's raw-function-pointer constructor. TrackRecon.C
// used to go through a `TF1 pcfix_func` purely to get a callable out of it,
// calling only pcfix_func.Eval(z) at every site (npar `p` was never used,
// and nothing called .Draw()/.Integral()/anything else TF1-specific) --
// so the TF1 wrapper bought nothing. a1c2_zfix has the same scalar
// in/scalar out shape as a1c0_wirePos/a1c1_solve above instead.
inline double a1c2_zfix(double z)
{
  double result = z;
  double slope = 0.52;
  double z_grid[8] = {147.998, 101.946, 59.7634, 19.6965, -19.6965, -59.7634, -101.946, -147.998};
  for (int i = 0; i < 7; i++)
  {
    if (z <= z_grid[i] && z > z_grid[i + 1])
    {
      double zavg = (z_grid[i] + z_grid[i + 1]) * 0.5; // midpoint about which we pivot
      result = (z - zavg) / slope + zavg;
      break;
    }
  }
  return result;
}

#endif