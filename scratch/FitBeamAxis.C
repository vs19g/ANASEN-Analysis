// FitBeamAxis.C
// ---------------------------------------------------------------------------
// Single-pass beam-axis line fit. Replaces the iterative
// FitBeamAxis.C + iterate_beam.sh scheme, which diverges.
//
// WHY THE ITERATIVE SCHEME FAILS
// ------------------------------
// FitBeamAxis.C fits beamAxis_all_vertexX_vs_Z and feeds the residual back as a
// correction to BEAM_AXIS_X / BEAM_TILT_X. But that vertex comes from
// beamVertex(), which returns the track's point of closest approach to the
// ASSUMED beam line -- so the reconstructed vertex is pulled onto whatever axis
// the run was configured with, by construction. Its residual therefore measures
// acceptance asymmetry and geometry artifacts far more than it measures "how far
// off is my assumed axis". Feeding that back adds a roughly fixed bias each pass
// with no restoring force, so the parameters drift monotonically and never
// converge. Measured behaviour: step L2 grew 1.48 -> 2.65 over 6 passes at full
// step, and still grew (1.27 -> 1.50 over 10 passes) at relax=0.3 -- damping only
// scales the drift rate, which is the signature of a wrong-direction update
// rather than overshoot.
//
// WHAT THIS DOES INSTEAD
// ----------------------
// TrackRecon.C's fillBeamProfile() also histograms where each track CROSSES a set
// of fixed z-planes (beamAxis_<tag>_crossX_z**, crossY_z**). A track crossing is a
// property of the track alone -- it never references the assumed beam axis. Tracks
// originate on the true beam line and fan out with roughly uniform azimuth, so at
// each plane the crossing distribution centres on the true beam position at that z.
// Take a robust centre per plane, line-fit those centres against z, and the axis
// falls out directly. One pass, no iteration, no feedback.
//
// WHICH TRACKS FEED THE PLANES
// ----------------------------
// fillBeamProfile() only sends a track to the crossing planes when its PC z is a
// real measurement: a1c2 (cathode charge division), a2c0 (two wires, geometry
// fixes z), and a1c1 whose a1c1_cfrac_pcz() solve landed in band. Bare a1c1 and
// a1c0 are excluded -- their z is a Gaussian dither about the raw wire z, which
// puts them on the assumed axis by construction and would make this fit return
// whatever axis the run was configured with. The vertexX_vs_Z and pocaDist
// diagnostics in the same folder are still filled for every topology, so their
// entry counts will not match the crossing planes'. Do not fit those.
//
// CAVEAT -- azimuthal symmetry. The per-plane centre is only unbiased if
// acceptance is uniform in phi. Dead channels break that (this analysis already
// excludes QQQ wedge 48/49, ring 63, and several anode wires). The median is used
// rather than the mean to blunt the resulting asymmetric tails, but a strongly
// lopsided acceptance will still pull the centre. If the fitted axis looks
// implausible, check the crossX/crossY distributions at a few planes for
// visible asymmetry before trusting the number.
//
// USAGE:
//   root -l -b -q 'scratch/FitBeamAxisCrossings.C("Output_27Al/Output_27Al.root")'
//   root -l -b -q 'scratch/FitBeamAxisCrossings.C("out.root", -440, 40, "all")'
//   root -l -b -q 'scratch/FitBeamAxisCrossings.C("out.root", -440, 40, "elastic")'
//
// Arguments:
//   filename  TrackRecon.C output ROOT file
//   zmin/zmax only planes whose z falls in this window are used in the line fit.
//             Defaults span every plane; narrow it to exclude planes you have
//             reason to distrust, not to chase a better chi2.
//   tag       which branch's crossings to use. "all" pools everything; the
//             per-branch tags are "reaction_<channel>_<det>" (e.g.
//             "reaction_m27Alax_qqq", "..._sx3") and "elastic_<DET>" ("elastic_QQQ",
//             "elastic_SX3"). Fitting the qqq and sx3 tags separately is the check
//             that a tilt is the beam and not one branch's acceptance: same beam,
//             different geometry and phi coverage, so they must agree.
//             Which tags exist at all depends on run configuration -- on a dataset
//             where doPCSX3ClusterAnalysis / doPCQQQClusterAnalysis /
//             process_alpha_proton_scattering are false and OUT_DIR is not
//             Output_p, the sx3a1c2, qqqa1c2, apCoinc and elastic_* tags are never
//             filled and "all" is exactly the reaction_* tags pooled.
//
// NOTE on chi2/ndf: medianAndError reports a STATISTICS-ONLY uncertainty
// (~8 um at a 4e5-entry plane) while observed plane-to-plane scatter in a
// well-behaved region is ~0.4 mm, i.e. ~50x larger. chi2/ndf is therefore
// inflated by a missing systematic term and will not approach 1 even for a
// correct axis. Judge the fit on whether the medians lie on a straight line,
// not on chi2/ndf alone.
//
// Prints the fitted BEAM_AXIS_X/Y and BEAM_TILT_X/Y as ready-to-paste export
// lines. These are ABSOLUTE values, not increments -- set them directly in
// run_27Al.sh; do NOT add them to the current values.
// ---------------------------------------------------------------------------

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Must match fillBeamProfile() in TrackRecon.C. The range follows the physical
// vertex acceptance (z_entrance = -454 mm up to the +100 mm cut), not a
// symmetric window about 0, and stops well below the QQQ plane at z = 105 --
// see the comment in fillBeamProfile() for why planes near a detector are
// unusable.
static const double kZLo = -440.0;
static const double kZHi = 40.0;
static const int kNSlice = 16;

static TH1 *findTH1(TDirectory *dir, const std::string &name)
{
  if (TObject *o = dir->Get(name.c_str()))
    if (o && o->InheritsFrom(TH1::Class()) && !o->InheritsFrom("TH2"))
      return static_cast<TH1 *>(o);
  TIter next(dir->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next()))
  {
    TObject *obj = key->ReadObj();
    if (!obj)
      continue;
    if (obj->InheritsFrom(TDirectory::Class()))
    {
      if (TH1 *h = findTH1(static_cast<TDirectory *>(obj), name))
        return h;
    }
    else if (obj->InheritsFrom(TH1::Class()) && !obj->InheritsFrom("TH2") && name == obj->GetName())
      return static_cast<TH1 *>(obj);
  }
  return nullptr;
}

// Robust centre of a crossing distribution: the median. Chosen over the mean
// because the crossing distribution has long tails (tracks originating far from
// this plane) that are not necessarily symmetric once dead channels bite.
// Uncertainty on the median ~ 1.253 * sigma / sqrt(N), with sigma taken from the
// interquartile range (IQR/1.349) so a few wild tails don't inflate it.
static bool medianAndError(TH1 *h, double &med, double &err, double &nEff, int minEntries = 500)
{
  if (!h)
    return false;
  double n = h->Integral();
  if (n < minEntries)
    return false;
  double probs[3] = {0.25, 0.5, 0.75};
  double q[3] = {0, 0, 0};
  h->GetQuantiles(3, q, probs);
  med = q[1];
  double sigma = (q[2] - q[0]) / 1.349; // IQR -> Gaussian-equivalent sigma
  if (sigma <= 0.0)
    return false;
  err = 1.253 * sigma / TMath::Sqrt(n);
  nEff = n;
  return true;
}

static bool fitOneProjection(TFile *f, const std::string &base, const char *coord,
                             double zmin, double zmax, double z0,
                             double &axis0, double &tilt, double &axis0_err, double &tilt_err,
                             double &chi2ndf, int &nUsed)
{
  std::vector<double> zs, cs, ces, zes;
  std::cout << "  plane z (mm)   median (mm)      +/-       entries\n";
  for (int k = 0; k < kNSlice; k++)
  {
    double sliceW = (kZHi - kZLo) / kNSlice;
    double zPlane = kZLo + (k + 0.5) * sliceW;
    if (zPlane < zmin || zPlane > zmax)
      continue;
    char kbuf[8];
    snprintf(kbuf, sizeof(kbuf), "%02d", k);
    TH1 *h = findTH1(f, base + "cross" + coord + "_z" + kbuf);
    double med = 0, err = 0, nEff = 0;
    if (!medianAndError(h, med, err, nEff))
    {
      std::cout << "  " << zPlane << "\t(skipped -- missing or too few entries)\n";
      continue;
    }
    std::cout << "  " << zPlane << "\t\t" << med << "\t" << err << "\t" << (long)nEff << "\n";
    zs.push_back(zPlane);
    zes.push_back(0.0);
    cs.push_back(med);
    ces.push_back(err);
  }
  nUsed = static_cast<int>(zs.size());
  if (nUsed < 3)
  {
    std::cout << "  [" << coord << "] only " << nUsed << " usable planes -- need >= 3 for a line fit\n";
    return false;
  }

  TGraphErrors g(nUsed, zs.data(), cs.data(), zes.data(), ces.data());
  TF1 fL("fL", Form("[0]+[1]*(x-%.8g)", z0), zmin, zmax);
  fL.SetParameters(cs[nUsed / 2], 0.0);
  if (g.Fit(&fL, "QRN") != 0)
  {
    std::cout << "  [" << coord << "] line fit failed\n";
    return false;
  }
  axis0 = fL.GetParameter(0);
  tilt = fL.GetParameter(1);
  axis0_err = fL.GetParError(0);
  tilt_err = fL.GetParError(1);
  chi2ndf = (fL.GetNDF() > 0) ? fL.GetChisquare() / fL.GetNDF() : -1.0;
  return true;
}

void FitBeamAxis(const char *filename, double zmin = -440.0, double zmax = 40.0,
                          const char *tag = "all")
{
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  double z0 = 0.0;
  if (const char *s = std::getenv("BEAM_AXIS_Z0"))
    z0 = std::atof(s);

  std::string base = std::string("beamAxis_") + tag + "_";

  std::cout << "\n============ Beam-axis fit (track crossings) ============\n";
  std::cout << "file:  " << filename << "\n";
  std::cout << "tag:   " << tag << "   (histogram prefix: " << base << ")\n";
  std::cout << "range: " << zmin << " <= z <= " << zmax << " mm\n";
  std::cout << "z0:    " << z0 << " mm  (BEAM_AXIS_Z0)\n\n";

  double bx = 0, tx = 0, bxe = 0, txe = 0, chi2x = -1;
  double by = 0, ty = 0, bye = 0, tye = 0, chi2y = -1;
  int nx = 0, ny = 0;

  std::cout << "--- X crossings ---\n";
  bool okx = fitOneProjection(f, base, "X", zmin, zmax, z0, bx, tx, bxe, txe, chi2x, nx);
  std::cout << "\n--- Y crossings ---\n";
  bool oky = fitOneProjection(f, base, "Y", zmin, zmax, z0, by, ty, bye, tye, chi2y, ny);

  std::cout << "\n--- fitted beam line ---\n";
  if (okx)
    std::cout << "  x(z0) = " << bx << " +/- " << bxe << " mm,  dx/dz = " << tx << " +/- " << txe
              << "  (" << tx * 1000.0 << " mrad, chi2/ndf = " << chi2x << ", " << nx << " planes)\n";
  if (oky)
    std::cout << "  y(z0) = " << by << " +/- " << bye << " mm,  dy/dz = " << ty << " +/- " << tye
              << "  (" << ty * 1000.0 << " mrad, chi2/ndf = " << chi2y << ", " << ny << " planes)\n";

  if (!okx || !oky)
  {
    std::cout << "\nOne or both projections failed -- not emitting exports.\n";
    std::cout << "=========================================================\n\n";
    f->Close();
    return;
  }

  std::cout << "\n--- set these directly (ABSOLUTE values, do not add to current) ---\n";
  std::cout << "export BEAM_AXIS_X=" << bx << "\n";
  std::cout << "export BEAM_AXIS_Y=" << by << "\n";
  std::cout << "export BEAM_AXIS_Z0=" << z0 << "\n";
  std::cout << "export BEAM_TILT_X=" << tx << "\n";
  std::cout << "export BEAM_TILT_Y=" << ty << "\n";
  std::cout << "=========================================================\n\n";

  f->Close();
}