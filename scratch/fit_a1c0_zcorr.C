// fit_a1c0_zcorr.C
// ---------------------------------------------------------------------------
// Linear-regression fit of the anode-only (A1C0) PC-z against the reference z,
// to extract the per-detector scale + offset corrections used by a1c1_zcorr()
// in TrackRecon.C.
//
// WORKFLOW
//   1. Set the 4 correction params to 0 (the default in TrackRecon.C now) and
//      run all SOURCE runs so the A1C0 benchmark histograms hold the RAW,
//      uncorrected anode-only z vs the reference z.
//   2. hadd the per-run outputs into one file (e.g. Output_av/all.root).
//   3. root -l -b -q 'fit_a1c0_zcorr.C("Output_av/all.root")'
//
// MODEL
//   The benchmark TH2s are filled as (x = reference/truth z, y = z_a1c0 raw).
//   Profiling y over x and fitting   y = A*x + B   gives the linear trend.
//   a1c1_zcorr does   z_corr = z_a1c0*(1-scale) - off,  and we want z_corr = x,
//   so:   scale = 1 - 1/A      off = B/A
//   (Equivalent to OLS on the scatter; the profile gives clean per-bin means
//    and avoids the per-event spread destabilising the fit.)
//
// CAVEAT
//   Smaller residual width alone is NOT proof of a correct scale: as scale->1
//   the estimate collapses to a constant. The fit below pins the slope, so the
//   corrected A1C0-vs-ref slope should come out ~1 by construction. Always
//   re-run with the printed params and confirm the corrected slope is ~1.
// ---------------------------------------------------------------------------

#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TKey.h>
#include <TDirectory.h>
#include <iostream>
#include <string>
#include <vector>

// Recursively search a directory tree for the first TH2 whose name matches.
static TH2 *findTH2(TDirectory *dir, const std::string &name)
{
  // direct hit in this directory
  if (TObject *o = dir->Get(name.c_str()))
    if (o->InheritsFrom(TH2::Class()))
      return static_cast<TH2 *>(o);

  TIter next(dir->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next()))
  {
    TObject *obj = key->ReadObj();
    if (obj->InheritsFrom(TDirectory::Class()))
    {
      if (TH2 *h = findTH2(static_cast<TDirectory *>(obj), name))
        return h;
    }
    else if (obj->InheritsFrom(TH2::Class()) && name == obj->GetName())
    {
      return static_cast<TH2 *>(obj);
    }
  }
  return nullptr;
}

// Try a list of candidate histogram names, fit the first one found.
// Returns true and sets scale/off on success.
static bool fitOne(TFile *f, const char *label,
                   const std::vector<std::string> &candidates,
                   double fitlo, double fithi,
                   double &scale, double &off)
{
  TH2 *h2 = nullptr;
  std::string used;
  for (const auto &nm : candidates)
  {
    h2 = findTH2(f, nm);
    if (h2)
    {
      used = nm;
      break;
    }
  }
  if (!h2)
  {
    std::cout << "[" << label << "] none of the candidate histograms were found:\n";
    for (const auto &nm : candidates)
      std::cout << "    " << nm << "\n";
    return false;
  }

  // ProfileX: mean of y (z_a1c0) in each x (reference) bin.
  TProfile *prof = h2->ProfileX((used + "_pfx").c_str());
  if (prof->GetEntries() < 1)
  {
    std::cout << "[" << label << "] profile from '" << used << "' is empty.\n";
    return false;
  }

  TF1 fpol("fpol", "pol1", fitlo, fithi);
  prof->Fit(&fpol, "QR"); // Quiet, restrict to range
  double B = fpol.GetParameter(0); // intercept
  double A = fpol.GetParameter(1); // slope
  double Aerr = fpol.GetParError(1);

  if (A == 0.0)
  {
    std::cout << "[" << label << "] degenerate slope (A=0), cannot invert.\n";
    return false;
  }

  scale = 1.0 - 1.0 / A;
  off = B / A;

  std::cout << "------------------------------------------------------------\n";
  std::cout << "[" << label << "]  hist: " << used << "\n";
  std::cout << "    entries (2D)      : " << (Long64_t)h2->GetEntries() << "\n";
  std::cout << "    fit range         : [" << fitlo << ", " << fithi << "] mm\n";
  std::cout << "    z_a1c0 = A*z_ref + B\n";
  std::cout << "      A (slope)       : " << A << "  +/- " << Aerr << "\n";
  std::cout << "      B (intercept)   : " << B << " mm\n";
  std::cout << "    => scale = 1-1/A  : " << scale << "\n";
  std::cout << "       off   = B/A    : " << off << " mm\n";
  return true;
}

void fit_a1c0_zcorr(const char *filename,
                    double fitlo = -150.0, double fithi = 150.0)
{
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::cout << "\n==================== A1C0 z-correction fit ====================\n";
  std::cout << "file: " << filename << "\n";

  double scale_qqq = 0, off_qqq = 0, scale_sx3 = 0, off_sx3 = 0;

  // QQQ: only a _vs_ref benchmark is filled (A1C2 crossover is the truth).
  bool okQ = fitOne(f, "QQQ",
                    {"Benchmark_QQQ_PCZ_A1C0_Hyb_vs_ref",
                     "Benchmark_QQQ_PCZ_A1C0_vs_ref",
                     "Benchmark_QQQ_PCZ_A1C0_Hyb_vs_qqqpczguess"},
                    fitlo, fithi, scale_qqq, off_qqq);

  // SX3: prefer the source-run geometric guess (fixed-vertex truth), fall
  // back to the A1C2 reference.
  bool okS = fitOne(f, "SX3",
                    {"Benchmark_SX3_PCZ_A1C0_Hyb_vs_sx3pczguess",
                     "Benchmark_SX3_PCZ_A1C0_Hyb_vs_ref",
                     "Benchmark_SX3_PCZ_A1C0_vs_sx3pczguess"},
                    fitlo, fithi, scale_sx3, off_sx3);

  std::cout << "\n==================== paste-ready output ====================\n";
  if (okQ || okS)
  {
    std::cout << "// --- TrackRecon.C Begin() dataset block ---\n";
    if (okQ)
    {
      std::cout << "double a1c1_z_scale_qqq = " << scale_qqq << ";\n";
      std::cout << "double a1c1_z_off_qqq   = " << off_qqq << ";\n";
    }
    if (okS)
    {
      std::cout << "double a1c1_z_scale_sx3 = " << scale_sx3 << ";\n";
      std::cout << "double a1c1_z_off_sx3   = " << off_sx3 << ";\n";
    }
    std::cout << "\n// --- or as environment overrides ---\n";
    if (okQ)
    {
      std::cout << "export A1C1_Z_SCALE_QQQ=" << scale_qqq << "\n";
      std::cout << "export A1C1_Z_OFF_QQQ=" << off_qqq << "\n";
    }
    if (okS)
    {
      std::cout << "export A1C1_Z_SCALE_SX3=" << scale_sx3 << "\n";
      std::cout << "export A1C1_Z_OFF_SX3=" << off_sx3 << "\n";
    }
  }
  else
  {
    std::cout << "No histograms fitted -- check the file and run with source data.\n";
  }
  std::cout << "============================================================\n\n";

  f->Close();
}
