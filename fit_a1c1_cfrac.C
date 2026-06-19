// fit_a1c1_cfrac.C
// =============================================================================
// Offline per-cell "gain match" for the A1C1 linear centre-fold model.
//
// The raw cfrac = cpMax/(apSum+cpMax) pedestal floats with the arbitrary
// anode/cathode gain, so the model cfrac = cfmin[cell] + k[cell]*fold has to be
// re-fit per dataset. This macro does that fit OFFLINE from an already-processed
// results file (no reprocessing of raw data), replacing the in-line autocal that
// used to live in TrackRecon.C.
//
// It reads the reference histograms
//     Benchmark_QQQ_A1C1_cfrac_vs_ref   (X = true z = pcz_ref, Y = cfrac)
//     Benchmark_SX3_A1C1_cfrac_vs_ref
// filled from A1C2 events (where the true z is known from the crossover), folds
// each cell about its centre (fold = |z - z_centre| / halfcell, in [0,1]) and
// does a weighted least-squares fit of cfrac vs fold per cathode cell. QQQ and
// SX3 are combined (same proportional-counter cells).
//
// Output: a per-cell table, paste-ready C++ array literals for TrackRecon.C
// (a1c1_cfmin_<DATASET>[7] / a1c1_k_<DATASET>[7]), and optionally a .dat file.
//
// Usage (compile once with + then call):
//   root -l -b -q 'fit_a1c1_cfrac.C+("Output_17F/Output_17F.root")'   // build
//   MAIN band:  fit_a1c1_cfrac("Output_17F/Output_17F.root","17F",0.25,1.05)
//   LOW band:   fit_a1c1_cfrac("Output_17F/Output_17F.root","17F",0.00,0.25,50,false,"2")
//   27Al:       fit_a1c1_cfrac("Output_27Al/output_27Al.root","27Al",0.10,1.05)
// Args: cfrac_min, cfrac_max (band window), min_counts (per-cell fit threshold),
//       write_dat (dump a1c1_cfrac_<DATASET>.dat), band ("" main, "2" low band ->
//       prints a1c1_cfmin2_<DATASET>), all_cells (pool all 7 cells into one fit
//       and broadcast -- robust when per-cell stats are thin / pedestals uniform).
//   ALL-CELLS: fit_a1c1_cfrac("Output_cal/17F.root","17F",0.25,1.05,50,false,"",true)
// =============================================================================

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TString.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {
// cathode wire grid; cell i spans zg[i] (high) .. zg[i+1] (low)
const double zg[8] = {147.998, 101.946, 59.7634, 19.6965,
                      -19.6965, -59.7634, -101.946, -147.998};

int cellOf(double z) {
  for (int i = 0; i < 7; ++i)
    if (z <= zg[i] && z > zg[i + 1]) return i;
  return -1;
}

// recursively search a directory tree for a TH2 by name (folders vary by group)
TH2 *findTH2(TDirectory *dir, const std::string &name) {
  if (TH2 *h = dynamic_cast<TH2 *>(dir->Get(name.c_str()))) return h;
  TIter next(dir->GetListOfKeys());
  TKey *key;
  while ((key = (TKey *)next())) {
    TObject *obj = key->ReadObj();
    if (TDirectory *sub = dynamic_cast<TDirectory *>(obj)) {
      if (TH2 *h = findTH2(sub, name)) return h;
    }
  }
  return nullptr;
}
} // namespace

// cfrac_min/cfrac_max bracket the band to fit: [0.25, 1.05] for the MAIN band,
// [0.0, 0.25] for the 17F low (incomplete-integration) band. band="" prints the
// a1c1_cfmin_<dataset> arrays; band="2" prints a1c1_cfmin2_<dataset> (the low band).
void fit_a1c1_cfrac(const char *filename, const char *dataset = "DATA",
                    double cfrac_min = 0.0, double cfrac_max = 1.05,
                    long min_counts = 50, bool write_dat = false,
                    const char *band = "", bool all_cells = false) {
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "fit_a1c1_cfrac: cannot open " << filename << std::endl;
    return;
  }

  const char *names[2] = {"Benchmark_QQQ_A1C1_cfrac_vs_ref",
                          "Benchmark_SX3_A1C1_cfrac_vs_ref"};

  // weighted OLS accumulators per cell:  cfrac = cfmin + k*fold
  double n[7] = {0}, sf[7] = {0}, sff[7] = {0}, sc[7] = {0}, sfc[7] = {0};

  for (int s = 0; s < 2; ++s) {
    TH2 *h = findTH2(f, names[s]);
    if (!h) {
      std::cout << "  (note: " << names[s] << " not found -- skipping)" << std::endl;
      continue;
    }
    std::cout << "Using " << names[s] << "  (" << h->GetEntries() << " entries)" << std::endl;
    const int nx = h->GetNbinsX(), ny = h->GetNbinsY();
    for (int ix = 1; ix <= nx; ++ix) {
      const double z = h->GetXaxis()->GetBinCenter(ix);
      const int c = cellOf(z);
      if (c < 0) continue;
      const double zc = 0.5 * (zg[c] + zg[c + 1]);
      const double half = 0.5 * (zg[c] - zg[c + 1]);
      if (half <= 0.0) continue;
      const double fold = TMath::Abs(z - zc) / half; // [0,1]
      for (int iy = 1; iy <= ny; ++iy) {
        const double w = h->GetBinContent(ix, iy);
        if (w <= 0.0) continue;
        const double cfrac = h->GetYaxis()->GetBinCenter(iy);
        if (cfrac < cfrac_min || cfrac > cfrac_max) continue;
        n[c] += w;
        sf[c] += w * fold;
        sff[c] += w * fold * fold;
        sc[c] += w * cfrac;
        sfc[c] += w * fold * cfrac;
      }
    }
  }

  // per-cell weighted least-squares; keep a fallback for starved cells
  const double cfmin_fallback = 0.40, k_fallback = 0.075;
  double cfmin[7], k[7];

  if (all_cells) {
    // pool every cell into one fit. fold is already cell-relative (|z-zc|/half),
    // so pooling the (fold, cfrac) pairs is valid -- robust when per-cell stats
    // are thin and the pedestals are ~uniform. Same value written to all 7 cells.
    double N = 0, SF = 0, SFF = 0, SC = 0, SFC = 0;
    for (int c = 0; c < 7; ++c) { N += n[c]; SF += sf[c]; SFF += sff[c]; SC += sc[c]; SFC += sfc[c]; }
    double cf = cfmin_fallback, kk = k_fallback;
    const char *status = "FALLBACK (starved)";
    if (N > min_counts) {
      const double denom = N * SFF - SF * SF;
      if (TMath::Abs(denom) > 1e-9) {
        const double kfit = (N * SFC - SF * SC) / denom;
        const double cffit = (SC - kfit * SF) / N;
        if (kfit > 0.01) { cf = cffit; kk = kfit; status = "fit"; }
        else status = "FALLBACK (flat k)";
      }
    }
    for (int c = 0; c < 7; ++c) { cfmin[c] = cf; k[c] = kk; }
    printf("\n ALL-CELLS pooled fit: cfmin=%.6f  k=%.6f  n(eff)=%.0f   %s\n", cf, kk, N, status);
  } else {
    std::cout << "\n cell    cfmin         k         n(eff)   status" << std::endl;
    std::cout << " ----  ---------   ---------   --------   ------" << std::endl;
    for (int c = 0; c < 7; ++c) {
      double cf = cfmin_fallback, kk = k_fallback;
      const char *status = "FALLBACK (starved)";
      if (n[c] > min_counts) {
        const double denom = n[c] * sff[c] - sf[c] * sf[c];
        if (TMath::Abs(denom) > 1e-9) {
          const double kfit = (n[c] * sfc[c] - sf[c] * sc[c]) / denom;
          const double cffit = (sc[c] - kfit * sf[c]) / n[c];
          if (kfit > 0.01) { // reject degenerate/flat fit
            cf = cffit;
            kk = kfit;
            status = "fit";
          } else {
            status = "FALLBACK (flat k)";
          }
        }
      }
      cfmin[c] = cf;
      k[c] = kk;
      printf("  %d    %9.6f   %9.6f   %8.0f   %s\n", c, cf, kk, n[c], status);
    }
  }

  // paste-ready literals for TrackRecon.C (band="" -> main set, "2" -> low band)
  printf("\n// ---- paste into TrackRecon.C (cfrac in [%.2f, %.2f]) ----\n", cfrac_min, cfrac_max);
  printf("static const double a1c1_cfmin%s_%s[7] = {", band, dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", cfmin[c], c < 6 ? ", " : "};\n");
  printf("static const double a1c1_k%s_%s[7]     = {", band, dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", k[c], c < 6 ? ", " : "};\n");

  if (write_dat) {
    const std::string fn = std::string("a1c1_cfrac_") + dataset + ".dat";
    std::ofstream out(fn);
    out << "# cell  cfmin  k    (offline fit cfrac = cfmin + k*fold; dataset=" << dataset << ")\n";
    for (int c = 0; c < 7; ++c) out << c << "  " << cfmin[c] << "  " << k[c] << "\n";
    out.close();
    std::cout << "\nwrote " << fn << std::endl;
  }

  f->Close();
}

// =============================================================================
// REFERENCE-FREE gain match.
//
// Does NOT use the A1C2 true-z reference. Reads the per-cell A1C1 cfrac
// histogram Benchmark_{QQQ,SX3}_trueA1C1_cfrac_vs_cell (cell assignment is pure
// geometry) and recovers the constants from the EDGES of each cell's cfrac
// distribution: within a cell cfrac = cfmin + k*fold with fold in [0,1], so the
// low percentile ~ cfmin (fold->0, cell centre) and the high percentile ~
// cfmin+k (fold->1, at the wire). Robust percentiles (plo/phi) are used instead
// of min/max.
//
// This aligns the per-cell pedestals (the actual gain match) WITHOUT a reference,
// but it cannot fix the absolute z scale -- it assumes the linear model and that
// each cell is illuminated across its length. If shared_k=true, a single
// geometric k (median of the per-cell spans) is used for every cell and only the
// per-cell cfmin pedestal is taken from the data -- the most robust mode.
//
// Usage:
//   root -l -b -q 'fit_a1c1_cfrac.C+("Output_17F/Output_17F.root")'   // compile
//   root -l -b -q 'fit_a1c1_cfrac_reffree("Output_17F/Output_17F.root","17F")'
// =============================================================================
void fit_a1c1_cfrac_reffree(const char *filename, const char *dataset = "DATA",
                            double plo = 0.05, double phi = 0.95,
                            bool shared_k = false, double min_counts = 50,
                            bool write_dat = false) {
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "fit_a1c1_cfrac_reffree: cannot open " << filename << std::endl;
    return;
  }

  const char *names[2] = {"Benchmark_QQQ_trueA1C1_cfrac_vs_cell",
                          "Benchmark_SX3_trueA1C1_cfrac_vs_cell"};

  // sum the QQQ + SX3 per-cell cfrac histograms into one (cell x cfrac)
  TH2 *acc = nullptr;
  for (int s = 0; s < 2; ++s) {
    TH2 *h = findTH2(f, names[s]);
    if (!h) {
      std::cout << "  (note: " << names[s] << " not found -- skipping)" << std::endl;
      continue;
    }
    std::cout << "Using " << names[s] << "  (" << h->GetEntries() << " entries)" << std::endl;
    if (!acc) { acc = (TH2 *)h->Clone("acc_cfrac_vs_cell"); acc->SetDirectory(nullptr); }
    else acc->Add(h);
  }
  if (!acc) { std::cerr << "no cfrac_vs_cell histograms found" << std::endl; f->Close(); return; }

  const double cfmin_fallback = 0.40, k_fallback = 0.075;
  double cfmin[7], k[7], span[7];
  bool ok[7];

  std::cout << "\n cell    cfmin       cfmin+k        k        n     status" << std::endl;
  std::cout << " ----  ---------   ---------   ---------   ------   ------" << std::endl;
  for (int c = 0; c < 7; ++c) {
    TH1D *proj = acc->ProjectionY(Form("proj_c%d", c), c + 1, c + 1);
    const double n = proj->Integral();
    double lo = cfmin_fallback, hi = cfmin_fallback + k_fallback;
    ok[c] = false;
    if (n > min_counts) {
      double q[2], p[2] = {plo, phi};
      proj->GetQuantiles(2, q, p);
      lo = q[0];
      hi = q[1];
      ok[c] = (hi - lo > 0.01);
    }
    cfmin[c] = ok[c] ? lo : cfmin_fallback;
    span[c] = ok[c] ? (hi - lo) : k_fallback;
    k[c] = span[c];
    printf("  %d    %9.6f   %9.6f   %9.6f   %6.0f   %s\n", c, cfmin[c],
           cfmin[c] + span[c], span[c], n, ok[c] ? "fit" : "FALLBACK");
    delete proj;
  }

  if (shared_k) {
    // median of the well-determined per-cell spans -> one geometric k for all cells
    double v[7];
    int m = 0;
    for (int c = 0; c < 7; ++c) if (ok[c]) v[m++] = span[c];
    double kg = k_fallback;
    if (m > 0) {
      std::sort(v, v + m);
      kg = (m % 2) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
    }
    for (int c = 0; c < 7; ++c) k[c] = kg;
    std::cout << "\nshared geometric k = " << kg << " (median of " << m << " fitted cells)" << std::endl;
  }

  printf("\n// ---- paste into TrackRecon.C (reference-free gain match) ----\n");
  printf("static const double a1c1_cfmin_%s[7] = {", dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", cfmin[c], c < 6 ? ", " : "};\n");
  printf("static const double a1c1_k_%s[7]     = {", dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", k[c], c < 6 ? ", " : "};\n");

  if (write_dat) {
    const std::string fn = std::string("a1c1_cfrac_") + dataset + ".dat";
    std::ofstream out(fn);
    out << "# cell  cfmin  k   (reference-free gain match, percentiles " << plo << "/" << phi
        << "; dataset=" << dataset << ")\n";
    for (int c = 0; c < 7; ++c) out << c << "  " << cfmin[c] << "  " << k[c] << "\n";
    out.close();
    std::cout << "\nwrote " << fn << std::endl;
  }

  delete acc;
  f->Close();
}