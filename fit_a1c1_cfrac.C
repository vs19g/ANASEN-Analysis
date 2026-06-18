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

namespace {
// cathode wire grid; cell i spans zg[i] (high) .. zg[i+1] (low)
const double zg[8] = {147.998, 101.946, 59.7634, 19.6965,
                      -19.6965, -59.7634, -101.946, -147.998};

int cellOf(double z) {
  for (int i = 0; i < 7; ++i)
    if (z <= zg[i] && z > zg[i + 1]) return i;
  return -1;
}

// recursively search a directory tree for a TH2 by name
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

// =============================================================================
// GLOBAL WITH REFERENCE (Least Squares Fit)
// =============================================================================
void fit_a1c1_cfrac(const char *filename, const char *dataset = "DATA",
                           double cfrac_min = 0.0, double cfrac_max = 1.05,
                           long min_counts = 50, bool write_dat = false,
                           const char *band = "") {
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return;
  }

  const char *names[2] = {"Benchmark_QQQ_A1C1_cfrac_vs_ref",
                          "Benchmark_SX3_A1C1_cfrac_vs_ref"};

  // Global weighted OLS accumulators across ALL cells
  double n_tot = 0, sf_tot = 0, sff_tot = 0, sc_tot = 0, sfc_tot = 0;

  for (int s = 0; s < 2; ++s) {
    TH2 *h = findTH2(f, names[s]);
    if (!h) continue;
    
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
        
        // Accumulate globally
        n_tot += w;
        sf_tot += w * fold;
        sff_tot += w * fold * fold;
        sc_tot += w * cfrac;
        sfc_tot += w * fold * cfrac;
      }
    }
  }

  // Global least-squares calculation
  double cfmin_global = 0.40, k_global = 0.075;
  const char *status = "FALLBACK";

  if (n_tot > min_counts) {
    const double denom = n_tot * sff_tot - sf_tot * sf_tot;
    if (TMath::Abs(denom) > 1e-9) {
      const double kfit = (n_tot * sfc_tot - sf_tot * sc_tot) / denom;
      const double cffit = (sc_tot - kfit * sf_tot) / n_tot;
      if (kfit > 0.01) { // reject degenerate/flat fit
        cfmin_global = cffit;
        k_global = kfit;
        status = "FIT";
      }
    }
  }

  printf("\n=== GLOBAL FIT RESULT ===\n");
  printf(" cfmin: %9.6f  |  k: %9.6f  |  N: %8.0f  | Status: %s\n", 
         cfmin_global, k_global, n_tot, status);

  // Print 7-element array with the identical global value for backward compatibility
  printf("\n// ---- paste into TrackRecon.C (Global Fit) ----\n");
  printf("static const double a1c1_cfmin%s_%s[7] = {", band, dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", cfmin_global, c < 6 ? ", " : "};\n");
  
  printf("static const double a1c1_k%s_%s[7]     = {", band, dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", k_global, c < 6 ? ", " : "};\n");

  if (write_dat) {
    const std::string fn = std::string("a1c1_cfrac_global_") + dataset + ".dat";
    std::ofstream out(fn);
    out << "# GLOBAL cfmin k (offline fit cfrac = cfmin + k*fold; dataset=" << dataset << ")\n";
    for (int c = 0; c < 7; ++c) out << c << "  " << cfmin_global << "  " << k_global << "\n";
    out.close();
  }
  
  f->Close();
}

// =============================================================================
// GLOBAL REFERENCE-FREE (Percentile Fit)
// =============================================================================
void fit_a1c1_cfrac_reffree_global(const char *filename, const char *dataset = "DATA",
                                   double plo = 0.05, double phi = 0.95,
                                   double min_counts = 50, bool write_dat = false) {
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) return;

  const char *names[2] = {"Benchmark_QQQ_trueA1C1_cfrac_vs_cell",
                          "Benchmark_SX3_trueA1C1_cfrac_vs_cell"};

  TH2 *acc = nullptr;
  for (int s = 0; s < 2; ++s) {
    TH2 *h = findTH2(f, names[s]);
    if (!h) continue;
    if (!acc) { acc = (TH2 *)h->Clone("acc_cfrac_vs_cell"); acc->SetDirectory(nullptr); }
    else acc->Add(h);
  }
  if (!acc) { std::cerr << "no cfrac_vs_cell histograms found\n"; f->Close(); return; }

  // Project the ENTIRE detector onto Y to get a global distribution
  TH1D *proj = acc->ProjectionY("proj_global");
  const double n_tot = proj->Integral();

  double cfmin_global = 0.40, k_global = 0.075;
  const char *status = "FALLBACK";

  if (n_tot > min_counts) {
    double q[2], p[2] = {plo, phi};
    proj->GetQuantiles(2, q, p);
    double lo = q[0];
    double hi = q[1];
    
    if (hi - lo > 0.01) {
      cfmin_global = lo;
      k_global = hi - lo;
      status = "FIT";
    }
  }

  printf("\n=== GLOBAL REFERENCE-FREE FIT ===\n");
  printf(" cfmin: %9.6f  |  k: %9.6f  |  N: %8.0f  | Status: %s\n", 
         cfmin_global, k_global, n_tot, status);

  // Print arrays
  printf("\n// ---- paste into TrackRecon.C (Global Ref-Free) ----\n");
  printf("static const double a1c1_cfmin_%s[7] = {", dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", cfmin_global, c < 6 ? ", " : "};\n");
  
  printf("static const double a1c1_k_%s[7]     = {", dataset);
  for (int c = 0; c < 7; ++c) printf("%.6f%s", k_global, c < 6 ? ", " : "};\n");

  delete proj;
  delete acc;
  f->Close();
}