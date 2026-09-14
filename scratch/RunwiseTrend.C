// =============================================================================
// RunwiseTrend.C
//
// Track any TrackRecon.C histogram-derived quantity as a function of run
// number, to separate a DRIFTING condition (gas pressure, beam energy, beam
// steering, detector gain) from a FIXED analysis systematic.
//
// No reprocessing needed. run_27Al.sh / run_17F.sh already write one file per
// run -- "<OUT_DIR>/results_run###.root" -- and only hadd them afterwards, so
// the per-run files are sitting next to the merged one. (Note they are wiped
// by the "rm -f ${OUT_DIR}/*.root" at the top of the next run, so they reflect
// the most recent processing pass only.)
//
// Not to be confused with RunTimeSummary.C, which reads the *mapped* files for
// a single raw timing histogram. This one reads TrackRecon output.
//
// Usage:
//   // does the Ex-vs-vertexZ slope drift across the experiment?
//   root -l -b -q 'scratch/RunwiseTrend.C("Output_27Al",
//        "m27Alax_VertexReconZ_vs_Ex_p_a1c2fix_qqq","slope",-300,0)'
//
//   // mean Ex per run
//   root -l -b -q 'scratch/RunwiseTrend.C("Output_27Al",
//        "m27Alax_Ex_vs_theta_p_a1c2fix_qqq","mean")'
//
//   // beam position per run, from one crossing plane
//   root -l -b -q 'scratch/RunwiseTrend.C("Output_27Al",
//        "beamAxis_all_crossY_z007","mean")'
//
// modes: "mean"    <y> of a TH2 (or <x> of a TH1), error on the mean
//        "slope"   ProfileX of a TH2, fitted pol1 over [xLo,xHi] -- the slope
//        "rms"     spread rather than centre
//        "entries" statistics per run, as a sanity/normalisation check
//
// The headline output is the fit to a CONSTANT. chi2/ndf near 1 means the
// quantity is stable run to run and any effect is a fixed systematic;
// chi2/ndf >> 1 means something genuinely changed during the experiment.
// Compare the scatter of the points against the typical error bar -- that
// ratio is the same "is this real" question the beam-axis plane medians pose,
// and it is answered the same way.
// =============================================================================

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cctype>

// HistPlotter buries objects in per-tier folders, so search recursively.
static TH1 *findHist(TDirectory *d, const TString &name)
{
  if (TObject *o = d->Get(name.Data()))
    if (o->InheritsFrom(TH1::Class()))
      return (TH1 *)o;
  TIter next(d->GetListOfKeys());
  while (TKey *k = (TKey *)next())
  {
    TObject *o = k->ReadObj();
    if (o && o->InheritsFrom(TDirectory::Class()))
      if (TH1 *h = findHist((TDirectory *)o, name))
        return h;
  }
  return nullptr;
}

// results_run024.root -> 24 ; -1 if the name does not match
static int runNumberOf(const TString &fname)
{
  Ssiz_t i = fname.Index("results_run");
  if (i == kNPOS)
    return -1;
  TString digits;
  for (Ssiz_t j = i + 11; j < fname.Length() && isdigit(fname[j]); ++j)
    digits += fname[j];
  return digits.Length() ? digits.Atoi() : -1;
}

void RunwiseTrend(const char *dir,
                  const char *histName,
                  const char *mode = "mean",
                  double xLo = 0.0, double xHi = 0.0,
                  double minEntries = 200.0,
                  const char *outPng = "")
{
  TString m(mode);

  // --- collect and sort the per-run files ---
  std::vector<std::pair<int, TString>> files;
  TSystemDirectory sd("sd", dir);
  if (TList *l = sd.GetListOfFiles())
  {
    TIter next(l);
    while (TSystemFile *sf = (TSystemFile *)next())
    {
      TString fn = sf->GetName();
      if (sf->IsDirectory() || !fn.EndsWith(".root"))
        continue;
      int r = runNumberOf(fn);
      if (r >= 0)
        files.push_back({r, TString(dir) + "/" + fn});
    }
  }
  std::sort(files.begin(), files.end());

  if (files.empty())
  {
    std::cerr << "ERROR: no results_run*.root under " << dir << "\n"
              << "       (the next run_*.sh wipes them -- they only exist for\n"
              << "        the most recent processing pass)\n";
    return;
  }

  std::cout << "\n=========== " << histName << "   [" << m << "] ===========\n"
            << dir << "   " << files.size() << " runs\n\n"
            << "   run     value       +/-        N\n";

  TGraphErrors *g = new TGraphErrors();
  int np = 0;
  std::vector<double> vals, errs;

  for (auto &fr : files)
  {
    TFile *f = TFile::Open(fr.second, "READ");
    if (!f || f->IsZombie())
      continue;
    TH1 *h = findHist(f, histName);
    if (!h) { f->Close(); continue; }

    double v = 0.0, e = 0.0, n = h->GetEntries();
    bool ok = (n >= minEntries);

    if (ok && m == "entries") { v = n; e = std::sqrt(n); }
    else if (ok && m == "slope")
    {
      TH2 *h2 = dynamic_cast<TH2 *>(h);
      if (!h2) { std::cerr << "  (slope needs a TH2)\n"; f->Close(); return; }
      TProfile *p = h2->ProfileX(Form("p_%d", fr.first));
      p->SetDirectory(0);
      TF1 fit("fit", "pol1", (xLo == xHi ? p->GetXaxis()->GetXmin() : xLo),
                             (xLo == xHi ? p->GetXaxis()->GetXmax() : xHi));
      if (p->Fit(&fit, "QNR") == 0) { v = fit.GetParameter(1); e = fit.GetParError(1); }
      else ok = false;
      delete p;
    }
    else if (ok)
    {
      TH2 *h2 = dynamic_cast<TH2 *>(h);
      TH1 *proj = h2 ? h2->ProjectionY(Form("py_%d", fr.first)) : h;
      if (m == "rms") { v = proj->GetRMS(); e = proj->GetRMSError(); }
      else            { v = proj->GetMean(); e = proj->GetMeanError(); }
      if (h2) delete proj;
    }

    if (ok)
    {
      printf("  %4d  %10.5f  %9.5f  %9.0f\n", fr.first, v, e, n);
      g->SetPoint(np, fr.first, v);
      g->SetPointError(np, 0.0, e);
      vals.push_back(v); errs.push_back(e);
      ++np;
    }
    else
      printf("  %4d  (skipped -- %.0f entries)\n", fr.first, n);

    f->Close();
  }

  if (np < 2) { std::cout << "\nonly " << np << " usable runs.\n"; return; }

  // --- is the run-to-run scatter real, or just statistics? ---
  TF1 c("c", "pol0");
  g->Fit(&c, "QN");
  double mean = c.GetParameter(0), emean = c.GetParError(0);
  double chi2 = c.GetChisquare();
  int ndf = c.GetNDF();

  double sum = 0, sum2 = 0, esum = 0;
  for (int i = 0; i < np; ++i) { sum += vals[i]; sum2 += vals[i] * vals[i]; esum += errs[i]; }
  double rms = std::sqrt(std::max(0.0, sum2 / np - (sum / np) * (sum / np)));
  double etyp = esum / np;

  std::cout << "\n  constant fit : " << mean << " +/- " << emean << "\n"
            << "  chi2/ndf     : " << (ndf > 0 ? chi2 / ndf : 0.0) << "  over " << np << " runs\n"
            << "  scatter      : RMS " << rms << "  vs typical error " << etyp
            << "   (ratio " << (etyp > 0 ? rms / etyp : 0.0) << ")\n\n"
            << (ndf > 0 && chi2 / ndf > 2.0
                    ? "  => scatter exceeds statistics: something changed DURING the\n"
                      "     experiment (gas pressure, beam energy, steering, gains).\n"
                      "     A single run-averaged constant will not describe it.\n"
                    : "  => consistent with a constant: no run-to-run drift at this\n"
                      "     precision, so any effect here is a FIXED systematic\n"
                      "     (tables, geometry, beamE0, z_entrance) rather than a drift.\n");

  if (TString(outPng).Length())
  {
    TCanvas *cv = new TCanvas("cv", "", 1200, 700);
    g->SetTitle(Form("%s [%s];run number;%s", histName, mode, mode));
    g->SetMarkerStyle(20);
    g->Draw("AP");
    cv->SaveAs(outPng);
    std::cout << "  saved: " << outPng << "\n";
  }
}