// fit_xy_offset.C
// ---------------------------------------------------------------------------
// Tests the "beam/source XY offset (tilted rod)" hypothesis by fitting the
// phi-dependence of the A1C2 vertex-z residual and the A1C0 z residual.
//
//   z_reco - z_true = -(d/tan(theta)) * cos(phi - phi_d)
//
// where (d, phi_d) is the transverse source offset. So a cos(phi - phi_d)
// modulation of the residual is the signature of an XY offset; a flat line
// means no offset. The PHASE phi_d must be the SAME for SX3 and QQQ if the
// offset is a real, shared source/beam displacement (the falsifiable test).
// The AMPLITUDE scales as 1/tan(theta), so QQQ (forward) > SX3 (transverse)
// for the same physical d.
//
// USAGE  (run ONE source position at a time -- combining z-positions and
//         per-run offsets washes out the modulation):
//   root -l -b -q 'fit_xy_offset.C("Output_a/results_run018.root")'
//
// Reads the Diag_XYoffset histograms filled by TrackRecon.C:
//   Diag_{SX3,QQQ}_A1C2_vtxZ_resid_vs_phi   (x=phi deg, y=vtxZ - source_vertex)
//   Diag_{SX3,QQQ}_A1C0_zresid_vs_phi        (x=phi deg, y=pcz_a1c0 - pcz_ref)
// ---------------------------------------------------------------------------

#include <TFile.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TMath.h>
#include <iostream>
#include <string>
#include <vector>

static TH2 *findTH2(TDirectory *dir, const std::string &name)
{
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
      return static_cast<TH2 *>(obj);
  }
  return nullptr;
}

// Fit A + B*cos(x*pi/180 - C) to the ProfileX of a (phi, residual) TH2.
// Returns true and reports the offset constant A, amplitude B, phase C(deg).
static bool fitCos(TFile *f, const char *label, const std::string &hname,
                   double &A, double &B, double &Cdeg)
{
  TH2 *h2 = findTH2(f, hname);
  if (!h2)
  {
    std::cout << "[" << label << "] histogram not found: " << hname << "\n";
    return false;
  }
  TProfile *prof = h2->ProfileX((hname + "_pfx").c_str());
  if (prof->GetEntries() < 10)
  {
    std::cout << "[" << label << "] too few entries to fit (" << prof->GetEntries() << ")\n";
    return false;
  }

  // A + B*cos(phi - C), phi in radians (x stored in degrees)
  TF1 fc("fc", "[0]+[1]*cos(x*TMath::DegToRad()-[2])", -180, 180);
  fc.SetParameters(prof->GetMean(2), prof->GetRMS(2), 0.0);
  fc.SetParLimits(1, 0.0, 1e4);   // amplitude >= 0
  prof->Fit(&fc, "QR");

  A = fc.GetParameter(0);
  B = fc.GetParameter(1);
  Cdeg = fc.GetParameter(2) * TMath::RadToDeg();
  // fold phase into (-180,180]
  while (Cdeg > 180.0) Cdeg -= 360.0;
  while (Cdeg <= -180.0) Cdeg += 360.0;

  double chi2ndf = (fc.GetNDF() > 0) ? fc.GetChisquare() / fc.GetNDF() : -1.0;
  std::cout << "------------------------------------------------------------\n";
  std::cout << "[" << label << "]  " << hname << "\n";
  std::cout << "    constant  A : " << A << " mm  (mean residual)\n";
  std::cout << "    amplitude B : " << B << " mm  +/- " << fc.GetParError(1) << "\n";
  std::cout << "    phase phi_d : " << Cdeg << " deg  (offset direction)\n";
  std::cout << "    chi2/ndf    : " << chi2ndf << "\n";
  return true;
}

void fit_xy_offset(const char *filename)
{
  TFile *f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "Cannot open file: " << filename << std::endl;
    return;
  }

  std::cout << "\n================= XY-offset (rod-tilt) test =================\n";
  std::cout << "file: " << filename << "\n";
  std::cout << "Model: residual = A - (d/tan(theta))*cos(phi - phi_d)\n";
  std::cout << "  -> amplitude B ~ d/tan(theta), phase phi_d = offset azimuth.\n";
  std::cout << "  -> phi_d must MATCH between SX3 and QQQ for a real offset.\n";

  double A, B, Cdeg;
  double phiSX3_vtx = 0, phiQQQ_vtx = 0;
  bool okSvtx = false, okQvtx = false;

  std::cout << "\n--- A1C2 vertex-z residual (cleanest: theta-projection) ---\n";
  if ((okSvtx = fitCos(f, "SX3 vtxZ", "Diag_SX3_A1C2_vtxZ_resid_vs_phi", A, B, Cdeg)))
    phiSX3_vtx = Cdeg;
  if ((okQvtx = fitCos(f, "QQQ vtxZ", "Diag_QQQ_A1C2_vtxZ_resid_vs_phi", A, B, Cdeg)))
    phiQQQ_vtx = Cdeg;

  std::cout << "\n--- A1C0 z residual (anode-only minus A1C2 ref) ---\n";
  fitCos(f, "SX3 A1C0", "Diag_SX3_A1C0_zresid_vs_phi", A, B, Cdeg);
  fitCos(f, "QQQ A1C0", "Diag_QQQ_A1C0_zresid_vs_phi", A, B, Cdeg);

  std::cout << "\n==================== verdict ====================\n";
  if (okSvtx && okQvtx)
  {
    double dphi = phiSX3_vtx - phiQQQ_vtx;
    while (dphi > 180.0) dphi -= 360.0;
    while (dphi <= -180.0) dphi += 360.0;
    std::cout << "vertex-z phase:  SX3 = " << phiSX3_vtx
              << " deg,  QQQ = " << phiQQQ_vtx << " deg\n";
    std::cout << "phase difference = " << dphi << " deg\n";
    if (TMath::Abs(dphi) < 30.0)
      std::cout << "=> phases AGREE -> consistent with a real shared XY offset"
                   " (rod tilt). Offset direction ~ " << phiSX3_vtx << " deg.\n";
    else
      std::cout << "=> phases DISAGREE -> not a single shared XY offset;"
                   " likely per-detector calibration instead.\n";
  }
  else
  {
    std::cout << "Could not fit both detectors -- check the run has A1C2 stats.\n";
  }
  std::cout << "=================================================\n\n";

  f->Close();
}
