// Aggregates every run's raw PC energy-calibration points (written by
// TrackRecon.C's Terminate() into pc_calib_raw/points_*.dat, one file per run
// so parallel jobs never collide) into a single per-wire linear fit, pooling
// BOTH the 17F and 27Al datasets together -- the PC's wire-level gain isn't
// expected to differ meaningfully between the two campaigns, so calibration
// statistics from either dataset's alpha-source or proton-scattering runs are
// combined into one fit rather than kept separate.
//
// Run non-interactively from the top-level ANASEN-Analysis directory once all
// desired calibration runs have completed:
//   root -q -l -b -e '.L pccal/fit_pc_energy_calibration.C' -e 'fit_pc_energy_calibration()'
//
// Writes pc_energy_calibration.dat (wire slope intercept, 48 rows), which
// TrackRecon.C's Begin() loads automatically on the next run to populate
// PC_Events_calibrated. Re-run this any time pc_calib_raw/ has new points --
// it always re-reads and re-fits everything currently on disk, so it's safe
// to call repeatedly as more calibration runs accumulate.
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <vector>

void fit_pc_energy_calibration()
{
  std::vector<std::pair<double, double>> pts[48]; // [wire] -> (ADC, dE_gas MeV)

  TSystemDirectory dir("pc_calib_raw", "pc_calib_raw");
  TList *files = dir.GetListOfFiles();
  if (!files)
  {
    std::cerr << "fit_pc_energy_calibration: pc_calib_raw/ not found or empty -- "
              << "run TrackRecon.C with PC energy calibration enabled first." << std::endl;
    return;
  }

  int nFiles = 0;
  TIter next(files);
  TSystemFile *f;
  while ((f = (TSystemFile *)next()))
  {
    TString name = f->GetName();
    if (f->IsDirectory() || !name.BeginsWith("points_") || !name.EndsWith(".dat"))
      continue;
    std::ifstream infile(std::string("pc_calib_raw/") + name.Data());
    if (!infile.is_open())
      continue;
    int wire;
    double adc, dE_gas;
    while (infile >> wire >> adc >> dE_gas)
      if (wire >= 0 && wire < 48)
        pts[wire].push_back({adc, dE_gas});
    ++nFiles;
  }
  std::cout << "fit_pc_energy_calibration: read " << nFiles << " run file(s) from pc_calib_raw/" << std::endl;

  std::ofstream outfile("pc_energy_calibration.dat");
  outfile << std::scientific << std::setprecision(6);
  for (int wire = 0; wire < 48; ++wire)
  {
    double slope = 1.0, intercept = 0.0;
    bool ok = pts[wire].size() >= 2;
    if (ok)
    {
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      for (const auto &p : pts[wire])
      {
        sx += p.first;
        sy += p.second;
        sxx += p.first * p.first;
        sxy += p.first * p.second;
      }
      double n = static_cast<double>(pts[wire].size());
      double denom = n * sxx - sx * sx;
      ok = std::isfinite(denom) && std::abs(denom) > 1e-12;
      if (ok)
      {
        slope = (n * sxy - sx * sy) / denom;
        intercept = (sy - slope * sx) / n;
      }
    }
    if (!ok)
      std::cerr << "fit_pc_energy_calibration: wire " << wire << " has too few points (" << pts[wire].size()
                << ") to fit -- writing identity (slope=1, intercept=0)" << std::endl;
    outfile << wire << " " << slope << " " << intercept << "\n";
  }
  outfile.close();
  std::cout << "fit_pc_energy_calibration: wrote pc_energy_calibration.dat" << std::endl;
}
