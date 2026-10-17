// Step 1 in ROOT: yields out of ROOT, into a plain CSV histogram. The ROOT twin of extract_yields.py:
// same config, same output (<out_dir>/yields_<det>.csv), same format. Use whichever is handier.
//
//   root -l -b -q 'extract_yields.C+("examples/27Al_ap_runA.toml")'
//   root -l -b -q 'extract_yields.C+("examples/27Al_ap_runA.toml", "out/elsewhere")'   // another output folder
//                                                                                       // (relative to where root runs)
// THE CSV FORMAT (cross_section.py reads only this)
//   1D:  <var>_lo, <var>_hi, counts, variance
//   2D:  <var>_lo, <var>_hi, theta_cm_lo, theta_cm_hi, counts, variance        (angular mode)
// <var> is z (mm), e_beam (MeV) or ecm (MeV). One row per bin; under/overflow are not written.
//
// CONFIG: the [extract] section, `dets`, `out_dir`, `title` and [det_vars.<det>] of the study's TOML, exactly as for
// extract_yields.py (keys listed there). Two differences:
//   - `cut` is handed to TTree::Draw. Write it with and / or (e.g. "det == 1 and e_si > 0.5"): they become && / ||
//     here and stay valid pandas syntax for the Python version.
//   - A tree histogram puts an entry exactly on the last upper edge in the overflow (numpy puts it in the last bin).
// The TOML reader below handles only what these configs use: "strings", numbers, one-line [lists] and
// { start = .., stop = .., step = .. } tables. Anything else stops with an error.

#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace xs {

using Section = std::map<std::string, std::string>;     // key -> raw value text
using Config = std::map<std::string, Section>;          // section name ("" = top level) -> keys

// ============================================================ minimal TOML reading

std::string Trim(const std::string& s)
{
   const auto a = s.find_first_not_of(" \t\r\n");
   if (a == std::string::npos) return "";
   const auto b = s.find_last_not_of(" \t\r\n");
   return s.substr(a, b - a + 1);
}

std::string StripComment(const std::string& line)
{
   bool inString = false;
   for (size_t i = 0; i < line.size(); ++i) {
      if (line[i] == '"') inString = !inString;
      if (line[i] == '#' && !inString) return line.substr(0, i);
   }
   return line;
}

Config ReadToml(const std::string& file)
{
   std::ifstream in(file);
   if (!in) throw std::runtime_error("cannot open config " + file);
   Config config;
   std::string section, line;
   while (std::getline(in, line)) {
      line = Trim(StripComment(line));
      if (line.empty()) continue;
      if (line.rfind("[[", 0) == 0) { section = line; continue; }          // array tables (e.g. [[reference]]): not used here
      if (line[0] == '[') { section = Trim(line.substr(1, line.size() - 2)); continue; }
      const auto eq = line.find('=');
      if (eq == std::string::npos) throw std::runtime_error("cannot read config line: " + line);
      const std::string value = Trim(line.substr(eq + 1));
      if ((value[0] == '[' && value.back() != ']') || (value[0] == '{' && value.back() != '}'))
         throw std::runtime_error("lists and tables must be on one line: " + line);
      config[section][Trim(line.substr(0, eq))] = value;
   }
   return config;
}

std::string Unquote(const std::string& v)
{
   if (v.size() >= 2 && v.front() == '"' && v.back() == '"') return v.substr(1, v.size() - 2);
   return v;
}

// "[a, b]" or "{ k = v, ... }" -> its items, split on commas outside quotes
std::vector<std::string> SplitItems(const std::string& v)
{
   std::vector<std::string> items;
   std::string inner = Trim(v.substr(1, v.size() - 2)), item;
   bool inString = false;
   for (char c : inner) {
      if (c == '"') inString = !inString;
      if (c == ',' && !inString) { items.push_back(Trim(item)); item.clear(); }
      else item += c;
   }
   if (!Trim(item).empty()) items.push_back(Trim(item));
   return items;
}

// Bin edges from [e0, e1, ...] or { start = .., stop = .., step = .. } (stop included), computed as numpy does
std::vector<double> Edges(const std::string& v)
{
   std::vector<double> edges;
   if (v.front() == '[') {
      for (const auto& item : SplitItems(v)) edges.push_back(std::stod(item));
      return edges;
   }
   std::map<std::string, double> table;
   for (const auto& item : SplitItems(v)) {
      const auto eq = item.find('=');
      table[Trim(item.substr(0, eq))] = std::stod(Trim(item.substr(eq + 1)));
   }
   const double start = table.at("start"), stop = table.at("stop"), step = table.at("step");
   const long nBins = std::lround((stop - start) / step);
   for (long i = 0; i <= nBins; ++i) edges.push_back(start + step * i);
   return edges;
}

// The [extract] section of one detector: every string with {det} and the [det_vars.<det>] keys replaced
Section ExtractSection(const Config& config, const std::string& det)
{
   std::map<std::string, std::string> replacements = {{"det", det}};
   const auto vars = config.find("det_vars." + det);
   if (vars != config.end())
      for (const auto& kv : vars->second) replacements[kv.first] = Unquote(kv.second);
   const auto ex = config.find("extract");
   if (ex == config.end()) throw std::runtime_error("the config has no [extract] section");
   Section out;
   for (const auto& kv : ex->second) {
      TString value = kv.second.c_str();
      for (const auto& r : replacements) value.ReplaceAll(("{" + r.first + "}").c_str(), r.second.c_str());
      out[kv.first] = value.Data();
   }
   return out;
}

std::string Get(const Section& s, const std::string& key, const std::string& fallback = "\x01")
{
   const auto it = s.find(key);
   if (it != s.end()) return Unquote(it->second);
   if (fallback == "\x01") throw std::runtime_error("[extract] needs the key " + key);
   return fallback;
}

// ============================================================ the yields table

struct Row {
   double lo, hi, thetaLo, thetaHi, counts, variance;
};

// Sum of weights squared of a bin (global index); a histogram without Sumw2 counts each entry once
double Variance(const TH1* h, int bin)
{
   return h->GetSumw2N() > 0 ? h->GetSumw2()->At(bin) : h->GetBinContent(bin);
}

// The histogram at "folder/name". A name that itself contains '/' (e.g. "dE_(z/100)=1") is not reachable with Get(),
// so each split into an existing folder + a key of that folder is also tried. Not found: the error names the
// histograms whose path matches when case is ignored (e.g. "elastic+misc_qqq_a" for "elastic+misc_QQQ_a").
TH1* FindHistogram(TFile* file, const std::string& path)
{
   if (auto* h = dynamic_cast<TH1*>(file->Get(path.c_str()))) return h;
   // every split into folder + name: at each '/', from the last one back, and finally no folder at all
   std::vector<std::pair<std::string, std::string>> splits;
   for (size_t i = path.size(); i-- > 0;)
      if (path[i] == '/') splits.emplace_back(path.substr(0, i), path.substr(i + 1));
   splits.emplace_back("", path);
   for (const auto& [folder, name] : splits) {
      TDirectory* dir = folder.empty() ? static_cast<TDirectory*>(file) : file->GetDirectory(folder.c_str());
      if (!dir) continue;
      if (TKey* key = dir->GetKey(name.c_str()))
         if (auto* h = dynamic_cast<TH1*>(key->ReadObj())) return h;
   }
   std::string close;
   TString lower = path.c_str();
   lower.ToLower();
   std::function<void(TDirectory*, const std::string&)> walk = [&](TDirectory* dir, const std::string& prefix) {
      for (auto* obj : *dir->GetListOfKeys()) {
         auto* key = static_cast<TKey*>(obj);
         const std::string full = prefix + key->GetName();
         if (TString(key->GetClassName()).BeginsWith("TDirectory"))
            walk(dir->GetDirectory(key->GetName()), full + "/");
         else if (TString(full.c_str()).EqualTo(lower, TString::kIgnoreCase))
            close += " '" + full + "'";
      }
   };
   walk(file, "");
   throw std::runtime_error("no histogram " + path + " in " + file->GetName() + (close.empty() ? "" : "; same path ignoring case:" + close) +
                            ". Paths are case-sensitive; for upper-case names per detector use [det_vars], e.g. DET = \"QQQ\"");
}

std::vector<Row> FromHistogram(const Section& ex, const std::string& path, bool& twoD)
{
   TFile* file = TFile::Open(path.c_str());
   if (!file || file->IsZombie()) throw std::runtime_error("cannot open " + path);
   const std::string name = Get(ex, "hist");
   TH1* h = FindHistogram(file, name);

   std::vector<Row> rows;
   twoD = false;
   if (h->GetDimension() == 1) {
      for (int i = 1; i <= h->GetNbinsX(); ++i)
         rows.push_back({h->GetXaxis()->GetBinLowEdge(i), h->GetXaxis()->GetBinUpEdge(i), 0, 0,
                         h->GetBinContent(i), Variance(h, i)});
      file->Close();
      return rows;
   }

   // TH2: `axis` holds the variable; the other axis is theta_cm (angle_axis given) or is summed
   const bool varOnX = Get(ex, "axis", "x") == "x";
   TH2* h2 = static_cast<TH2*>(h);
   TAxis* varAxis = varOnX ? h2->GetXaxis() : h2->GetYaxis();
   TAxis* otherAxis = varOnX ? h2->GetYaxis() : h2->GetXaxis();
   const int nVar = varAxis->GetNbins(), nOther = otherAxis->GetNbins();
   auto global = [&](int iVar, int iOther) { return varOnX ? h2->GetBin(iVar, iOther) : h2->GetBin(iOther, iVar); };

   if (ex.count("angle_axis")) {                       // keep it 2D: one row per (variable bin, theta bin)
      twoD = true;
      for (int i = 1; i <= nVar; ++i)
         for (int j = 1; j <= nOther; ++j)
            rows.push_back({varAxis->GetBinLowEdge(i), varAxis->GetBinUpEdge(i), otherAxis->GetBinLowEdge(j),
                            otherAxis->GetBinUpEdge(j), h2->GetBinContent(global(i, j)), Variance(h2, global(i, j))});
      file->Close();
      return rows;
   }

   // project on the variable: the other axis summed whole (with under/overflow), or the fraction inside other_range
   std::vector<double> weight(nOther + 2, 1.0);
   if (ex.count("other_range")) {
      const auto range = Edges(ex.at("other_range"));
      weight.front() = weight.back() = 0.0;
      for (int j = 1; j <= nOther; ++j) {
         const double lo = otherAxis->GetBinLowEdge(j), hi = otherAxis->GetBinUpEdge(j);
         weight[j] = std::max(0.0, std::min(hi, range[1]) - std::max(lo, range[0])) / (hi - lo);
      }
   }
   for (int i = 1; i <= nVar; ++i) {
      double counts = 0, variance = 0;
      for (int j = 0; j <= nOther + 1; ++j) {
         counts += weight[j] * h2->GetBinContent(global(i, j));
         variance += weight[j] * weight[j] * Variance(h2, global(i, j));
      }
      rows.push_back({varAxis->GetBinLowEdge(i), varAxis->GetBinUpEdge(i), 0, 0, counts, variance});
   }
   file->Close();
   return rows;
}

std::vector<Row> FromTree(const Section& ex, const std::vector<std::string>& files, bool& twoD)
{
   TChain chain(Get(ex, "tree").c_str());
   for (const auto& f : files) chain.Add(f.c_str());        // TChain::Add takes wildcards in the file name
   if (chain.GetEntries() <= 0)                              // a missing tree or file must not become a table of zeros
      throw std::runtime_error("no entries in tree " + Get(ex, "tree") + " of " + files.front() + (files.size() > 1 ? " ..." : ""));

   // The cut, with and / or as && / ||; a weight branch multiplies it
   TString cut = Get(ex, "cut", "").c_str();
   cut.ReplaceAll(" and ", " && ").ReplaceAll(" or ", " || ");
   if (cut.IsNull()) cut = "1";
   TString selection = ex.count("weight_branch") ? TString::Format("(%s)*(%s)", Get(ex, "weight_branch").c_str(), cut.Data())
                                                 : cut;

   const auto edges = Edges(ex.at("bins"));
   const TString branch = Get(ex, "branch").c_str();
   gROOT->cd();
   std::vector<Row> rows;
   twoD = ex.count("angle_branch") > 0;
   if (!twoD) {
      TH1D h("h_extract_yields", "", edges.size() - 1, edges.data());
      h.Sumw2();
      const Long64_t n = chain.Draw(branch + ">>h_extract_yields", selection, "goff");
      printf("  %d file(s); %lld entries after the cut\n", chain.GetListOfFiles()->GetEntries(), n);
      for (int i = 1; i <= h.GetNbinsX(); ++i)
         rows.push_back({edges[i - 1], edges[i], 0, 0, h.GetBinContent(i), Variance(&h, i)});
      return rows;
   }
   const auto angleEdges = Edges(ex.at("angle_bins"));
   TH2D h("h_extract_yields2", "", edges.size() - 1, edges.data(), angleEdges.size() - 1, angleEdges.data());
   h.Sumw2();
   const Long64_t n = chain.Draw(TString(Get(ex, "angle_branch").c_str()) + ":" + branch + ">>h_extract_yields2",
                                 selection, "goff");                 // "y:x" in TTree::Draw
   printf("  %d file(s); %lld entries after the cut\n", chain.GetListOfFiles()->GetEntries(), n);
   for (size_t i = 1; i < edges.size(); ++i)
      for (size_t j = 1; j < angleEdges.size(); ++j)
         rows.push_back({edges[i - 1], edges[i], angleEdges[j - 1], angleEdges[j], h.GetBinContent(h.GetBin(i, j)),
                         Variance(&h, h.GetBin(i, j))});
   return rows;
}

void WriteCsv(const std::vector<Row>& rows, bool twoD, const std::string& var, const std::vector<std::string>& header,
              const std::string& outFile)
{
   FILE* f = fopen(outFile.c_str(), "w");
   if (!f) throw std::runtime_error("cannot write " + outFile);
   for (const auto& line : header) fprintf(f, "# %s\n", line.c_str());
   fprintf(f, twoD ? "%s_lo,%s_hi,theta_cm_lo,theta_cm_hi,counts,variance\n" : "%s_lo,%s_hi,counts,variance\n",
           var.c_str(), var.c_str());
   for (const auto& r : rows) {
      if (twoD) fprintf(f, "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", r.lo, r.hi, r.thetaLo, r.thetaHi, r.counts, r.variance);
      else      fprintf(f, "%.6g,%.6g,%.6g,%.6g\n", r.lo, r.hi, r.counts, r.variance);
   }
   fclose(f);
}

} // namespace xs

void extract_yields(const char* configFile, const char* outDirOverride = "")
{
   using namespace xs;
   try {
      const Config config = ReadToml(configFile);
      const std::string configDir = gSystem->GetDirName(configFile).Data();
      auto fullPath = [&](const std::string& p) -> std::string {
         if (p.rfind("~/", 0) == 0 && getenv("HOME")) return std::string(getenv("HOME")) + p.substr(1);  // "~/..." = home
         return p.empty() || p[0] == '/' ? p : configDir + "/" + p;
      };

      const Section& top = config.count("") ? config.at("") : Section();
      std::vector<std::string> dets = {""};
      if (top.count("dets")) {
         dets.clear();
         for (const auto& d : SplitItems(top.at("dets"))) dets.push_back(Unquote(d));
      }
      const std::string outDir = *outDirOverride ? std::string(outDirOverride)
                                                 : fullPath(top.count("out_dir") ? Unquote(top.at("out_dir")) : "out");
      gSystem->mkdir(outDir.c_str(), kTRUE);
      const std::string title = top.count("title") ? Unquote(top.at("title")) : "";

      for (const auto& det : dets) {
         const Section ex = ExtractSection(config, det);
         const std::string source = Get(ex, "source"), var = Get(ex, "variable");
         if (source != "hist" && source != "tree")
            throw std::runtime_error("[extract] source must be \"hist\" or \"tree\", not " + source +
                                     (source == "csv" ? " (event-list CSVs: use extract_yields.py)" : ""));
         std::vector<std::string> files;                      // `file` may be one path or a list
         if (ex.at("file").front() == '[') for (const auto& f : SplitItems(ex.at("file"))) files.push_back(fullPath(Unquote(f)));
         else files.push_back(fullPath(Get(ex, "file")));
         const std::string fileText = ex.at("file").front() == '[' ? ex.at("file") : Get(ex, "file");
         printf("--- %s: %s : %s\n", det.empty() ? "data" : det.c_str(), fileText.c_str(),
                Get(ex, source == "hist" ? "hist" : "tree").c_str());

         bool twoD = false;
         std::vector<Row> rows;
         std::string from;
         if (source == "hist") {
            rows = FromHistogram(ex, files.front(), twoD);
            from = "histogram " + Get(ex, "hist");
         } else if (source == "tree") {
            rows = FromTree(ex, files, twoD);
            from = "tree " + Get(ex, "tree") + ", branch " + Get(ex, "branch");
            if (!Get(ex, "cut", "").empty()) from += ", cut: " + Get(ex, "cut");
         } else {
            throw std::runtime_error("[extract] source must be \"hist\" or \"tree\", not " + source +
                                     (source == "csv" ? " (event-list CSVs: use extract_yields.py)" : ""));
         }

         TString detUpper = det.c_str();
         detUpper.ToUpper();
         const std::string outFile = outDir + "/" + (det.empty() ? "yields.csv" : "yields_" + det + ".csv");
         WriteCsv(rows, twoD, var,
                  {title + " " + detUpper.Data() + ": yields, binned in " + var, "from " + fileText + " : " + from,
                   "written by extract_yields.C"},
                  outFile);
         double total = 0;
         for (const auto& r : rows) total += r.counts;
         printf("  %.0f counts in %zu bins -> %s\n", total, rows.size(), outFile.c_str());
      }
   } catch (const std::exception& e) {
      fprintf(stderr, "extract_yields.C: %s\n", e.what());
      gSystem->Exit(1);                                   // a failed extraction must not look like a finished one
   }
}
