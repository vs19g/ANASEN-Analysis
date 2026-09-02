// overlay_2d.C
//
// Overlays several 2D histograms from one ROOT file onto a single canvas,
// one fill colour per histogram, with a dashed y = x reference line and a
// legend. Originally built around "scat"/marker-coloured points; now uses
// "box" so each histogram's colour comes from SetFillColor (which "colz"
// does not respect when several histograms are drawn on the same pad).
//
// Which histograms to draw is now HARDCODED at the top of overlay_2d(),
// built from the selectors below rather than typed out as a full CSV
// string each time. Histogram names follow TrackRecon.C's convention:
//
//   <reaction>_ETrack_vs_EKin<state><ejtag>_<topology>_<detector>
//
// e.g. "m27Alax_ETrack_vs_EKin2235keV_p_a1c2fix_sx3". Edit reaction /
// ejtag / topology / detector / states / labels below to change what's
// plotted -- everything else in the file stays the same.
//
// Other notes:
//   - rebinX / rebinY: applied to every cloned histogram before drawing,
//     to smooth out the speckled look of narrow source binning.
//   - minz: floor applied via SetMinimum() so bins with only one or two
//     entries (noise at the edges of each locus) don't draw at all.
//   - xlo[] / xhi[]: hardcoded per-state display windows (MeV), same
//     order as states[]. Applied by zeroing each histogram's bins outside
//     its own [xlo[i], xhi[i)) window before drawing -- NOT via
//     SetRangeUser, because SetRangeUser on the first-drawn histogram
//     would also shrink the pad's FRAME to that histogram's own narrow
//     window (the first Draw() call is what sets the frame's axis
//     limits), clipping every other state out of the picture entirely
//     rather than just cropping its own content. Zeroing bins keeps the
//     frame at its full native range while still restricting each
//     histogram's drawn content to its own tier.
//     THIS IS CURRENTLY A NO-OP either way: TrackRecon.C still fills each
//     of the six histograms exclusively within its own tier at fill time
//     (the if/else chain at TrackRecon.C:4015-4041), so every bin outside
//     a given histogram's own window is already zero -- there's nothing
//     there to zero out. This starts doing real work only once
//     TrackRecon.C fills all six histograms unconditionally for every
//     event (dropping that if/else so every proton event contributes to
//     all six, each under its own Ex assumption); only then does this
//     become the thing separating the states, rather than a redundant
//     no-op.
//   - The diagonal line's range is read from the first histogram's own
//     axis limits (not hardcoded), with enough sample points (SetNpx) to
//     reach both corners cleanly.
//
// Call, now just: rootFile, axis labels, and the drawing knobs --
//
//   root -l -b -q 'scratch/overlay_2d.C("Output_27Al/output_27Al.root", "Tracked Beam Energy (MeV)", "Kinematic Energy (MeV)", 2, 2, 2.0)'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TF1.h"
#include <iostream>
#include <vector>

void overlay_2d(TString rootFile, TString xAxisLabel, TString yAxisLabel,
                int rebinX = 8, int rebinY = 8, double minz = 5.0)
{
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    // ---- Hardcoded selection: edit these to choose what to plot ----
    TString reaction = "m27Alax"; // reaction tag, e.g. "m27Alax", "m17Fax"
    TString ejtag = "p";          // ejectile: "p" (proton) or "a" (alpha)
    TString topology = "a2c0";    // e.g. "a1c0", "a1c1", "a1c2fix", "a2c0", "a1c1c2"
    TString detector = "sx3";     // "sx3" or "qqq"

    // Excited states to overlay, in the same order as labels[]/xlo[]/xhi[]
    // below. Must match TrackRecon.C's EKin<state> suffixes exactly --
    // note 3498's tag is lowercase "kev" (TrackRecon.C:4028), unlike every
    // other state's "keV". 3774 keV is deliberately left out (too close
    // to 3498, spin-parity predicts it's weaker anyway).
    TString states[] = {"GS", "2235keV", "3498kev", "4809keV", "5614keV", "6550keV"};
    TString labels[] = {"Ground State", "2.235 MeV", "3.498 MeV",
                        "4.809 MeV", "5.614 MeV", "6.550 MeV"};

    // Per-state display window (MeV), same order as states[] above.
    // See the header comment for why these are currently a no-op.
    double xlo[] = {0, 4, 12, 24, 36, 42};
    double xhi[] = {4, 12, 24, 36, 42, 90};

    const size_t nStates = sizeof(states) / sizeof(states[0]);
    if (sizeof(labels) / sizeof(labels[0]) != nStates ||
        sizeof(xlo) / sizeof(xlo[0]) != nStates ||
        sizeof(xhi) / sizeof(xhi[0]) != nStates)
    {
        std::cerr << "overlay_2d: states[]/labels[]/xlo[]/xhi[] must all be "
                     "the same length -- check your edits.\n";
        return;
    }

    // Standard ANASEN color sequence
    int colors[] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1};
    std::vector<TH2 *> hists;
    std::vector<TString> usedLabels;

    TFile *f = TFile::Open(rootFile, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "overlay_2d: could not open " << rootFile << "\n";
        return;
    }

    for (size_t i = 0; i < nStates; i++)
    {
        TString hName = "ETrackvsKin_assumed/" + reaction + "_ETrack_vs_EKin" + states[i] + "_" + ejtag + "_" + topology + "_" + detector;
        TH2 *h = (TH2 *)f->Get(hName);
        if (h)
        {
            TH2 *clone = (TH2 *)h->Clone(Form("h_%zu", i));
            clone->SetDirectory(0); // Detach from file
            hists.push_back(clone);
            usedLabels.push_back(labels[i]);
        }
        else
        {
            std::cerr << "Warning: Could not find " << hName << "\n";
        }
    }
    f->Close();

    if (hists.empty())
        return;

    // Apply the per-state display gate by zeroing bins outside [xlo[i],
    // xhi[i]) directly, rather than SetRangeUser -- see header comment.
    for (size_t i = 0; i < hists.size(); i++)
    {
        TAxis *ax = hists[i]->GetXaxis();
        for (int bx = 1; bx <= hists[i]->GetNbinsX(); bx++)
        {
            double xc = ax->GetBinCenter(bx);
            if (xc < xlo[i] || xc >= xhi[i])
                for (int by = 1; by <= hists[i]->GetNbinsY(); by++)
                    hists[i]->SetBinContent(bx, by, 0);
        }
        hists[i]->RebinX(rebinX);
        hists[i]->RebinY(rebinY);
        hists[i]->SetMinimum(minz);
    }

    TCanvas *c = new TCanvas("c", "", 1800, 1800);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);

    // Turn on the X and Y grid lines
    c->SetGridx(1);
    c->SetGridy(1);

    TLegend *leg = new TLegend(0.65, 0.75, 0.9, 0.9);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001); // Solid background so grid doesn't bleed through
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.03);

    for (size_t i = 0; i < hists.size(); i++)
    {
        hists[i]->SetFillColor(colors[i % 5]);
        if (i == 0)
        {
            hists[i]->SetTitle("");
            hists[i]->GetXaxis()->SetTitle(xAxisLabel);
            hists[i]->GetYaxis()->SetTitle(yAxisLabel);
            hists[i]->GetXaxis()->CenterTitle();
            hists[i]->GetYaxis()->CenterTitle();

            hists[i]->Draw("box");

            // Draw the y = x diagonal dashed line across the plot's own
            // native range (read off the histogram, not hardcoded) with
            // enough sample points to reach both corners cleanly.
            double xMin = hists[i]->GetXaxis()->GetXmin();
            double xMax = hists[i]->GetXaxis()->GetXmax();
            TF1 *diag = new TF1("diag", "x", xMin, xMax);
            diag->SetNpx(1000);
            diag->SetLineStyle(2);
            diag->SetLineColor(kRed);
            diag->Draw("same");
        }
        else
        {
            hists[i]->Draw("box same");
        }

        leg->AddEntry(hists[i], usedLabels[i].Data(), "f"); // "f" -- box uses fill swatches
    }
    leg->Draw();

    c->SaveAs("scratch/kinematic_states_overlay.png");
    std::cout << "Saved: scratch/kinematic_states_overlay.png\n";
}
