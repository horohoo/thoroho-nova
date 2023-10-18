///
/// \brief Demo script to create and plot simple spectrum for the neutrino magnetic moment analysis
/// \author R. Kralik rkralik@fnal.gov
///
// To run interactively: cafe -bq -l NFiles --numagmomentana -ss DemoSpectrum.cxx
//
// To run on the grid: submit_cafana.py -n NJobs --numagmomentana -ss
//                                      -r development -o /path/to/outDir
//                                      DemoSpectrum.cxx
//

// CAFAna includes
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"
#include "CAFAna/Analysis/Exposures.h"

// NuMM includes
#include "NuMagMomentAna/Vars/NuoneVars.h"
#include "NuMagMomentAna/Vars/NuoneWeights.h"
#include "NuMagMomentAna/Cuts/NuoneCuts.h"
#include "NuMagMomentAna/Cuts/SinglesCVNCuts.h"
#include "NuMagMomentAna/Cuts/EventCVNCuts.h"
#include "NuMagMomentAna/Vars/NuoneHistAxis.h"

// ROOT includes
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"

void DemoSpectrum(){
  // Definitions used:
  /// flat sum decaf ND definition used by the singlesCVN (Athula and Barnali)
  /// analysis.
  /// includes following cuts: 1) true neutrino vertex is contained
  ///                          2) contains at least 1 kalman track or fuxxyK prong
  ///                          3) at least 20 hits in each slice
  ///                          4) at least 4 continuous planes with hits
  ///                          5) ND containment cuts
  const std::string nd_decaf = "prod_flatsumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1";

  ana::SpectrumLoader NDLoader(nd_decaf);
  NDLoader.SetSpillCut(ana::kStandardSpillCuts);

  /// Nuone signal events for full EventCVN event selection using enhanced sample
  ana::Spectrum sFullEventCVNSig(NDLoader,
                                 nuone::RecoETh2Axis,
                                 nuone::kEventCVNSig && nuone::EventCVNCutFlow[12].cut,
                                 ana::kNoShift, ana::kXSecCVWgt2020 * ana::kPPFXFluxCVWgt);

  /// Nuone background events for full EventCVN event selection
  ana::Spectrum sFullEventCVNBkg(NDLoader,
                                 nuone::RecoETh2Axis,
                                 nuone::kEventCVNBkg && nuone::EventCVNCutFlow[12].cut,
                                 ana::kNoShift, ana::kXSecCVWgt2020 * ana::kPPFXFluxCVWgt);

  // Populate spectra from MC loaders
  NDLoader.Go();

  TFile* fOut = TFile::Open("DemoSpectrumOutput.root", "recreate");
  sFullEventCVNSig.SaveTo(fOut,"sFullEventCVNSig");
  sFullEventCVNBkg.SaveTo(fOut,"sFullEventCVNBkg");

  // Take out the histograms from the spectra
  TH1D* hFullEventCVNSig = sFullEventCVNSig.ToTH1(ana::kProd5p1NDFHCPOT);
  TH1D* hFullEventCVNBkg = sFullEventCVNBkg.ToTH1(ana::kProd5p1NDFHCPOT);

  // Make a comparison plot
  TCanvas c("c","c",800,600);
  hFullEventCVNSig->SetTitle("Signal-background comparison for EventCVN nu-on-e selection");
  hFullEventCVNSig->SetLineColor(kRed+1);
  hFullEventCVNSig->Draw("hist");
  hFullEventCVNBkg->SetLineColor(kBlue-4);
  hFullEventCVNBkg->Draw("hist same");

  TLegend l(0.65,0.75,0.89,0.89);
  l.AddEntry(hFullEventCVNSig,"EventCVN signal","l");
  l.AddEntry(hFullEventCVNBkg,"EventCVN background","l");
  l.Draw("same");

  c.SaveAs("DemoSpectrumPlot.pdf");
  c.Clear();

  delete hFullEventCVNSig;
  delete hFullEventCVNBkg;
}
