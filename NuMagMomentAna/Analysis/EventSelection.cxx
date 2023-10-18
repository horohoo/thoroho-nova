///
/// \brief Loading and saving spectra for event selection for the neutrino magnetic moment analysis
/// \author R. Kralik rkralik@fnal.gov
///
// To run interactively: cafe -bq -l NFiles --numagmomentana -ss DemoSpectrum.cxx
//
// To run on the grid: submit_cafana.py -n NJobs --numagmomentana -ss
//                                      -r development -o /path/to/outDir
//                                      EventSelection.cxx
//

// CAFAna includes
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/SpillCuts.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"

// NuMM includes
#include "NuMagMomentAna/Vars/NuoneVars.h"
#include "NuMagMomentAna/Vars/NuoneWeights.h"
#include "NuMagMomentAna/Cuts/NuoneCuts.h"
#include "NuMagMomentAna/Cuts/SinglesCVNCuts.h"
#include "NuMagMomentAna/Cuts/EventCVNCuts.h"
#include "NuMagMomentAna/Vars/NuoneHistAxis.h"

// ROOT includes
#include "TFile.h"

void EventSelection(){
  // Definitions used:

  /// full ND CAF definition used by the eventCVN (Wenjie and Yiwen) analysis
  const std::string nd_caf = "prod_caf_R20-11-25-prod5.1reco.a_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1";

  /// flat sum decaf ND definition used by the singlesCVN (Athula and Barnali)
  /// analysis.
  /// includes following cuts: 1) true neutrino vertex is contained
  ///                          2) contains at least 1 kalman track or fuxxyK prong
  ///                          3) at least 20 hits in each slice
  ///                          4) at least 4 continuous planes with hits
  ///                          5) ND containment cuts
  const std::string nd_decaf = "prod_flatsumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1";

  /// special neutrino-on-electron enhanced overlay sample produced by Wenjie Wu
  /// (docdb:55189)
  /// to be used only for the signal (nu-on-e) events as it doesn't contain many
  /// background events!
  const std::string nuone_overlay = "prod_caf_R20-11-25-prod5.1reco.g_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_v1_nuone_overlay";

  ana::SpectrumLoader loader(nd_decaf);
  //SpectrumLoader loader(nuone_overlay);
  loader.SetSpillCut(ana::kStandardSpillCuts);

  const ana::Weight kReweight = ana::kXSecCVWgt2020 * ana::kPPFXFluxCVWgt;

  ana::Spectrum *spectrumSinglesCVN[nuone::kSinglesCVNNAxes][nuone::kSinglesCVNNChns][nuone::kSinglesCVNNSels];
  ana::Spectrum *spectrumSinglesCVNNuMM[nuone::kSinglesCVNNAxes][nuone::kSinglesCVNNChns][nuone::kSinglesCVNNSels];

  // creating spectra for the singlesCVN analysis event selection replication
  for(Int_t iAxis = 0; iAxis < nuone::kSinglesCVNNAxes; ++iAxis){
    for(Int_t iSig = 0; iSig < nuone::kSinglesCVNNChns; ++iSig){
      for(Int_t iSel = 0; iSel < nuone::kSinglesCVNNSels; ++iSel){
        spectrumSinglesCVN[iAxis][iSig][iSel] = new ana::Spectrum(loader,
                                                                  nuone::SinglesCVNAxis[iAxis].axis, 
                                                                  nuone::SinglesCVNChannels[iSig].cut &&
                                                                  nuone::SinglesCVNCutFlow[iSel].cut, 
                                                                  ana::kNoShift, kReweight);
        spectrumSinglesCVNNuMM[iAxis][iSig][iSel] = new ana::Spectrum(loader,
                                                                      nuone::SinglesCVNAxis[iAxis].axis,
                                                                      nuone::SinglesCVNChannels[iSig].cut &&
                                                                      nuone::SinglesCVNCutFlow[iSel].cut,
                                                                      ana::kNoShift, kReweight*nuone::kNuMM);
      }
    }
  }

  // Populate spectra from MC loaders
  loader.Go();

  TFile* fOut = TFile::Open("EventSelection_fullCutFlow_noSysts.root", "recreate");
  for(Int_t iAxis = 0; iAxis < nuone::kSinglesCVNNAxes; ++iAxis){
    TDirectory* topDir = fOut->mkdir(nuone::SinglesCVNAxis[iAxis].name.c_str());
    topDir->cd();
    for(Int_t iSig = 0; iSig < nuone::kSinglesCVNNChns; ++iSig){
      TDirectory* Dir = topDir->mkdir(nuone::SinglesCVNChannels[iSig].name.c_str());
      Dir->cd();
      for(Int_t iSel = 0; iSel < nuone::kSinglesCVNNSels; ++iSel){
        spectrumSinglesCVN[iAxis][iSig][iSel]->SaveTo(Dir, nuone::SinglesCVNCutFlow[iSel].name.c_str());
        spectrumSinglesCVNNuMM[iAxis][iSig][iSel]->SaveTo(Dir, (nuone::SinglesCVNCutFlow[iSel].name+"NuMM").c_str());
      }
    }
  }

  for(Int_t iAxis = 0; iAxis < nuone::kSinglesCVNNAxes; ++iAxis){
    for(Int_t iSig = 0; iSig < nuone::kSinglesCVNNChns; ++iSig){
      for(Int_t iSel = 0; iSel < nuone::kSinglesCVNNSels; ++iSel){
        delete spectrumSinglesCVN[iAxis][iSig][iSel];
        delete spectrumSinglesCVNNuMM[iAxis][iSig][iSel];
      }
    }
  }
}
