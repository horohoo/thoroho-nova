#pragma once

///
/// Cuts used for the EventCVN (Wenjie's and Yiwen's) nu-on-e analysis
/// More details can be found at doc-db:55336
///

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Vars/Vars.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "NuMagMomentAna/Vars/NuoneVars.h"
#include "NuMagMomentAna/Cuts/NuoneCuts.h"
#include "NDAna/Classifiers/MuonID.h"
#include "NDAna/Classifiers/NueID.h"
#include "TVector3.h"

namespace nuone
{
  // Fiducial regions
  const TVector3 kEventCVNFiducialMin(-185, -175, 95);
  const TVector3 kEventCVNFiducialMax(175, 175, 1095);

  extern const ana::NuTruthCut kEventCVNFiducialCut_NT;

  // Signal definitions
  // TO DO: try different signal definitions
  const ana::Cut kEventCVNSig = ana::CutFromNuTruthCut(kEventCVNFiducialCut_NT) &&
                                kIsNuone;
  const ana::Cut kEventCVNBkg = !ana::CutFromNuTruthCut(kEventCVNFiducialCut_NT) ||
                                kIsNotNuone;

  // TO DO: try adding the ND decaf cuts to see the difference between using different definitions

  // remove obvious numu cc interactions
  const ana::Cut kLongestProngCut = ana::kLongestProng < 800 &&
                                    ana::kLongestProng >= 0;
  const ana::Cut kNPlaneCut       = nuone::kNPlane < 120;
  const ana::Cut kNCellCut        = nuone::kNCell < 600;

  // fiducial cut
  extern const ana::Cut kEventCVNFiducialCut;

  // containment cut
  // suppress backgrounds induced by neutrino interactions in the rock
  // (mostly upstream of the ND)
  const ana::Cut kminx = nuone::kMinX > -190;
  const ana::Cut kmaxx = nuone::kMaxX < 180;
  const ana::Cut kminy = nuone::kMinY > -180;
  const ana::Cut kmaxy = nuone::kMaxY < 190;
  const ana::Cut kminz = nuone::kMinZ > 105;
  const ana::Cut kmaxz = nuone::kMaxZ < 1275;
  const ana::Cut kEventCVNContainmentCut = kminx && kmaxx &&
                                           kminy && kmaxy &&
                                           kminz && kmaxz;

  // single particle requirement
  // the topology of signal event requires one EM shower with no other
  // particles in the final state
  const ana::Cut kShwEFracCut = nuone::kShwEFrac > 0.8;
  const ana::Cut kVtxECut     = nuone::kEvtx < 0.02;
  const ana::Cut kGapCut      = nuone::kGap < 20;

  // shower energy
  // exclude low energy events which are hard to be distinguished
  const ana::Cut kEventCVNShowerECut = nuone::kShowerCalE > 0.5 &&
                                       nuone::kShowerCalE < 5.0;

  // nu-on-e classifiers
  // separating signals from substantial backgrounds using CNN technique
  const ana::Cut kNuoneIDCut = nuone::kNuoneID > 0.73;
  const ana::Cut kEpi0IDCut  = nuone::kEpi0ID  > 0.92;

  // Etheta2
  // the scattered electron is very forward going
  const ana::Cut kEventCVNETheta2Cut = nuone::kETheta2 < 0.005;
  // selecting the sideband region
  const ana::Cut kEventCVNETheta2Cut_sb = nuone::kETheta2 >= 0.005 &&
                                          nuone::kETheta2 < 0.04;

  const ana::Cut kPreselection =
    kLongestProngCut && kNPlaneCut && kNCellCut && kEventCVNFiducialCut &&
    kEventCVNContainmentCut && kShwEFracCut && kVtxECut && kGapCut && kEventCVNShowerECut;

  const ana::Cut kEventCVNSigPreselection    = kEventCVNSig && kPreselection;
  const ana::Cut kEventCVNBkgPreselection    = kEventCVNBkg && kPreselection;

  const int kEventCVNNSels = 18;
  const int kEventCVNNChns = 7;

  const SelDef EventCVNCutFlow[kEventCVNNSels] = {
    {"slicing",      ana::kNoCut},
    {"longestProng", kLongestProngCut},
    {"nplane",       kLongestProngCut && kNPlaneCut},
    {"ncell",        kLongestProngCut && kNPlaneCut && kNCellCut},
    {"fiducial",     kLongestProngCut && kNPlaneCut && kNCellCut && kEventCVNFiducialCut},
    {"containment",  kLongestProngCut && kNPlaneCut && kNCellCut && kEventCVNFiducialCut &&
                     kEventCVNContainmentCut},
    {"showerEFrac",  kLongestProngCut && kNPlaneCut && kNCellCut && kEventCVNFiducialCut &&
                     kEventCVNContainmentCut && kShwEFracCut},
    {"vtxE",         kLongestProngCut && kNPlaneCut && kNCellCut && kEventCVNFiducialCut &&
                     kEventCVNContainmentCut && kShwEFracCut && kVtxECut},
    {"gap",          kLongestProngCut && kNPlaneCut && kNCellCut && kEventCVNFiducialCut &&
                     kEventCVNContainmentCut && kShwEFracCut && kVtxECut && kGapCut},
    {"showerE",      kPreselection},
    {"nuoneid",      kPreselection && kNuoneIDCut},
    {"epi0id",       kPreselection && kNuoneIDCut && kEpi0IDCut},
    {"etheta2",      kPreselection && kNuoneIDCut && kEpi0IDCut && kEventCVNETheta2Cut},
    {"nuoneid_excep",kPreselection && kEpi0IDCut && kEventCVNETheta2Cut},
    {"epi0id_excep", kPreselection && kNuoneIDCut && kEventCVNETheta2Cut},
    {"pid_excep",    kPreselection && kEventCVNETheta2Cut},
    {"pid_excep_sb", kPreselection && kEventCVNETheta2Cut_sb},
    {"etheta2_sb",   kPreselection && kNuoneIDCut && kEpi0IDCut && kEventCVNETheta2Cut_sb},
  };

  const SelDef EventCVNChannels[kEventCVNNChns]  = {
    {"all",        ana::kNoCut},
    {"nuonefid",   kEventCVNSig},
    {"bkg",        kEventCVNBkg},
    {"bkg_nuecc",  kEventCVNBkg && ana::kIsBeamNue},
    {"bkg_numucc", kEventCVNBkg && ana::kIsNumuCC},
    {"bkg_nc",     kEventCVNBkg && ana::kIsNC},
    {"bkg_other",  kEventCVNBkg && !ana::kIsBeamNue && !ana::kIsNumuCC && !ana::kIsNC},
  };
} // namespace
