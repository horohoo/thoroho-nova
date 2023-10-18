#pragma once

///
/// Cuts used for the SinglesCVN (Athula's and Barnali's) nu-on-e analysis
/// More details can be found at doc-db:53583
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
  const TVector3 kSinglesCVNFiducialMin(-130, -140, 150);
  const TVector3 kSinglesCVNFiducialMax(150, 140, 800);

  // True fiducial cut
  extern const ana::NuTruthCut kSinglesCVNFiducialCut_NT;

  // Signal definitions
  // TO DO: try different signal definitions
  const ana::Cut kSinglesCVNSig = ana::CutFromNuTruthCut(kSinglesCVNFiducialCut_NT) &&
                                  nuone::kIsNuone;
  const ana::Cut kSinglesCVNBkg = !ana::CutFromNuTruthCut(kSinglesCVNFiducialCut_NT) ||
                                  nuone::kIsNotNuone;

  // Quality cut
  const ana::Cut kSinglesCVNQualityCut = ana::kNHit > 20;

  // Reco fiducial and containment cut
  extern const ana::Cut kSinglesCVNFiducialCut;
  extern const ana::Cut kSinglesCVNContainmentCut;

  // Prong selection
  const ana::Cut kSingleProngCut = nuone::kNPng == 1;
  const ana::Cut kNo2DProngCut = nuone::kNPng == 1 && nuone::kNPng2D == 0;
  const ana::Cut kSinglesCVNShowerECut = nuone::kShowerCalE < 4.1;
  const ana::Cut kProng3DVtxEvol10Cut = nuone::kProng3DvertexEnergyVol10 > 0.0 &&
                                        nuone::kProng3DvertexEnergyVol10 < 0.03;

  // Particle ID
  const ana::Cut kProngCVNEleCut = nuone::kProngCVNEle > 0.89;
  const ana::Cut kMuonIDCut  = ana::muonid_classifier::kMuonID < -0.1;
  const ana::Cut k5labelEleIDCut = nuone::k5labelEleID > 0.5;
  const ana::Cut kNueIDCut = ana::nueid_classifier::kNueID > -0.05;

  // Hadronic energy cut
  const ana::Cut kCVNHadECut = nuone::kCVNhadE < 0.035;

  // Cut excluding non-forward going particles
  const ana::Cut kSinglesCVNETh2Cut = nuone::kETheta2 < 0.01;

  const int kSinglesCVNNSels = 10;
  const int kSinglesCVNNChns = 7;

  const SelDef SinglesCVNCutFlow[kSinglesCVNNSels] = {
    {"slicing",     ana::kNoCut},
    {"SingleProng", kSingleProngCut},
    {"Fiducial",    kSingleProngCut && kSinglesCVNFiducialCut},
    {"Containment", kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut},
    {"ProngCVN",    kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut &&
                    kProngCVNEleCut},
    {"CVNHadE",     kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut &&
                    kProngCVNEleCut && kCVNHadECut},
    {"ETh2",        kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut &&
                    kProngCVNEleCut && kCVNHadECut && kProng3DVtxEvol10Cut && kSinglesCVNETh2Cut},
    {"ShwCalE",     kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut &&
                    kProngCVNEleCut && kCVNHadECut && kProng3DVtxEvol10Cut && kSinglesCVNETh2Cut &&
                    kSinglesCVNShowerECut},
    {"5lableEleID", kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut &&
                    kProngCVNEleCut && kCVNHadECut && kProng3DVtxEvol10Cut && kSinglesCVNETh2Cut &&
                    kSinglesCVNShowerECut && k5labelEleIDCut},
    {"NueID",       kSingleProngCut && kSinglesCVNFiducialCut && kSinglesCVNContainmentCut &&
                    kProngCVNEleCut && kCVNHadECut && kProng3DVtxEvol10Cut && kSinglesCVNETh2Cut &&
                    kSinglesCVNShowerECut && k5labelEleIDCut && kNueIDCut}
  };

  const SelDef SinglesCVNChannels[kSinglesCVNNChns]  = {
    {"all",        ana::kNoCut},
    {"nuonefid",   kSinglesCVNSig},
    {"bkg",        kSinglesCVNBkg},
    {"bkg_nuecc",  kSinglesCVNBkg && ana::kIsBeamNue},
    {"bkg_numucc", kSinglesCVNBkg && ana::kIsNumuCC},
    {"bkg_nc",     kSinglesCVNBkg && ana::kIsNC},
    {"bkg_other",  kSinglesCVNBkg && !ana::kIsBeamNue && !ana::kIsNumuCC && !ana::kIsNC},
  };
} // namespace
