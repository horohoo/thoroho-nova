#pragma once

///
/// Cuts used for ldm analysis

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Vars/Vars.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "NuMagMomentAna/Vars/NuoneVars.h"
#include "NuMagMomentAna/Cuts/NuoneCuts.h"

#include "NDAna/Classifiers/MuonID.h"
#include "NDAna/Classifiers/NueID.h"
#include "TVector3.h"

namespace ldmone
{

  struct SelDef
  {
//    std::string suffix;
    std::string label;
    float efficiency;
    const ana::Cut cut;
  };

  // Defining nu-on-e signal with interaction mode
  // For the enhanced nu-on-e sample the mode was changed to 10005
  // and only contains the neutrino-on-electron elastic scattering.
  // Interaction mode 5 is generally electron scattering interaction,
  // which can also be inelastic and produce muons, so these should
  // probably be excluded.
//  const ana::Cut kIsNuone    = ana::kMode == 10005 && nuone::kPrimPDG == 11;  ///< Nu-on-e signal definition
//  const ana::Cut kIsNotNuone = !(ana::kMode == 5     && nuone::kPrimPDG == 11);      ///< Nu-on-e background definition
  const ana::Cut kIsNuone    = nuone::kintType == 1098;
  const ana::Cut kIsNotNuone = nuone::kintType != 1098;



  //Pre-selection
  const ana::Cut kLongestProngCut = ana::kLongestProng < 800 && ana::kLongestProng >= 0;
  const ana::Cut kNPlaneCut       = nuone::kNPlane < 120;
  const ana::Cut kNCellCut        = nuone::kNCell < 600;
  
  const ana::Cut kPreselection =  kLongestProngCut && kNPlaneCut && kNCellCut;

  // fiducial cut
  const ana::Cut kminvtxX = nuone::kVtxX > -185;
  const ana::Cut kmaxvtxX = nuone::kVtxX < 175;
  const ana::Cut kminvtxY = nuone::kVtxY > -175;
  const ana::Cut kmaxvtxY = nuone::kVtxY < 175;
  const ana::Cut kminvtxZ = nuone::kVtxZ > 95;
  const ana::Cut kmaxvtxZ = nuone::kVtxZ < 1095;

  const ana::Cut kFiducialization = kminvtxX && kmaxvtxX
                                 && kminvtxY && kmaxvtxY 
                                 && kminvtxZ && kmaxvtxZ;

  // containment cut
  // suppress backgrounds induced by neutrino interactions in the rock
  // (mostly upstream of the ND)
  const ana::Cut kminx = nuone::kMinX > -190;
  const ana::Cut kmaxx = nuone::kMaxX < 180;
  const ana::Cut kminy = nuone::kMinY > -180;
  const ana::Cut kmaxy = nuone::kMaxY < 190;
  const ana::Cut kminz = nuone::kMinZ > 105;
  const ana::Cut kmaxz = nuone::kMaxZ < 1275;
  
  const ana::Cut kContainment = kminx && kmaxx
                             && kminy && kmaxy
                             && kminz && kmaxz;

  // single particle requirement
  // the topology of signal event requires one EM shower with no other particles in the final state
  const ana::Cut kShwEFracCut = nuone::kShwEFrac > 0.8;
  const ana::Cut kVtxECut     = nuone::kEvtx < 0.02;
  const ana::Cut kGapCut      = nuone::kGap < 20;

  const ana::Cut kSingleParticle = kShwEFracCut && kVtxECut && kGapCut;

  // shower energy
  // exclude low energy events which are hard to be distinguished
  const ana::Cut kPrimaryShower = nuone::kShowerCalE > 0.5 && nuone::kShowerCalE < 5.0;
  
  // nu-on-e classifiers
  // separating signals from substantial backgrounds using CNN technique
  const ana::Cut kNuoneIDCut = nuone::kNuoneID > 0.73;
  const ana::Cut kEpi0IDCut  = nuone::kEpi0ID  > 0.92;

  const ana::Cut kClassifier = kNuoneIDCut && kEpi0IDCut;
    
  const SelDef SingleElecEventCVNCutFlow[6] = {
    {"nocut",           0.5977,  ana::kNoCut},
    {"Pre-selection",   0.5975, kPreselection},
    {"Fiducialization", 0.4221, kPreselection && kFiducialization && kContainment},
    {"SingleParticle",  0.3576, kPreselection && kFiducialization && kContainment && kSingleParticle},
    {"PriEnShower",     0.2400, kPreselection && kFiducialization && kContainment && kSingleParticle && kPrimaryShower},
    {"Classifiers",     0.1611, kPreselection && kFiducialization && kContainment && kSingleParticle && kPrimaryShower && kClassifier}
    };

} // namespace
