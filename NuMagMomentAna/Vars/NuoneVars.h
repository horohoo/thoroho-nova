#pragma once

///
/// Variables useful for the neutrino magnetic moment analysis
///
/// Combined from:
/// Wenjie's and Yiwen's work (doc-db:55336):
/// https://github.com/novaexperiment/novasoft/tree/feature/wenjiewu_nuone_ana/CAFAna/Vars
/// which uses event level CVNs (eventCVN);
///
/// and Athula's and Barnali's work (doc-db:53583):
/// /nova/app/users/bbrahma/nue/test/AthulaCodes/
/// which uses single-particle trained CVNs (singlesCVN).
///
////////////////////////////////////////////////////////////////

#include "CAFAna/Core/Var.h"
#include "CAFAna/Vars/Vars.h"
#include "StandardRecord/Proxy/FwdDeclare.h"

namespace nuone
{
  // Useful truth information from rec.mc.nu
  extern const ana::Var kNuVtxX;
  extern const ana::Var kNuVtxY;
  extern const ana::Var kNuVtxZ;
  extern const ana::Var kNuIntType; ///< interaction type - depreciated
  extern const ana::Var kNuPdgorig;
  extern const ana::Var kNuTgtZ;
  extern const ana::Var kNuTgtA;
  extern const ana::Var kNuX;
  extern const ana::Var kNuY;
  extern const ana::Var kPrimPDG;
  extern const ana::Var kMotherPrimPDG;
  extern const ana::Var kNuSlcEff; ///< slice efficiency for all events
  extern const ana::Var kNuonESlcEff_OV; ///< slice efficiency for overlay nuone events
  extern const ana::Var kNuonESlcEff; ///< slice efficiency for nuone events
  
  extern const ana::Var kintType;
  extern const ana::Var kisVtxCont;

  // Reconstructed variables:
  extern const ana::Var kVtxX;
  extern const ana::Var kVtxY;
  extern const ana::Var kVtxZ;

  extern const ana::Var kNPng; ///< number of prongs
  extern const ana::Var kNPng2D; ///< number of 2D prongs

  extern const ana::Var kShowLen; ///< length of 1st shower
  extern const ana::Var kProngLen; ///< length of 1st prong

  extern const ana::Var kNPlane; ///< number of planes
  extern const ana::Var kNCell; ///< number of cells

  // Extreme positions of reconstructed showers
  // Used for containment cuts
  extern const ana::Var kMinX;
  extern const ana::Var kMaxX;
  extern const ana::Var kMinY;
  extern const ana::Var kMaxY;
  extern const ana::Var kMinZ;
  extern const ana::Var kMaxZ;

  // shower start positions
  extern const ana::Var ShowStartRecoX;
  extern const ana::Var ShowStartRecoY;
  extern const ana::Var ShowStartRecoZ;

  // single particle
  extern const ana::Var kShwEFrac; ///< fraction of energy of leading shower out of total energy of slice
  extern const ana::Var kEvtx; ///< Energy of slice in vertex region
  extern const ana::Var kGap; ///< gap from vertex to start of shower

  // energy of the primary shower
  extern const ana::Var kShowerCalE; ///< energy based on summed calibrated deposited charge [GeV]
  extern const ana::Var kShowerE; ///< reconstructed shower energy [GeV]

  // CVN variables from the SRNuonEResult class
  // Developed by Wenjie and Yiwen (doc-db:49767)
  extern const ana::Var kNuoneID; ///< likelihood of nu-on-e event
  extern const ana::Var kEpi0ID; ///< likelihood primary prong is a pi0

  // Useful variables for the nu-on-e analyses

  /// electron energy multiplied by electron angle from the beam direction squared
  /// This variable is used to effectively distinguish nu-on-e from CC events
  /// as nu-on-e events have a forward going angle dependence
  extern const ana::Var kETheta2;
  extern const ana::Var kTheta; ///< reconstructed angle between electron and beam directions
  extern const ana::Var kTrueElectronTheta; ///< true angle between electron and beam directions
  extern const ana::Var kTrueElectronE; ///< true electron energy
  /// true electron energy multiplied by electron angle squared
  const ana::Var kTrueElectronEThetaSq = kTrueElectronE *
                                    kTrueElectronTheta *
                                    kTrueElectronTheta;

  /// vertex energy (GeV) calculated using all 3D prong hits for 10cm^3 around vertex
  /// Developed by Chatura Kuruppu (doc-db:47585) as part of NDRecoVertexObj class,
  /// used by Athula (singlesCVN) for event selection
  extern const ana::Var kProng3DvertexEnergyVol10;

  // maxima of CVN PID outputs over all fuzzyK prongs
  extern const ana::Var kProngCVNEle; ///< maximum electron likelihood
  extern const ana::Var kProngCVNMuon; ///< maximum muon likelihood
  extern const ana::Var kProngCVNPhoton; ///< maximum photon likelihood
  extern const ana::Var kProngCVNPhotonPlusEle; ///< sum of maximum photon and electron likelihoods

  /// hadronic energy calculated as total calorimetric energy - electromagnetic calE
  /// uses energy based on summed calibrated deposited charge [GeV] from fuzzyK
  /// and CVN scores from SRCVNParticleResult also from fuzzyK
  extern const ana::Var kCVNhadE;

  /// single particle-trained prong CVN, Electron/Photon/Proton/Pion/Muon (doc-db:51176)
  extern const ana::Var k5labelEleID;

  // truth information from mc.cosmic for dark matter analysis
  extern const ana::Var kTrueElTheta;
  extern const ana::Var kTrueElE;
  const ana::Var kTrueElEThta2 = kTrueElE * kTrueElTheta * kTrueElTheta;
} // end of namespace
