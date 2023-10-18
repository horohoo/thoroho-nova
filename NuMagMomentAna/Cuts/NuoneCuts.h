#pragma once

///
/// General cuts useful for the neutrino magnetic moment analysis
///

#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Vars/Vars.h"
#include "NuMagMomentAna/Vars/NuoneVars.h"

namespace nuone
{
  // Structure holding selection cut flows with names
  struct SelDef{
    std::string name;
    ana::Cut cut;
  };

  // Defining nu-on-e signal with interaction mode
  // For the enhanced nu-on-e sample the mode was changed to 10005
  // and only contains the neutrino-on-electron elastic scattering.
  // Interaction mode 5 is generally electron scattering interaction,
  // which can also be inelastic and produce muons, so these should
  // probably be excluded.
  const ana::Cut kIsNuone    = ana::kMode == 10005;  ///< Nu-on-e signal definition
  const ana::Cut kIsNotNuone = ana::kMode != 5;      ///< Nu-on-e background definition

  /// Cut on any muons in the final state (NOTE: is this enough?)
  const ana::Cut kHasElectronInFinState([](const caf::SRProxy* sr)
                                        {
                                          if(sr->mc.nnu==0) return false;
                                          int nprims = sr->mc.nu[0].prim.size();
                                          for(int iprim=0; iprim<nprims; ++iprim){
                                            if(abs(sr->mc.nu[0].prim[iprim].pdg)==11){
                                              return true;
                                            }
                                          }
                                          return false;
                                        });

  // interaction types - DEPRECIATED - use IntModes instead!
  // TO DO: use cuts from CAFAna/Cuts/TruthCuts.h
  const ana::Cut kIsCCQE =  nuone::kNuIntType == 1001;
  const ana::Cut kIsNCQE =  nuone::kNuIntType == 1002;
  const ana::Cut kIsResNCNuNeutronPi0 =  nuone::kNuIntType == 1008; ///< resonant neutral current, \f$\nu n \to \nu n \pi^0\f$
  const ana::Cut kIsCCCOH =  nuone::kNuIntType == 1097; ///< charged current coherent pion
  const ana::Cut kIsNCCOH =  nuone::kNuIntType == 1096; ///< neutral current coherent
  const ana::Cut kIsNuoneType = nuone::kNuIntType == 1098;

} // namespace
