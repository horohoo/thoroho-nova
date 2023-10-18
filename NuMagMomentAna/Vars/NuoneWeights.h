#pragma once

///
/// Event weights useful for the neutrino magnetic moment analysis
///

#include "CAFAna/Core/Weight.h"
#include "StandardRecord/Proxy/FwdDeclare.h"

namespace nuone
{
  /// Neutrino magnetic moment weight written by Robert Kralik (rkralik@fnal.gov)
  /// used as a proxy for simulating the neutrino magnetic moment events on top
  /// of true neutrino-on-electron events.
  /// Based on the equation in docdb:51208 slide 9.
  /// Uses true information on electron ID and electron and neutrino energies
  extern const ana::NuTruthWeight kNuMM_NT;
  const ana::Weight kNuMM = WeightFromNuTruthWeight(kNuMM_NT, 1);
} // namespace
