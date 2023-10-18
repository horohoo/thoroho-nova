#pragma once

///
/// Histogram axes useful for the neutrino magnetic moment analysis
///

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/HistAxis.h"
#include "NuMagMomentAna/Vars/NuoneVars.h"

namespace nuone{
  struct AxisDef{
    std::string name;
    ana::HistAxis axis;
  };

  // axis copied from Athula's (single particle CVN) event selection slides (doc-db:53583):
  const ana::HistAxis showerEAxis = ana::HistAxis("Reconstructed shower E (GeV)", ana::Binning::Simple(50,0,5), kShowerE);
  const ana::HistAxis showerCalEAxis = ana::HistAxis("Reconstructed shower calorimetric E (GeV)", ana::Binning::Simple(50,0,5), kShowerCalE);
  const ana::HistAxis trueElectronEAxis = ana::HistAxis("True electron E (GeV)", ana::Binning::Simple(50,0,5), kTrueElectronE);
  const ana::HistAxis TrueETh2Axis = ana::HistAxis("True electron E #dot #theta^{2}", ana::Binning::Simple(20,0,0.01), kTrueElectronEThetaSq);
  const ana::HistAxis RecoETh2Axis = ana::HistAxis("Reconstructed electron E #dot #theta^{2}", ana::Binning::Simple(20,0,0.01), kETheta2);
  const ana::HistAxis cvnPngAxis = ana::HistAxis("Prong CVNe", ana::Binning::Simple(21,0,1.1), kProngCVNEle);
  const ana::HistAxis cvnHadEAxis = ana::HistAxis("CVNhadE [GeV]", ana::Binning::Simple(5,0,0.012), kCVNhadE);

  // group of axes use by Athula (single particle CVN)
  const int kSinglesCVNNAxes = 7;
  const AxisDef SinglesCVNAxis[kSinglesCVNNAxes] = {
    {"shwE",     showerEAxis},
    {"shwCalE",  showerCalEAxis},
    {"trueElE",  trueElectronEAxis},
    {"TrueETh2", TrueETh2Axis},
    {"RecoETh2", RecoETh2Axis},
    {"cvnPng",   cvnPngAxis},
    {"cvnHadE",  cvnHadEAxis}
  };
}
