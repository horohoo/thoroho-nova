#pragma once


#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Binning.h"

#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"
#include "NuMagMomentAna/Vars/NuoneVars.h"

#include <cassert>
#include <map>

#include "TH1.h"
#include "TMath.h"

namespace ana
{
  const Binning kTrueEBins = Binning::Simple(500, 0, 15);
  /// Transition probability for any one channel as a function of energy
  class NDSingleElectronWeight: public Ratio
  {
  public:
    NDSingleElectronWeight(osc::IOscCalc* calc);
    virtual ~NDSingleElectronWeight();

    NDSingleElectronWeight(const NDSingleElectronWeight& rhs) = default;
    NDSingleElectronWeight& operator=(const NDSingleElectronWeight& rhs) = default;

    TH1D* ToTH1(bool title = false) const;
  };
}

