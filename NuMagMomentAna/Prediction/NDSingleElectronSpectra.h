/////////////////////////////////////////////////////////////////////////////
// \file    NDSingleElectronSpectra.h 
// \brief   
// \author  Muve (wmu@fnal.gov)
/////////////////////////////////////////////////////////////////////////////

#pragma once

#include "CAFAna/Core/ReweightableSpectrum.h"

#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/OscCalcFwdDeclare.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"

#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Utilities.h"
#include "CAFAna/Core/FwdDeclare.h"

#include "CAFAna/Core/Stan.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include "OscLib/IOscCalc.h"

#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"
#include "NuMagMomentAna/Prediction/NDSingleElectronWeight.h"

#include "TDirectory.h"
#include "TH2.h"
#include "TObjString.h"

#include <cassert>
#include <memory>
#include <string>
#include "TMD5.h"

class TH2;
class TH2D;


namespace ana
{
  class Binning;

  /// Spectrum with true information, allowing it to be reweighted
  class NDSingleElectronSpectra: public ReweightableSpectrum
  {
  public:
    friend class SpectrumLoaderBase;
    friend class SpectrumLoader;
    friend class NullLoader;

    NDSingleElectronSpectra(SpectrumLoaderBase& loader,
                           const HistAxis& recoaxis,
                           const HistAxis& trueaxis,
                           const Cut& cut = kNoCut,
                           const SystShifts& shift = kNoShift,
                           const Weight& wei = kUnweighted);

    /// The only valid thing to do with such a spectrum is to assign something else into it.
    static NDSingleElectronSpectra Uninitialized(){return NDSingleElectronSpectra();}

    ~NDSingleElectronSpectra();

    /// Copy constructor
    NDSingleElectronSpectra(const NDSingleElectronSpectra& rhs);
    NDSingleElectronSpectra(NDSingleElectronSpectra&& rhs);
    /// Assignment operator
    NDSingleElectronSpectra& operator=(const NDSingleElectronSpectra& rhs);
    NDSingleElectronSpectra& operator=(NDSingleElectronSpectra&& rhs);

    Spectrum Reweight(osc::IOscCalc* calc);
    void Scale(double x) {fMat *= x;}
    void AjustPOT(double pot) {fPOT = pot;}
    void AjustLiveTime(double ltime) {fLivetime = ltime;}
    double GetPOT(){return fPOT;}
    
    void SaveTo(TDirectory* dir, const std::string& name) const;
    static std::unique_ptr<NDSingleElectronSpectra> LoadFrom(TDirectory* dir, const std::string& name);

  protected:
    
    /// Constructor for Uninitialized()
    NDSingleElectronSpectra()
    {
    }    
  };
}
