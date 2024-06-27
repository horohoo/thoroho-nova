#pragma once

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/LoadFromRegistry.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Decomp/IDecomp.h"
#include "CAFAna/Extrap/IExtrap.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Vars/Vars.h" // for kMeanTime
//#include "CAFAna/Weights/PPFXWeights.h"
//#include "CAFAna/Weights/XsecTunes.h"


#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjString.h"

#include "NuMagMomentAna/Vars/FitVarsSingleElectron.h"
#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"
#include "NuMagMomentAna/Prediction/NDSingleElectronSpectra.h"
//#include "NuMagMomentAna/Prediction/SpectrumMirror.h"

namespace ana
{
  /// Take the output of an extrapolation and oscillate it as required
  class NDPredictionSingleElectron: public IPrediction
  {
  public:    
    NDPredictionSingleElectron(std::vector<std::unique_ptr<NDSingleElectronSpectra>>& spectra);

    //--------------------------------------------------------------------
    static NDPredictionSingleElectron NDPredictionSingleElectron_c(SpectrumLoaderBase& signalloaders,
             SpectrumLoaderBase& ibkgloaders,
             SpectrumLoaderBase& bkgloaders,
	     SpectrumLoaderBase& mecloaders,
             const HistAxis& axis,
             const Cut& cutSignal,
             const Cut& cutIBkg,
             const Cut& cutBkg,
	     const Cut& cutMEC,
             const SystShifts& shiftSignal = kNoShift,
             const SystShifts& shiftIBkg   = kNoShift,
             const SystShifts& shiftBkg    = kNoShift,
    	     const SystShifts& shiftMEC = kNoShift,
             const Weight& weightIBkg      = kUnweighted,
	     const Weight& weightBkg       = kUnweighted,
	     const Weight& weightMEC       = kUnweighted);

    NDPredictionSingleElectron(SpectrumLoaderBase& signalloaders,
             SpectrumLoaderBase& ibkgloaders,
             SpectrumLoaderBase& bkgloaders,
	     SpectrumLoaderBase& mecloaders,
             const HistAxis& axis,
             const Cut& cutSignal,
             const Cut& cutIBkg,
             const Cut& cutBkg,
	     const Cut& cutMEC,
             const SystShifts& shiftSignal = kNoShift,
             const SystShifts& shiftIBkg   = kNoShift,
             const SystShifts& shiftBkg    = kNoShift,
	     const SystShifts& shiftMEC    = kNoShift,
             const Weight& weightIBkg      = kUnweighted,
	     const Weight& weightBkg       = kUnweighted,
	     const Weight& weightMEC       = kUnweighted);
             

    virtual ~NDPredictionSingleElectron() {};

    NDPredictionSingleElectron(const NDPredictionSingleElectron&) = delete;
    NDPredictionSingleElectron& operator=(const NDPredictionSingleElectron&) = delete;
    NDPredictionSingleElectron(NDPredictionSingleElectron&&) = default;
    NDPredictionSingleElectron& operator=(NDPredictionSingleElectron&&) = default;

    //------------------------------------------------------------------------

    // un-hide inherited method stubs so we don't get warnings from the compiler
    using IPrediction::Predict;
    using IPrediction::PredictComponent;
    using IPrediction::PredictSyst;

    virtual Spectrum Predict(osc::IOscCalc* calc) const override;
    virtual Spectrum PredictComponent(osc::IOscCalc* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override;

    virtual Spectrum FakeData(osc::IOscCalc* calc, double POT);
    void SavePrediction(TDirectory* dir , const std::string& name);
    static std::unique_ptr<NDPredictionSingleElectron> LoadFrom(TDirectory* dir , const std::string& name);
    void AjustPOT(double pot);
    void SetPOT(double pot);
    double GetPOT();
    double GetSignalNum(osc::IOscCalc* calc);


  protected:
    NDPredictionSingleElectron() {};

    std::unique_ptr<NDSingleElectronSpectra> fSignalSpectra;
    std::unique_ptr<NDSingleElectronSpectra> fIrreducibleBkgSpectra;
    std::unique_ptr<NDSingleElectronSpectra> fBkgSpectra;
    std::unique_ptr<NDSingleElectronSpectra> fMECSpectra;
    
    //---------------------------------------
    double fPOT;
//    SpectrumLoaderBase* fSpecLoad;

  }; // end of NDPredictionSingleElectron
}
