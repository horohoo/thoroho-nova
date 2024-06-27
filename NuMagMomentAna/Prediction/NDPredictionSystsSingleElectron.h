#pragma once

#include "CAFAna/Prediction/IPrediction.h"

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/LoadFromRegistry.h"
#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/Registry.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoaderBase.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Core/ThreadLocal.h"
#include "CAFAna/Core/Utilities.h"

#include "CAFAna/Decomp/IDecomp.h"
#include "CAFAna/Extrap/IExtrap.h"

//#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Systs/BeamSysts.h"
#include "CAFAna/Systs/XSecSystLists.h"

#include "CAFAna/Vars/Vars.h"

#include "OscLib/IOscCalc.h"

#include "Utilities/func/MathUtil.h"
#include "Utilities/func/Stan.h"
#include "Utilities/func/StanUtils.h"

#include <malloc.h>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TObjString.h"
#include "TVectorD.h"

#include "NuMagMomentAna/Vars/FitVarsSingleElectron.h"
#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"
#include "NuMagMomentAna/Prediction/NDSingleElectronSpectra.h"
#include "NuMagMomentAna/Prediction/NDPredictionSingleElectron.h"
#include "NuMagMomentAna/Systs/LDMSysts.h"


namespace ana
{
  class Loaders;

  class NDPredictionSystsSingleElectron: public IPrediction
  {
  public:
    enum EMode_t
    {
        kCombineSigns, kSplitBySign
    };
    enum SpectraType
    {
        kSpectraSignal=0, kSpectraNuone=1, kSpectraOther=2, kNSpectraTypes=3
    };

    struct Coeffs
    {
        Coeffs(double _a, double _b, double _c, double _d)
        : a(_a), b(_b), c(_c), d(_d)
        {
        }
        double a, b, c, d;
    };

    // ISyst    What systematics we should be capable of interpolating
    // Others   From which NDPredictionSingleElectrons (IPredictions) are constructed
    NDPredictionSystsSingleElectron(SpectrumLoaderBase& signalloaders,
                                    SpectrumLoaderBase& ibkgloaders,
                                    SpectrumLoaderBase& bkgloaders,
				    SpectrumLoaderBase& mecloaders,//added in
                                    const HistAxis& axis,
                                    const Cut& cutSignal,
                                    const Cut& cutIBkg,
                                    const Cut& cutBkg,
				    const Cut& cutMEC,
                                    const Weight& weightIBkg = kUnweighted,
                                    const Weight& weightBkg = kUnweighted,
				    const Weight& weightMEC = kUnweighted);
    virtual ~NDPredictionSystsSingleElectron();

    // un-hide inherited method stubs so we don't get warnings from the compiler
    using IPrediction::Predict;
    using IPrediction::PredictComponent;
    using IPrediction::PredictSyst;
    using IPrediction::PredictComponentSyst;

    Spectrum Predict(osc::IOscCalc* calc) const override;
    Spectrum PredictSyst(osc::IOscCalc* calc, const SystShifts& shift) const override;
    Spectrum PredictComponent(osc::IOscCalc* calc,
                              Flavors::Flavors_t flav,
                              Current::Current_t curr = Current::kBoth,
                              Sign::Sign_t sign       = Sign::kBoth
                              ) const override;
    Spectrum PredictComponentSyst(osc::IOscCalc* calc,
                                  const SystShifts& shift,
                                  Flavors::Flavors_t flav,
                                  Current::Current_t curr = Current::kBoth,
                                  Sign::Sign_t sign       = Sign::kBoth
                                  ) const override;
    NDPredictionSystsSingleElectron() :  fBinning(Spectrum::Uninitialized()) {;};
    
    Spectrum ShiftedComponent(osc::IOscCalc* calc, const SystShifts& shift, Flavors::Flavors_t flav, SpectraType type) const;
    Spectrum ShiftSpectrum(const Spectrum& s, SpectraType type, const SystShifts& shift) const;
    void ShiftBins(unsigned int N, double* arr, SpectraType type, const SystShifts& shift) const;

    // Find coefficients describing the ratios from the component  and the set of shifts
    std::vector<std::vector<Coeffs>> FitComponent(osc::IOscCalc* calc, const std::vector<double>& shifts, const std::vector<NDPredictionSingleElectron*>& preds, Flavors::Flavors_t flav, const std::string& systName) const;
    std::vector<std::vector<Coeffs>> FitRatios(const std::vector<double>& shifts, const std::vector<Eigen::ArrayXd>& ratios) const;

    Spectrum GetPredictSyst(osc::IOscCalc* calc, const ISyst* syst, double shift);

    std::vector<const ISyst*> GetAllSysts() const;

    // After calling this, DebugPlots won't work fully and SaveTo won't work at all.
    void MinimizeMemory(osc::IOscCalc* calc);

    void SaveAs(osc::IOscCalc* calc, TDirectory* dir, const std::string& name);
    static std::unique_ptr<NDPredictionSystsSingleElectron> LoadFrom(TDirectory* dir, const std::string& name);    
        
    static void LoadFromBody(TDirectory* dir, NDPredictionSystsSingleElectron* ret, std::vector<const ISyst*> veto = {});
    void SetOscSeed(osc::IOscCalc* oscSeed);

    void DebugPlot(osc::IOscCalc* calc, const ISyst* syst, Flavors::Flavors_t flav) const;
    void DebugPlotColz(osc::IOscCalc* calc, const ISyst* syst, Flavors::Flavors_t flav) const;
    
    void AjustPOT(double pot);
    double GetPOT();
//    double fPOT;

  protected:
    void InitFits(osc::IOscCalc* calc) const;
    void SetDefaultNDSingleElectronSysts() const;
    
    struct ShiftedPreds
    {
        std::string systName;                            // What systematic we're interpolating
        std::vector<double> shifts;                      // Shift values sampled
        std::vector<NDPredictionSingleElectron*> preds;  // Predictions with systematic shifts from -3 to 3 sigma
        
        int nCoeffs;                                     // Faster than calling size()
        
        std::vector<std::vector<std::vector<Coeffs>>> fits;       // Indices: [type][histogram bin][shift bin]
        std::vector<std::vector<std::vector<Coeffs>>> fitsRemap;  // Indices: [type][shift bin][histogram bin].
        void FillRemaps();                                        // Fill fitsRemap from fits
        
        double Stride() const {return shifts.size() > 1 ? shifts[1]-shifts[0] : 1;}
    };
    
    std::unique_ptr<NDPredictionSingleElectron> fPredNom;              // The nominal prediction
    mutable Spectrum fBinning;                                         // Dummy spectrum to provide binning
    mutable std::unordered_map<const ISyst*, ShiftedPreds> fPreds;

    mutable std::vector<const ISyst*> fNDSignalSysts;
    mutable std::vector<const ISyst*> fNDIBkgSysts;
    mutable std::vector<const ISyst*> fNDBkgSysts;        
    mutable std::vector<const ISyst*> fNDMECSysts;
  };
}
