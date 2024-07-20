#ifndef OSC_OSCCALCULATORSINGLEELECTRON_H
#define OSC_OSCCALCULATORSINGLEELECTRON_H

/////////////////////////////////////////////////////////////////////////////
// \file    OscCalcSingleElectron.h 
// \brief   Head file for signle electron spectra scaling
//          For LDM analysis, two paramters: the mass, and the scaling factor
//          For other single electron analysis: scaling factor
// \author  Muve (wmu@fnal.gov)
/////////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalc.h"
#include "TMD5.h"

#include "TH1.h"
#include "TFile.h"


#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>

//#include "NuMagMomentAna/Prediction/NDSingleElectronWeight.h"


namespace osc
{
  namespace label
  {
      enum SpectraLabel_t
      {
          kSignal                = 1<<0,
          kIrreducibleBackground = 1<<1,
          kBackground            = 1<<2,
	  kMEC                   = 1<<3
      };
  }
  
  class OscCalcSingleElectron: public IOscCalcAdjustable
  {
  public:
    using IOscCalcAdjustable::P;
    OscCalcSingleElectron();
    virtual ~OscCalcSingleElectron();

//    IOscCalcAdjustable* Copy() const override;
    double P(int spectralable, int dmmass, double E) override;
    IOscCalcAdjustable* Copy() const override;
    
    virtual void SetL(double L) override {fL = L;}
    virtual void SetRho(double rho) override {fRho = rho;}
    virtual void SetDmsq21(const double& dmsq21) override {fDmsq21 = dmsq21;}
    virtual void SetDmsq32(const double& dmsq32) override {fDmsq32 = dmsq32;}
    virtual void SetTh12(const double& th12) override {fTh12 = th12;};
    virtual void SetTh13(const double& th13) override {fTh13 = th13;};
    virtual void SetTh23(const double& th23) override {fTh23 = th23;};
    virtual void SetdCP(const double& dCP) override {fdCP = dCP;};

    void SetAna        (bool        isDM   )  {fIsDM         = isDM;   } 
    void SetDMMass     (double      mass   )  {fDMMass       = mass;   }
    void SetSigScale   (double      scale  )  {fSigScale     = scale;  }
    void SetIBkgScale  (double      scale  )  {fIBkgScale    = scale;  }
    void SetBkgScale   (double      scale  )  {fBkgScale     = scale;  }
    void SetMECScale   (double      scale  )  {fMECScale     = scale;  }
    void SetDMFile     (std::string dmfile )  {fDMFile       = dmfile; }
    void SetSpectra    (int         spectra)  {fSpectralable = spectra;}
    
    bool   GetAna       () const { return fIsDM;        }
    double GetDMMass    () const { return fDMMass;      }
    double GetSigScale  () const { return fSigScale;    }
    double GetIBkgScale () const { return fIBkgScale;   }
    double GetBkgScale  () const { return fBkgScale;    }
    int  GetSpectra     () const { return fSpectralable;}

//    TMD5* GetParamsHash() const override;
    
  protected:    
    double fSigScale;
    double fIBkgScale;
    double fBkgScale;
    double fMECScale;
    double fDMMass;
    
    int fSpectralable;  //SpectraLabel: 0-signal, 1-irreducible background (nu-on-e), 2-other background
    
    bool fIsDM;
    float fDMEnergyBinw;
    std::string fDMFile;
    std::map<int, TH1F*> fMapHists;
    std::vector<int> fDMMassPoints = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 450};  //dark matter mass in GeV
    
    void LoadDMHists();
  };
} // namespace
#endif
