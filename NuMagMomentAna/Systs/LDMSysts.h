#pragma once

#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include "NuMagMomentAna/Vars/NuoneVars.h"

#include <iostream>

#include "TFile.h"
#include "TH1.h"


namespace ana
{
  //----------------------------------------------------------------------  
  /// ND calibration, light level and cherenkov systematics:
  class SystFromHist: public ISyst
  {
  public:
    SystFromHist(const std::string &fname, const std::string &shortname, const std::string &latexname);
    ~SystFromHist();
  
    void Shift(double sigma, caf::SRProxy* sr, double &weight) const override;  
    double WeightFor(double sigma,    double energy) const;
  
  private:
    std::string fFileName;  
    void LoadHists(const std::string &shortname) const;  
    int GetOscChannel(caf::SRProxy* sr) const; // get oscillation channel for current event
    mutable std::vector< std::pair<int, TH1D*> > fHists; // Store hist for one Syst, store a pair of a sigma shift value and a histogram evaluated at that shift
  };

  extern const SystFromHist kNDldmCalibSyst;
  extern const SystFromHist kNDldmLightSyst;
  extern const SystFromHist kNDldmCherSyst;

  extern const SystFromHist kNDNuoneCalibSyst;
  extern const SystFromHist kNDNuoneLightSyst;
  extern const SystFromHist kNDNuoneCherSyst;

  extern const SystFromHist kNDNumiCalibSyst;
  extern const SystFromHist kNDNumiLightSyst;
  extern const SystFromHist kNDNumiCherSyst;

  //----------------------------------------------------------------------
  
  // The systematic to ignore pipeup effect for single electron in ND
  class NDPileupEffectSyst: public ISyst
  {
  public:
    NDPileupEffectSyst() : ISyst("pileup", "Pileup Effect") {}
  
    void Shift(double sigma, caf::SRProxy* sr, double& weight) const override;
  };
  
  extern const NDPileupEffectSyst kNDPileupEffectSyst;

  //----------------------------------------------------------------------

  // A dummy systematic (if needed)
  class DummySyst: public ISyst
  {
  public:
  DummySyst() : ISyst("dummy", "dummy syst") {}

    void Shift (double sigma, caf::SRProxy* sr, double& weight) const override;
  };

  extern const DummySyst kDummySyst;

  //---------------------------------------------------------------------

  // The systematic for LDM flux
  class LDMFluxSyst: public ISyst
  {
  public:
    LDMFluxSyst() : ISyst("flux", "LDM Flux") {}

    void Shift (double sigma, caf::SRProxy* sr, double& weight) const override;
  };

  extern const LDMFluxSyst kLDMFluxSyst;

} // namespace
