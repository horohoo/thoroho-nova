/////////////////////////////////////////////////////////////////////////////
// \file    FitVarSingleElectron.h
// \brief   
// \author  Muve (wmu@fnal.gov)
/////////////////////////////////////////////////////////////////////////////
#pragma once

#include "CAFAna/Core/IFitVar.h"
#include "OscLib/IOscCalc.h"

#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"

namespace ana
{
  //----------------------------------------------------------------------
  class FitDMMassSingleElectron: public IFitVar
  {
    public:
      FitDMMassSingleElectron() : IFitVar("dmmass", "Dark Matter Mass [GeV]") {};

      double GetValue(const osc::IOscCalcAdjustable* osc) const override;
      void SetValue(osc::IOscCalcAdjustable* osc, double val) const override;
  };
  extern const FitDMMassSingleElectron kFitDMMassSingleElectron;

  //----------------------------------------------------------------------
  class FitSigScalingSingleElectron: public IConstrainedFitVar
  {
    public:
      FitSigScalingSingleElectron() : IConstrainedFitVar("scaling_signal", "Signal") {};

      double GetValue(const osc::IOscCalcAdjustable* osc) const override;
      void SetValue(osc::IOscCalcAdjustable* osc, double val) const override;

      double LowLimit() const override {return 0.0;}
      double HighLimit() const override {return 100.0;}
  };
  extern const FitSigScalingSingleElectron kFitSigScalingSingleElectron;

  //----------------------------------------------------------------------
  class FitIBkgScalingSingleElectron: public IConstrainedFitVar
  {
    public:
      FitIBkgScalingSingleElectron() : IConstrainedFitVar("scaling_ibkg", "Irreducible Background") {};

      double GetValue(const osc::IOscCalcAdjustable* osc) const override;
      void SetValue(osc::IOscCalcAdjustable* osc, double val) const override;

      double LowLimit() const override {return 0.8;}
      double HighLimit() const override {return 1.2;}
  };
  extern const FitIBkgScalingSingleElectron kFitIBkgScalingSingleElectron;
  
  //----------------------------------------------------------------------
  class FitBkgScalingSingleElectron: public IConstrainedFitVar
  {
    public:
      FitBkgScalingSingleElectron() : IConstrainedFitVar("scaling_bkg", "Background") {};

      double GetValue(const osc::IOscCalcAdjustable* osc) const override;
      void SetValue(osc::IOscCalcAdjustable* osc, double val) const override;

      double LowLimit() const override {return 0.8;}
      double HighLimit() const override {return 1.2;}
  };
  extern const FitBkgScalingSingleElectron kFitBkgScalingSingleElectron;
} // namespace
