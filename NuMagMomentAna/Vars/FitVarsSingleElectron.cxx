/////////////////////////////////////////////////////////////////////////////
// \file    FitVarSingleElectron.cxx
// \brief   
// \author  Muve (wmu@fnal.gov)
/////////////////////////////////////////////////////////////////////////////

#include "NuMagMomentAna/Vars/FitVarsSingleElectron.h"


namespace ana
{
  // declared as 'extern' in FitVarsSingleElectron.cxx
  const FitDMMassSingleElectron      kFitDMMassSingleElectron;
  const FitSigScalingSingleElectron  kFitSigScalingSingleElectron;
  const FitIBkgScalingSingleElectron kFitIBkgScalingSingleElectron;
  const FitBkgScalingSingleElectron  kFitBkgScalingSingleElectron;
  
  //---------------------------------------------------------------------------
  double FitDMMassSingleElectron::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    const osc::OscCalcSingleElectron* selcalc = dynamic_cast<const osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return 0.0;
    }
    
    double dmmass = selcalc->GetDMMass();
    return dmmass;
  }

  //---------------------------------------------------------------------------
  void FitDMMassSingleElectron::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return;
    }
    
    selcalc->SetDMMass(val);
  }

  //---------------------------------------------------------------------------
  double FitSigScalingSingleElectron::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    const osc::OscCalcSingleElectron* selcalc = dynamic_cast<const osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return 0.0;
    }
    
    double scale = selcalc->GetSigScale();
    return scale;
  }

  //---------------------------------------------------------------------------
  void FitSigScalingSingleElectron::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return;
    }
    selcalc->SetSigScale(val);
  }

  //---------------------------------------------------------------------------
  double FitIBkgScalingSingleElectron::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    const osc::OscCalcSingleElectron* selcalc = dynamic_cast<const osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return 0.0;
    }
    
    double scale = selcalc->GetIBkgScale();
    return scale;
  }
  
  //---------------------------------------------------------------------------
  void FitIBkgScalingSingleElectron::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return;
    }
    selcalc->SetIBkgScale(val);
  }
  
  //---------------------------------------------------------------------------
  double FitBkgScalingSingleElectron::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    const osc::OscCalcSingleElectron* selcalc = dynamic_cast<const osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return 0.0;
    }
    
    double scale = selcalc->GetBkgScale();
    return scale;
  }
  
  //---------------------------------------------------------------------------
  void FitBkgScalingSingleElectron::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(osc);
    if (!selcalc)
    {
        std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
        return;
    }
    selcalc->SetBkgScale(val);
  }
} // namespace

