/////////////////////////////////////////////////////////////////////////////
// \file    OscCalcSingleElectron.cxx 
// \brief   Source file for signle electron spectra scaling
// \author  Muve (wmu@fnal.gov)
/////////////////////////////////////////////////////////////////////////////
#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"

namespace osc
{
  //---------------------------------------------------------------------------
  OscCalcSingleElectron::OscCalcSingleElectron()
    : fSigScale(1.0),
      fIBkgScale(1.0),
      fBkgScale(1.0),
      fMECScale(1.0),
      fDMMass(10),
      fIsDM(true),
      fDMEnergyBinw(0.3),
      fDMFile("/exp/nova/app/users/thoroho/ldmanalysis/NuMagMomentAna/data/ldmspectra/ldmone.root")
  {
      LoadDMHists();
  }

  //---------------------------------------------------------------------------
  OscCalcSingleElectron::~OscCalcSingleElectron()
  {
  }

//  //---------------------------------------------------------------------------
//  IOscCalcAdjustable* OscCalcSingleElectron::Copy() const
//  {
//      return new OscCalcSingleElectron(*this);
//  }

  //---------------------------------------------------------------------------
  double OscCalcSingleElectron::P(int spectralable, int dmmass, double E)
  {
//      std::cout << "Spectra: " << spectralable << " mass: " << dmmass  << " E: " << E << "\n";
      if (spectralable == 0)
      {
          if (fIsDM == true)
          {
              int bin    = ceil(E/fDMEnergyBinw);        //E: true electron energy
//              std::cout << "Bin " << bin << "\n";
//              std::cout << "Scale " << fMapHists[dmmass]->GetBinContent(bin) << " fSigScale " << fSigScale << std::endl;
              return fSigScale * fMapHists[dmmass]->GetBinContent(bin);
          }
          else
          {
              return fSigScale; //do something if it is not a simple overall scaling factor.
          }
      }

      if (spectralable == 1)
      {
//          std::cout << "IBkgScale " << fIBkgScale << std::endl;
          return fIBkgScale;
      }

      if (spectralable == 2)
      {
//          std::cout << "BkgScale " << fIBkgScale << std::endl;
          return fBkgScale;
      }

      if (spectralable == 3)
      {
	  return fMECScale;
      }

      
      std::cout << "Wrong Scale " << fSigScale << std::endl;
      return fSigScale;
  }

  IOscCalcAdjustable* OscCalcSingleElectron::Copy() const
  {
      return new OscCalcSingleElectron(*this);
  }

//  //---------------------------------------------------------------------------
//  TMD5* OscCalcSingleElectron::GetParamsHash() const
//  {
//      TMD5* ret = new TMD5;
//      std::string txt = "SingleElectronAnalysis";
//      ret->Update((unsigned char*)txt.c_str(), txt.size());
//      std::vector<double> buf = GetState();
//      ret->Update((unsigned char*)&buf[0], sizeof(double)*buf.size());
//      ret->Final();
//      return ret;
//  }
  
  //---------------------------------------------------------------------------
  void OscCalcSingleElectron::LoadDMHists()
  {
      if (!fIsDM) return;
      if (!fMapHists.empty()) return;
      
      TFile fin(fDMFile.c_str(), "read");
      if (fin.IsZombie())
      {
          std::cout << "Warning: couldn't open " << fDMFile << ".  Crashing" <<  std::endl;
          abort();
      }
  
      TH1F *pmag_el = (TH1F*)fin.Get("pmag_el");
      fDMEnergyBinw = pmag_el->GetBinWidth(1);
      
      for (int i = 0; i < (int)fDMMassPoints.size(); i++)
      {
          TString hName = TString::Format("norm_%d", fDMMassPoints[i]);
          TH1F* htmp = (TH1F*)fin.Get(hName);
          if(!htmp)
          {
            std::cout << "Error: can't find " << hName << " histogram in file " << fDMFile << ".  Existing..." << std::endl;
            abort();
          }
          htmp->SetDirectory(0);
          fMapHists.insert(std::pair < int, TH1F* > (fDMMassPoints[i], htmp));
      }

      std::cout << "DMHist Loaded!" <<  std::endl;
  }
} // namespace
