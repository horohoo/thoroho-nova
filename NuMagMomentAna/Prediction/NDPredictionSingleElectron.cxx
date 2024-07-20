#include "NuMagMomentAna/Prediction/NDPredictionSingleElectron.h"

using std::make_unique;
using std::move;
using std::string;
using std::unique_ptr;
using std::vector;

using osc::IOscCalc;

namespace ana {

  REGISTER_LOADFROM("NDPredictionSingleElectron", IPrediction, NDPredictionSingleElectron);

  //----------------------------------------------------------------------
  NDPredictionSingleElectron NDPredictionSingleElectron::NDPredictionSingleElectron_c(
             SpectrumLoaderBase& signalloaders,
             SpectrumLoaderBase& ibkgloaders,
             SpectrumLoaderBase& bkgloaders,
	     SpectrumLoaderBase& mecloaders,
             const HistAxis& axis,
             const Cut& cutSignal,
             const Cut& cutIBkg,
             const Cut& cutBkg,
	     const Cut& cutMEC,
             const SystShifts& shiftSignal,
             const SystShifts& shiftIBkg,
             const SystShifts& shiftBkg,
	     const SystShifts& shiftMEC,
             const Weight& weightIBkg,
             const Weight& weightBkg,
	     const Weight& weightMEC)
  {
    NDPredictionSingleElectron extrap(signalloaders, ibkgloaders, bkgloaders, mecloaders, axis, cutSignal, cutIBkg, cutBkg, cutMEC, shiftSignal, shiftIBkg, shiftBkg, shiftMEC, weightIBkg, weightBkg, weightMEC);
      return extrap;
  }

  //----------------------------------------------------------------------
  NDPredictionSingleElectron::NDPredictionSingleElectron(vector<unique_ptr<NDSingleElectronSpectra>>& spectra)
  {
      ;
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSingleElectron::Predict(IOscCalc* calc) const
  {
      osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(calc);
      if (!selcalc)
      {
          std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;

          abort();
      }

      std::cout << "Setting signal spectra ...\n";
      selcalc->SetSpectra(0);
      Spectrum ret = fSignalSpectra->Reweight(selcalc);
//      SpectrumMirror* retmirror = (SpectrumMirror*)(&ret);
//      retmirror->fHist.ResetErrors();
      
      std::cout << "Setting irreducible background spectra ...\n";
      selcalc->SetSpectra(1);
      ret += fIrreducibleBkgSpectra->Reweight(selcalc);

      std::cout << "Setting nominal background spectra ...\n";
      selcalc->SetSpectra(2);
      ret += fBkgSpectra->Reweight(selcalc);

      std::cout << "Setting MEC spectra ...\n";
      selcalc->SetSpectra(3);
      ret += fMECSpectra->Reweight(selcalc);
      
      return ret;    
  }
 
  //----------------------------------------------------------------------
  Spectrum NDPredictionSingleElectron::PredictComponent(IOscCalc* calc,
              Flavors::Flavors_t flav,
              Current::Current_t curr,
              Sign::Sign_t sign) const
  {
      osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(calc);
      if (!selcalc)
      {
          std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
      }
      
      if (flav == Flavors::kNuEToNuE) //let's take it as signal
      {
          selcalc->SetSpectra(0);
//          std::cout << "singal scale: " << selcalc->GetSigScale() << std::endl;
          Spectrum ret = fSignalSpectra->Reweight(selcalc);
          return ret;
      }
      else if (flav == Flavors::kNuEToNuMu)  // irreducible background (nu-on-e)
      {
          selcalc->SetSpectra(1);
//          std::cout << "ibkg scale: " << selcalc->GetIBkgScale() << std::endl;
          Spectrum ret = fIrreducibleBkgSpectra->Reweight(selcalc);
          return ret;    
      }
      else if (flav == Flavors::kNuEToNuTau) // background
      {
          selcalc->SetSpectra(2);
          Spectrum ret = fBkgSpectra->Reweight(selcalc);

	  selcalc->SetSpectra(3);
	  ret += fMECSpectra->Reweight(selcalc);
          return ret;    
      }
      else if(flav == Flavors::kAll)
      {
          selcalc->SetSpectra(0);
          Spectrum ret = fSignalSpectra->Reweight(selcalc);
          
          selcalc->SetSpectra(1);
          ret += fIrreducibleBkgSpectra->Reweight(selcalc);
          
          selcalc->SetSpectra(2);
          ret += fBkgSpectra->Reweight(selcalc);
	  
	  selcalc->SetSpectra(3);
	  ret += fMECSpectra->Reweight(selcalc);
          
          return ret;
      }
      else
      {
          std::cout << "PredictComponent wrong position " << selcalc->GetSpectra() << ", Flavor: " << flav << std::endl;
          selcalc->SetSpectra(0);
          Spectrum ret = fSignalSpectra->Reweight(selcalc);
          return ret;
      }
      
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSingleElectron::FakeData(IOscCalc* calc, double POT)
  {
      return Predict(calc).FakeData(POT);
  }

  //----------------------------------------------------------------------
  NDPredictionSingleElectron::NDPredictionSingleElectron(
             SpectrumLoaderBase& signalloaders,
             SpectrumLoaderBase& ibkgloaders,
             SpectrumLoaderBase& bkgloaders,
	     SpectrumLoaderBase& mecloaders,
             const HistAxis& axis,
             const Cut& cutSignal,
             const Cut& cutIBkg,
             const Cut& cutBkg,
	     const Cut& cutMEC,
             const SystShifts& shiftSignal,
             const SystShifts& shiftIBkg,
             const SystShifts& shiftBkg,
	     const SystShifts& shiftMEC,
             const Weight& weightIBkg,
             const Weight& weightBkg,
	     const Weight& weightMEC)
   : fSignalSpectra(new NDSingleElectronSpectra(signalloaders, axis, 
                        HistAxis("Electron True Energy (GeV)", kTrueEBins, nuone::kTrueElE),
//                        HistAxis("Electron True Energy (GeV)", kTrueEBins, nuone::kTrueElectronE),
                        cutSignal, shiftSignal))
   , fIrreducibleBkgSpectra(new NDSingleElectronSpectra(ibkgloaders, axis,
                        HistAxis("Electron True Energy (GeV)", kTrueEBins, nuone::kTrueElectronE),
                        cutIBkg, shiftIBkg, weightIBkg))
   , fBkgSpectra(new NDSingleElectronSpectra(bkgloaders, axis,
                        HistAxis("Electron True Energy (GeV)", kTrueEBins, nuone::kTrueElectronE),
                        cutBkg, shiftBkg, weightBkg))
   , fMECSpectra(new NDSingleElectronSpectra(mecloaders, axis,
					     HistAxis("Electron True Energy(GeV)", kTrueEBins, nuone::kTrueElectronE),
					     cutMEC, shiftMEC, weightMEC))

  {
      ;
  }

  //----------------------------------------------------------------------
  void NDPredictionSingleElectron::AjustPOT(double pot)
  {
      std::cout << "Rest POT according to the simulated number of events in BdNMC and Novasoft...\n";
      fSignalSpectra->AjustPOT(1.0);
      double nSimu = fSignalSpectra->ToTH2(1)->Integral();
      fSignalSpectra->AjustPOT(pot*nSimu);
      std::cout << "Ajusted POT: " << fSignalSpectra->POT() << std::endl;
  }
  
  //----------------------------------------------------------------------
  void NDPredictionSingleElectron::SetPOT(double pot)
  {
      fSignalSpectra->AjustPOT(pot);
  }
  
  //----------------------------------------------------------------------
  double NDPredictionSingleElectron::GetPOT()
  {
      std::cout << "NDPredictionSingleElectron, Signal POT: " << fSignalSpectra->POT() << std::endl;
      return fSignalSpectra->POT();
  } 

  double NDPredictionSingleElectron:: GetSignalNum(IOscCalc* calc)
  {
      fSignalSpectra->AjustPOT(1.0);
      
      double nSimu = fSignalSpectra->ToTH2(1)->Integral();
      std::cout << "NDPredictionSingleElectron, simulated " << nSimu << " events before reweight.\n";
      
      osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(calc);
      if (!selcalc)
      {
          std::cout << "Input calculator was not of type OscCalcSingleElectron." << std::endl;
      }
      selcalc->SetSpectra(0);
      Spectrum spec = fSignalSpectra->Reweight(selcalc);
      nSimu = spec.ToTH1(1.0, kPOT, kBinContent)->Integral();
      
      std::cout << "NDPredictionSingleElectron, simulated " << nSimu << " events after reweight.\n";
      return nSimu;      
  }

  //----------------------------------------------------------------------
  void NDPredictionSingleElectron::SavePrediction(TDirectory* dir , const std::string& name)
  {
      TDirectory* tmp = gDirectory;
      dir = dir->mkdir(name.c_str()); // switch to subdir
      dir->cd();
  
      TObjString("NDPredictionSingleElectron").Write("type");      
      fSignalSpectra->SaveTo(dir->mkdir("signal"), "signal");
      fIrreducibleBkgSpectra->SaveTo(dir->mkdir("ibkg"), "ibkg");
      fBkgSpectra->SaveTo(dir->mkdir("bkg"), "bkg");
      fMECSpectra->SaveTo(dir->mkdir("mec"), "mec");

      dir->Write();
      delete dir;
      tmp->cd();
  }

  //----------------------------------------------------------------------
  unique_ptr<NDPredictionSingleElectron> NDPredictionSingleElectron::LoadFrom(TDirectory* dir, const std::string& name)
  {
//      std::cout << "NDPredictionSingleElectron Loading " << name << "...\n";
      dir = dir->GetDirectory(name.c_str()); // switch to subdir
      TObjString* tag = (TObjString*)dir->Get("type");
      
      assert(dir);
      assert(tag);

      unique_ptr<NDPredictionSingleElectron> ret(new NDPredictionSingleElectron);
      ret->fSignalSpectra = ana::LoadFrom<NDSingleElectronSpectra>(dir->GetDirectory("signal"), "signal");
      ret->fIrreducibleBkgSpectra = ana::LoadFrom<NDSingleElectronSpectra>(dir->GetDirectory("ibkg"), "ibkg");
      ret->fBkgSpectra = ana::LoadFrom<NDSingleElectronSpectra>(dir->GetDirectory("bkg"), "bkg");
      ret->fMECSpectra = ana::LoadFrom<NDSingleElectronSpectra>(dir->GetDirectory("mec"), "mec");

//      std::cout << "NDPredictionSingleElectron::LoadFrom, POT: " << ret->fSignalSpectra->GetPOT() <<std::endl;
      
      delete dir;
  
      return ret;
  }
}
