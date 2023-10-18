
/////////////////////////////////////////////////////////////////////////////
// \file    NDSingleElectronSpectra.cxx
// \brief   
// \author  Muve (wmu@fnal.gov)
/////////////////////////////////////////////////////////////////////////////

#include "NuMagMomentAna/Prediction/NDSingleElectronSpectra.h"

using namespace nuone;

namespace ana
{

  //----------------------------------------------------------------------
  NDSingleElectronSpectra:: NDSingleElectronSpectra(SpectrumLoaderBase& loader,
                         const HistAxis& recoaxis,
                         const HistAxis& trueaxis,
                         const Cut& cut,
                         const SystShifts& shift,
                         const Weight& wei)
    : ReweightableSpectrum(loader, recoaxis, trueaxis, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  NDSingleElectronSpectra::~NDSingleElectronSpectra()
  {
  }

  //----------------------------------------------------------------------
  NDSingleElectronSpectra::NDSingleElectronSpectra(const NDSingleElectronSpectra& rhs)
    : ReweightableSpectrum(rhs)
  {
      assert( rhs.fReferences.empty() ); // Copying with pending loads is unexpected
  }


  //----------------------------------------------------------------------
  NDSingleElectronSpectra::NDSingleElectronSpectra(NDSingleElectronSpectra&& rhs)
    : ReweightableSpectrum(rhs)
  {
      assert( rhs.fReferences.empty() ); // Copying with pending loads is unexpected
  }

  //----------------------------------------------------------------------
  NDSingleElectronSpectra& NDSingleElectronSpectra::operator=(const NDSingleElectronSpectra& rhs)
  {
      if(this == &rhs) return *this;
  
      ReweightableSpectrum::operator=(rhs);
  
      assert( rhs.fReferences.empty() ); // Copying with pending loads is unexpected
      assert( fReferences.empty() ); // Copying with pending loads is unexpected
  
      return *this;
  }

  //----------------------------------------------------------------------
  NDSingleElectronSpectra& NDSingleElectronSpectra::operator=(NDSingleElectronSpectra&& rhs)
  {
      if(this == &rhs) return *this;
  
      ReweightableSpectrum::operator=(rhs);
  
      assert( rhs.fReferences.empty() ); // Copying with pending loads is unexpected
      assert( fReferences.empty() ); // Copying with pending loads is unexpected
  
      return *this;
  }

  //----------------------------------------------------------------------
  Spectrum NDSingleElectronSpectra::Reweight(osc::IOscCalc* calc)
  {
      osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(calc);
//      std::cout << "NDSingleElectronSpectra::Reweight... \n";
      
//      int spectralable = selcalc->GetSpectra();

//      std::cout << "spectralable: " << spectralable << "\n";

//      if (spectralable == 0)
//      {
//          if (selcalc->GetAna() == true)
//          {
//              NDSingleElectronWeight weight(selcalc);
//              std::cout << "Reweight LDM \n";
//              return WeightedBy(weight);
//          }
//          else
//          {
//              //Scale(selcalc->P(spectralable, 0, 0)); // if it is NOT a simple overall scaling factor
//              double scale = selcalc->GetSigScale();
//              Scale(scale);
//          }
//      }
//      else if (spectralable == 1)
//      {
//          double scale = selcalc->GetIBkgScale();
//          std::cout << "Scale IBKG " << scale << std::endl;
//          Scale(scale);
//      }
//      else if (spectralable == 2)
//      {
//          double scale = selcalc->GetBkgScale();
//          std::cout << "Scale BKG " << scale << std::endl;
//          Scale(scale); 
//      }
//      std::cout <<"Signal scale: " << selcalc->GetSigScale() << std::endl;
//      std::cout <<"IBkg scale: " << selcalc->GetIBkgScale() << std::endl;
//      std::cout <<"Bkg scale: " << selcalc->GetBkgScale() << std::endl;

      NDSingleElectronWeight weight(selcalc);

//      const Eigen::VectorXd& vec = weight.GetEigen();
//      return Spectrum(Eigen::ArrayXd(vec.transpose() * fMat), fAxisX, fPOT, fLivetime);
//      
      return WeightedBy(weight);
//      return WeightingVariable();
  }

  //----------------------------------------------------------------------
  void NDSingleElectronSpectra::SaveTo(TDirectory* dir, const std::string& name) const
  {
      _SaveTo(dir, name, "NDSingleElectronSpectra");
  }

  //----------------------------------------------------------------------
  std::unique_ptr<NDSingleElectronSpectra> NDSingleElectronSpectra::LoadFrom(TDirectory* dir, const std::string& name)
  {
//      std::cout << "NDSingleElectronSpectra Loading " << name << "...\n";
      dir = dir->GetDirectory(name.c_str()); // switch to subdir
      assert(dir);
      
      DontAddDirectory guard;
      
      TObjString* tag = (TObjString*)dir->Get("type");
      assert(tag);
      assert(tag->GetString() == "NDSingleElectronSpectra");
      delete tag;
      
      TH2D* spect = (TH2D*)dir->Get("hist");
      assert(spect);
      TH1* hPot = (TH1*)dir->Get("pot");
      assert(hPot);
      TH1* hLivetime = (TH1*)dir->Get("livetime");
      assert(hLivetime);
      
      std::vector<std::string> labelsX;
      std::vector<Binning> binsX;

      std::vector<std::string> labelsY;
      std::vector<Binning> binsY;
      
      for(int i = 0; ; ++i)
      {
          const std::string nameX = TString::Format("bins%d", i).Data();
          const std::string nameY = TString::Format("binsy%d", i).Data();
          
          TDirectory* subdir = dir->GetDirectory(nameX.c_str());
          if(!subdir) break;
          delete subdir;
          
          binsX.push_back(*Binning::LoadFrom(dir, nameX));
          binsY.push_back(*Binning::LoadFrom(dir, nameY));
          
          TObjString* labelX = (TObjString*)dir->Get(TString::Format("label%d", i));
          TObjString* labelY = (TObjString*)dir->Get(TString::Format("labely%d", i));
          
          labelsX.push_back(labelX ? labelX->GetString().Data() : "");
          labelsY.push_back(labelY ? labelY->GetString().Data() : "");
//          std::cout << "labelX->GetString().Data(): " << labelX->GetString().Data() << std::endl;
          delete labelX;
          delete labelY;
      }
      
      delete dir;
      
      auto ret = std::make_unique<NDSingleElectronSpectra>(kNullLoader, HistAxis(labelsX, binsX), HistAxis(labelsY, binsY));
      
      // ROOT histogram storage is row-major, but Eigen is column-major by default
      typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen:: Dynamic, Eigen::RowMajor> MatRowMajor;
      ret->fMat = Eigen::Map<MatRowMajor>(spect->GetArray(), ret->fMat.rows(), ret->fMat.cols());
      
      delete spect;
      ret->AjustPOT(hPot->Integral(0, -1));
//      std::cout << name << ", pot: " << ret->GetPOT() << std::endl;
      ret->AjustLiveTime(hLivetime->Integral(0, -1));
      
      delete hPot;
      delete hLivetime;
      
      return ret;
  }
}
