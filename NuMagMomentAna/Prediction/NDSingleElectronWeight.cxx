#include "NuMagMomentAna/Prediction/NDSingleElectronWeight.h"


namespace ana
{
  
  //----------------------------------------------------------------------
  /// Helper for constructors
  Eigen::ArrayXd ToEigenSingleElectron(osc::IOscCalc* calc)
  {
      osc::OscCalcSingleElectron* selcalc = dynamic_cast<osc::OscCalcSingleElectron*>(calc);
            
      unsigned int N = kTrueEBins.NBins();
      
      // underflow and overflow
      Eigen::ArrayXd val(N+2);
      val[0] = 0;
      val[N+1] = 0;
      
      int mass;
      if (selcalc->GetAna() == true)
      {
          mass = static_cast<int>(selcalc->GetDMMass());
      }
      else
      {
          mass = 0;
      }
      
      for(unsigned int i = 1; i <= N; ++i)
      {
          double x0 = kTrueEBins.Edges()[i-1];
          double x1 = kTrueEBins.Edges()[i];
          double ElTrueE = (x0+x1)/2.0;
      
          if(ElTrueE > 0)
          {
              val[i] = selcalc->P(selcalc->GetSpectra(), mass, ElTrueE);
          }
      }
      return val;
  }

  //----------------------------------------------------------------------
  NDSingleElectronWeight::NDSingleElectronWeight(osc::IOscCalc* calc)
    : Ratio(Hist::Adopt(ToEigenSingleElectron(calc)),
            std::vector<Binning>(1, kTrueEBins),
            std::vector<std::string>(1, "Electron True Energy (GeV)"))
  {
      ;
  }

  
  //----------------------------------------------------------------------
  NDSingleElectronWeight::~NDSingleElectronWeight()
  {
      ;
  }

  //----------------------------------------------------------------------
  TH1D* NDSingleElectronWeight::ToTH1(bool title) const
  {
    TH1D* ret = Ratio::ToTH1();
    ret->GetYaxis()->SetTitle("Probability");

    if(title)
    {
      ret->SetTitle("Single Electron");
    }

    return ret;
  }
}

