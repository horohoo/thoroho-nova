#include "NuMagMomentAna/Prediction/NDPredictionSystsSingleElectron.h"

namespace ana
{
  REGISTER_LOADFROM("NDPredictionSystsSingleElectron", IPrediction, NDPredictionSystsSingleElectron);

  //----------------------------------------------------------------------
  NDPredictionSystsSingleElectron::NDPredictionSystsSingleElectron(SpectrumLoaderBase& signalloaders,
                                                                   SpectrumLoaderBase& ibkgloaders,
                                                                   SpectrumLoaderBase& bkgloaders,
								   SpectrumLoaderBase& mecloaders, //added in
                                                                   const HistAxis& axis,
                                                                   const Cut& cutSignal,
                                                                   const Cut& cutIBkg,
                                                                   const Cut& cutBkg,
								   const Cut& cutMEC,
                                                                   const Weight& weightIBkg,
                                                                   const Weight& weightBkg,
								   const Weight& weightMEC)
  : fBinning(Spectrum::Uninitialized())
  {
    NDPredictionSingleElectron* prednom = new NDPredictionSingleElectron(signalloaders, ibkgloaders, bkgloaders, mecloaders, 
									 axis, cutSignal, cutIBkg, cutBkg, cutMEC, 
									 kNoShift, kNoShift, kNoShift, kNoShift,
									 weightIBkg, weightBkg, weightMEC);
      fPredNom.reset(prednom);
      
      SetDefaultNDSingleElectronSysts();
      
      for(int i_syst = 0; i_syst < (int)fNDSignalSysts.size(); ++ i_syst)
      {
          ShiftedPreds sp;
          sp.systName = fNDSignalSysts[i_syst]->ShortName();
          sp.shifts = {-3, -2, -1, 0, +1, +2, +3};
          
          //std::cout << "NDPredictionSystsSingleElectron is processing " << sp.systName << "... \n";
          
          for(int sigma: sp.shifts)
          {
              SystShifts shiftSignal;
              shiftSignal.SetShift(fNDSignalSysts[i_syst], sigma);

              SystShifts shiftIBkg;
              shiftIBkg.SetShift(fNDIBkgSysts[i_syst], sigma);
              
              SystShifts shiftBkg;
              shiftBkg.SetShift(fNDBkgSysts[i_syst], sigma);

	      SystShifts shiftMEC;
	      shiftMEC.SetShift(fNDMECSysts[i_syst], sigma);

              NDPredictionSingleElectron* preds = new NDPredictionSingleElectron(signalloaders, ibkgloaders, bkgloaders, mecloaders, axis, cutSignal, cutIBkg, cutBkg, cutMEC, shiftSignal, shiftIBkg, shiftBkg, shiftMEC,  weightIBkg, weightBkg, weightMEC);
              sp.preds.emplace_back(preds);
          }
          fPreds.emplace(fNDSignalSysts[i_syst], std::move(sp));
	  if (fNDIBkgSysts[i_syst]->ShortName() != fNDSignalSysts[i_syst]->ShortName()) {
	  }
      }

      for(int i_syst = 0; i_syst < (int)fNDIBkgSysts.size(); ++ i_syst)
      {
	// loop through background systs to collect any that were dummy systs in the signal
	if (fNDIBkgSysts[i_syst]->ShortName() != fNDSignalSysts[i_syst]->ShortName())
	{
	  ShiftedPreds sp;
	  sp.systName = fNDIBkgSysts[i_syst]->ShortName();
	  sp.shifts = {-3, -2, -1, 0, +1, +2, +3};
	  
	  //std::cout << "NDPredictionSystsSingleElectron is processing " << sp.systName << "... \n";

	  for(int sigma: sp.shifts)
	    {
	      SystShifts shiftSignal;
	      shiftSignal.SetShift(fNDSignalSysts[i_syst], sigma);
	      
	      SystShifts shiftIBkg;
	      shiftIBkg.SetShift(fNDIBkgSysts[i_syst], sigma);
	      
	      SystShifts shiftBkg;
	      shiftBkg.SetShift(fNDBkgSysts[i_syst], sigma);
	      
	      SystShifts shiftMEC;
	      shiftMEC.SetShift(fNDMECSysts[i_syst], sigma);
	      
	      NDPredictionSingleElectron* preds = new NDPredictionSingleElectron(signalloaders, ibkgloaders, bkgloaders, mecloaders, axis, cutSignal, cutIBkg, cutBkg, cutMEC, shiftSignal, shiftIBkg, shiftBkg, shiftMEC,  weightIBkg, weightBkg, weightMEC);
	      sp.preds.emplace_back(preds);
	    }
	  fPreds.emplace(fNDIBkgSysts[i_syst], std::move(sp));
	}
      }

      std::cout << "NDPredictionSystsSingleElectron constructed. \n";
  }

  //----------------------------------------------------------------------
  NDPredictionSystsSingleElectron::~NDPredictionSystsSingleElectron()
  {
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::Predict(osc::IOscCalc* calc) const
  {
      return fPredNom->Predict(calc);
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::PredictComponent(osc::IOscCalc* calc, Flavors::Flavors_t flav, Current::Current_t curr, Sign::Sign_t sign) const
  {
      return fPredNom->PredictComponent(calc, flav, curr, sign);
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::PredictSyst(osc::IOscCalc* calc, const SystShifts& shift) const
  {
      InitFits(calc);
      return PredictComponentSyst(calc, shift, Flavors::kAll, Current::kBoth, Sign::kBoth);
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::PredictComponentSyst(osc::IOscCalc* calc, const SystShifts& shift, Flavors::Flavors_t flav, Current::Current_t curr, Sign::Sign_t sign) const
  {
      InitFits(calc);
      Spectrum ret = fBinning;
      ret.Clear();
      
      // Check that we're able to handle all the systs we were passed
      for(const ISyst* syst: shift.ActiveSysts())
      {
          if(fPreds.find(syst) == fPreds.end())
          {
              std::cout << "NDPredictionSystsSingleElectron is not set up to handle the systematic: " << syst->ShortName() << std::endl;
              abort();
          }
      }
      
      if (flav & Flavors::kNuEToNuE) //let's take it as signal
      {
          ret += ShiftedComponent(calc, shift, flav, kSpectraSignal);
      }
      else if (flav & Flavors::kNuEToNuMu)  // irreducible background (nu-on-e)
      {
          ret += ShiftedComponent(calc, shift, flav, kSpectraNuone);
      }
      else if (flav & Flavors::kNuEToNuTau) // background
      {
          ret += ShiftedComponent(calc, shift, flav, kSpectraOther);
      }
      else if(flav == Flavors::kAll)
      {
          ret += ShiftedComponent(calc, shift, Flavors::kNuEToNuE,   kSpectraSignal);
          ret += ShiftedComponent(calc, shift, Flavors::kNuEToNuMu,  kSpectraNuone);
          ret += ShiftedComponent(calc, shift, Flavors::kNuEToNuTau, kSpectraOther);
      }

      return ret;
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::SetDefaultNDSingleElectronSysts() const
  {
      std::vector<const ISyst*> XSectSys = getAllXsecSysts_2020_GSF();
      
      fNDSignalSysts.push_back(&kNDPileupEffectSyst);
      fNDSignalSysts.push_back(&kNDldmCalibSyst);
      fNDSignalSysts.push_back(&kNDldmLightSyst);
      fNDSignalSysts.push_back(&kNDldmCherSyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kDummySyst);
      fNDSignalSysts.push_back(&kLDMFluxSyst);
      for (long unsigned int i=0; i < XSectSys.size(); i++) {
	fNDSignalSysts.push_back(&kDummySyst);
      }

      fNDIBkgSysts.push_back(&kNDPileupEffectSyst);//&kDummySyst);
      fNDIBkgSysts.push_back(&kNDNuoneCalibSyst);
      fNDIBkgSysts.push_back(&kNDNuoneLightSyst);
      fNDIBkgSysts.push_back(&kNDNuoneCherSyst);
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(0)); 
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(1)); 
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(2)); 
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(3)); 
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(4)); 
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(5)); 
      fNDIBkgSysts.push_back(GetFluxPrincipalsND2020(6)); 
      fNDIBkgSysts.push_back(&kDummySyst);
      fNDIBkgSysts.insert(fNDIBkgSysts.end(), XSectSys.begin(), XSectSys.end());

      fNDBkgSysts.push_back(&kDummySyst);
      fNDBkgSysts.push_back(&kNDNumiCalibSyst);
      fNDBkgSysts.push_back(&kNDNumiLightSyst);
      fNDBkgSysts.push_back(&kNDNumiCherSyst);
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(0)); 
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(1)); 
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(2)); 
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(3)); 
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(4)); 
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(5)); 
      fNDBkgSysts.push_back(GetFluxPrincipalsND2020(6));
      fNDBkgSysts.push_back(&kDummySyst);
      fNDBkgSysts.insert(fNDBkgSysts.end(), XSectSys.begin(), XSectSys.end());

      fNDMECSysts.push_back(&kDummySyst);
      fNDMECSysts.push_back(&kNDNumiCalibSyst);
      fNDMECSysts.push_back(&kNDNumiLightSyst);
      fNDMECSysts.push_back(&kNDNumiCherSyst);
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(0));
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(1));
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(2));
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(3));
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(4));
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(5));
      fNDMECSysts.push_back(GetFluxPrincipalsND2020(6));
      fNDMECSysts.push_back(&kDummySyst);
      fNDMECSysts.insert(fNDMECSysts.end(), XSectSys.begin(), XSectSys.end());
  }
  
  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::InitFits(osc::IOscCalc* calc) const
  {
      if(fPreds.empty()) // No systs
      {
          return;
      }
      else if(!fPreds.begin()->second.fits.empty()) // Already initialized
      {
          return;
      }
      
      for(auto& it: fPreds)
      {
          ShiftedPreds& sp = it.second;
          
          sp.fits.resize(kNSpectraTypes);
          
          sp.fits[kSpectraSignal] = FitComponent(calc, sp.shifts, sp.preds, Flavors::kNuEToNuE,   sp.systName);  //use Flavors::kNuEToNuE for signal
          sp.fits[kSpectraNuone]  = FitComponent(calc, sp.shifts, sp.preds, Flavors::kNuEToNuMu,  sp.systName);
          sp.fits[kSpectraOther]  = FitComponent(calc, sp.shifts, sp.preds, Flavors::kNuEToNuTau, sp.systName);    
	  
          sp.nCoeffs = sp.fits[0][0].size();
          
          sp.FillRemaps(); // Copy the outputs into the remapped indexing order.
      }
      // Predict something, anything, so that we can know what binning to use
      fBinning = fPredNom->Predict(calc);
      fBinning.Clear();
  }
  
  // This function is used to find the ratio between the spectra predicted with and without systematic shifts
  //----------------------------------------------------------------------
  std::vector<std::vector<NDPredictionSystsSingleElectron::Coeffs>> NDPredictionSystsSingleElectron:: FitComponent(osc::IOscCalc* calc, const std::vector<double>& shifts, const std::vector<NDPredictionSingleElectron*>& preds, Flavors::Flavors_t flav, const std::string& systName) const
  {
      NDPredictionSingleElectron* pNom = 0;
      for(unsigned int i = 0; i < shifts.size(); ++i)
      {
          if(shifts[i] == 0)
          {
              pNom = preds[i];
          }
      }
      assert(pNom);
  
      // Do it this way rather than via fPredNom so that systematics evaluated relative to some alternate nominal can work.
      const Spectrum nom = pNom->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth);
  
      std::vector<Eigen::ArrayXd> ratios;
      ratios.reserve(preds.size());
      for(auto& p: preds)
      {
          ratios.emplace_back(Ratio(p->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth), nom).GetEigen());

          
          // Make sure none of the ratio values is crazy (larger than 500)
          Eigen::ArrayXd& r = ratios.back();
          for(int i = 0; i < r.size(); ++i)
          {
              if(r[i] > 500)
              {
                  std::cout << "NDPredictionSystsSingleElectron: WARNING, ratio in bin " << i << " for " 
                            << shifts[&p-&preds.front()] << " sigma shift of "
                            << systName  << " is " << r[i] << ". Ignoring." << std::endl;
                  r[i] = 1;
              }
          }
      }
      return FitRatios(shifts, ratios);
  }

  //----------------------------------------------------------------------
  // This function is used to do a cubic interpolation. For each adjacent set of four points we determine coefficients for a cubic which will be the curve between the center two.
  // We constrain the function to match the two center points and to have the right mean gradient at them.
  // This causes this patch to match smoothly with the next one along. The resulting function is continuous, first and second differentiable.
  // At the ends of the range we fit a quadratic instead with only one constraint on the slope.
  // The coordinate conventions are that point y1 sits at x=0 and y2 at x=1. The matrices are simply the inverses of writing out the constraints expressed above.
  std::vector<std::vector<NDPredictionSystsSingleElectron::Coeffs>> NDPredictionSystsSingleElectron::FitRatios(const std::vector<double>& shifts, const std::vector<Eigen::ArrayXd>& ratios) const
  {
      if(ratios.size() < 2)
      {
          std::cout << "NDPredictionSystsSingleElectron::FitRatios(): ratios.size() = " << ratios.size() << " - how did that happen?" << std::endl;
          abort();
      }
      assert(shifts.size() == ratios.size());
      
      std::vector<std::vector<Coeffs>> ret;
      const int binMax = ratios[0].size();
      for(int binIdx = 0; binIdx < binMax; ++binIdx)
      {
          ret.push_back({});      
          // Special-case for linear interpolation
          if(ratios.size() == 2)
          {
              const double y0 = ratios[0][binIdx];
              const double y1 = ratios[1][binIdx];
              
              ret.back().emplace_back(0, 0, y1-y0, y0);
              continue;
          }
          
          {
          const double y1 = ratios[0][binIdx];
          const double y2 = ratios[1][binIdx];
          const double y3 = ratios[2][binIdx];
          const double v[3] = {y1, y2, (y3-y1)/2};
          const double m[9] = { 1, -1,  1,
                               -2,  2, -1,
                               1,  0,  0};
          const TVectorD res = TMatrixD(3, 3, m) * TVectorD(3, v);
          ret.back().emplace_back(0, res(0), res(1), res(2));
          }
          
          // We're assuming here that the shifts are separated by exactly 1 sigma.
          for(unsigned int shiftIdx = 1; shiftIdx < ratios.size()-2; ++shiftIdx)
          {
              const double y0 = ratios[shiftIdx-1][binIdx];
              const double y1 = ratios[shiftIdx  ][binIdx];
              const double y2 = ratios[shiftIdx+1][binIdx];
              const double y3 = ratios[shiftIdx+2][binIdx];
              
              const double v[4] = {y1, y2, (y2-y0)/2, (y3-y1)/2};
              const double m[16] = { 2, -2,  1,  1,
                                    -3,  3, -2, -1,
                                     0,  0,  1,  0,
                                     1,  0,  0,  0};
              const TVectorD res = TMatrixD(4, 4, m) * TVectorD(4, v);
              ret.back().emplace_back(res(0), res(1), res(2), res(3));
          }

          {
          const int N = ratios.size()-3;
          const double y0 = ratios[N  ][binIdx];
          const double y1 = ratios[N+1][binIdx];
          const double y2 = ratios[N+2][binIdx];
          const double v[3] = {y1, y2, (y2-y0)/2};
          const double m[9] = {-1,  1, -1,
                                0,  0,  1,
                                1,  0,  0};
          const TVectorD res = TMatrixD(3, 3, m) * TVectorD(3, v);
          ret.back().emplace_back(0, res(0), res(1), res(2));
          }
      }
      
      double stride = -1;
      for(unsigned int i = 0; i < shifts.size()-1; ++i)
      {
          const double newStride = shifts[i+1]-shifts[i];
          assert((stride < 0 || fabs(stride-newStride) < 1e-3) && "Variably-spaced syst templates are unsupported");
          stride = newStride;
      }
      
      if(stride != 1) // Need to rescale all the coefficients
      {
          std::cout << "NDPredictionSystsSingleElectron FitRatios stride = " << stride << std::endl;
          for(std::vector<Coeffs>& cs: ret)
          {
              for(Coeffs& c: cs)
              {
                  c = Coeffs(c.a/(stride*stride*stride), c.b/(stride*stride), c.c/stride, c.d);
              }
          }
      }
      
      return ret;
  }

  //----------------------------------------------------------------------
  // Copy fits into the remapped indexing order...
  void NDPredictionSystsSingleElectron::ShiftedPreds::FillRemaps() 
  {
      // Allocate the space
      fitsRemap.resize(fits.size());
      for(auto& it: fitsRemap)
      {
          it.resize(fits[0][0].size());
          for(auto& it2: it)
          {
              it2.resize(fits[0].size(), Coeffs(0, 0, 0, 0));
          }
      }
      
      // Copy with the transposed indices
      for(unsigned int i = 0; i < fitsRemap.size(); ++i)
      {
          for(unsigned int j = 0; j < fitsRemap[i].size(); ++j)
          {
              for(unsigned int k = 0; k < fitsRemap[i][j].size(); ++k)
              {
                  fitsRemap[i][j][k] = fits[i][k][j];
              }
          }
      }
  }
  
  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::ShiftedComponent(osc::IOscCalc* calc, const SystShifts& shift, Flavors::Flavors_t flav, SpectraType type) const
  {      
      const Spectrum nom = fPredNom->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth);
      return ShiftSpectrum(nom, type, shift);
  }

  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::ShiftSpectrum(const Spectrum& s, SpectraType type, const SystShifts& shift) const
  {
      Eigen::ArrayXd vec = s.GetEigen(s.POT());

      ShiftBins(vec.size(), vec.data(), type, shift);

      return Spectrum(std::move(vec), HistAxis(s.GetLabels(), s.GetBinnings()), s.POT(), s.Livetime());
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::ShiftBins(unsigned int N, double* arr, SpectraType type, const SystShifts& shift) const
  {
      double corr[N];
      for(unsigned int i = 0; i < N; ++i)
      {
          corr[i] = 1;
      }

      for(auto& it: fPreds)
      {
          const ISyst* syst = it.first;
          const ShiftedPreds& sp = it.second;
          
          double x = shift.GetShift<double>(syst);
          
          int shiftBin = (util::GetValAs<double>(x) - sp.shifts[0])/sp.Stride();
          shiftBin = std::max(0, shiftBin);
          shiftBin = std::min(shiftBin, sp.nCoeffs - 1);
          
          const Coeffs* fits = &sp.fitsRemap[type][shiftBin].front();
          
          x -= sp.shifts[shiftBin];
          
          const double x_cube = util::cube(x);
          const double x_sqr = util::sqr(x);
          
          for(unsigned int n = 0; n < N; ++n)
          {
              // Uncomment to debug crashes in this function
              // assert(type < fits.size());
              // assert(n < sp.fits[type].size());
              // assert(shiftBin < int(sp.fits[type][n].size()));
              const Coeffs& f = fits[n];
              corr[n] *= f.a*x_cube + f.b*x_sqr + f.c*x + f.d;
	      //if (type == 2) {
	      //if (n == 2 || n == 6 || n == 9 || n == 16) {
	      //  std::cout << n << " " << corr[n] << " " << f.a << " " << f.b << " " << f.c << " " << f.d << " " << x << std::endl;
	      //}
	      //}
          }
      }

      for(unsigned int n = 0; n < N; ++n)
      {
          // arr[n] *= (corr[n] > 0.) ? corr[n] : 0.; // std::max() doesn't work with stan::math::var.
          arr[n] *= std::max(corr[n], 0.0);
      }
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::AjustPOT(double pot)
  {      
      fPredNom->SetPOT(pot);
      for(auto& it: fPreds)
      {
          ShiftedPreds& sp = it.second;
          for(auto& preds: sp.preds)
          {
              preds->SetPOT(pot);
          }
      }
  }
  
  //----------------------------------------------------------------------
  double NDPredictionSystsSingleElectron::GetPOT()
  {      
      double POT = fPredNom->GetPOT();
      
      std::cout << "Signal POT: " << POT << std::endl;

      return POT;
  }
  
  //----------------------------------------------------------------------
  std::vector<const ISyst*> NDPredictionSystsSingleElectron::GetAllSysts() const
  {
      std::vector<const ISyst*> allsysts;
      for (const auto &p : fPreds)
      {
          allsysts.push_back(p.first);
      }
      
      return allsysts;
  } 

  //----------------------------------------------------------------------
  Spectrum NDPredictionSystsSingleElectron::GetPredictSyst(osc::IOscCalc* calc, const ISyst* syst, double shift)
  {
    if (fPreds.empty()) // No systs
      {
	std::cout << "No systs" << std::endl;
	abort();
      }
    auto it = fPreds.find(syst);
    if (it == fPreds.end())
      {
	std::cout << "NDPredictionSystsSingleElectron is not set up to handle the systematic: " << syst->ShortName() << std::endl;
	abort();
      }
    
    ShiftedPreds& sp = it->second;
    long unsigned int i = 0;
    for (; i < sp.shifts.size(); i++)
      {
	if (sp.shifts[i] == shift)
	  {
	    break;
	  }
      }
    std::cout << "Spepctrum for " << sp.shifts[i] << " sigma shift" << std::endl;
    return sp.preds[i]->Predict(calc);
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::SetOscSeed(osc::IOscCalc* oscSeed)
  {
      for(auto& it: fPreds)
      {
          it.second.fits.clear();
      }
      InitFits(oscSeed);
  }


  void NDPredictionSystsSingleElectron::SaveAs(osc::IOscCalc* calc, TDirectory* dir, const std::string& name) 
  {
      InitFits(calc);
      TDirectory* tmp = gDirectory;
      dir = dir->mkdir(name.c_str()); // switch to subdir
      dir->cd();
      TObjString("NDPredictionSystsSingleElectron").Write("type");
      
      fPredNom->SavePrediction(dir, "pred_nom");
      
      for(auto& it: fPreds)
      {
          const ShiftedPreds& sp = it.second;
          
          for(unsigned int i = 0; i < sp.shifts.size(); ++i)
          {
              if(!sp.preds[i])
              {
                  std::cout << "Can't save a NDPredictionSystsSingleElectron after MinimizeMemory()" << std::endl;
                  abort();
              }
              
              sp.preds[i]->SavePrediction(dir, TString::Format("pred_%s_%+d", sp.systName.c_str(),int(sp.shifts[i])).Data());
          }
      }
      
      if(!fPreds.empty())
      {
          TH1F hSystNames("syst_names", ";Syst names", fPreds.size(), 0, fPreds.size());
          int binIdx = 1;
          for(auto& it: fPreds)
          {
              hSystNames.GetXaxis()->SetBinLabel(binIdx++, it.second.systName.c_str());
          }
          dir->cd();
          hSystNames.Write("syst_names");
      }
      
      dir->Write();
      delete dir;
      
      tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<NDPredictionSystsSingleElectron> NDPredictionSystsSingleElectron::LoadFrom(TDirectory* dir, const std::string& name)
  {
//      std::cout << "NDPredictionSystsSingleElectron::LoadFrom: " << name << std::endl;
      dir = dir->GetDirectory(name.c_str()); // switch to subdir
      assert(dir);
      TObjString* tag = (TObjString*)dir->Get("type");
      assert(tag);
      
      std::unique_ptr<NDPredictionSystsSingleElectron> ret(new NDPredictionSystsSingleElectron);
      LoadFromBody(dir, ret.get());
//      
//      NDPredictionSystsSingleElectron* ret(new NDPredictionSystsSingleElectron);
//      LoadFromBody(dir, ret);
      
      delete dir;
      
      return ret;
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::LoadFromBody(TDirectory* dir, NDPredictionSystsSingleElectron* ret, std::vector<const ISyst*> veto)
  {
      ret->fPredNom = ana::LoadFrom<NDPredictionSingleElectron>(dir, "pred_nom");
//      std::cout << "NDPredictionSystsSingleElectron::LoadFromBody...\n";
//      ret->GetPOT();
      
      TH1* hSystNames = (TH1*)dir->Get("syst_names");
      if(hSystNames)
      {
          for(int systIdx = 0; systIdx < hSystNames->GetNbinsX(); ++systIdx)
          {
              ShiftedPreds sp;
              sp.systName = hSystNames->GetXaxis()->GetBinLabel(systIdx+1);
              
              const ISyst* syst = Registry<ISyst>::ShortNameToPtr(sp.systName, true);
              
              if(std::find(veto.begin(), veto.end(), syst) != veto.end()) continue;
              
              for(int shift = -3; shift <= +3; ++shift)
              {
                  const std::string subname = TString::Format("pred_%s_%+d", sp.systName.c_str(), shift).Data();
                  TDirectory *preddir = dir->GetDirectory(subname.c_str());
                  if(!preddir) continue; // Can happen for genie systs
                  delete preddir;
                  
                  sp.shifts.push_back(shift);
                  sp.preds.emplace_back(ana::LoadFrom<NDPredictionSingleElectron>(dir, subname).release());
              } // end for shift
              
              ret->fPreds.emplace(syst, std::move(sp));
          }
      }
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::MinimizeMemory(osc::IOscCalc* calc)
  {
      InitFits(calc);
      
      std::set<NDPredictionSingleElectron*> todel;
      for(auto& it: fPreds)
      {
          std::vector<NDPredictionSingleElectron*>& preds = it.second.preds;
          for(unsigned int i = 0; i < preds.size(); ++i)
          {
              if(preds[i] != fPredNom.get())
              {
                  todel.insert(preds[i]);
                  preds[i] = 0;
              }
          }
      }
      
      for(auto* p: todel) delete p;
      
      // We probably just freed up a lot of memory, but malloc by default hangs on to all of it as cache.
      malloc_trim(0);
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::DebugPlot(osc::IOscCalc* calc, const ISyst* syst, Flavors::Flavors_t flav) const
  {
      InitFits(calc);
      
      auto it = fPreds.find(syst);
      if(it == fPreds.end())
      {
          std::cout << "NDPredictionSystsSingleElectron::DebugPlots(): "
          << syst->ShortName() << " not found" << std::endl;
          return;
      }
      
      std::unique_ptr<TH1> nom(fPredNom->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth).ToTH1(18e20));
      const int nbins = nom->GetNbinsX();
      
      TGraph* curves[nbins];
      TGraph* points[nbins];
      
      for(int i = 0; i <= 80; ++i)
      {
          const double x = .1*i-4;
          const SystShifts ss(it->first, x);
          std::unique_ptr<TH1> h(PredictComponentSyst(calc, ss, flav, Current::kBoth, Sign::kBoth).ToTH1(18e20));
          
          for(int bin = 0; bin < nbins; ++bin)
          {
              if(i == 0)
              {
                  curves[bin] = new TGraph;
                  points[bin] = new TGraph;
              }
              
              const double ratio = h->GetBinContent(bin+1)/nom->GetBinContent(bin+1);
              
              if(!std::isnan(ratio)) curves[bin]->SetPoint(curves[bin]->GetN(), x, ratio);
              else curves[bin]->SetPoint(curves[bin]->GetN(), x, 1);
          }
      } 
      
      // As elswhere, to allow BirksC etc that need a different nominal to plot
      // right.
      IPrediction* pNom = 0;
      for(unsigned int shiftIdx = 0; shiftIdx < it->second.shifts.size(); ++shiftIdx)
      {
          if(it->second.shifts[shiftIdx] == 0) pNom = it->second.preds[shiftIdx];;
      }
      if(pNom)
      { // if not, probably MinimizeMemory() was called
          std::unique_ptr<TH1> hnom(pNom->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth).ToTH1(18e20));
          
          for(unsigned int shiftIdx = 0; shiftIdx < it->second.shifts.size(); ++shiftIdx)
          {
              if(!it->second.preds[shiftIdx]) continue; // Probably MinimizeMemory()
              std::unique_ptr<TH1> h(it->second.preds[shiftIdx]->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth).ToTH1(18e20));
              
              for(int bin = 0; bin < nbins; ++bin)
              {
                  const double ratio = h->GetBinContent(bin+1)/hnom->GetBinContent(bin+1);
                  if(!std::isnan(ratio)) points[bin]->SetPoint(points[bin]->GetN(), it->second.shifts[shiftIdx], ratio);
                  else points[bin]->SetPoint(points[bin]->GetN(), it->second.shifts[shiftIdx], 1);
              }
          }
      }
      
      
      int nx = int(sqrt(nbins));
      int ny = int(sqrt(nbins));
      if(nx*ny < nbins) ++nx;
      if(nx*ny < nbins) ++ny;
      
      TCanvas* c = new TCanvas;
      c->Divide(nx, ny);
      
      for(int bin = 0; bin < nbins; ++bin)
      {
          c->cd(bin+1);
          (new TH2F("",
          TString::Format("%s %g < %s < %g;Shift;Ratio",
          it->second.systName.c_str(),
          nom->GetXaxis()->GetBinLowEdge(bin+1),
          nom->GetXaxis()->GetTitle(),
          nom->GetXaxis()->GetBinUpEdge(bin+1)),
          100, -4, +4, 100, .5, 1.5))->Draw();
          curves[bin]->Draw("l same");
          points[bin]->SetMarkerStyle(kFullDotMedium);
          points[bin]->Draw("p same");
      }

      c->SaveAs(Form("Debug_flav_%d.root", int(flav)));
      
      c->cd(0);
  }

  //----------------------------------------------------------------------
  void NDPredictionSystsSingleElectron::DebugPlotColz(osc::IOscCalc* calc, const ISyst* syst, Flavors::Flavors_t flav) const
  {
      InitFits(calc);
      
      std::unique_ptr<TH1> nom(fPredNom->PredictComponent(calc, flav, Current::kBoth, Sign::kBoth).ToTH1(18e20));
      const int nbins = nom->GetNbinsX();
      
      TH2* h2 = new TH2F("", (syst->LatexName()+";;#sigma").c_str(),nbins, nom->GetXaxis()->GetXmin(), nom->GetXaxis()->GetXmax(), 80, -4, +4);
      h2->GetXaxis()->SetTitle(nom->GetXaxis()->GetTitle());
      
      for(int i = 1; i <= 80; ++i)
      {
          const double y = h2->GetYaxis()->GetBinCenter(i);
          const SystShifts ss(syst, y);
          std::unique_ptr<TH1> h(PredictComponentSyst(calc, ss, flav, Current::kBoth, Sign::kBoth).ToTH1(18e20));
          
          for(int bin = 0; bin < nbins; ++bin)
          {
              const double ratio = h->GetBinContent(bin+1)/nom->GetBinContent(bin+1);
              
              if(!isnan(ratio) && !isinf(ratio))
              h2->Fill(h2->GetXaxis()->GetBinCenter(bin), y, ratio);
          }
      }
      
      h2->Draw("colz");
      h2->SetMinimum(0.5);
      h2->SetMaximum(1.5);
      h2->SaveAs(Form("DebugColz_flav_%d.root", int(flav)));
  }
}
