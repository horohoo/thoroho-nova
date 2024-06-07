#include "const.h"

using namespace ana;

void MakeMaps(std::vector <const IFitVar*> profVars,
              std::map<const IFitVar*, TGraph*> &profVarsMap,
              std::vector <const ISyst* > profSysts,
              std::map<const ISyst*, TGraph*> &profSystsMap);
void SaveMaps(TDirectory * dir,
              std::vector <const IFitVar*> profVars,
              std::vector <TString > profnames,
              std::map<const IFitVar*, TGraph*> &profVarsMap,
              std::vector <const ISyst* > profSysts,
              std::map<const ISyst*, TGraph*> &profSystsMap);
TGraph* DMSlice(const IExperiment* expt,
                osc::IOscCalcAdjustable* calc, const IFitVar* v,
                int nbinsx, double minx, double maxx, double minchi,
                std::vector<const IFitVar*> profVars,
                std::vector<const ISyst*> profSysts,
                const SeedList& seedPts,
                const std::vector<SystShifts>& systSeedPts,
                std::map<const IFitVar*, TGraph*>& profVarsMap,
                std::map<const ISyst*, TGraph*>& profSystsMap,
                double &upperLimit);
TGraph* DMProfile(const IExperiment* expt,
                  osc::IOscCalcAdjustable* calc, const IFitVar* v,
                  int nbinsx, double minx, double maxx,
                  double input_minchi,
                  const std::vector<const IFitVar*>& profVars,
                  const std::vector<const ISyst*>& profSysts,
                  const SeedList& seedPts,
                  const std::vector<SystShifts>& systSeedPts,
                  std::map<const IFitVar*, TGraph*>& profVarsMap,
                  std::map<const ISyst*, TGraph*>& profSystsMap);
TH1D* PullTerm(const SystShifts & shifts, bool sortName);
TH1F* LoadLdmNum();


void run_fit_spectra(int dmmass = 10, TString options="nd_flux_pileup_xsec", int iCuts = 5, bool isFHC = true, bool mock = false)
{
    TStopwatch sw;
    sw.Start();

    bool statonly  = options.Contains("statonly");
    bool ndsys     = options.Contains("nd");
    bool fluxsys   = options.Contains("flux");
    bool pileupsys = options.Contains("pileup");
    bool xsecsys   = options.Contains("xsec");
    bool fullsys   = ndsys && fluxsys && pileupsys  && xsecsys;
    assert (!(statonly && (ndsys || fluxsys || pileupsys || xsecsys)));

    if (statonly)
    {
        std::cout << "Fit with statistical uncertainties only...\n";
    }
    else if (fullsys)
    {
        std::cout << "Fit with statistical and full systematical uncertainties only...\n";
    }
    else if (ndsys)
    {
        std::cout << "Fit with statistical and ND detector systematical uncertainties ...\n";
    }
    else if (fluxsys)
    {
        std::cout << "Fit with statistical and ND flux systematical uncertainties ...\n";
    }
    else if (pileupsys)
    {
        std::cout << "Fit with statistical and pileup systematical uncertainties ...\n";
    }
    else if (xsecsys)
    {
        std::cout << "Fit with statistical and cross section systematical uncertainties ...\n";
    }

    TH1F* ldmNum = LoadLdmNum();

    //The calculator will be used to generate MC prediction
    osc::OscCalcSingleElectron calc;

    //Systematics                                                        
    std::vector<const ISyst*> systs;
    std::vector <SystShifts> seedShifts = {};
    SystShifts auxShifts = SystShifts::Nominal();  // A container to hold the best shifts of the Systs 

    if (ndsys)
    {
	std::cout << "Adding ND detector systematics... \n";
        systs.push_back(&kNDldmCalibSyst);
        systs.push_back(&kNDldmLightSyst);
        systs.push_back(&kNDldmCherSyst);

        for (double systshift:{-2, +2})
	{
            SystShifts Shiftscali (&kNDldmCalibSyst, systshift);
            seedShifts.emplace_back(std::move(Shiftscali));
            SystShifts Shiftsllevl (&kNDldmLightSyst, systshift);
            seedShifts.emplace_back(std::move(Shiftsllevl));
            SystShifts Shiftscher (&kNDldmCherSyst, systshift);
            seedShifts.emplace_back(std::move(Shiftscher));
	}
    }

    if (fluxsys)
    {
	std::cout << "Adding ND flux systematics... \n";
        systs.push_back(GetFluxPrincipalsND2020(0));
        systs.push_back(GetFluxPrincipalsND2020(1));
        systs.push_back(GetFluxPrincipalsND2020(2));
        systs.push_back(GetFluxPrincipalsND2020(3));
        systs.push_back(GetFluxPrincipalsND2020(4));
        systs.push_back(GetFluxPrincipalsND2020(5));
        systs.push_back(GetFluxPrincipalsND2020(6));

        for (double systshift:{-2, +2})
	{
            SystShifts ShiftFPP0 (GetFluxPrincipalsND2020(0), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP0));
            SystShifts ShiftFPP1 (GetFluxPrincipalsND2020(1), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP1));
            SystShifts ShiftFPP2 (GetFluxPrincipalsND2020(2), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP2));
            SystShifts ShiftFPP3 (GetFluxPrincipalsND2020(3), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP3));
            SystShifts ShiftFPP4 (GetFluxPrincipalsND2020(4), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP4));
            SystShifts ShiftFPP5 (GetFluxPrincipalsND2020(5), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP5));
            SystShifts ShiftFPP6 (GetFluxPrincipalsND2020(6), systshift);
            seedShifts.emplace_back(std::move(ShiftFPP6));
	}
    }

    if (pileupsys)
    {
	std::cout << "Adding pileup systematics... \n";
        systs.push_back(&kNDPileupEffectSyst);

        for (double systshift:{-2, +2})
	{
            SystShifts tempShifts (&kNDPileupEffectSyst, systshift);
            seedShifts.emplace_back(std::move(tempShifts));
	}
    }

    if (xsecsys)
    {
	std::cout << "Adding cross section systematics... \n";
	std::vector<const ISyst*> XSectSys = getAllXsecSysts_2020_GSF();
        systs.insert(systs.end(), XSectSys.begin(), XSectSys.end());
    }

    TString outsuffix = options+Form("_sys_%s_%s_data_cut_%d.root", isFHC ? "fhc" : "rhs", mock ? "mock" : "fake", iCuts);
    std::string outsuffix_string(outsuffix.Data());
    std::string pred_outname = "preds_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root.root";
    std::cout << "pred_outname: " << pred_outname << std::endl;
    //std::string pred_outname = "preds_" + outsuffix_string + ".root";

    //Vars to be fitted
    
    std::vector <const IFitVar*> fitvars = {&kFitSigScalingSingleElectron,
                                            &kFitIBkgScalingSingleElectron,
                                            &kFitBkgScalingSingleElectron,
                                            &kFitMECScalingSingleElectron};
    
    //std::vector <const IFitVar*> fitvars = {&kFitSigScalingSingleElectron};

    //Seeds
    std::vector<double> dmscale_seeds   = {1e-20, 1e-5};
    std::vector<double> ibkgscale_seeds = {1.00, 1.01};//{1.04, 1.06};
    std::vector<double> bkgscale_seeds  = {1.00, 1.01};//{0.94, 0.96};
    std::vector<double> mecscale_seeds  = {1.00, 1.01};
    const SeedList& seedFitVars = SeedList({{&kFitSigScalingSingleElectron,  dmscale_seeds},
    	                                    {&kFitIBkgScalingSingleElectron, ibkgscale_seeds},
					    {&kFitBkgScalingSingleElectron,  bkgscale_seeds},
					    {&kFitMECScalingSingleElectron, mecscale_seeds}});


    //Cuts for event selection
    const Cut sel_cut   = ldmone::SingleElecEventCVNCutFlow[iCuts].cut;
    const Cut nuone_cut = nuone::kintType == 1098 && nuone::kisVtxCont == 1;
    const Cut mec_cut = MEC;
    const Cut nuecc_cut = kIsNueCC;
    const Cut numu_cut  = !nuone_cut && !mec_cut && !nuecc_cut;

    
    std::vector<double> DMMass;
    std::vector<double> ULimit;
    std::vector<double> UScale;
    std::vector<double> DMFit;
    std::vector<double> IBkgFit;
    std::vector<double> BkgFit;
    std::vector<double> MECFit;
    std::vector<double> Chi2Fit;

    std::cout << "ldmscale: "  << ldmScale << std::endl;
    double bkgscale;
    double nuonescale;
    double mecscale;
    std::string fnominal;
    if(isFHC)
      {
	std::cout << "Loading FHC samples ... \n" << std::endl;
        bkgscale   = fhcbkgscale;
        nuonescale = 0.93*fhcnuonescale;
	mecscale = fhcmecscale;

        fnominal.assign(ffhc);
      }
    else
      {
	std::cout << "Loading RHC samples ... \n" << std::endl;
        bkgscale   = rhcbkgscale;
        nuonescale = rhcnuonescale;

        fnominal.assign(frhc);
      }

    // construct prediction (needed to LoadFrom() later)
    SpectrumLoader loaderldm(fldm);
    SpectrumLoader loadernuone(fnuone);
    SpectrumLoader loadernominal(fnominal);
    SpectrumLoader loadermec(fmec);

    //std::vector<double> bin_edges = {0.0, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005, 0.0055, 0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.0095, 0.01, 0.0105, 0.011, 0.0115, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, 0.02};
    //const Binning bins  = Binning::Custom(bin_edges);
    const Binning bins  = Binning::Simple(20, 0, 0.02);
    const HistAxis etheta2Axis("E #theta^{2} (GeV Rad^{2})", bins, nuone::kETheta2);
    NDPredictionSystsSingleElectron pred(loaderldm, loadernuone, loadernominal, loadermec, etheta2Axis, sel_cut, sel_cut&&nuone_cut, sel_cut&&numu_cut, sel_cut&&mec_cut, ana::kPPFXFluxCVWgt*kradWt, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51);


    const double pot = 1.25e21;
    for(auto i_dminf : dmInf)
    {
      int i_dm = i_dminf.DmMassPoints;
      if (i_dm > dmmass)
      {
	break;
      }

      double npred = ldmNum->GetBinContent(ldmNum->FindBin(i_dm));
      std::cout << "\nDM mass: " << i_dm << " MeV, Number of predicted LDM events: " << npred << std::endl;

      // Load the prediction from file
      std::cout << "Loading prediction of mass " << i_dm << " MeV from file..." << std::endl;
      std::cout << "Using file " << "Prediction_DM_5MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root" << std::endl;
      TFile *predFile = new TFile(outDir+"Prediction_DM_5MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root");//outsuffix);
      predFile->cd();
      static std::unique_ptr<NDPredictionSystsSingleElectron> pred_ptr = pred.LoadFrom(predFile, pred_outname);

      double potAdjust = i_dminf.DMSimuPOT;
      pred_ptr->AjustPOT(potAdjust);
      std::cout << "Signal POT: " << pred_ptr->GetPOT() << std::endl;

      calc.SetAna(true);
      calc.SetBkgScale(bkgscale);
      calc.SetIBkgScale(nuonescale);
      calc.SetMECScale(mecscale);
      calc.SetSigScale(ldmScale);
      calc.SetDMMass(i_dm);
      calc.SetDMFile(dmFile);

      std::cout << "Prediction Loaded, making prediction... \n";
      Spectrum spred = pred_ptr->Predict(&calc);

      //Generate data                       
      Spectrum data = Spectrum::Uninitialized();
      if (mock)
        {
	  std::cout <<"Building Mock Data... \n";
	  data = spred.MockData(pot);
        }
      else
        {
	  std::cout <<"Building Fake Data... \n";
	  data = spred.FakeData(pot);
        }

      calc.SetBkgScale(1.0);
      calc.SetIBkgScale(1.0);
      calc.SetMECScale(1.0);

        std::cout << "\nFitting ...\n";
        const IExperiment* expt = new SingleSampleExperiment(&*pred_ptr, data);
        
        // Find the best fit points
        double minichi = 1E20;
        calc.SetSigScale(1e-30);
        
        MinuitFitter fitdm(expt, fitvars, systs);        
        
        auto thisminchi   = fitdm.Fit(&calc, auxShifts, seedFitVars, seedShifts, IFitter::kQuiet)->EvalMetricVal();
        double nominalVal = kFitBkgScalingSingleElectron.GetValue(&calc);
        double nuoneVal   = kFitIBkgScalingSingleElectron.GetValue(&calc);
	double mecVal     = kFitMECScalingSingleElectron.GetValue(&calc);
        double ldmVal     = kFitSigScalingSingleElectron.GetValue(&calc);
                
        std::cout << "Min Chi:" << thisminchi << std::endl;
        std::cout << "DM Y value  True: " << sqrt(ldmScale*epsilon4)*constY << ", Fitted scale: " << ldmVal << ", Y: "<< sqrt(ldmVal*epsilon4)*constY << std::endl;
        std::cout << "Nuone scale True: " << nuonescale << ", Fitted: " << nuoneVal << std::endl;
        std::cout << "Numi scale True: " << bkgscale << ", Fitted: " << nominalVal << std::endl;
	std::cout << "MEC scale True: " << mecscale << ", Fitted: " << mecVal << std::endl;
        
        //Save the details for each mass point
        TFile * fitFile = new TFile(outDir+Form("Fit_DM_%dMeV_", i_dm)+outsuffix, "recreate");
        fitFile->cd();
        
        //Systematic pulls
        if(ndsys || fluxsys || pileupsys)
        {
            std::cout << "Pull systematics ... \n";
            TH1D* shifts = PullTerm(auxShifts, true);
            TString str = "Best fit " ;
            for (auto &v:fitvars)
            {
                str += TString::Format(" %s=%.3f ",v->LatexName().c_str(),v->GetValue(&calc));
            }
            str+= TString::Format(" LL=%.6f", thisminchi);
            shifts->SetTitle(str);
            TCanvas* c_systs = new TCanvas ("systs", "systs", 800, 600);
            c_systs->SetBottomMargin(0.6);
            shifts->Draw("lp hist");
            gPad->Update();
            TLine *l = new TLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);
            l->Draw("same");
            c_systs->Write();
            shifts->Write();
            delete c_systs;
        }
        
        //Slicing
        std::cout << "Slicing ... \n";
        TVectorD v(1);
        v[0] = thisminchi;
        v.Write("minchi");
        
        int steps = 2000;
        const IFitVar * fitvar;
        double minX;
        double maxX;
        std::vector <const IFitVar * > profVars;
        std::vector <TString> profvarnames;
        std::map  <const IFitVar *, std::vector <double> >profseeds;
        std::vector <SystShifts> profsysseeds;  
        TString shortname;
        
        fitvar       = &kFitSigScalingSingleElectron;
        shortname    = "DM Y";
        minX         = 0.0;
        maxX         = i_dminf.DMExclC;
        profVars     = {&kFitBkgScalingSingleElectron, &kFitIBkgScalingSingleElectron, &kFitMECScalingSingleElectron};
        profvarnames = {"Bkg", "Nuone", "MEC"};
        profseeds    = {{&kFitIBkgScalingSingleElectron, {nuoneVal}},
                        {&kFitBkgScalingSingleElectron,  {nominalVal}},
			{&kFitMECScalingSingleElectron, {mecVal}}};
        profsysseeds.emplace_back(std::move(auxShifts)); 
        
        std::map<const IFitVar*, TGraph*> profVarsMap;
        std::map<const ISyst*,   TGraph*> profSystsMap;

        double upperLimit;
        
        MakeMaps (profVars, profVarsMap, systs, profSystsMap);
        TGraph* slice = DMSlice(expt, &calc, fitvar, steps, minX, maxX, thisminchi, profVars, systs,
                                profseeds, profsysseeds, profVarsMap, profSystsMap, upperLimit);

        if(i_dm > 5 && upperLimit < UScale.back() )
        {
            std::cout << "\nCannot find limit for DM mass: " << i_dm << "MeV, skip this point\n";
            fitFile->Write();
            fitFile->Close();
            delete fitFile;
            continue;
        }
        std::cout << "\nUpper Limit: " << sqrt(upperLimit*epsilon4)*constY << std::endl;
        DMMass.push_back(i_dm);
        DMFit.push_back(ldmVal);
        BkgFit.push_back(nominalVal);
        IBkgFit.push_back(nuoneVal);
	MECFit.push_back(mecVal);
        Chi2Fit.push_back(thisminchi);

        UScale.push_back(upperLimit);
        ULimit.push_back(sqrt(upperLimit*epsilon4)*constY);
        
        slice->Write(shortname);
        SaveMaps (fitFile, profVars, profvarnames, profVarsMap, systs, profSystsMap);
        
        fitFile->Write();
        fitFile->Close();
        delete fitFile;

	predFile->Close();
	delete predFile;
    }

    std::cout << "\nMaking confidence level ...\n";
    TFile* exFile = new TFile(outDir+"Exclusion_"+outsuffix, "recreate");
    exFile->cd();
    TGraph* gNOvA = new TGraph(DMMass.size(), &*DMMass.begin(), &*ULimit.begin());
    gNOvA->SetName(Form("gNOvA%s", isFHC ? "fhc" : "rhc"));
    gNOvA->SetTitle("; m_{#chi}(MeV/c^{2}); Y=#epsilon^{2}#alpha^{'}(m_{#chi}/m_{V})^{4}");
    gNOvA->Write();

    std::cout << "\nPloting sample spectra ...\n";
    TH1* hData = NULL;
    Spectrum data = Spectrum::Uninitialized();    
    for(int i = 0; i < DMMass.size(); i ++)
    {
        int i_dm = (int)DMMass[i];
        calc.SetAna(true);
        calc.SetDMFile(dmFile);
        calc.SetBkgScale(bkgscale);
        calc.SetIBkgScale(nuonescale);
	calc.SetMECScale(mecscale);
        calc.SetSigScale(ldmScale);
        calc.SetDMMass(DMMass[i]);

	TFile *predFile = new TFile(outDir+"Prediction_DM_5MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root");// + outsuffix);
	predFile->cd();
	static std::unique_ptr<NDPredictionSystsSingleElectron> pred_ptr = pred.LoadFrom(predFile, pred_outname);
        
	exFile->cd();

        Spectrum spred = pred_ptr->Predict(&calc);
        int binnum;
        
        if(i == 0)
        {        
            if (mock)
            {
                data = spred.MockData(pot);
            }
            else
            {
                data = spred.FakeData(pot);
            }
            
            hData = data.ToTH1(pot, 2, 1);
            hData->SetName("hData");
            hData->SetMarkerStyle(33);
            hData->Sumw2();
            hData->Write();
        }

        THStack * hStackBkg = new THStack(Form("hStackBkg_%dMeV", i_dm), Form("hStackBkg_%dMeV", i_dm));
        TH1* hFit   = NULL;
        TH1* hIBkg = NULL;
        TH1* hBkg  = NULL;
	TH1* hMEC = NULL;
        
        calc.SetBkgScale(BkgFit[i]);
        calc.SetIBkgScale(IBkgFit[i]);
	calc.SetMECScale(MECFit[i]);
        calc.SetSigScale(DMFit[i]);
        
        hFit = pred_ptr->Predict(&calc).ToTH1(pot, kBlack, 1);
        hFit->SetName(Form("hFit_%dMeV", i_dm));
        binnum = hFit->GetNbinsX();
        for(int bin = 1; bin <= binnum ; bin++)
        {
            hFit->SetBinError(bin, sqrt(hFit->GetBinContent(bin)));
        }
        hFit->Sumw2();

        hIBkg = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuMu, Current::kCC, Sign::kNu).ToTH1(pot, 16, 4);
        hIBkg->SetName(Form("hNuone_%dMeV", i_dm));
        hIBkg->SetFillColor(17);
        
        hBkg = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuTau, Current::kCC, Sign::kNu).ToTH1(pot, 13, 5);
        hBkg->SetName(Form("hNumi_%dMeV", i_dm));
        hBkg->SetFillColor(15);

	hMEC = pred_ptr->PredictComponent(&calc, Flavors::kNuMuToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 13, 5);
	hMEC->SetName(Form("hMEC_%dMeV", i_dm));
	hMEC->SetFillColor(15);
        
        hStackBkg->Add(hBkg);
	hStackBkg->Add(hMEC);
        hStackBkg->Add(hIBkg);
        
        calc.SetSigScale(UScale[i]);
        std::cout << "Plot spectrum for DM mass: " << i_dm << "MeV, Y: " << ULimit[i] << ", scale: " << UScale[i]
                  << ", mass: " << dmInf[i].DmMassPoints  << ", effective POT: " << dmInf[i].DMSimuPOT*(bdnmcY/ULimit[i]) << std::endl;

        TH1* hExDM = NULL;
        hExDM = pred_ptr->Predict(&calc).ToTH1(pot, i+3, 1);
        hExDM->SetName(Form("hExDM_%dMeV", i_dm));
        binnum = hExDM->GetNbinsX();
        for(int bin = 1; bin <= binnum ; bin++)
        {
            hExDM->SetBinError(bin, sqrt(hExDM->GetBinContent(bin)));
        }
        hExDM->Sumw2();
        
        TH1* hSig = NULL;
        hSig = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, i+3, 1);
        hSig->SetName(Form("hSig_%dMeV", i_dm));
        hSig->SetFillColor(i+3);
        binnum = hSig->GetNbinsX();
        for(int bin = 1; bin <= binnum ; bin++)
        {
            hSig->SetBinError(bin, sqrt(hSig->GetBinContent(bin)));
        }
        
        THStack * hStackSig = new THStack(Form("hStackSig_%dMeV", i_dm), Form("hStackSig_%dMeV", i_dm));
        hStackSig->Add(hBkg);
	hStackSig->Add(hMEC);
        hStackSig->Add(hIBkg);
        hStackSig->Add(hSig);
                
        TH1* hRatioData = new TH1D(Form("hRatioData_%dMeV", i_dm), Form("hRatioData_%dMeV", i_dm), 20, 0, 0.02);
        hRatioData->Divide(hData, hFit, 1.0, 1.0, "B");
        hRatioData->SetMarkerStyle(33);
        hRatioData->SetMarkerColor(2);
        hRatioData->SetLineColor(2);
        hRatioData->SetLineStyle(1);
        hRatioData->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
        hRatioData->GetYaxis()->SetTitle("Ratio to Fit Result");
        hRatioData->GetYaxis()->CenterTitle(true);

        TH1* hRatioExDM = new TH1D(Form("hRatioExDM_%dMeV", i_dm), Form("hRatioExDM_%dMeV", i_dm), 20, 0, 0.02);        
        hRatioExDM->Divide(hExDM, hFit, 1.0, 1.0, "B");
        hRatioExDM->SetLineColor(i+3);
        hRatioExDM->SetLineStyle(1);
        hRatioExDM->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
        hRatioExDM->GetYaxis()->SetTitle("Ratio to Fit Result");
        hRatioExDM->GetYaxis()->CenterTitle(true);

        gROOT->Reset();
        gROOT->ForceStyle();
        TCanvas* cSpectra = new TCanvas (Form("Spectra_%dMeV", i_dm), Form("Spectra_%dMeV", i_dm), 800, 800);
        cSpectra->SetBottomMargin(0.15);
        cSpectra->SetLeftMargin(0.15);
        gStyle->SetPadBorderMode(0);
        gStyle->SetFrameBorderMode(0);
        float small = 1e-5;

        cSpectra->Divide(1, 2, small, small);

        cSpectra->cd(1);
        gPad->SetBottomMargin(small);
        hStackSig->Draw("hist");
        hData->Draw("same P E");
        hFit->Draw("same hist");
        hStackSig->GetYaxis()->SetTitle(Form("Events (%.2e POT)", pot));
        hStackSig->GetYaxis()->CenterTitle(true);
        gStyle->SetOptTitle(1);
        gPad->Modified();
        gPad->Update();
        
        TLegend legend(0.5, 0.35, 0.89, 0.88);  
        legend.AddEntry(hData, Form("%s Data", mock ? "Mock" : "Fake"));
        legend.AddEntry(hFit, Form("Fit Results, LL: %.3e", Chi2Fit[i]));
        legend.AddEntry(hIBkg, Form("#nu-e, True: %.3f, Fit: %.3f", nuonescale, IBkgFit[i]));
	legend.AddEntry(hMEC, Form("MEC, True: %.f, Fit: %.f", mecscale, MECFit[i]));
        legend.AddEntry(hBkg,  Form("Others, True: %.3f, Fit: %.3f", bkgscale, BkgFit[i]));
        legend.AddEntry(hSig,  Form("#chi-e Mass: %dMeV, Y: %.3e", i_dm, ULimit[i]));
        legend.SetBorderSize(0);
        legend.SetTextFont(42);
        legend.Draw("same");
        gPad->Update();
        
        cSpectra->cd(2);
        gPad->SetTopMargin(small);
        gPad->SetTickx();
        hRatioExDM->Draw("P E");
        hRatioExDM->GetYaxis()->SetRangeUser(0.4, 1.6);
        hRatioData->Draw("same");
        gStyle->SetOptTitle(0);
        gPad->Modified();
        gPad->Update();
        TLine *l1 = new TLine(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
        l1->Draw("same");
                
        cSpectra->Write();
        
        hFit->Write();
        hBkg->Write();
        hIBkg->Write();
	hMEC->Write();
        
        hExDM->Write();
        hSig->Write();

        hRatioData->Write();
        hRatioExDM->Write();
        
        hStackBkg->Write();
        hStackSig->Write();
    }
    
    exFile->Write();
    exFile->Close();
    delete exFile;    
    
    sw.Stop();
    std::cout << "\nAll done in: ";
    sw.Print();
    return;
}

void MakeMaps(std::vector <const IFitVar*> profVars,
              std::map<const IFitVar*, TGraph*> &profVarsMap,
              std::vector <const ISyst* > profSysts,
              std::map<const ISyst*, TGraph*> &profSystsMap)
{
    for (const IFitVar * var : profVars)
    {
        profVarsMap.insert(std::pair<const IFitVar*, TGraph*> (var, new TGraph()));
    }
                
    for (const ISyst* syst : profSysts)
    {
        profSystsMap.insert(std::pair<const ISyst*,   TGraph*> (syst, new TGraph()));
    }
}

void SaveMaps(TDirectory * dir,
              std::vector <const IFitVar*> profVars,
              std::vector <TString > profnames,
              std::map<const IFitVar*, TGraph*> &profVarsMap,
              std::vector <const ISyst* > profSysts,
              std::map<const ISyst*, TGraph*> &profSystsMap)
{
    TDirectory *tmp = gDirectory;
    dir->cd();
    
    for (int i = 0; i < (int) profVars.size(); ++i)
    {
        profVarsMap[profVars[i]]->Write(("ProfVar-" + profnames[i]));
    }

    for (int i = 0; i < (int) profSysts.size(); ++i)
    {
        profSystsMap[profSysts[i]]->Write("ProfSysts");
    }
      
    tmp->cd();
}

TGraph* DMSlice(const IExperiment* expt,
                osc::IOscCalcAdjustable* calc, const IFitVar* v,
                int nbinsx, double minx, double maxx, double minchi,
                std::vector<const IFitVar*> profVars,
                std::vector<const ISyst*> profSysts,
                const SeedList& seedPts,
                const std::vector<SystShifts>& systSeedPts,
                std::map<const IFitVar*, TGraph*>& profVarsMap,
                std::map<const ISyst*, TGraph*>& profSystsMap,
                double &upperLimit)
{
    double bestfit = minx > 0.0 ? sqrt(minx*2*epsilon4)*constY: 0.0;
    double maxchi;
    double maxY;
    if(minx < 0.0)
    {
        nbinsx /= 10;
    }
    
    TGraph* ret = DMProfile(expt, calc,
                            v, nbinsx, minx, maxx,
                            minchi, profVars, profSysts, seedPts, systSeedPts,
                            profVarsMap, profSystsMap);
    for(int i = 0; i < ret->GetN(); ++i)
    {
        double x, y;
        ret->GetPoint(i, x, y);
        ret->SetPoint(i, sqrt(x*epsilon4)*constY, y > 0 ? y : 0);
        if(y < Chi2for90CL)
        {
            upperLimit = x;
        }
        if(i == ret->GetN()-1)
        {
            std::cout << "Max Chi2 = " << y << ", x = " << x << ", iteration = " << i << std::endl;
            maxchi = y;
            maxY   = sqrt(x*epsilon4)*constY;
        }
    }

    TCanvas* c_chisquare = new TCanvas ("chisquare", "chisquare", 800, 600);
    c_chisquare->SetBottomMargin(0.15);
    c_chisquare->SetLeftMargin(0.12);
    ret->SetTitle("#chi^{2} vs. Y value; Y=#epsilon^{2}#alpha_{D}(m_{#chi}/m_{V})^{4}; #Delta#chi^{2}");
    ret->SetMarkerStyle(1);
    ret->SetMarkerColor(1);
    ret->SetLineColor(1);
    ret->SetLineStyle(1);
    ret->SetLineWidth(2);
    
    ret->Draw("AC");
    ret->GetXaxis()->SetLimits(bestfit>0.0?0.0:bestfit*1.2, maxY);
    ret->SetMinimum(0);
    ret->SetMaximum(Chi2for90CL*1.2);   

    gPad->Update();   
    
    TLine* lChiSquare = new TLine(gPad->GetUxmin(), Chi2for90CL, gPad->GetUxmax(), Chi2for90CL);
    TLine* lLimit     = new TLine(sqrt(upperLimit*epsilon4)*constY, gPad->GetUymin(), sqrt(upperLimit*epsilon4)*constY, gPad->GetUymax());
    TLine* lBestFit   = new TLine(bestfit, gPad->GetUymin(), bestfit, gPad->GetUymax());
        
    lChiSquare->SetLineWidth(2);
    lChiSquare->SetLineStyle(3);
    lChiSquare->SetLineColor(8);
    
    lLimit->SetLineWidth(2);
    lLimit->SetLineStyle(7);
    lLimit->SetLineColor(3);
    
    lBestFit->SetLineWidth(2);
    lBestFit->SetLineStyle(7);
    lBestFit->SetLineColor(6);
    
    lChiSquare->Draw("same");
    lLimit->Draw("same");
    lBestFit->Draw("same");

    TLatex* tlimit = new TLatex(sqrt(upperLimit*epsilon4)*constY*1.01, 2.5, Form(".9 CL U-Limit Y = %6.2e", sqrt(upperLimit*epsilon4)*constY));
    tlimit->SetTextColor(12);
    tlimit->SetTextSize(0.05);
    tlimit->SetTextAngle(-90);

    TLatex* tbestfit = new TLatex(bestfit*1.01, 2.5, Form("Best Fit Y = %6.2e", bestfit));
    tbestfit->SetTextColor(12);
    tbestfit->SetTextSize(0.05);
    tbestfit->SetTextAngle(-90);

    TLatex* tchisquare = new TLatex(sqrt(upperLimit*epsilon4)*constY/3.0, Chi2for90CL+0.1, "#chi^{2} = 2.70554");
    tchisquare->SetTextColor(12);
    tchisquare->SetTextSize(0.05);

    tlimit->Draw("same");
    tbestfit->Draw("same");
    tchisquare->Draw("same");

    c_chisquare->Write();
        
    delete c_chisquare;

    for(const IFitVar* var: profVars)
    {
        for(int i = 0; i < profVarsMap[var]->GetN(); ++i)
        {
            double x, y;
            profVarsMap[var]->GetPoint(i, x, y);
            profVarsMap[var]->SetPoint(i, sqrt(x*epsilon4)*constY, y);
        }
    }
    for(const ISyst* s: profSysts)
    {
        for(int i = 0; i < profSystsMap[s]->GetN(); ++i)
        {
            double x, y;
            profSystsMap[s]->GetPoint(i, x, y);
            profSystsMap[s]->SetPoint(i, sqrt(x*epsilon4)*constY, y);
        }
    }
    
    return ret;
}

TGraph* DMProfile(const IExperiment* expt,
                  osc::IOscCalcAdjustable* calc, const IFitVar* v,
                  int nbinsx, double minx, double maxx,
                  double input_minchi,
                  const std::vector<const IFitVar*>& profVars,
                  const std::vector<const ISyst*>& profSysts,
                  const SeedList& seedPts,
                  const std::vector<SystShifts>& systSeedPts,
                  std::map<const IFitVar*, TGraph*>& profVarsMap,
                  std::map<const ISyst*, TGraph*>& profSystsMap)
{
    for(auto it: profVarsMap) delete it.second;
    for(auto it: profSystsMap) delete it.second;
    profVarsMap.clear();
    profSystsMap.clear();
    for(const IFitVar* v: profVars) profVarsMap[v] = new TGraph;
    for(const ISyst* s: profSysts) profSystsMap[s] = new TGraph;

    TGraph* ret = new TGraph;

    // Save the values of the fit vars as they were in the seed so we can put them back to that value every iteration.
    std::vector<double> seedValues;
    for(const IFitVar* v: profVars)
    {
        seedValues.push_back(v->GetValue(calc));
    }

    MinuitFitter fit(expt, profVars, profSysts);

    int count_over = 0;

    for(int n = 0; n <= nbinsx; ++n)
    {
        const double x = minx + (maxx-minx)*double(n)/nbinsx;
        v->SetValue(calc, x);
    
        // Put oscillation values back to their seed position each iteration
        for(unsigned int i = 0; i < seedValues.size(); ++i)
        {
            profVars[i]->SetValue( calc, seedValues[i] );
        }
        SystShifts systshift = SystShifts::Nominal();
        
        // Zero-subtract the result graph
        const double chi = fit.Fit(calc, systshift, seedPts, systSeedPts, MinuitFitter::kQuiet)->EvalMetricVal() - input_minchi;
        
        ret->SetPoint(ret->GetN(), x, chi);
    
        for(const IFitVar* var: profVars)
        {
            profVarsMap[var]->SetPoint(n, x, var->GetValue(calc));
        }
        for(const ISyst* s: profSysts)
        {
            profSystsMap[s]->SetPoint(n, x, systshift.GetShift(s));
        }        
    
        if(chi > Chi2for90CL)
        {            
            count_over += 1;
            if(count_over > 500)
            {
                break;
            }
        }
    }
  
    return ret;
}

TH1D* PullTerm(const SystShifts & shifts, bool sortName)
{
    auto systs = shifts.ActiveSysts();
    int nsysts = systs.size();
    
    TH1D* h = new TH1D ("hPullTerm",";; Pull (N#sigma)", nsysts, -0.5, nsysts + 0.5);
    
    for (int systIdx = 0; systIdx  < nsysts; ++systIdx)
    {
        double shiftThis = shifts.GetShift(systs[systIdx]);
        h->SetBinContent( systIdx + 1, shiftThis);
        
        TString tempName = systs[systIdx]->ShortName();
        h->GetXaxis()->SetBinLabel( systIdx + 1, systs[systIdx]->ShortName().c_str());
    }
    
    h->SetMarkerStyle(kFullCircle);
    h->SetTitleOffset(3.0);
    h->GetXaxis()->SetNdivisions(nsysts,kFALSE);
    h->GetXaxis()->SetLabelSize(0.065);
    h->GetXaxis()->LabelsOption("v");
    h->GetYaxis()->SetRangeUser(-1,1);
    return h;
}

TH1F* LoadLdmNum()
{
    // For LDM POT calculation
    TFile* bdnmcFile = TFile::Open("/exp/nova/app/users/thoroho/ldmanalysis/NuMagMomentAna/data/ldmspectra/ldmprediction.root", "READ");
    if(!bdnmcFile)
    {
        std::cout << "\nBdNMC prediction file does not exist, exit." << std::endl;
        abort();
    }
    TH1F* ldmNum = (TH1F*)bdnmcFile->Get("ldmnum");
    if(!ldmNum)
    {
        std::cout << "Cannot open ldmNum, exit." << std::endl;
        abort();
    }

    return ldmNum;
}
