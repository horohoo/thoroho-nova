#include "const.h"

using namespace ana;

TH1F* LoadLdmNum();

void load_prediction(int dmmass = 200, TString options="nd_flux_pileup_xsec", int iCuts = 5, bool isFHC = true, bool mock = false)
{
  std::cout << "event started" << std::endl;
    TStopwatch sw;
    sw.Start();

    TString outsuffix = options+Form("_sys_%s_%s_data_cut_%d.root", isFHC ? "fhc" : "rhc", mock ? "mock" : "fake", iCuts);
    std::string outsuffix_string(outsuffix.Data());

    std::string pred_outname = "preds_" + outsuffix_string + ".root";
    
    std::cout << "Loading BdNMC files... \n";
    TH1F* ldmNum = LoadLdmNum();
    
    //The calculator will be used to genreate MC prediction
    osc::OscCalcSingleElectron calc; 

                                            
    //Cuts for event selection
    const Cut sel_cut   = ldmone::SingleElecEventCVNCutFlow[iCuts].cut;
    const Cut nuone_cut = nuone::kintType == 1098 && nuone::kisVtxCont == 1;
    const Cut mec_cut = MEC && kIsNueCC;
    const Cut numu_cut  = !nuone_cut && !mec_cut;

    std::vector<double> DMMass;
    std::vector<double> ULimit;
    std::vector<double> UScale;
    std::vector<double> DMFit;
    std::vector<double> IBkgFit;
    std::vector<double> BkgFit;
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
        nuonescale = fhcnuonescale;
	mecscale   = fhcmecscale;
        
        fnominal.assign(ffhc);
    }
    else
    {
        std::cout << "Loading RHC samples ... \n" << std::endl;
        bkgscale   = rhcbkgscale;
        nuonescale = rhcnuonescale;
        
        fnominal.assign(frhc);
    }

    //Construct prediction
    SpectrumLoader loaderldm(fldm);
    SpectrumLoader loadernuone(fnuone);
    SpectrumLoader loadernominal(fnominal);
    SpectrumLoader loadermec(fmec);
    
    //bins and axis for the analysis
    //??how to fit a sub-range of the specra?
    const Binning bins  = Binning::Simple(20, 0, 0.02);
    const HistAxis etheta2Axis("E #theta^{2} (GeV Rad^{2})", bins, nuone::kETheta2);
    NDPredictionSystsSingleElectron pred(loaderldm, loadernuone, loadernominal, loadermec, etheta2Axis, sel_cut, sel_cut&&nuone_cut, sel_cut&&numu_cut, sel_cut&&mec_cut, ana::kPPFXFluxCVWgt*kradWt, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51);
    
    std::cout << "\nLoading dark matter ... \n";
    loaderldm.Go();
    std::cout << "\nLoading nu-on-e ... \n";
    loadernuone.Go();
    std::cout << "\nLoading nominal ...  \n";
    loadernominal.Go();
    std::cout <<"\nLoading special MEC file ... \n";
    loadermec.Go();

    std::cout << "\nFinished loading all samples ... \n";

    //process at each DM mass point
    const double pot = 1.5e21;//2.5e21
    for(auto i_dminf : dmInf)
    {
        int i_dm = i_dminf.DmMassPoints;
        if (i_dm > dmmass)
        {
            break;
        }

	std::cout << "\nProcessing DM mass point " << i_dm << " MeV" << std::endl;
        
        double npred  = ldmNum->GetBinContent(ldmNum->FindBin(i_dm));
        std::cout << "\nDM mass: " << i_dm << " MeV, Number of predicted LDM events: " << npred << std::endl;
        
        double potAjust = i_dminf.DMSimuPOT;
        std::cout << "\nPOT ajust to: " << potAjust << std::endl;
        pred.AjustPOT(potAjust); //Adjust the signal POT according to the BdNMC simulation before making prediction
        
        double potset = pred.GetPOT();

        calc.SetAna(true);
        calc.SetBkgScale(bkgscale);
        calc.SetIBkgScale(nuonescale);
	calc.SetMECScale(mecscale);
        calc.SetSigScale(ldmScale);
        calc.SetDMMass(i_dm);
        calc.SetDMFile(dmFile);

        std::cout << "Prediction Loaded, making prediction... \n";
        Spectrum spred = pred.Predict(&calc);
        
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
	
	TFile *fitFile = new TFile(outDir+Form("Prediction_DM_%dMeV_", i_dm) + outsuffix, "recreate"); 
	fitFile->cd();
	pred.SaveAs(&calc, fitFile, pred_outname);

	fitFile->Close();
	delete fitFile;
    }

    sw.Stop();
    std::cout << "\nAll done in: ";
    sw.Print();

    return;
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
