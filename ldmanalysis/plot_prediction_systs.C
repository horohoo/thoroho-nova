#include "const.h"

using namespace ana;

TH1F* LoadLdmNum();

int ColorShift(int shift);

int StyleShift(int shift);

void plot_prediction_systs(int dmmass = 200, TString options="nd_flux_pileup_xsec", int iCuts = 5, bool isFHC = true, bool mock = false)
{
    TStopwatch sw;
    sw.Start();

    TH1F* ldmNum = LoadLdmNum();

    std::string pred_outname = "preds_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root.root";

    //The calculator will be used to generate MC prediction
    osc::OscCalcSingleElectron calc;


    std::vector<int> shifts = {+3, +2, +1, 0, -1, -2, -3};
    std::vector<const ISyst*> systs;
    
    systs.push_back(&kNDldmCalibSyst);
    systs.push_back(&kNDldmLightSyst);
    systs.push_back(&kNDldmCherSyst);

    systs.push_back(GetFluxPrincipalsND2020(0));
    systs.push_back(GetFluxPrincipalsND2020(1));
    systs.push_back(GetFluxPrincipalsND2020(2));
    systs.push_back(GetFluxPrincipalsND2020(3));
    systs.push_back(GetFluxPrincipalsND2020(4));
    systs.push_back(GetFluxPrincipalsND2020(5));
    systs.push_back(GetFluxPrincipalsND2020(6));
    systs.push_back(&kLDMFluxSyst);
    
    systs.push_back(&kNDPileupEffectSyst);

    std::vector<const ISyst*> XSectSys = getAllXsecSysts_2020_GSF();
    systs.insert(systs.end(), XSectSys.begin(), XSectSys.end());


    //Cuts for event selection
    const Cut sel_cut   = ldmone::SingleElecEventCVNCutFlow[iCuts].cut;
    const Cut nuone_cut = nuone::kintType == 1098 && nuone::kisVtxCont == 1;
    const Cut mec_cut = MEC;
    const Cut nuecc_cut = kIsNueCC;
    const Cut numu_cut  = !nuone_cut && !mec_cut && !nuecc_cut;

    // construct prediction (needed to LoadFrom() later)
    SpectrumLoader loaderldm(fldm);
    SpectrumLoader loadernuone(fnuone);
    SpectrumLoader loadernominal(ffhc);
    SpectrumLoader loadermec(fmec);

    const Binning bins  = Binning::Simple(20, 0, 0.02);
    const HistAxis etheta2Axis("E #theta^{2} (GeV Rad^{2})", bins, nuone::kETheta2);
    NDPredictionSystsSingleElectron pred(loaderldm, loadernuone, loadernominal, loadermec, etheta2Axis, sel_cut, sel_cut&&nuone_cut, sel_cut&&numu_cut, sel_cut&&mec_cut, ana::kPPFXFluxCVWgt*kradWt, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51);

    const double pot = 1.25e21;

    TFile *predFile = new TFile(outDir+"Prediction_DM_5MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root");
    predFile->cd();
    static std::unique_ptr<NDPredictionSystsSingleElectron> pred_ptr = pred.LoadFrom((predFile, pred_outname);
    
    TFile* exFile = new TFile(pdfDir+"PredictionSysts.root", "recreate");
    exFile->cd();
    
    calc.SetAna(true);
    calc.SetDMFile(dmFile);
    calc.SetBkgScale(1.0);
    calc.SetIBkgScale(1.0);
    calc.SetMECScale(1.0);
    calc.SetSigScale(1.0);

    for (const ISyst* syst : systs) {

      TCanvas* cBkg = new TCanvas (Form("Bkg_%s", syst->ShortName().c_str()), Form("Bkg_%s", syst->ShortName().c_str()), 800, 600);
      TCanvas* cSig1 = new TCanvas (Form("Sig_1MeV_%s", syst->ShortName().c_str()), Form("Sig_1MeV_%s", syst->ShortName()), 800, 600);
      TCanvas* cSig10 = new TCanvas (Form("Sig_10MeV__%s", syst->ShortName().c_str()), Form("Sig_10MeV__%s", syst->ShortName().c_str()), 800, 600);
      TCanvas* cSig20 = new TCanvas (Form("Sig_20MeV_%s", syst->ShortName().c_str()), Form("Sig_20MeV_%s", syst->ShortName().c_str()), 800, 600);
      TCanvas* cSig100 = new TCanvas (Form("Sig_100MeV_%s", syst->ShortName().c_str()), Form("Sig_100MeV_%s", syst->ShortName().c_str()), 800, 600);

      for (int shift : shifts) {

	SystShifts systshift (syst, shift);
	
	TH1* hIBkg = NULL;//new TH1F();
	TH1* hBkg  = NULL;//new TH1F();
	TH1* hMEC = NULL;//new TH1F();

        hIBkg = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuMu, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
        hBkg = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuTau, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
	hMEC = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuMuToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
        
	hBkg->Add(hMEC);
        hBkg->Add(hIBkg);

	cBkg->cd();

	hBkg->SetLineColor(ColorShift(shift));
	hBkg->SetLineStyle(StyleShift(shift));
	hBkg->SetLineWidth(4);

	if (shift == 3) {
	  hBkg->Draw("histo");
	  hBkg->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
	}
	else {
	  hBkg->Draw("same hist");
	}

	cSig1->cd();
        calc.SetSigScale(1);
        TH1* hSig1 = NULL;//new TH1F();
        hSig1 = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
	hSig1->SetLineColor(ColorShift(shift));
	hSig1->SetLineStyle(StyleShift(shift));
	hSig1->SetLineWidth(4);
	if (shift == 3) {
	  hSig1->Draw("histo");
	  hSig1->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
	}
	else {
	  hSig1->Draw("same hist");
	}

	cSig10->cd();
        calc.SetSigScale(10);
        TH1* hSig10 = NULL;
        hSig10 = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
        hSig10->SetLineColor(ColorShift(shift));
        hSig10->SetLineStyle(StyleShift(shift));
        hSig10->SetLineWidth(4);
        if (shift == 3) {
          hSig10->Draw("histo");
	  hSig10->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
	}
        else {
          hSig10->Draw("same hist");
	}

	cSig20->cd();
        calc.SetSigScale(20);
        TH1* hSig20 = NULL;
        hSig20 = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
        hSig20->SetLineColor(ColorShift(shift));
        hSig20->SetLineStyle(StyleShift(shift));
        hSig20->SetLineWidth(4);
        if (shift == 3) {
          hSig20->Draw("histo");
          hSig20->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
        }
        else {
          hSig20->Draw("same hist");
        }

	cSig100->cd();
        calc.SetSigScale(100);
        TH1* hSig100 = NULL;
        hSig100 = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
        hSig100->SetLineColor(ColorShift(shift));
        hSig100->SetLineStyle(StyleShift(shift));
        hSig100->SetLineWidth(4);
        if (shift == 3) {
          hSig100->Draw("histo");
          hSig100->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
        }
        else {
          hSig100->Draw("same hist");
        }
      }
      cBkg->Write();
      cSig1->Write();
      cSig10->Write();
      cSig20->Write();
      cSig100->Write();

      cBkg->SaveAs(pdfDir+Form("Bkg_%s.root", syst->ShortName().c_str()));
      cSig1->SaveAs(pdfDir+Form("Sig_1MeV_%s.root", syst->ShortName().c_str()));
      cSig10->SaveAs(pdfDir+Form("Sig_10MeV_%s.root", syst->ShortName().c_str()));
      cSig20->SaveAs(pdfDir+Form("Sig_20MeV__%s.root", syst->ShortName().c_str()));
      cSig100->SaveAs(pdfDir+Form("Sig_100MeV_%s.root", syst->ShortName().c_str()));

    }
    
    exFile->Write();
    exFile->Close();
    delete exFile;    
    
    sw.Stop();
    std::cout << "\nAll done in: ";
    sw.Print();
    return;
}

int ColorShift(int shift)
{
  int color;
  if (shift == 0) // nominal spectra
    {
      color = 1;
    }
  if (shift == 1 || shift == -1) // systematic shift by 1-sigma
    {
      color = 417;
    }
  if (shift == 2 || shift == -2) // systematic shift by 2-sigma
    {
      color = 600;
    }
  if (shift == 3 || shift == -3) // systematic shift by 3-sigma
    {
      color = 881;
    }
  return color;
}

int StyleShift(int shift)
{
  int style;
  if (shift == 0) // nominal spectra                                                                  
    {
      style = 1;
    }
  else
    {
      style = 2;
    }// shifted by systematic
  return style;
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
