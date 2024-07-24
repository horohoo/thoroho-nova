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


    std::vector<int> shifts = {+1, -1};
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
    const Cut mec_cut = MEC && kIsNueCC;
    const Cut numu_cut  = !nuone_cut && !mec_cut;

    // construct prediction (needed to LoadFrom() later)
    SpectrumLoader loaderldm(fldm);
    SpectrumLoader loadernuone(fnuone);
    SpectrumLoader loadernominal(ffhc);
    SpectrumLoader loadermec(fmec);

    const Binning bins  = Binning::Simple(20, 0, 0.02);
    const HistAxis etheta2Axis("E #theta^{2} (GeV Rad^{2})", bins, nuone::kETheta2);
    NDPredictionSystsSingleElectron pred(loaderldm, loadernuone, loadernominal, loadermec, etheta2Axis, sel_cut, sel_cut&&nuone_cut, sel_cut&&numu_cut, sel_cut&&mec_cut, ana::kPPFXFluxCVWgt*kradWt, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51, ana::kPPFXFluxCVWgt*ana::kXSecCVWgt2020GSFProd51);

    const double pot = 1.25e21;

    TFile *predFile = new TFile(outDir+"Prediction_DM_1MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root");
    predFile->cd();
    static std::unique_ptr<NDPredictionSystsSingleElectron> pred_ptr = pred.LoadFrom(predFile, pred_outname);

    TFile *predFile10 = new TFile(outDir+"Prediction_DM_10MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root");
    predFile10->cd();
    static std::unique_ptr<NDPredictionSystsSingleElectron> pred_ptr10 = pred.LoadFrom(predFile10, pred_outname);

    TFile *predFile100 = new TFile(outDir+"Prediction_DM_100MeV_nd_flux_pileup_xsec_sys_fhc_fake_data_cut_5.root");
    predFile100->cd();
    static std::unique_ptr<NDPredictionSystsSingleElectron> pred_ptr100 = pred.LoadFrom(predFile100, pred_outname);
    
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
      cBkg->SetBottomMargin(1e-5);
      cBkg->SetLeftMargin(0.15);
      gStyle->SetPadBorderMode(0);
      gStyle->SetFrameBorderMode(0);
      cBkg->Divide(1, 2, 1e-5, 1e-5);


      TCanvas* cSig1 = new TCanvas (Form("Sig_1MeV_%s", syst->ShortName().c_str()), Form("Sig_1MeV_%s", syst->ShortName()), 800, 600);
      cSig1->SetBottomMargin(1e-5);
      cSig1->SetLeftMargin(0.15);
      gStyle->SetPadBorderMode(0);
      gStyle->SetFrameBorderMode(0);
      cSig1->Divide(1, 2, 1e-5, 1e-5);


      TCanvas* cSig10 = new TCanvas (Form("Sig_10MeV__%s", syst->ShortName().c_str()), Form("Sig_10MeV__%s", syst->ShortName().c_str()), 800, 600);
      cSig10->SetBottomMargin(1e-5);
      cSig10->SetLeftMargin(0.15);
      gStyle->SetPadBorderMode(0);
      gStyle->SetFrameBorderMode(0);
      cSig10->Divide(1, 2, 1e-5, 1e-5);


      TCanvas* cSig100 = new TCanvas (Form("Sig_100MeV_%s", syst->ShortName().c_str()), Form("Sig_100MeV_%s", syst->ShortName().c_str()), 800, 600);
      cSig100->SetBottomMargin(1e-5);
      cSig100->SetLeftMargin(0.15);
      gStyle->SetPadBorderMode(0);
      gStyle->SetFrameBorderMode(0);
      cSig100->Divide(1, 2, 1e-5, 1e-5);

      SystShifts systshiftnom (syst, 0);

      // calculate nominal spectrum for ratio plots
      TH1* hIBkg_nom = NULL;
      TH1* hBkg_nom = NULL;
      TH1* hMEC_nom = NULL;

      hIBkg_nom = pred_ptr->PredictComponentSyst(&calc, systshiftnom, Flavors::kNuEToNuMu, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
      hBkg_nom = pred_ptr->PredictComponentSyst(&calc, systshiftnom, Flavors::kNuEToNuTau, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
      //hMEC_nom = pred_ptr->PredictComponentSyst(&calc, systshiftnom, Flavors::kNuMuToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
      
      TH1* hIBkg_nosyst = NULL;
      TH1* hBkg_nosyst = NULL;
      hBkg_nosyst = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuTau, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
      hIBkg_nosyst = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuMu, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);

      //hBkg_nom->Add(hMEC_nom);
      //hBkg_nom->Add(hIBkg_nom);

      cBkg->cd(1);
      gPad->SetBottomMargin(1e-5);
      gStyle->SetOptTitle(1);
      gPad->Modified();
      gPad->Update();

      hBkg_nom->SetLineColor(ColorShift(0));
      hBkg_nom->SetLineStyle(StyleShift(0));
      hBkg_nom->SetLineWidth(4);

      hBkg_nom->Draw("histo");

      TLegend lBkg(0.68, 0.68, 0.98, 0.98);
      lBkg.SetHeader(Form("Bkg_%s", syst->ShortName().c_str()));
      lBkg.AddEntry(hBkg_nom, "Nominal");

      cBkg->cd(2);
      gPad->SetTopMargin(1e-5);
      gPad->SetTickx();
      gStyle->SetOptTitle(0);
      gPad->Modified();
      gPad->Update();

      
      calc.SetDMMass(1);
      calc.SetSigScale(1e-5);
      TH1* hSig1_nom = NULL;
      hSig1_nom = pred_ptr->PredictComponentSyst(&calc, systshiftnom, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 2);

      TH1* hSig1_nosyst = NULL;
      hSig1_nosyst = pred_ptr->PredictComponent(&calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 2);

      std::cout << "Bkg->PredictComponentSyst(): " << hBkg_nom->Integral() << std::endl;
      std::cout << "Bkg->PredictComponent(): " << hBkg_nosyst->Integral() << std::endl;
      std::cout << "Nuone->PredictComponentSyst(): " << hIBkg_nom->Integral() << std::endl;
      std::cout << "Nuone->PredictComponent(): " << hIBkg_nosyst->Integral() << std::endl;
      std::cout << "Sig->PredictComponentSyst(): " << hSig1_nom->Integral() << std::endl;
      std::cout << "Sig->PredictComponent(): " << hSig1_nosyst->Integral() << std::endl;


      cSig1->cd(1);
      gPad->SetBottomMargin(1e-5);
      gStyle->SetOptTitle(1);
      gPad->Modified();
      gPad->Update();

      hSig1_nom->SetLineColor(ColorShift(0));
      hSig1_nom->SetLineStyle(StyleShift(0));
      hSig1_nom->SetLineWidth(4);

      hSig1_nom->Draw("histo");

      TLegend lSig1(0.68, 0.68, 0.98, 0.98);
      lSig1.SetHeader(Form("DM_1MeV_%s", syst->ShortName().c_str()));
      lSig1.AddEntry(hSig1_nom, "Nominal");

      cSig1->cd(2);
      gPad->SetTopMargin(1e-5);
      gPad->SetTickx();
      gStyle->SetOptTitle(0);
      gPad->Modified();
      gPad->Update();

      calc.SetDMMass(10);
      calc.SetSigScale(1e-5);
      TH1* hSig10_nom = NULL;
      hSig10_nom = pred_ptr10->PredictComponentSyst(&calc, systshiftnom, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 2);

      cSig10->cd(1);
      gPad->SetBottomMargin(1e-5);
      gStyle->SetOptTitle(1);
      gPad->Modified();
      gPad->Update();

      hSig10_nom->SetLineColor(ColorShift(0));
      hSig10_nom->SetLineStyle(StyleShift(0));
      hSig10_nom->SetLineWidth(4);

      hSig10_nom->Draw("histo");

      TLegend lSig10(0.68, 0.68, 0.98, 0.98);
      lSig10.SetHeader(Form("DM_10MeV_%s", syst->ShortName().c_str()));
      lSig10.AddEntry(hSig10_nom,"Nominal");

      cSig10->cd(2);
      gPad->SetTopMargin(1e-5);
      gPad->SetTickx();
      gStyle->SetOptTitle(0);
      gPad->Modified();
      gPad->Update();

      calc.SetDMMass(100);
      calc.SetSigScale(1.0);
      TH1* hSig100_nom = NULL;
      hSig100_nom = pred_ptr100->PredictComponentSyst(&calc, systshiftnom, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 2);

      cSig100->cd(1);
      gPad->SetBottomMargin(1e-5);
      gStyle->SetOptTitle(1);
      gPad->Modified();
      gPad->Update();

      hSig100_nom->SetLineColor(ColorShift(0));
      hSig100_nom->SetLineStyle(StyleShift(0));
      hSig100_nom->SetLineWidth(4);

      hSig100_nom->Draw("histo");

      TLegend lSig100(0.68, 0.68, 0.98, 0.98);
      lSig100.SetHeader(Form("DM_100MeV_%s", syst->ShortName().c_str()));
      lSig100.AddEntry(hSig100_nom,"Nominal");

      cSig100->cd(2);
      gPad->SetTopMargin(1e-5);
      gPad->SetTickx();
      gStyle->SetOptTitle(0);
      gPad->Modified();
      gPad->Update();

      TH1* hRatioBkgUp = (TH1*)hBkg_nom->Clone("hRatioBkg_+1sigma");
      TH1* hRatioBkgDown = (TH1*)hBkg_nom->Clone("hRatioBkg_-1sigma");

      TH1* hRatioSig1Up = (TH1*)hSig1_nom->Clone("hRatioSig1MeV_+1sigma");
      TH1* hRatioSig1Down = (TH1*)hSig1_nom->Clone("hRatioSig1MeV_-1sigma");

      TH1* hRatioSig10Up = (TH1*)hSig10_nom->Clone("hRatioSig10MeV_+1sigma");
      TH1* hRatioSig10Down = (TH1*)hSig10_nom->Clone("hRatioSig10MeV_-1sigma");

      TH1* hRatioSig100Up = (TH1*)hSig100_nom->Clone("hRatioSig100MeV_+1sigma");
      TH1* hRatioSig100Down = (TH1*)hSig100_nom->Clone("hRatioSig100MeV_-1sigma");

      double maxBkgUp;
      double maxBkgDown;
      
      double maxSig1Up;
      double maxSig1Down;

      double maxSig10Up;
      double maxSig10Down;

      double maxSig100Up;
      double maxSig100Down;

      std::cout << syst->ShortName() << std::endl;

      std::cout << "Nominal background: " << hBkg_nom->Integral() << std::endl;
      //std::cout << "Nominal signal 1 MeV: " << hSig1_nom->Integral() << std::endl;
      std::cout << "Nominal signal 10 MeV: " << hSig10_nom->Integral() <<std::endl;
      //std::cout << "Nominal signal 100 MeV: " << hSig100_nom->Integral() <<std::endl;


      for (int shift : shifts) {

	SystShifts systshift (syst, shift);
	
	TH1* hIBkg = NULL;
	TH1* hBkg  = NULL;
	TH1* hMEC = NULL;

        hIBkg = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuMu, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
        hBkg = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuTau, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
	hMEC = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuMuToNuE, Current::kCC, Sign::kNu).ToTH1(pot, 1, 5);
        
	hBkg->Add(hMEC);
        hBkg->Add(hIBkg);

	cBkg->cd(1);

	hBkg->SetLineColor(ColorShift(shift));
	hBkg->SetLineStyle(StyleShift(shift));
	hBkg->SetLineWidth(4);

	hBkg->Draw("same hist");

	if (shift == 1) {
	  std::cout << "+1 sigma background: " << hBkg->Integral() << std::endl;
	  maxBkgUp = hBkg->GetMaximum();

	  hRatioBkgUp->Divide(hBkg);
	  hRatioBkgUp->SetLineColor(ColorShift(shift));
	  hRatioBkgUp->SetLineWidth(4);

	  lBkg.AddEntry(hBkg, "+1#sigma");
	}

	else {	  
	  std::cout << "-1 sigma background: " << hBkg->Integral() << std::endl;
	  maxBkgDown = hBkg->GetMaximum();

	  hRatioBkgDown->Divide(hBkg);
	  hRatioBkgDown->SetLineColor(ColorShift(shift));
	  hRatioBkgDown->SetLineStyle(StyleShift(shift));
	  hRatioBkgDown->SetLineWidth(4);

	  lBkg.AddEntry(hBkg, "-1#sigma");
	}

	
	cSig1->cd(1);
        calc.SetDMMass(1);
	calc.SetSigScale(1e-5);
        TH1* hSig1 = NULL;
        hSig1 = pred_ptr->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
	hSig1->SetLineColor(ColorShift(shift));
	hSig1->SetLineStyle(StyleShift(shift));
	hSig1->SetLineWidth(4);
	hSig1->Draw("same hist");

	if (shift == 1) {
	  //std::cout << "+1 sigma signal 1 MeV: " << hSig1->Integral() << std::endl;
	  hRatioSig1Up->Divide(hSig1);
	  hRatioSig1Up->SetLineColor(ColorShift(shift));
	  hRatioSig1Up->SetLineWidth(4);

	  maxSig1Up = hSig1->GetMaximum();

	  lSig1.AddEntry(hSig1, "+1#sigma");

	}
	else {
	  //std::cout << "-1 sigma signal 1 MeV: " << hSig1->Integral() << std::endl;
	  hRatioSig1Down->Divide(hSig1);
          hRatioSig1Down->SetLineColor(ColorShift(shift));
	  hRatioSig1Down->SetLineStyle(StyleShift(shift));
          hRatioSig1Down->SetLineWidth(4);

          maxSig1Down = hSig1->GetMaximum();

	  lSig1.AddEntry(hSig1, "-1#sigma");
	}

	cSig10->cd(1);
        calc.SetDMMass(10);
	calc.SetSigScale(1e-5);
        TH1* hSig10 = NULL;
        hSig10 = pred_ptr10->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
        hSig10->SetLineColor(ColorShift(shift));
        hSig10->SetLineStyle(StyleShift(shift));
        hSig10->SetLineWidth(4);
	hSig10->Draw("same hist");

	if (shift == 1) {
	  std::cout << "+1 sigma signal 10 MeV: " << hSig10->Integral() << std::endl;
          hRatioSig10Up->Divide(hSig10);
          hRatioSig10Up->SetLineColor(ColorShift(shift));
          hRatioSig10Up->SetLineWidth(4);

          maxSig10Up = hSig10->GetMaximum();

	  lSig10.AddEntry(hSig10,"+1#sigma");

        }
	else {
	  std::cout << "-1 sigma signal 10 MeV: " << hSig10->Integral() << std::endl;
          hRatioSig10Down->Divide(hSig10);
          hRatioSig10Down->SetLineColor(ColorShift(shift));
	  hRatioSig10Down->SetLineStyle(StyleShift(shift));
          hRatioSig10Down->SetLineWidth(4);

          maxSig10Down = hSig10->GetMaximum();

	  lSig10.AddEntry(hSig10,"-1#sigma");
	}


	cSig100->cd(1);
        calc.SetDMMass(100);
	calc.SetSigScale(1.0);
        TH1* hSig100 = NULL;
        hSig100 = pred_ptr100->PredictComponentSyst(&calc, systshift, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(pot, ColorShift(shift), 2);
        hSig100->SetLineColor(ColorShift(shift));
        hSig100->SetLineStyle(StyleShift(shift));
        hSig100->SetLineWidth(4);
	hSig100->Draw("same hist");

	if (shift == 1) {
	  //std::cout << "+1 sigma signal 100 MeV: " << hSig100->Integral() << std::endl;
          hRatioSig100Up->Divide(hSig100);
          hRatioSig100Up->SetLineColor(ColorShift(shift));
          hRatioSig100Up->SetLineWidth(4);

          maxSig100Up = hSig100->GetMaximum();
	  
	  lSig100.AddEntry(hSig100,"+1#sigma");
        }
	else {
	  //std::cout << "-1 sigma signal 100 MeV: " << hSig100->Integral() << std::endl;
          hRatioSig100Down->Divide(hSig100);
          hRatioSig100Down->SetLineColor(ColorShift(shift));
	  hRatioSig100Down->SetLineStyle(StyleShift(shift));
          hRatioSig100Down->SetLineWidth(4);

          maxSig100Down = hSig100->GetMaximum();

	  lSig100.AddEntry(hSig100,"-1#sigma");
	}

      }


      cBkg->cd(1);
      hBkg_nom->GetYaxis()->SetRangeUser(0., 1.05*std::max(maxBkgUp, maxBkgDown));

      lBkg.SetBorderSize(1);
      lBkg.SetFillColor(0);
      lBkg.SetTextFont(42);
      lBkg.SetTextSize(0.05);
      lBkg.Draw("same");

      cBkg->Update();


      cBkg->cd(2);
      hRatioBkgUp->GetYaxis()->SetRangeUser(0.99*std::min(hRatioBkgUp->GetMinimum(), hRatioBkgDown->GetMinimum()), 1.01*std::max(hRatioBkgUp->GetMaximum(), hRatioBkgDown->GetMaximum()));
      hRatioBkgUp->Draw("hist");
      hRatioBkgDown->Draw("same hist");
      hRatioBkgUp->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
      hRatioBkgUp->GetYaxis()->SetTitle("Ratio");

      cSig1->cd(1);
      hSig1_nom->GetYaxis()->SetRangeUser(0., 1.05*std::max(maxSig1Up, maxSig1Down));

      lSig1.SetBorderSize(1);
      lSig1.SetFillColor(0);
      lSig1.SetTextFont(42);
      lSig1.SetTextSize(0.05);
      lSig1.Draw("same");

      cSig1->Update();

      cSig1->cd(2);
      hRatioSig1Up->GetYaxis()->SetRangeUser(0.99*std::min(hRatioSig1Up->GetMinimum(), hRatioSig1Down->GetMinimum()), 1.01*std::max(hRatioSig1Up->GetMaximum(), hRatioSig1Down->GetMaximum()));
      hRatioSig1Up->Draw("hist");
      hRatioSig1Down->Draw("same hist");
      hRatioSig1Up->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
      hRatioSig1Up->GetYaxis()->SetTitle("Ratio");

      cSig10->cd(1);
      hSig10_nom->GetYaxis()->SetRangeUser(0., 1.05*std::max(maxSig10Up, maxSig10Down));

      lSig10.SetBorderSize(1);
      lSig10.SetFillColor(0);
      lSig10.SetTextFont(42);
      lSig10.SetTextSize(0.05);
      lSig10.Draw("same");

      cSig10->Update();

      cSig10->cd(2);
      hRatioSig10Up->GetYaxis()->SetRangeUser(0.99*std::min(hRatioSig10Up->GetMinimum(), hRatioSig10Down->GetMinimum()), 1.01*std::max(hRatioSig10Up->GetMaximum(), hRatioSig10Down->GetMaximum()));
      hRatioSig10Up->Draw("hist");
      hRatioSig10Down->Draw("same hist");
      hRatioSig10Up->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
      hRatioSig10Up->GetYaxis()->SetTitle("Ratio");

      cSig100->cd(1);
      hSig100_nom->GetYaxis()->SetRangeUser(0., 1.05*std::max(maxSig100Up, maxSig100Down));

      lSig100.SetBorderSize(1);
      lSig100.SetFillColor(0);
      lSig100.SetTextFont(42);
      lSig100.SetTextSize(0.05);
      lSig100.Draw("same");

      cSig100->Update();

      cSig100->cd(2);
      hRatioSig100Up->GetYaxis()->SetRangeUser(0.99*std::min(hRatioSig100Up->GetMinimum(), hRatioSig100Down->GetMinimum()), 1.01*std::max(hRatioSig100Up->GetMaximum(), hRatioSig100Down->GetMaximum()));
      hRatioSig100Up->Draw("hist");
      hRatioSig100Down->Draw("same hist");
      hRatioSig100Up->GetXaxis()->SetTitle("E#theta^{2}(GeV Rad^{2})");
      hRatioSig100Up->GetYaxis()->SetTitle("Ratio");
      

      cBkg->Write();
      cSig1->Write();
      cSig10->Write();
      cSig100->Write();

      cBkg->SaveAs(pdfDir+Form("Bkg_%s.pdf", syst->ShortName().c_str()));
      cSig1->SaveAs(pdfDir+Form("Sig_1MeV_%s.pdf", syst->ShortName().c_str()));
      cSig10->SaveAs(pdfDir+Form("Sig_10MeV_%s.pdf", syst->ShortName().c_str()));
      cSig100->SaveAs(pdfDir+Form("Sig_100MeV_%s.pdf", syst->ShortName().c_str()));

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
  if (shift == 1) // systematic shift by 1-sigma
    {
      color = 600;
    }
  if (shift == -1)
    {
      color = 632;
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
  else // shifted by systematic
    {
      style = 2;
    }
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
