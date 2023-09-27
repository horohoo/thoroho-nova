#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/Cuts.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Fit/Fit.h"
#include "CAFAna/Systs/XSecSystLists.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Vars/Vars.h"
#include "CAFAna/Weights/PPFXWeights.h"
#include "CAFAna/Weights/XsecTunes.h"

#include "Utilities/rootlogon.C"

#include "StandardRecord/Proxy/SRProxy.h"


#include "TAxis.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TText.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVectorD.h"
#include "TVector3.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "iostream"
#include "string"
#include "dirent.h"


#include "NuMagMomentAna/Cuts/ldmCuts.h"
#include "NuMagMomentAna/OscCalc/OscCalcSingleElectron.h"
#include "NuMagMomentAna/Prediction/NDPredictionSingleElectron.h"
#include "NuMagMomentAna/Prediction/NDPredictionSystsSingleElectron.h"
#include "NuMagMomentAna/Systs/LDMSysts.h"
#include "NuMagMomentAna/Vars/NuoneVars.h"
#include "NuMagMomentAna/Vars/FitVarsSingleElectron.h"

//degree of freedom: 1, probability (confidence level): 0.9, chi-square critical value: 2.70554
const double Chi2for90CL = 2.70554;

//parameters in the BdNMC simulation for DM production
const float detEfficiency = 40.2/9.0;
const float potbdnmc      = 1e21;    
const float alpha         = 0.5;
const float mchimvratio   = 1.0/3.0;
const float epsilon4      = 1e-12;
const float constY        = alpha*pow(mchimvratio, 4)/detEfficiency;
const float bdnmcY        = sqrt(epsilon4)*constY;

//initial parameters for "oscillation calculator" to produce fake data
double ldmScale       = 0.00;

double fhcbkgscale    = 0.95;
double fhcnuonescale  = 1.05;

double rhcbkgscale    = 1.05;
double rhcnuonescale  = 0.90;

const std::string dmFile = "/nova/app/users/thoroho/ldmanalysis/NuMagMomentAna/data/ldmspectra/ldmone.root";

const double DMSampleNum = 1.25524e+06;
struct DmInf
{
    int DmMassPoints;
    int DMSimuNum;
    double DMSimuPOT;
    double DMExclC;
};

const std::vector<DmInf> dmInf = {{5,   73666,  2.03434e+17, 1.5e-6}, {6,   85346,  3.61850e+17, 2.1e-6}, {7,   94751,  5.94086e+17, 2.8e-6}
                                 ,{8,   102437, 8.92630e+17, 3.6e-6}, {9,   108635, 1.33293e+18, 4.8e-6}, {10,  113729, 1.87307e+18, 6.4e-6}
                                 ,{20,  138171, 2.37604e+19, 6.6e-5}, {30,  146560, 2.38331e+20, 7.1e-4}, {40,  149881, 9.44323e+21, 5.6e-2}
                                 ,{50,  151135, 3.62749e+22, 1.2e-1}, {60,  149583, 7.73356e+22, 2.6e-1}, {70,  149292, 1.55033e+23, 5.3e-1}
                                 ,{80,  148692, 3.13427e+23, 1.1e-0}, {90,  137096, 5.79043e+23, 2.2e-0}, {100, 127174, 9.74147e+23, 4.0e-0}
                                 ,{200, 142212, 7.46202e+25, 200}, {300, 144485, 1.18906e+26, 270}, {400, 145540, 1.57898e+27, 300}
                                 ,{450, 145831, 9.39244e+27, 300}
                                 };


//MC samples, including the signal (dm-on-e), the irreduicible background (nu-on-e), and other nominal background
const std::string fldm   = "/nova/ana/users/wmu/mcsample/dm_caf/*";
const std::string fnuone = "/nova/ana/users/wmu/mcsample/nuone_caf/*"; //1.34928e+23 POT
const std::string ffhc   = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_fhc_nova_v08_full_ndphysics_contain_v1"; //5.54495e+21 POT
const std::string frhc   = "prod_sumdecaf_development_nd_genie_N1810j0211a_nonswap_rhc_nova_v08_full_ndphysics_contain_v1"; //5.06177e+21 POT

//Test MC samples
//const std::string fldm   = "/pnfs/nova/persistent/users/wmu/test/ldmone_caf/*";
//const std::string fnuone = "/pnfs/nova/persistent/users/wmu/test/nuone_caf/*";
//const std::string ffhc   = "/pnfs/nova/persistent/users/wmu/test/nominal_caf/*";
//const std::string frhc   = "/pnfs/nova/persistent/users/wmu/test/nominal_caf/*";

TString preDir = "./prediction/";
TString fitDir = "./fits/";
TString excDir = "./exclusion/";

TString outDir = "./outputs/";
