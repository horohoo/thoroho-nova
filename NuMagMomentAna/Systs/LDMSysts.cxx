#include "NuMagMomentAna/Systs/LDMSysts.h"

namespace ana
{
  const NDPileupEffectSyst kNDPileupEffectSyst;

  const LDMFluxSyst kLDMFluxSyst;

  const DummySyst kDummySyst;
  
  const std::string kSystDir = FindPackageDir("NuMagMomentAna")+"/data/systs";
  
  const std::string kDmSystsFile    = kSystDir+"/nd_systs_dm.root";
  const std::string kNuoneSystsFile = kSystDir+"/nd_systs_nuone.root";
  const std::string kNumiSystsFile  = kSystDir+"/nd_systs_numi.root";

  const SystFromHist kNDldmCalibSyst   (kDmSystsFile,    "calib_dm",  "Calibration DM");
  const SystFromHist kNDldmLightSyst   (kDmSystsFile,    "ll_dm",     "Light Level DM");
  const SystFromHist kNDldmCherSyst    (kDmSystsFile,    "ckv_dm",    "Cherenkov DM");

  const SystFromHist kNDNuoneCalibSyst (kNuoneSystsFile, "calib_nuone",  "Calibration Nuone");
  const SystFromHist kNDNuoneLightSyst (kNuoneSystsFile, "ll_nuone",     "Light Level Nuone");
  const SystFromHist kNDNuoneCherSyst  (kNuoneSystsFile, "ckv_nuone",    "Cherenkov  Nuone");

  const SystFromHist kNDNumiCalibSyst  (kNumiSystsFile, "calib_numi", "Calibration Numi");
  const SystFromHist kNDNumiLightSyst  (kNumiSystsFile, "ll_numi",     "Light Level Numi");
  const SystFromHist kNDNumiCherSyst   (kNumiSystsFile, "ckv_numi",    "Cherenkov Numi");

  //----------------------------------------------------------------------
  SystFromHist::SystFromHist(const std::string &fname, const std::string &shortname, const std::string &latexname)
    : ISyst(shortname, latexname),
      fFileName(fname)
  {
      LoadHists(shortname);
  }
  
  //----------------------------------------------------------------------
  void SystFromHist::LoadHists(const std::string &shortname) const
  {
      if(!fHists.empty()) return;
      TFile* fin = TFile::Open(fFileName.c_str(), "read");
      if (fin->IsZombie()||!fin)
      {
          std::cout << "Warning: couldn't open " << fFileName << ".  Crashing" <<  std::endl;
          abort();
      }
      
      TDirectory* dir = (TDirectory*)fin->Get("systs");
      if (!dir)
      {
          std::cout << "Warning: couldn't open systs directory. Crashing" <<  std::endl;
          abort();
      }
           
      const std::vector<int> sigmas = {-1, +1};
      for (int i = 0; i < (int)sigmas.size(); i++)
      {
          int strend = shortname.find_first_of('_');
          std::string hName = i==0 ? Form("%sdown", shortname.substr(0, strend).c_str())  // -1 sigma
                                   : Form("%sup",   shortname.substr(0, strend).c_str()); // +1 sigma
                                   
          TH1D* h = (TH1D*)dir->Get(hName.c_str());
          if(!h)
          {
              std::cout << "Error: can't find necessary " << hName << " histogram in file " << fFileName << ".  Crashing... \n";
              abort();
          }
          h->SetDirectory(0);
          fHists.emplace_back(sigmas[i], h);
      }  
  }
  
  //----------------------------------------------------------------------
  SystFromHist::~SystFromHist()
  {
  }
  
  //----------------------------------------------------------------------
  double SystFromHist::WeightFor(double sigma, double etheta2) const
  {
      const int bin = fHists[0].second->FindBin(etheta2);
      double ret = 1.0;
      if(sigma < 0)
      {
          ret = -1 * sigma * fHists[0].second->GetBinContent(bin);  // -1 sigma
      }
      else if (sigma > 0)
      {
          ret = sigma * fHists[1].second->GetBinContent(bin);  // +1 sigma
      }
      
      return std::max(0., ret);
  }
  
  //----------------------------------------------------------------------
  void SystFromHist::Shift(double sigma, caf::SRProxy* sr, double& weight) const
  {
      if(!sr->vtx.elastic.IsValid || sr->vtx.elastic.fuzzyk.nshwlid < 1 || sr->vtx.elastic.fuzzyk.npng < 1)
      {
          weight *= 1.0;
      }
      else
      {
          float shwE = sr->vtx.elastic.fuzzyk.png[0].shwlid.calE;
          TVector3 prongDir = (TVector3)sr->vtx.elastic.fuzzyk.png[0].dir;
          TVector3 beamDir = ana::NuMIBeamDirection(sr->hdr.det);
          float ct = float(prongDir.Dot(beamDir));
    
          double etheta2 = shwE*acos(ct)*acos(ct);
    
          weight *= WeightFor(sigma, etheta2);
      }
  }

  //----------------------------------------------------------------------
  void NDPileupEffectSyst:: Shift(double sigma, caf::SRProxy* /*sr*/, double& weight) const
  {
    weight *= 1+.106*sigma;
  }

  //----------------------------------------------------------------------
  void LDMFluxSyst:: Shift(double sigma, caf::SRProxy* /*sr*/, double& weight) const
  {
    weight *= 1+0.06*sigma;
  }

  void DummySyst:: Shift(double sigma, caf::SRProxy* /*sr*/, double& weight) const
  {
    weight *= 1+0.0*sigma;
  }

}
