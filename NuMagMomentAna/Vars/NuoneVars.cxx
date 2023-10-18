#include "NuMagMomentAna/Vars/NuoneVars.h"
#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"

namespace nuone
{
  const ana::Var kNuVtxX([](const caf::SRProxy* sr)
                         {
                           if(sr->mc.nnu==0) return -2000.f;;
                           return float(sr->mc.nu[0].vtx.x);
                         });
  const ana::Var kNuVtxY([](const caf::SRProxy* sr)
                         {
                           if(sr->mc.nnu==0) return -2000.f;;
                           return float(sr->mc.nu[0].vtx.y);
                         });
  const ana::Var kNuVtxZ([](const caf::SRProxy* sr)
                         {
                           if(sr->mc.nnu==0) return 3000.f;;
                           return float(sr->mc.nu[0].vtx.z);
                         });
  const ana::Var kNuIntType([](const caf::SRProxy* sr)
                            {
                              if(sr->mc.nnu==0) return 950;
                              return int(sr->mc.nu[0].inttype);
                            });
  const ana::Var kNuPdgorig([](const caf::SRProxy* sr)
                            {
                              if(sr->mc.nnu==0) return -1;
                              return int(sr->mc.nu[0].pdgorig);
                            });
  const ana::Var kNuTgtZ([](const caf::SRProxy* sr)
                         {
                           if(sr->mc.nnu==0) return -1;
                           return int(sr->mc.nu[0].tgtZ);
                         });
  const ana::Var kNuTgtA([](const caf::SRProxy* sr)
                         {
                           if(sr->mc.nnu==0) return -1;
                           return int(sr->mc.nu[0].tgtA);
                         });
  const ana::Var kNuX([](const caf::SRProxy* sr)
                      {
                        if(sr->mc.nnu==0) return -2.0f;
                        return float(sr->mc.nu[0].x);
                      });
  const ana::Var kNuY([](const caf::SRProxy* sr)
                      {
                        if(sr->mc.nnu==0) return -2.0f;
                        return float(sr->mc.nu[0].y);
                      });
  const ana::Var kPrimPDG([](const caf::SRProxy* sr)
                          {
                            if(!sr->vtx.elastic.IsValid) return 0;
                            if(sr->vtx.elastic.fuzzyk.npng < 1) return 0;
                            return int(sr->vtx.elastic.fuzzyk.png[0].truth.pdg);
                          });
  const ana::Var kMotherPrimPDG([](const caf::SRProxy* sr)
                                {
                                  if(!sr->vtx.elastic.IsValid) return 0;
                                  if(sr->vtx.elastic.fuzzyk.npng < 1) return 0;
                                  return int(sr->vtx.elastic.fuzzyk.png[0].truth.motherpdg);
                                });
  const ana::Var kNuSlcEff([](const caf::SRProxy* sr)
                           {
                             if(sr->mc.nnu==0) return -2.0f;
                             return float(sr->mc.nu[0].eff);
                           });
  const ana::Var kNuonESlcEff_OV([](const caf::SRProxy* sr)
                                 {
                                   for(int i=0; i<(int)sr->mc.nallnus; i++)
                                     if(sr->mc.allnus[i].mode > 9999 ||
                                        !sr->spill.ismc)
                                       return float(sr->mc.allnus[i].eff);
                                   return -5.f;
                                 });
  const ana::Var kNuonESlcEff([](const caf::SRProxy* sr)
                              {
                                for(int i=0; i<(int)sr->mc.nallnus; i++)
                                  if(sr->mc.allnus[i].inttype == 1098)
                                    return float(sr->mc.allnus[i].eff);
                                return -5.f;
                              });
                              

  extern const ana::Var kintType([](const caf::SRProxy* sr)
                              {return (sr->mc.nnu == 0) ? -1 : int(sr->mc.nu[0].inttype);});
                              
  extern const ana::Var kisVtxCont([](const caf::SRProxy* sr)
                              {return (sr->mc.nnu == 0) ? -1 : int(sr->mc.nu[0].isvtxcont);}); 

                              
  const ana::Var kVtxX([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -999.f;
                         return float(sr->vtx.elastic.vtx.x);
                       });
  const ana::Var kVtxY([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -999.f;
                         return float(sr->vtx.elastic.vtx.y);
                       });
  const ana::Var kVtxZ([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -999.f;
                         return float(sr->vtx.elastic.vtx.z);
                       });
  const ana::Var kNPng([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -1;
                         return int(sr->vtx.elastic.fuzzyk.npng);
                       });
  const ana::Var kNPng2D([](const caf::SRProxy* sr)
                         {
                           if(!sr->vtx.elastic.IsValid) return -1;
                           return int(sr->vtx.elastic.fuzzyk.npng2d);
                         });
  const ana::Var kShowLen([](const caf::SRProxy* sr)
                          {
                            if(!sr->vtx.elastic.IsValid) return -1000.f;
                            if(sr->vtx.elastic.fuzzyk.npng < 1) return -500.f;
                            return float(sr->vtx.elastic.fuzzyk.png[0].shwlid.len);
                          });
  const ana::Var kProngLen([](const caf::SRProxy* sr)
                           {
                             if(!sr->vtx.elastic.IsValid) return -1000.f;
                             if(sr->vtx.elastic.fuzzyk.npng < 1) return -500.f;
                             return float(sr->vtx.elastic.fuzzyk.png[0].len);
                           });
  const ana::Var kNPlane([](const caf::SRProxy* sr)
                         {
                           if(!sr->vtx.elastic.IsValid) return -10;
                           if(sr->vtx.elastic.fuzzyk.npng==0) return -5;
                           auto idx = sr->vtx.elastic.fuzzyk.longestidx;
                           return int(sr->vtx.elastic.fuzzyk.png[idx].nplane);
                         });
  const ana::Var kNCell([](const caf::SRProxy* sr)
                        {
                          if(!sr->vtx.elastic.IsValid) return -40;
                          if(sr->vtx.elastic.fuzzyk.npng == 0) return -20;
                          int ncell = 0;
                          for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                            ncell += sr->vtx.elastic.fuzzyk.png[i].shwlid.nhit;
                          }
                          return ncell;
                        });
  const ana::Var kMinX([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return 490.f;
                         if(sr->vtx.elastic.fuzzyk.npng == 0) return 490.f;
                         float _minx = 490.f;
                         for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                           auto startX = sr->vtx.elastic.fuzzyk.png[i].shwlid.start.x;
                           auto stopX = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop.x;
                           float _f_min_x = std::min(startX,stopX);
                           _minx = std::min(_minx,_f_min_x);
                         }
                         return _minx;
                       });
  const ana::Var kMaxX([](const caf::SRProxy* sr)
                       {
                         if (!sr->vtx.elastic.IsValid) return -490.f;
                         if (sr->vtx.elastic.fuzzyk.npng == 0) return -490.f;
                         float _maxx = -490.f;
                         for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                           auto startX = sr->vtx.elastic.fuzzyk.png[i].shwlid.start.x;
                           auto stopX = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop.x;
                           float _f_max_x = std::max(startX,stopX);
                           _maxx = std::max(_maxx,_f_max_x);
                         }
                         return _maxx;
                       });
  const ana::Var kMinY([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return 490.f;
                         if(sr->vtx.elastic.fuzzyk.npng == 0) return 490.f;
                         float _miny = 490.f;
                         for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                           auto startY = sr->vtx.elastic.fuzzyk.png[i].shwlid.start.y;
                           auto stopY = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop.y;
                           float _f_min_y = std::min(startY,stopY);
                           _miny = std::min(_miny,_f_min_y);
                         }
                         return _miny;
                       });
  const ana::Var kMaxY([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -490.f;
                         if(sr->vtx.elastic.fuzzyk.npng == 0) return -490.f;
                         float _maxy = -490.;
                         for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                           auto startY = sr->vtx.elastic.fuzzyk.png[i].shwlid.start.y;
                           auto stopY = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop.y;
                           float _f_max_y = std::max(startY,stopY);
                           _maxy = std::max(_maxy,_f_max_y);
                         }
                         return _maxy;
                       });
  const ana::Var kMinZ([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return 1790.f;
                         if(sr->vtx.elastic.fuzzyk.npng == 0) return 1790.f;
                         float _minz = 1790.;
                         for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                           auto startZ = sr->vtx.elastic.fuzzyk.png[i].shwlid.start.z;
                           auto stopZ = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop.z;
                           float _f_min_z = std::min(startZ,stopZ);
                           _minz = std::min(_minz,_f_min_z);
                         }
                         return _minz;
                       });
  const ana::Var kMaxZ([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -90.f;
                         if(sr->vtx.elastic.fuzzyk.npng == 0) return -90.f;
                         float _maxz = -90.;
                         for(unsigned int i=0; i<sr->vtx.elastic.fuzzyk.npng; ++i){
                           auto startZ = sr->vtx.elastic.fuzzyk.png[i].shwlid.start.z;
                           auto stopZ = sr->vtx.elastic.fuzzyk.png[i].shwlid.stop.z;
                           float _f_max_z = std::max(startZ,stopZ);
                           _maxz = std::max(_maxz,_f_max_z);
                         }
                         return _maxz;
                       });
  const ana::Var ShowStartRecoX([](const caf::SRProxy* sr)->float
                                {
                                  if(sr->vtx.elastic.IsValid == false) return -90.f;
                                  TVector3 start = sr->vtx.elastic.fuzzyk.png[0].shwlid.start;
                                  return (float)start.X();
                                });
  const ana::Var ShowStartRecoY([](const caf::SRProxy* sr)->float
                                {
                                  if(sr->vtx.elastic.IsValid == false) return -90.f;
                                  TVector3 start = sr->vtx.elastic.fuzzyk.png[0].shwlid.start;
                                  return (float)start.Y();
                                });
  const ana::Var ShowStartRecoZ([](const caf::SRProxy* sr)->float
                                {
                                  if(sr->vtx.elastic.IsValid == false) return -90.f;
                                  TVector3 start = sr->vtx.elastic.fuzzyk.png[0].shwlid.start;
                                  return (float)start.Z();
                                });
  const ana::Var kShwEFrac([](const caf::SRProxy* sr)
                           {
                             if(!sr->vtx.elastic.IsValid) return -0.1f;
                             if(sr->vtx.elastic.fuzzyk.nshwlid < 1) return -0.1f;
                             return float(sr->vtx.elastic.fuzzyk.png[0].shwlid.lid.shwEFrac);
                           });
  const ana::Var kEvtx([](const caf::SRProxy* sr)
                       {
                         if(!sr->vtx.elastic.IsValid) return -1.f;
                         if(sr->vtx.elastic.fuzzyk.nshwlid < 1) return -1.f;
                         return float(sr->vtx.elastic.fuzzyk.png[0].shwlid.lid.vtxgev);
                       });
  const ana::Var kGap([](const caf::SRProxy* sr)
                      {
                        if(!sr->vtx.elastic.IsValid) return -5.f;
                        if(sr->vtx.elastic.fuzzyk.nshwlid < 1) return -5.f;
                        return float(sr->vtx.elastic.fuzzyk.png[0].shwlid.lid.gap);
                      });
  const ana::Var kShowerCalE([](const caf::SRProxy* sr)
                             {
                               if(!sr->vtx.elastic.IsValid) return -2.0f;
                               if(sr->vtx.elastic.fuzzyk.npng < 1) return -2.0f;
                               return float(sr->vtx.elastic.fuzzyk.png[0].shwlid.calE);
                             });
  const ana::Var kShowerE([](const caf::SRProxy* sr)
                          {
                            if(!sr->vtx.elastic.IsValid) return -2.0f;
                            if(sr->vtx.elastic.fuzzyk.npng < 1) return -2.0f;
                            return float(sr->vtx.elastic.fuzzyk.png[0].shwlid.shwE);
                          });
  
  const ana::Var kNuoneID = SIMPLEVAR(sel.nuone.nuoneid);
  const ana::Var kEpi0ID = SIMPLEVAR(sel.nuone.epi0nuoneid);
  
  const ana::Var kETheta2([](const caf::SRProxy* sr)
                          {
                            if(!sr->vtx.elastic.IsValid) return -0.1f;
                            if(sr->vtx.elastic.fuzzyk.nshwlid < 1) return -0.1f;
                            if(sr->vtx.elastic.fuzzyk.npng < 1) return -0.1f;
                            float shwE = sr->vtx.elastic.fuzzyk.png[0].shwlid.calE;
                            TVector3 prongDir = (TVector3)sr->vtx.elastic.fuzzyk.png[0].dir;
                            TVector3 beamDir = ana::NuMIBeamDirection(sr->hdr.det);
                            float ct = float(prongDir.Dot(beamDir));
                            if(fabs(ct)>1) return -0.1f;
                            return shwE*acos(ct)*acos(ct);
                          });
  const ana::Var kTheta([](const caf::SRProxy* sr)
                        {
                          if(!sr->vtx.elastic.IsValid) return -0.5f;
                          if(sr->vtx.elastic.fuzzyk.nshwlid < 1) return -0.5f;
                          if(sr->vtx.elastic.fuzzyk.npng < 1) return -0.5f;
                          TVector3 prongDir = (TVector3)sr->vtx.elastic.fuzzyk.png[0].dir;
                          TVector3 beamDir = ana::NuMIBeamDirection(sr->hdr.det);
                          float ct = float(prongDir.Dot(beamDir));
                          if(fabs(ct)>1) return -0.5f;
                          return acos(ct);
                        });
  const ana::Var kTrueElectronTheta([](const caf::SRProxy* sr)
                                    {
                                      if(sr->mc.nnu==0) return -0.5f;
                                      float angle = -0.5f;
                                      int nprims = sr->mc.nu[0].prim.size();
                                      for(int iprim=0; iprim<nprims; ++iprim){
                                        if(abs(sr->mc.nu[0].prim[iprim].pdg)==11){
                                          TVector3 eDir = (TVector3)sr->mc.nu[0].prim[iprim].p.Vect().Unit();
                                          TVector3 beamDir = ana::NuMIBeamDirection(sr->hdr.det);
                                          float ct = float(eDir.Dot(beamDir));
                                          if(fabs(ct)>1) angle = -0.5f;
                                          angle = acos(ct);
                                        }
                                      }
                                      return angle;
                                    });
  const ana::Var kTrueElectronE([](const caf::SRProxy* sr)
                                {
                                  if(sr->mc.nnu==0) return -0.5f;
                                  float eE = -0.5f;
                                  int nprims = sr->mc.nu[0].prim.size();
                                  for(int iprim=0; iprim<nprims; ++iprim){
                                    if(abs(sr->mc.nu[0].prim[iprim].pdg)==11){
                                      eE = sr->mc.nu[0].prim[iprim].p.E;
                                    }
                                  }
                                  return eE;
                                });
  const ana::Var kProng3DvertexEnergyVol10([](const caf::SRProxy* sr)
                                           {
                                             if(!sr->vtx.elastic.IsValid) return -1000.f;
                                             if(sr->vtx.elastic.fuzzyk.npng<1) return -1000.f;
                                             return float(sr->vtx.elastic.prong3dvertexenergyvolume10);
                                           });
  const ana::Var kProngCVNEle([](const caf::SRProxy* sr)
                              {
                                if(!sr->vtx.elastic.IsValid) return -1000.f;
                                if(sr->vtx.elastic.fuzzyk.npng<1) return -1000.f;
                                int nProngs = sr->vtx.elastic.fuzzyk.npng;
                                float maxProngCVNe = -5.;
                                for(int iProng=0; iProng<nProngs; iProng++){
                                  if(sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.electronid>maxProngCVNe){
                                    maxProngCVNe = sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.electronid;
                                  }
                                }
                                return maxProngCVNe;
                              });
  const ana::Var kProngCVNMuon([](const caf::SRProxy* sr)
                               {
                                 if(!sr->vtx.elastic.IsValid) return -1000.f;
                                 if(sr->vtx.elastic.fuzzyk.npng<1) return -1000.f;
                                 int nProngs = sr->vtx.elastic.fuzzyk.npng;
                                 float maxVal = -5.;
                                 for(int iProng=0; iProng<nProngs; iProng++){
                                   if(sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.muonid>maxVal){
                                     maxVal = sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.muonid;
                                   }
                                 }
                                 return maxVal;
                               });
  const ana::Var kProngCVNPhoton([](const caf::SRProxy* sr)
                                 {
                                   if(!sr->vtx.elastic.IsValid) return -1000.f;
                                   if(sr->vtx.elastic.fuzzyk.npng<1) return -1000.f;
                                   int nProngs = sr->vtx.elastic.fuzzyk.npng;
                                   float maxVal = -5.;
                                   for(int iProng=0; iProng<nProngs; iProng++){
                                     if(sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.photonid>maxVal){
                                       maxVal = sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.photonid;
                                     }
                                   }
                                   return maxVal;
                                 });
  const ana::Var kProngCVNPhotonPlusEle([](const caf::SRProxy* sr)
                                        {
                                          if(!sr->vtx.elastic.IsValid) return -1000.f;
                                          if(sr->vtx.elastic.fuzzyk.npng<1) return -1000.f;
                                          int nProngs = sr->vtx.elastic.fuzzyk.npng;
                                          float maxVal = -5.;
                                          float maxVal1 = -5;
                                          float maxVal2 = -5;
                                          for(int iProng=0; iProng<nProngs; iProng++){
                                            if(sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.photonid>maxVal1){
                                              maxVal1 = sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.photonid;
                                            }
                                          }
                                          for(int iProng=0; iProng<nProngs; iProng++){
                                            if(sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.electronid>maxVal2){
                                              maxVal2 = sr->vtx.elastic.fuzzyk.png[iProng].cvnpart.electronid;
                                            }
                                          }
                                          maxVal = maxVal1+maxVal2;
                                          return maxVal;
                                        });
  const ana::Var kCVNhadE([](const caf::SRProxy* sr)
                          {
                            if(!sr->vtx.elastic.IsValid) return -1.0;
                            if(sr->vtx.elastic.fuzzyk.npng<1) return -1.0;
                            double CVNem_CalE = 0;
                            for(const caf::SRFuzzyKProngProxy& png: sr->vtx.elastic.fuzzyk.png){
                              double png_CalE = png.shwlid.calE;
                              double emPID = ((double)png.cvnpart.photonid +
                                              (double)png.cvnpart.pizeroid +
                                              (double)png.cvnpart.electronid);
                              double haPID = ((double)png.cvnpart.protonid +
                                              (double)png.cvnpart.pionid +
                                              (double)png.cvnpart.neutronid +
                                              (double)png.cvnpart.otherid +
                                              (double)png.cvnpart.muonid );
                              if(emPID<=0) continue; // for no cvn scores
                              if(emPID>=haPID) CVNem_CalE += png_CalE;
                            } // png
                            double kCVNemE = CVNem_CalE * ana::CalibrationBugCorrectionFactor(sr->hdr);
                            const double calE = sr->slc.calE * ana::CalibrationBugCorrectionFactor(sr->hdr);
                            return std::max(calE-kCVNemE,0.);
                          });
  const ana::Var k5labelEleID([](const caf::SRProxy* sr)
                              {
                                if(!sr->vtx.elastic.IsValid) return -1000.f;
                                if(sr->vtx.elastic.fuzzyk.npng<1) return -1000.f;
                                int nProngs = sr->vtx.elastic.fuzzyk.npng;
                                float max5labEleID = -5.;
                                for(int iProng=0; iProng<nProngs; iProng++){
                                  if(sr->vtx.elastic.fuzzyk.png[iProng].spprongcvnpart5label.electronid>max5labEleID){
                                    max5labEleID = sr->vtx.elastic.fuzzyk.png[iProng].spprongcvnpart5label.electronid;
                                  }
                                }
                                return max5labEleID;
                              });


  const ana::Var kTrueElTheta([](const caf::SRProxy* sr)
                                    {
                                      return -(sr->mc.cosmic[0].zenith);
                                    });
  const ana::Var kTrueElE([](const caf::SRProxy* sr)
                          {
                            if(!sr->vtx.elastic.IsValid) return -2.0f;
                            if(sr->vtx.elastic.fuzzyk.npng < 1) return -2.0f;
                            return float(sr->vtx.elastic.fuzzyk.png[0].truth.p.E);
                          });
}
