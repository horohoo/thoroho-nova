#include "CAFAna/Core/Cut.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "NuMagMomentAna/Cuts/SinglesCVNCuts.h"

namespace nuone
{
  const ana::Cut kSinglesCVNFiducialCut([](const caf::SRProxy* sr){
      if(!sr->vtx.elastic.IsValid) return false;
      return (sr->vtx.elastic.vtx.X() < kSinglesCVNFiducialMax.X() &&
              sr->vtx.elastic.vtx.X() > kSinglesCVNFiducialMin.X() &&
              sr->vtx.elastic.vtx.Y() < kSinglesCVNFiducialMax.Y() &&
              sr->vtx.elastic.vtx.Y() > kSinglesCVNFiducialMin.Y() &&
              sr->vtx.elastic.vtx.Z() < kSinglesCVNFiducialMax.Z() &&
              sr->vtx.elastic.vtx.Z() > kSinglesCVNFiducialMin.Z());
      return true;
    });

  const ana::NuTruthCut kSinglesCVNFiducialCut_NT([](const caf::SRNeutrinoProxy* sr){
      return (sr->vtx.X() < kSinglesCVNFiducialMax.X() &&
              sr->vtx.X() > kSinglesCVNFiducialMin.X() &&
              sr->vtx.Y() < kSinglesCVNFiducialMax.Y() &&
              sr->vtx.Y() > kSinglesCVNFiducialMin.Y() &&
              sr->vtx.Z() < kSinglesCVNFiducialMax.Z() &&
              sr->vtx.Z() > kSinglesCVNFiducialMin.Z());
    });

  const ana::Cut kSinglesCVNContainmentCut([](const caf::SRProxy* sr){
      if(sr->vtx.elastic.fuzzyk.nshwlid < 1) return false;
      TVector3 start = sr->vtx.elastic.fuzzyk.png[0].shwlid.start;
      TVector3 stop  = sr->vtx.elastic.fuzzyk.png[0].shwlid.stop;

      //Make sure none of the event goes into the muon catcher
      for(uint ix = 0; ix < sr->vtx.elastic.fuzzyk.nshwlid; ix++){
        TVector3 start_muoncatcher =
          sr->vtx.elastic.fuzzyk.png[ix].shwlid.start;
        TVector3 stop_muoncatcher  =
          sr->vtx.elastic.fuzzyk.png[ix].shwlid.stop;
        if (std::max(start_muoncatcher.Z(),
                     stop_muoncatcher.Z()) >  1250.0) return false;
      }

      if(sr->sel.nuecosrej.distallpngtop < 50) return false;
      if(sr->sel.nuecosrej.distallpngbottom < 30) return false;
      if(sr->sel.nuecosrej.distallpngeast < 50) return false;
      if(sr->sel.nuecosrej.distallpngwest < 30) return false;
      if(sr->sel.nuecosrej.distallpngfront < 150) return false;
      return true;
    });
}
