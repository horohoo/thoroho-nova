#include "CAFAna/Core/Cut.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "NuMagMomentAna/Cuts/EventCVNCuts.h"

namespace nuone
{
  const ana::Cut kEventCVNFiducialCut([](const caf::SRProxy* sr){
      if(!sr->vtx.elastic.IsValid) return false;
      return (sr->vtx.elastic.vtx.X() < kEventCVNFiducialMax.X() &&
              sr->vtx.elastic.vtx.X() > kEventCVNFiducialMin.X() &&
              sr->vtx.elastic.vtx.Y() < kEventCVNFiducialMax.Y() &&
              sr->vtx.elastic.vtx.Y() > kEventCVNFiducialMin.Y() &&
              sr->vtx.elastic.vtx.Z() < kEventCVNFiducialMax.Z() &&
              sr->vtx.elastic.vtx.Z() > kEventCVNFiducialMin.Z());
      return true;
    });

  const ana::NuTruthCut kEventCVNFiducialCut_NT([](const caf::SRNeutrinoProxy* sr){
      return (sr->vtx.X() < kEventCVNFiducialMax.X() &&
              sr->vtx.X() > kEventCVNFiducialMin.X() &&
              sr->vtx.Y() < kEventCVNFiducialMax.Y() &&
              sr->vtx.Y() > kEventCVNFiducialMin.Y() &&
              sr->vtx.Z() < kEventCVNFiducialMax.Z() &&
              sr->vtx.Z() > kEventCVNFiducialMin.Z());
    });
}
