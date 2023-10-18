#include "CAFAna/Core/Spectrum.h"

namespace ana
{
  class NDSpectrum: public Spectrum
  {
  public:
    friend class NDPredictionSingleElectron;
   
    Hist GetHist() {return fHist;};
  };
}
