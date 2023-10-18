#include "NuMagMomentAna/Vars/NuoneWeights.h"
#include "StandardRecord/Proxy/SRProxy.h"

namespace nuone
{
  const ana::NuTruthWeight kNuMM_NT([](const caf::SRNeutrinoProxy* nu)
    {
      double MMWeight = 1.0;
      // do nothing if it's not a true nu-on-e events
      if(nu->inttype != 1098) return 1.0;

      // get true energy of the electron
      double ElE = 0.;
      for(int i = 0; i < (int)nu->prim.size(); ++i){
        if( nu->prim[i].pdg == 11 ){
          ElE = nu->prim[i].p.E;
          break;
        }
      }
      double NuE = nu->E; // true neutrino energy

      double m_e = 0.51099895*pow(10,-3); // electron mass [GeV]
      double T_e = ElE - m_e; // electron kinetic energy
      if(T_e > (2*pow(NuE,2)/(2*NuE + m_e))) return 1.0; // upper limit for T_e
      double PI = 3.14159265; // hardcode PI for easy calculations

      // calculate the constant in front of the standard model (SM) term
      double G_F = 1.166379; // Fermi constant from PDG
      double SMConst = pow(G_F,2)*pow(10,-10)*m_e / (2*PI);

      // calculate the standard model interaction constants based on
      // type fo particle and chirality
      double g_V = 2*0.2313 - 0.5;
      double g_A = -0.5;
      if(nu->pdg == +12){
        g_V = 2*0.2313 + 0.5;
        g_A = +0.5;
      }else if(nu->pdg == -12) g_V = 2*0.2313 +0.5;
      else if(nu->pdg == -14) g_A = +0.5;

      // calculate the standard model cross section
      double SMXSec = SMConst*(pow(g_V+g_A,2)+pow(g_V-g_A,2)*pow(1-(T_e/NuE),2)+
                               (pow(g_A,2)-pow(g_V,2))*m_e*T_e/pow(NuE,2));

      // hardcode neutrino magnetic moment (in units of bohr magneton)
      // TO DO: make this an input variable!
      double NuMM = pow(10,-9);

      double alpha = 0.00729735257; // fine structure constant
      // calculate constant in front of the neutrino magnetic moment (nuMM) term
      double NuMMConst = PI*pow(alpha,2) / pow(m_e,2);
      // calculate the neutrino magnetic moment cross section
      double NuMMXSec = NuMMConst*((1/T_e)-(1/NuE))*pow(NuMM,2);

      // final weight is a multiplication term to the SM cross section
      MMWeight = 1. + NuMMXSec/SMXSec;
      return MMWeight;
    });
}
