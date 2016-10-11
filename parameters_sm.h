//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.0.1, 2014-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MATRIXELEMENT_Parameters_sm_H
#define MATRIXELEMENT_Parameters_sm_H

#include <complex>

class Parameters_sm {
  public:
    // Model parameters independent of aS
    double cw, sw, ee, MZ, WZ;
    // Model couplings independent of aS
    std::complex<double> GC_1, GC_2, GC_3, GC_50, GC_51, GC_58, GC_59;
    
    // Set parameters that are unchanged during the run
    void Set(double Gf, double aEWM1, double mz, double wz);
    
    double alpha;
    double MW = 80.419002;
    double WidthW = 2.085;
    
    double MH = 125.9;
    double WidthH = 4e-3;
    
    double MUQuark = 6.983e-2;
    double MDQuark = 6.983e-2;
    double MCQuark = 1.2;
    double MSQuark = 0.15;
    double MTQuark = 173.07;
    double MBquark = 4.6;
    
    double MElectron = 0.510998928e-3;
    double MMuon = 0.1056583715;
    double MTau = 1.77682;

};

struct Parameters_alphaS {
  // couplings
  std::complex<double> GC_11;
  // Model parameters dependent on aS
    double sqrt__aS, G, G__exp__2, aS;
    
    void Set(double alphaS);
  
};



#endif // Parameters_sm_H

