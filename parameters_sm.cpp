//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.0.1, 2014-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "parameters_sm.h"

using namespace std;

void Parameters_sm::Set(double Gf, double aEWM1, double mz, double wz) {
    WZ = wz;
    MZ = mz;
    double MZ__exp__2 = pow(MZ, 2.);
    double aEW = 1. / aEWM1;
    alpha = aEW;
    ee = 2. * sqrt(aEW) * sqrt(M_PI);
    double MW__exp__2 = pow(MW, 2.);
    double sw2 = 1. - MW__exp__2 / MZ__exp__2;
    cw = sqrt(1. - sw2);
    sw = sqrt(sw2);

    GC_1 = std::complex<double>(0.0, -ee*(1./ 3.));
    GC_2 = std::complex<double>(0.0, ee*(2./ 3.));
    GC_3 = std::complex<double>(0.0,-ee);
    GC_50 = std::complex<double>(0.0, -(cw * ee) / (2. * sw));
    GC_51 = std::complex<double>(0.0, (cw * ee) / (2. * sw));
    GC_58 = std::complex<double>(0.0, -(ee * sw) / (6. * cw));
    GC_59 = std::complex<double>(0.0, (ee * sw) / (2. * cw));
}

void Parameters_alphaS::Set(double alphaS) {
  aS = alphaS;
  sqrt__aS = sqrt(aS);
    G = 2. * sqrt__aS * sqrt(M_PI);
    G__exp__2 = pow(G, 2.);
    GC_11 = std::complex<double>(0.0, G);
}


