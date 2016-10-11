//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.14, 2013-11-27
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef uux_mupmum_H
#define uux_mupmum_H

#include <complex>

#include "phasespace.h"
#include "parameters_sm.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > mu+ mu- WEIGHTED=4
// Process: c c~ > mu+ mu- WEIGHTED=4
//--------------------------------------------------------------------------

/**
@brief matrix element for u u~ -> mu+ mu-
*/
class ME_uux_mupmum {
  public:
    // Constructor.
    ME_uux_mupmum() {
    }

    ~ME_uux_mupmum() {
    }

    // Calculate flavour-independent parts of cross section.
    double Calculate(const Phasespace::Phasespace &ps, int perm[],
                     const Parameters_sm &param,
                     const Parameters_alphaS &param_aS);

    double ColorCorrelated(int id1, int id2, int i, int j);

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 4;
    static const int nprocesses = 2;

  private:
    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[],
                                 const Parameters_sm &par,
                                 const Parameters_alphaS &param_aS);
    static const int nwavefuncs = 6;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 2;
    std::complex<double> amp[namplitudes];
    double matrix_uux_mupmum();

    // Store the matrix element value from sigmaKin
    double matrix_element;

    // Color flows, used when selecting color
    double jamp2;

    // vector with external particle masses
    double mME[4] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];
};

#endif // PP_MUPMUM_Sigma_sm_uux_mupmum_H
