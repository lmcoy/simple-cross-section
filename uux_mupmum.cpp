//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.14, 2013-11-27
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "uux_mupmum.h"

#include <complex>

#include "aloha/aloha.h"

using namespace std;

namespace {

void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
              double M3, double W3, complex<double> V3[]) {
    complex<double> cI = complex<double>(0., 1.);
    double P3[4];
    complex<double> denom;
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    denom = COUP / (pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) -
                    pow(P3[3], 2) - M3 * (M3 - cI * W3));
    V3[2] = denom * -cI *
            (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]);
    V3[3] = denom * -cI *
            (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3] * F2[4]);
    V3[4] = denom * -cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) +
                           cI * (F1[3] * F2[4] + F1[4] * F2[3]));
    V3[5] = denom * -cI *
            (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5] * F2[3]);
}

} // namespace

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: u u~ > mu+ mu- WEIGHTED=4
// Process: c c~ > mu+ mu- WEIGHTED=4

/**
@brief Evaluate |M|^2

You can get the result with a call of MatrixElement.

@param ps The phasespace which is used to evaluate the matrix element.
*/
double ME_uux_mupmum::Calculate(const Phasespace::Phasespace &ps, int perm[],
                                const Parameters_sm &param,
                                const Parameters_alphaS &param_aS) {
    for (int i = 0; i < 4; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }
    // Reset color flows
    jamp2 = 0.;

    // Local variables and constants
    const int ncomb = 16;
    static bool goodhel[ncomb] = { ncomb * false };
    static int ntry = 0, sum_hel = 0, ngood = 0;
    static int igood[ncomb];
    static int jhel;
    double t;
    // Helicities for the process
    static const int helicities[ncomb][nexternal] = { { -1, -1, -1, -1 },
                                                      { -1, -1, -1, 1 },
                                                      { -1, -1, 1, -1 },
                                                      { -1, -1, 1, 1 },
                                                      { -1, 1, -1, -1 },
                                                      { -1, 1, -1, 1 },
                                                      { -1, 1, 1, -1 },
                                                      { -1, 1, 1, 1 },
                                                      { 1, -1, -1, -1 },
                                                      { 1, -1, -1, 1 },
                                                      { 1, -1, 1, -1 },
                                                      { 1, -1, 1, 1 },
                                                      { 1, 1, -1, -1 },
                                                      { 1, 1, -1, 1 },
                                                      { 1, 1, 1, -1 },
                                                      { 1, 1, 1, 1 } };
    // Denominators: spins, colors and identical particles
    const int denominators = 36;

    ntry = ntry + 1;

    matrix_element = 0.0;

    if (sum_hel == 0 || ntry < 10) {
        // Calculate the matrix element for all helicities
        for (int ihel = 0; ihel < ncomb; ihel++) {
            if (goodhel[ihel] || ntry < 2) {
                calculate_wavefunctions(perm, helicities[ihel], param,
                                        param_aS);
                t = matrix_uux_mupmum();
                double tsum = 0;

                matrix_element += t;
                tsum += t;

                // Store which helicities give non-zero result
                if (tsum != 0. && !goodhel[ihel]) {
                    goodhel[ihel] = true;
                    ngood++;
                    igood[ngood] = ihel;
                }
            }
        }
        jhel = 0;
        sum_hel = min(sum_hel, ngood);
    } else {
        // Only use the "good" helicities
        for (int j = 0; j < sum_hel; j++) {
            jhel++;
            if (jhel >= ngood)
                jhel = 0;
            double hwgt = double(ngood) / double(sum_hel);
            int ihel = igood[jhel];
            calculate_wavefunctions(perm, helicities[ihel], param, param_aS);
            t = matrix_uux_mupmum();

            matrix_element += t * hwgt;
        }
    }

    matrix_element /= denominators;

    return matrix_element;
}

/**
@brief Returns the color correlated born amplitude

@param id1, id2 type of the initial state partons
@param i, j number of the particle which radiates a parton
*/
/*double ME_uux_mupmum::ColorCorrelated(int id1, int id2, int i, int j) {
    double me = MatrixElement(id1, id2);
    if (i == 0 && j == 1) {
        return (4.0 / 3.0) * me;
    }
    if (i == 1 && j == 0) {
        return (4.0 / 3.0) * me;
    }
    return 0.0;
}*/
//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void ME_uux_mupmum::calculate_wavefunctions(const int perm[], const int hel[],
                                            const Parameters_sm &pars,
                                            const Parameters_alphaS &param_aS) {
    // Calculate all wavefunctions
    ixxxxx(momenta[perm[0]], mME[0], hel[0], +1, w[0]);
    oxxxxx(momenta[perm[1]], mME[1], hel[1], -1, w[1]);
    ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]);
    FFV1P0_3(w[0], w[1], pars.GC_2, 0.0, 0.0, w[4]);
    FFV2_5_3(w[0], w[1], pars.GC_51, pars.GC_58, pars.MZ, pars.WZ, w[5]);

    // Calculate all amplitudes
    // Amplitude(s) for diagram number 0
    FFV1_0(w[2], w[3], w[4], pars.GC_3, amp[0]);
    FFV2_4_0(w[2], w[3], w[5], pars.GC_50, pars.GC_59, amp[1]);
}

double ME_uux_mupmum::matrix_uux_mupmum() {
    // Local variables
    std::complex<double> ztemp;
    std::complex<double> jamp;
    // The color matrix;
    static const double cf = 3;

    // Calculate color flows
    jamp = +amp[0] + amp[1];

    // Sum and square the color flows to get the matrix element
    double matrix = 0;

    ztemp = cf * jamp;
    matrix = matrix + real(ztemp * conj(jamp));

    // Store the leading color flows for choice of color
    jamp2 += real(jamp * conj(jamp));

    return matrix;
}

