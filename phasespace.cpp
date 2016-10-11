#include "phasespace.h"

namespace Phasespace {

namespace {

/// gamma returns gamma for a Lorentz boost for parton CMS -> lab frame.
double gamma(double x1, double x2) {
    return 0.5 * (x1 + x2) / sqrt(x1 * x2);
}

/// betagamma returns beta*gamma for a Lorentz boost for parton CMS -> lab
/// frame.
double betagamma(double x1, double x2) {
    return -0.5 * (x1 - x2) / sqrt(x1 * x2);
}

/// boost_z boosts the 4-momentm \p p in z direction.
void boost_z(double gamma, double betagamma, Math::FourMomentum *p) {
    const double E = p->E();
    const double pz = p->PZ();
    p->SetE(E * gamma - betagamma * pz);
    p->SetPZ(gamma * pz - E * betagamma);
}
} // end namespace

// ----------------------------------------------------------------------------
// Implementation of phasespace function.
//
// Note: They have to be implemented in the header because the compiler has to
// now the functions. If they would be in a cpp file, the compiler wouldn't know
// them if the program is compiled to a lib.
// ----------------------------------------------------------------------------

/**
SetToCMSFromLab sets this phasespace to ps_lab in the partonic CMS frame.

@param ps_lab phase space in lab frame
*/
void Phasespace::SetToCMSFromLab(Phasespace const *const ps_lab) {
    const double x1 = ps_lab->X1;
    const double x2 = ps_lab->X2;
    const double g = gamma(x1, x2);
    const double bg = -betagamma(x1, x2);

    N = ps_lab->N;
    X1 = x1;
    X2 = x2;
    S = ps_lab->S;
    Masses = ps_lab->Masses;
    Momenta = ps_lab->Momenta;
    for (int i = 0; i < N+2; i++) {
        boost_z(g, bg, &Momenta[i]);
    }
}

/**
SetToLabFromCMS sets this phasespace to ps_cms in the lab frame.

@param ps_cms phase space in partonic CMS.
*/
void Phasespace::SetToLabFromCMS(Phasespace const *const ps_cms) {
    const double x1 = ps_cms->X1;
    const double x2 = ps_cms->X2;
    const double g = gamma(x1, x2);
    const double bg = betagamma(x1, x2);

    N = ps_cms->N;
    X1 = x1;
    X2 = x2;
    S = ps_cms->S;
    Masses = ps_cms->Masses;
    Momenta = ps_cms->Momenta;
    for (int i = 0; i < N+2; i++) {
        boost_z(g, bg, &Momenta[i]);
    }
}



/**
@brief Generate 2 particle phase space.

A two particle phase space is created from the input parameters
-# x1 - momentum fraction of the 1st inital state parton
-# x2 - momentum fraction of the 2nd inital state parton
-# y  - cos(theta) of the first final state particle
-# phi - azimuthal angle of the 1st final state particle

These parameters have to be in the above order in the array x.
All parameters have to be from 0 to 1, i.e. the value for phi is multiplied
internally with 2pi to obtain the angle. (y is converted internally as well).

@param[out] ps the phase space is written to ps.
@param n number of parameters has to be 4.
@param x n parameters.
*/
void TwoParticleGenerator::operator()(Phasespace *ps, int n,
                                      double const *const x, double S, int m,
                                      double *masses) const {
    assert(n == 4 && "need four parameters");
    const double x1 = x[0];
    const double x2 = x[1];
    assert(x1 > 0.0 && x1 <= 1.0 && "x1 not in range (0,1]");
    assert(x2 > 0.0 && x2 <= 1.0 && "x2 not in range (0,1]");
    const double CosTheta = 2.0 * x[2] - 1.0;
    const double phi = 2.0 * M_PI * x[3];
    assert(fabs(CosTheta) <= 1.0 && "cos theta not in range [-1.0, 1.0]");
    // assert(phi >= 0.0 && phi < 2.0 * M_PI &&
    //        "azimuth angle not in range [0, 2pi)");

    ps->N = 2;
    ps->S = S;
    const double s = x1 * x2 * S;
    const double sqrts = sqrt(s);

    ps->X1 = x1;
    ps->X2 = x2;

    // initial state particles
    const double p = sqrts / 2.0;
    ps->Momenta[0].Set(p, 0.0, 0.0, p);
    ps->Momenta[1].Set(p, 0.0, 0.0, -p);

    ps->Masses[0] = masses[0];
    ps->Masses[1] = masses[1];
    const double m1 = ps->Masses[0];
    const double m2 = ps->Masses[1];
    // final state particles
    const double m1_2 = m1 * m1;
    const double m2_2 = m2 * m2;
    const double m1_4 = m1_2 * m1_2;
    const double m2_4 = m2_2 * m2_2;
    const double p_final = sqrt(m1_4 - 2.0 * m1_2 * m2_2 + m2_4 -
                                2.0 * m1_2 * s - 2.0 * m2_2 * s + s * s) /
                           (2.0 * sqrts);
    ps->Momenta[2].SetFromPMCosThetaPhi(p_final, m1, CosTheta, phi);
    const double px = -ps->Momenta[2].PX();
    const double py = -ps->Momenta[2].PY();
    const double pz = -ps->Momenta[2].PZ();
    ps->Momenta[3].Set(sqrt(p_final * p_final + m2_2), px, py, pz);

    ps->Jacobian = 1.0/(8.0*M_PI);
}

} // namespace Phasespace

void GenBornPhasespace(Phasespace::Phasespace *ps_out, const int ndim,
                       const double x[], const Data *userdata) {
    double SqrtS = userdata->SqrtS;
    double y_b = x[2];
    Phasespace::TwoParticleGenerator gen;
    double xt[] = { x[0], x[1], y_b, x[6] };
    double masses[2] = { 0.0 }; // massive case not tested!
    gen(ps_out, 4, xt, SqrtS * SqrtS, 2, masses);
}