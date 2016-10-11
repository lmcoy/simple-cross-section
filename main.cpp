#include <cmath>
#include <cassert>

#include <mpi.h>

#include "LHAPDF/LHAPDF.h"

#include "phasespace.h"
#include "vegas.h"
#include "rnd.h"
#include "data.h"
#include "fourmomentum.h"
#include "parameters_sm.h"
#include "uux_mupmum.h"

const double CrossSectionToPb = 3.8937944929e8;

bool applycuts(const Phasespace::Phasespace &ps, const Data & data) {
  // cut on inv. mass of final state particles
  const auto & momenta = ps.Momenta;
  Math::FourMomentum mll = momenta[2].Plus(momenta[3]);
  if (mll.Dot(mll) < data.mllmin * data.mllmin) {
        return false;
  }
  // implement cuts
  return true;
}

double ME2(const Phasespace::Phasespace &ps, const Parameters_sm & p, const Parameters_alphaS &p_as) {
  /* 
  // Example for phase space usage
  auto q = ps.Momenta[0].Plus(ps.Momenta[1]);
  double s = q.Dot(q);
  
  auto tmp = ps.Momenta[0].Minus(ps.Momenta[2]);
  double t = tmp.Dot(tmp);
  
  tmp = ps.Momenta[0].Minus(ps.Momenta[3]);
  double u = tmp.Dot(tmp);
  */
  
  // madgraph matrix element for u ux -> mu+ mu-
  ME_uux_mupmum me;
  
  
  int perm[] = {0, 1, 2, 3};
  
  // for ux u -> mu+ mu-
  //int perm[] = {1, 0, 2, 3};
  
  double m = me.Calculate(ps, perm, p, p_as);
  
  return m;
}

// x1, x2, x3 -- for NLO real radiation, can be ignored here
// wgt -- VEGAS weight
int XSecBorn(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *out, Data *params) {
  
    Phasespace::Phasespace ps_lab;
    ps_lab.SetToLabFromCMS(&ps);
    bool cuts = applycuts(ps_lab,*params);
    if (!cuts) {
      *out = 0.0;
      return 0;
    }
    
    double f1 = LHAPDF::xfx(ps.X1, 91.188, 2)/ps.X1; // LHAPDF
    double f2 = LHAPDF::xfx(ps.X2, 91.188, -2)/ps.X2; // LHAPDF
    double lumi = f1 * f2;
    
    Parameters_sm params_sm;
    params_sm.Set(1.16639E-005, 132.507, 91.1876, 2.4414039999999999);
    Parameters_alphaS params_as;
    params_as.Set(LHAPDF::alphasPDF(91.188)); 
    
    double me2 = ME2(ps,params_sm, params_as); // squared matrix element
	
    double s = ps.X1 * ps.X2 * ps.S;
    double dsigma = ps.Jacobian / (2.0 * s) * me2;

    double result_born = dsigma * CrossSectionToPb;
	
    *out = lumi * result_born;
    return 0;
}

typedef int (*function)(const Phasespace::Phasespace &, double, double, double,
                        double, double *, Data *);

double clip(double x) {
    if (x >= 1.0 && x < 1.0 + 1e-9) {
        return 1.0 - 1e-12;
    }
    if (x <= 0.0 && x > -1e-9) {
        return 1e-12;
    }
    return x;
}

// integrand for 2 -> Z -> 2
template <function func>
static int integrand_resonance(int n, const double *x, const double wgt, void *params,
                      int threadid, double *out) {
    constexpr double M = 91.188; // mass of the resonance 
    constexpr double G = 2.45; // width
    Data *userdata = (Data *)params;
    // implement mll cut in Breit-Wigner trafo.
    // can be set to 0.0 for small cut values (mll < ~100).
    double Mll = userdata->mllmin;
    double S = userdata->SqrtS * userdata->SqrtS;
    double Delta =
        atan((S - M * M) / (M * G)) - atan((Mll * Mll - M * M) / (M * G));
    double s = x[0];
    double z = x[1];
    double ts = Delta * s + atan((Mll * Mll - M * M) / (M * G));
    double tau = clip((G * tan(ts) + M) * M / S);
    assert(tau >= 0.0 && tau <= 1.0);
    double xx[n];
    double x1 = (1.0 - tau) * z + tau;
    xx[0] = clip(x1);
    xx[1] = clip(tau / x1);

    for (int i = 2; i < n; i++) {
        xx[i] = clip(x[i]);
    }

    double tmp = tau * S - M * M;
    double inv_breit_wigner = (M * M * G * G + tmp * tmp) / (M * G * S);
    double jac = Delta * (1.0 - tau) / x1 * inv_breit_wigner;

    Phasespace::Phasespace ps;
    GenBornPhasespace(&ps, n, xx, userdata);
    double integrand = 0.0;
    int ret = func(ps, xx[3], xx[4], xx[5], jac * wgt, &integrand, userdata);

    *out = integrand * jac;
    return ret;
}


template <function func>
static int integrand_simple(int n, const double *x, const double wgt, void *params,
                      int threadid, double *out) {
    
    Data *userdata = (Data *)params;
    double xx[n];

    for (int i = 0; i < n; i++) {
        xx[i] = clip(x[i]);
    }

    Phasespace::Phasespace ps;
    GenBornPhasespace(&ps, n, xx, userdata);
    double integrand = 0.0;
    int ret = func(ps, xx[3], xx[4], xx[5], wgt, &integrand, userdata);

    *out = integrand;
    return ret;
}

int main(int argc, char *argv[]) {
  MPI_Init(NULL, NULL);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // init pdfs
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet("cteq6ll", LHAPDF::LHPDF, 0);

  const int NDIM = 4; // 4d integral x_1, x_2, and two angles
  VegasState *state = vegas_new(NDIM, 50);
  vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd, Random::free_rnd);
  
  vegas_integrand integrand_ptr = 0;
  integrand_ptr = integrand_resonance<XSecBorn>;
  //integrand_ptr = integrand_simple<XSecBorn>;

  Data userdata;
  userdata.SqrtS = 14000.0;
  // only important for integrand_resonance.
  // This implents a cut on the inv. mass of the two final state particles.
  userdata.mllmin = 50.0; 

  // -----------------------------------------
  // VEGAS grid setup 
  // -----------------------------------------
  int verbosity = 1;

  int seed1 = 1; // random seed
  vegas_seed_random(state, seed1);

  // integration result
  double integral = 0.0, error = 0.0, prob = 0.0;
  // number of vegas iterations and integrand evaluations per interation
  int N_setup_points = 10000;
  int N_setup_iterations = 5;
  if (rank == 0) {
    printf("*******************************************************************\n");
    printf("*                                                                 *\n");
    printf("* setup VEGAS grid                                                *\n");
    printf("*                                                                 *\n");
    printf("*******************************************************************\n");
  }
  vegas_integrate(state, seed1, integrand_ptr, (void *)&userdata, 1,
                    N_setup_points,N_setup_iterations, verbosity, &integral, &error, &prob);
  // values in integral, error etc. should be discarded.

  // -----------------------------------------
  // VEGAS integration 
  // -----------------------------------------
  vegas_reset_int(state); // only keep grid from previous integration
  vegas_disable_update_grid(state); // don't update grid during integration

  int seed2 = 2; // seed for VEGAS integration
  vegas_seed_random(state, seed2);
  int N_points = 2000000;
  int N_iterations = 500;
  if (rank == 0) {
    printf("*******************************************************************\n");
    printf("*                                                                 *\n");
    printf("* VEGAS integration                                               *\n");
    printf("*                                                                 *\n");
    printf("*******************************************************************\n");
  }
  vegas_integrate(state, seed2, integrand_ptr, (void *)&userdata, 1,
                    N_points, N_iterations,
                    verbosity, &integral, &error, &prob);
  
  if (rank == 0) {
    printf("*******************************************************************\n");
    printf("*                                                                 *\n");
    printf("* result                                                          *\n");
    printf("*                                                                 *\n");
    printf("*******************************************************************\n");
    printf("\n");
    printf("I = %g +- %g (chi^2/ndf = %g)\n", integral, error, prob);
    printf("MadGraph: 376.814 +- 0.086 pb\n");
  }
  
  MPI_Finalize();
  
  return 0;
}
