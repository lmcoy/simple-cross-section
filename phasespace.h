#ifndef PHASESPACE_H_
#define PHASESPACE_H_

#include "fourmomentum.h"
#include "data.h"

namespace Phasespace {

const int MAXMOM = 6;

class Phasespace {
  public:
    void SetToCMSFromLab(Phasespace const *const ps_lab);
    void SetToLabFromCMS(Phasespace const *const ps_cms);

    double X1; ///< momentum fraction of the 1st initial state particles
    double X2; ///< momentum fraction of the 2nd initial state particles
    double S;  ///< CMS energy of the protons
    std::array<double, MAXMOM-2> Masses; ///< masses of the final state particles
    std::array<Math::FourMomentum, MAXMOM> Momenta;
    double Jacobian; ///< jacobian determinante
    int N;
};

/**
Generator for a two particle phase space.
*/
class TwoParticleGenerator {
  public:
    void operator()(Phasespace *ps, int n, double const *const x, double S,
                    int m, double *masses) const;
};

} // namespace Phasespace

// Generate two particle phase space
//
// ps_out  --  phase space will be written to ps_out
// dim -- number of dimensions (4 for two particle ps)
// x -- random numbers in (0,1)
// userdata -- addition variables
void GenBornPhasespace(Phasespace::Phasespace *ps_out, const int ndim,
                       const double x[], const Data *userdata);

#endif