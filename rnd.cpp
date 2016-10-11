#include "rnd.h"

namespace Random {

RNG::RNG(int seed) {
        r = new std::mt19937(seed);
        dis = new std::uniform_real_distribution<>(0.0, 1.0);
    }
    RNG::~RNG() {
        delete r;
        delete dis;
    }

    void RNG::Seed(int seed) {
        r->seed(seed);
    }

    double RNG::Random() {
        return (*dis)(*r);
    }



void *init_rnd(void *args, int seed) { return new RNG(seed); }

void free_rnd(void *r) {
    RNG *rng = (RNG *)r;
    delete rng;
}

double get_rnd(void *r) {
    RNG *rng = (RNG *)r;
    return (*rng->dis)(*rng->r);
}

} // end namespace random
