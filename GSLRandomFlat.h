#ifndef GSL_RANDOM_FLAT_H
#define GSL_RANDOM_FLAT_H

#include "GSLRandomBase.h"

#include "gsl/gsl_randist.h"

/**
  Random numbers in a flat, or uniform distribution

  The class constructor is given a seed and a lower and upper
  bound value for the uniform distribution.  The random numbers
  that result will be a uniform distribution in the range

<pre>
    lower <= randVal < upper
</pre>

 */

class GSLRandomFlat : public GSLRandomBase {
private:
    double lower_, upper_;

public:

    GSLRandomFlat(int seedVal,
            double lower,
            double upper) :
    GSLRandomBase(seedVal),
    lower_(lower),
    upper_(upper) {
        ;
    }

    ~GSLRandomFlat() {
        gsl_rng_free(rStatePtr_);
    }

    double nextRandVal() {
        return gsl_ran_flat(state(), lower_, upper_);
    }
};

#endif
