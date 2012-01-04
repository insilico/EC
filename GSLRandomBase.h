#ifndef GSL_RANDOM_BASE_H
#define GSL_RANDOM_BASE_H

//
// GNU Scientific Library includes
//
#include "gsl/gsl_rng.h"

/**
   A base class for GNU Scientific Library (GSL) random number
   functions.  The setup, initialization and clean-up is the same for
   all GSL random number functions.  This class abstracts away these
   details, placing the stup and initialization in the class
   constructor and the clean-up in the class destructor.  The class
   constructor is passed a seed value for the random number generator.

   A class that provides access to one or more GSL random number
   functions should be derived from this class.  This class must
   provide an implementation for the nextRandVal() pure virtual
   function.  The nextRandVal will call the specific random
   number function (for example gsl_ran_ugaussian() for Gaussian
   distribution or gsl_ran_flat() for a flat random number 
   distribution).

   This class uses the default random number generator.  At least on
   Windows XP using the Visual C++ 6.0 compiler the type definitions
   for the random functions (for example gsl_rng_mt19937 or
   gsl_rng_knuthran) would not link properly.  Perhaps they are not
   properly exported from the pre-built library.

   I decided to use the GSL because is is supported on all major
   platforms (UNIX, Linux and Windows) and provides high quality
   pseudo-random number generation support.  The standard POSIX rand()
   function is notorious for its poor quality.  While the random()
   function on UNIX provides better pseudo-random number quality, but
   is still not as good as functions like MT19937.
 */

class GSLRandomBase {
private:
    GSLRandomBase(const GSLRandomBase &rhs);

protected:

    gsl_rng *state() {
        return rStatePtr_;
    }
    gsl_rng *rStatePtr_;

public:

    GSLRandomBase(int seedVal) {
        const gsl_rng_type *T;

        // The gsl_rng_env_setup() function returns a pointer to the
        // type gsl_rng_mt19937, which is associated with the MT19937
        // random number generator.
        T = gsl_rng_env_setup();

        // Allocate a random number state structure
        rStatePtr_ = gsl_rng_alloc(T);

        // set the seed
        gsl_rng_set(rStatePtr_, seedVal);
    } // GSLRandomBase constructor

    virtual ~GSLRandomBase() {
    } // GSLRandomBase destructor

    virtual double nextRandVal() = 0;

};

#endif
