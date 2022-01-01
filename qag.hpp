#pragma once
#include"quadpack.hpp"
/* DQAG - Approximation to definite integral. (From QUADPACK)
 *
 *  Calls DQAGE with appropriate parameters assigned.
 *
 * PARAMETERS:
 *
 *    f() - double precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    b - upper limit of integration.
 *
 *    epsabs - absolute accuracy requested.
 *
 *    epsrel - relative accuracy requested.
 *
 *    irule - integration rule to be used as follows:
 *        irule = 1 -- G_K 7-15
 *        irule = 2 -- G_K 10-21
 *        irule = 3 -- G_K 15-31
 *        irule = 4 -- G_K 20-41
 *        irule = 5 -- G_K 25-51
 *        irule = 6 -- G_K 30-61
 */

namespace quadpack {
    template<typename T_fun, realtype T_real>
    T_real Quadpack<T_fun, T_real>::qag(T_fun f, T_real user_data[], T_real a, T_real b, T_real epsabs, T_real epsrel, int irule, T_real* abserr, int* neval, int* ier)
    {
        T_real result;
        int last;

        result = qage(f, user_data, a, b, epsabs, epsrel, irule, abserr, neval, ier, &last);

        return result;
    }
}