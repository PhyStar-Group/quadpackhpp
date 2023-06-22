#pragma once
#include"quadpack.hpp"

/*  DQAWC - Computation of Cauchy principal value
 *
 *  PARAMETERS:
 *
 *      f() -   double precision function defining the integrand.
 *
 *      a   -   lower limit of integration
 *
 *      b   -   upper limit of integration
 *
 *      c   -   parameter in the weight function
 *
 *      epsabs  -   absolute accuracy requested
 *
 *      epsrel  -   relative accuracy requested
 *
 *      abserr  -   estimate of the modulus of the absolute error
 *
 *      neval   -   number of function evaluations
 *
 *      ier     -   error code
 */
namespace quadpack {
    template<typename T_fun, typename T_param, realtype T_real>
    T_real Quadpack<T_fun, T_param, T_real>::qawc(T_fun f, T_param user_data[], T_real a, T_real b, T_real c, T_real epsabs,
        T_real epsrel, T_real* abserr, int* neval, int* ier)
    {
        T_real result;

        result = qawce(f, user_data, a, b, c, epsabs, epsrel, abserr, neval, ier);
        return result;
    }
}