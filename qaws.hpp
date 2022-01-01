#pragma once
#include"quadpack.hpp"

/*  DQAWS - Approximation to integral with algebraic and/or logarithmic
 *          singularities.
 *
 *  PARAMETERS:
 *
 *      f() - double precision function to be integrated.
 *
 *      a   - double lower limit of integration.
 *
 *      b   - upper limit of integration.
 *
 *      alfa - parameter in the weight function.
 *
 *      beta - parameter in the weight function.
 *
 *      wgtfunc - indicates which weight function is to be used.
 *                  = 1:    (x-a)^alfa * (b-x)^beta
 *                  = 2:    (x-a)^alfa * (b-x)^beta * log(x-a)
 *                  = 3:    (x-a)^alfa * (b-x)^beta * log(b-x)
 *                  = 4:    (x-a)^alfa * (b-x)^beta * log(x-a) * log(b-x)
 *
 *      epsabs  - absolute accuracy requested.
 *
 *      epsrel  - relative accuracy requested.
 *
 */
namespace quadpack {
    template<typename T_fun, realtype T_real>
    T_real Quadpack<T_fun, T_real>::qaws(T_fun f, T_real user_data[], T_real a, T_real b, T_real alfa, T_real beta, int wgtfunc,
        T_real epsabs, T_real epsrel, T_real* abserr, int* neval, int* ier)
    {
        T_real result;

        result = qawse(f, user_data, a, b, alfa, beta, wgtfunc, epsabs, epsrel, abserr,
            neval, ier);
        return result;
    }
}