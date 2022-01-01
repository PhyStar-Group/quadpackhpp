#pragma once
#include"quadpack.hpp"

/* DQAWF - Approximation to Fourier integral. (From QUADPACK)
 *
 *
 * PARAMETERS:
 *
 *    f() - T_real precision function to be integrated.
 *
 *    a - lower limit of integration.
 *
 *    omega - parameter in weight function.
 *
 *    sincos - indicates which weight function to use:
 *        sincos = COSINE (= 1) --- use cos(omega*x)
 *        sincos = SINE   (= 2) --- use sin(omega*x)
 *
 *    epsabs - absolute accuracy requested.
 *
 */
namespace quadpack {
    template<typename T_fun, realtype T_real>
    T_real Quadpack<T_fun, T_real>::qawo(T_fun f, T_real user_data[], T_real a, T_real b, T_real omega, int sincos, T_real epsabs,
        T_real epsrel, T_real* abserr, int* neval, int* ier)
    {
        T_real** chebmo, result;
        int i, momcom;

        if ((chebmo = (T_real**)calloc(MAXP1, sizeof(T_real*))) == NULL) {
            fprintf(stderr, "Out of memory in DQAWO!\n");
            exit(1);
        }
        for (i = 0; i < MAXP1; i++) {
            if ((chebmo[i] = (T_real*)calloc(25, sizeof(T_real))) == NULL) {
                fprintf(stderr, "Out of memory in DQAWO!\n");
                exit(1);
            }
        }

        momcom = 0;
        result = qfour(f, user_data, a, b, omega, sincos, epsabs, epsrel,
            1, MAXP1, abserr, neval, ier, &momcom, chebmo);
        for (i = 0; i < MAXP1; i++)
            free(chebmo[i]);
        free(chebmo);
        return result;
    }
}