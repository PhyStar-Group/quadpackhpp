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
    T_real Quadpack<T_fun, T_real>::qawf(T_fun f, T_real user_data[], T_real a, T_real omega, int sincos, T_real epsabs,
        T_real* abserr, int* neval, int* ier)
    {
        T_real** chebmo, erlst[50];
        T_real result, rslst[50];

        int ierlst[50], i;
        int limlst;

        if ((chebmo = (T_real**)calloc(MAXP1, sizeof(T_real*))) == NULL) {
            fprintf(stderr, "Out of memory in qawf!\n");
            exit(1);
        }
        for (i = 0; i < MAXP1; i++) {
            if ((chebmo[i] = (T_real*)calloc(25, sizeof(T_real))) == NULL) {
                fprintf(stderr, "Out of memory in qawf!\n");
                exit(1);
            }
        }
        *ier = 6;
        *neval = 0;
        result = 0.0;
        *abserr = 0.0;

        /* Dimensioning parameters.
         *    limlst - upper bound on number of cycles,
         *    MAXP1 - upper bound on the number of Chebyshev moments.
         */
        limlst = 50;

        /* Check validity of limlst and MAXP1. */
        if ((limlst < 3) || (MAXP1 < 1))
            goto _10;

        /* Prepare call for dqawfe. */
        result = qawfe(f, user_data, a, omega, sincos, epsabs, limlst, MAXP1,
            abserr, neval, ier, rslst, erlst, ierlst, chebmo);
    _10:
        for (i = 0; i < MAXP1; i++)
            free(chebmo[i]);
        free(chebmo);
        return result;
    }
}