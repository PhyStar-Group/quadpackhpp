#pragma once
#include"quadpack.hpp"


/* DQAWFE - Approximation to Fourier integral. (From QUADPACK)
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
 *    limlst - upper bound on the number of cycles.
 *
 *  limit - upper bound on the number of subintervals (XXX deleted XXX).
 *
 *    maxp1 - upper bound on the number of Chebyshev moments
 *        which can be stored.
 */
namespace quadpack {
    static constexpr auto P = 0.9;
    template<typename T_fun, typename T_param, realtype T_real>
        T_real Quadpack<T_fun, T_param, T_real>::qawfe(T_fun f, T_param user_data[], T_real a, T_real omega, int sincos, T_real epsabs, int limlst, int maxp1, T_real* abserr, int* neval, int* ier,
            T_real* rslst, T_real* erlst, int* ierlst, T_real** chebmo)
    {
        T_real abseps, correc, cycle, c1, c2, dl, drl;
        T_real ep, eps, epsa, errsum, fact, p1, psum[52], reseps;
        T_real res3la[3], result;

        int ktmin, l, ll, momcom, nev, nres, numrl2, lst;

        /* Test on validity of parameters. */
        result = 0.0;
        *abserr = 0.0;
        *neval = 0;
        *ier = 0;
        ll = 0;
        if (((sincos != COSINE) && (sincos != SINE)) || (epsabs <= 0.0) ||
            (limlst < 3)) *ier = 6;
        if (*ier == 6) goto _999;
        if (omega != 0.0) goto _10;

        /* Integration by DQAGI if omega is zero. */
        if (sincos == COSINE)
            result = qagi(f, user_data, 0.0, 1, epsabs, 0.0, abserr, neval, ier);
        rslst[0] = result;
        erlst[0] = *abserr;
        ierlst[0] = *ier;
        goto _999;

        /* Initialization. */
    _10:
        res3la[0] = 0.0;    /* res3la must be initialized to 0.0 */
        res3la[1] = 0.0;
        res3la[2] = 0.0;
        l = int(fabs(omega));
        dl = T_real(2 * l + 10.);
        cycle = dl * Pi / fabs(omega);
        *ier = 0;
        ktmin = 0;
        *neval = 0;
        numrl2 = -1;    /* used as array index. first use is after increment. */
        nres = 0;
        c1 = a;
        c2 = cycle + a;
        p1 = 1.0 - P;
        eps = epsabs;
        if (epsabs > (uflow / p1))
            eps = epsabs * p1;
        ep = eps;
        fact = 1.0;
        correc = 0.0;
        *abserr = 0.0;
        errsum = 0.0;

        /* Main Loop */
        for (lst = 0; lst < limlst; lst++) {

            /* Integrate over current subinterval. */
            /*    dla = lst;  This line is in the original code, but dla is unused. */
            epsa = eps * fact;
            rslst[lst] = qfour(f, user_data, c1, c2, omega, sincos, epsa, 0.0, lst + 1, maxp1,  // lst+1
                &erlst[lst], &nev, &ierlst[lst], &momcom, chebmo);
            *neval += nev;
            fact *= P;
            errsum += erlst[lst];
            drl = 50.0 * fabs(rslst[lst]);

            /* Test on accuracy with partial sum. */
            if (((errsum + drl) <= epsabs) && (lst >= 5))
                goto _80;
            correc = max(correc, erlst[lst]);
            if (ierlst[lst] != 0)
                eps = max(ep, correc * p1);
            if (ierlst[lst] != 0)
                *ier = 7;
            if ((*ier == 7) && ((errsum + drl) <= (correc * 10.0))
                && (lst > 4)) goto _80;
            numrl2++;
            if (lst > 0)
                goto _20;
            psum[0] = rslst[0];
            goto _40;
        _20:
            psum[numrl2] = psum[ll] + rslst[lst];

            if (lst == 1)
                goto _40;

            /* Test on maximum number of subintervals. */
            if (lst == limlst - 1)
                *ier = 8;

            /* Perform new extrapolation. */
            reseps = qext(&numrl2, psum, &abseps, res3la, &nres);

            /* Test whether extrapolated result is influenced by roundoff. */
            ktmin++;
            if ((ktmin >= 15) && (*abserr <= 0.001 * (errsum + drl)))
                *ier = 9;
            if ((abseps > *abserr) && (lst != 2))
                goto _30;
            *abserr = abseps;
            result = reseps;
            ktmin = 0;

            /* If ier is not 0, check whether direct result (partial sum) or
             * extrapolated result yields the best integral approximation.
             */
            if (((*abserr + 10.0 * correc) <= epsabs) || (*abserr <= epsabs) &&
                (10.0 * correc >= epsabs)) goto _60;
        _30:
            if ((*ier != 0) && (*ier != 7))
                goto _60;
        _40:
            ll = numrl2;
            c1 = c2;
            c2 += cycle;
            //        _50:
            ;
        }

        /* Set final result and error estimate. */
    _60:
        (*abserr) += (10.0 * correc);
        if (*ier == 0)
            goto _999;
        if ((result != 0.0) && (psum[numrl2] != 0.0))
            goto _70;
        if (*abserr > errsum)
            goto _80;
        if (psum[numrl2] == 0.0)
            goto _999;
    _70:
        if ((*abserr / fabs(result) > (errsum + drl) / fabs(psum[numrl2])))
            goto _80;
        if ((*ier >= 1) && (*ier != 7))
            (*abserr) += drl;
        goto _999;
    _80:
        result = psum[numrl2];
        *abserr = errsum + drl;
    _999:
        return result;
    }
}