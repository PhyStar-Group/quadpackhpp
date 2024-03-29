#pragma once
#include"quadpack.hpp"

/* DQAGE - Approximation to definite integral. (From QUADPACK)
 *
 *    Allows user's choice of Gauss-Kronrod integration rule.
 *
 * PARAMETERS:
 *
 *    f() - T_real precision function to be integrated.
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
 *
 *    limit - maximum number of subintervals.
 */

namespace quadpack {
    template<typename T_fun, typename T_param, realtype T_real>
    T_real Quadpack<T_fun, T_param, T_real>::qage(T_fun f, T_param user_data[], T_real a, T_real b, T_real epsabs, T_real epsrel, int irule, T_real* abserr, int* neval, int* ier, int* last)
    {
        T_real area, area1, area2, area12, a1, a2, b1, b2, c, defabs;
        T_real defab1, defab2, errbnd, errmax, error1, error2;
        T_real erro12, errsum, resabs, result;
        T_real alist[LIMIT], blist[LIMIT], rlist[LIMIT], elist[LIMIT];
        int iroff1, iroff2, k, keyf, maxerr, nrmax, iord[LIMIT], limit;

        limit = LIMIT - 1;
        *ier = 0;
        *neval = 0;
        *last = 0;
        result = 0.0;
        *abserr = 0.0;
        alist[0] = a;
        blist[0] = b;
        rlist[0] = 0.0;
        elist[0] = 0.0;
        iord[0] = 0;
        defabs = 0.0;
        resabs = 0.0;
        if ((epsabs < 0.0) && (epsrel < 0.0))
            *ier = 6;
        if (*ier == 6) return result;

        /* First approximation to the integral. */
        keyf = irule;
        if (irule <= 0) keyf = 1;
        if (irule >= 7) keyf = 6;
        c = keyf;
        *neval = 0;
        switch (keyf) {
        case 1:
            result = gk15(f, user_data, a, b, abserr, &defabs, &resabs);
            break;
        case 2:
            result = gk21(f, user_data, a, b, abserr, &defabs, &resabs);
            break;
        case 3:
            result = gk31(f, user_data, a, b, abserr, &defabs, &resabs);
            break;
        case 4:
            result = gk41(f, user_data, a, b, abserr, &defabs, &resabs);
            break;
        case 5:
            result = gk51(f, user_data, a, b, abserr, &defabs, &resabs);
            break;
        case 6:
            result = gk61(f, user_data, a, b, abserr, &defabs, &resabs);
            break;
        }
        *last = 0;
        rlist[0] = result;
        elist[0] = *abserr;
        iord[0] = 0;

        /* Test on accuracy. */
        errbnd = max(epsabs, epsrel * fabs(result));
        if ((*abserr <= 50.0 * epmach * defabs) && (*abserr > errbnd))
            *ier = 2;
        if (limit == 0) *ier = 1;
        if ((*ier != 0) || ((*abserr <= errbnd) && (*abserr != resabs)) ||
            (*abserr == 0.0)) goto _60;

        /* Initialization. */
        errmax = *abserr;
        maxerr = 0;
        area = result;
        errsum = *abserr;
        nrmax = 0;
        iroff1 = 0;
        iroff2 = 0;

        /* Main Loop. */
        for (*last = 1; *last <= limit; (*last)++) {
            /* Bisect the subinterval with the largest error estimate. */
            a1 = alist[maxerr];
            b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
            a2 = b1;
            b2 = blist[maxerr];
            switch (keyf) {
            case 1:
                area1 = gk15(f, user_data, a1, b1, &error1, &resabs, &defab1);
                area2 = gk15(f, user_data, a2, b2, &error2, &resabs, &defab2);
                break;
            case 2:
                area1 = gk21(f, user_data, a1, b1, &error1, &resabs, &defab1);
                area2 = gk21(f, user_data, a2, b2, &error2, &resabs, &defab2);
                break;
            case 3:
                area1 = gk31(f, user_data, a1, b1, &error1, &resabs, &defab1);
                area2 = gk31(f, user_data, a2, b2, &error2, &resabs, &defab2);
                break;
            case 4:
                area1 = gk41(f, user_data, a1, b1, &error1, &resabs, &defab1);
                area2 = gk41(f, user_data, a2, b2, &error2, &resabs, &defab2);
                break;
            case 5:
                area1 = gk51(f, user_data, a1, b1, &error1, &resabs, &defab1);
                area2 = gk51(f, user_data, a2, b2, &error2, &resabs, &defab2);
                break;
            case 6:
                area1 = gk61(f, user_data, a1, b1, &error1, &resabs, &defab1);
                area2 = gk61(f, user_data, a2, b2, &error2, &resabs, &defab2);
                break;
            }

            /* Improve previous approximations to integral and error,
                    and test for accuracy. */
            (*neval) += 1;
            area12 = area1 + area2;
            erro12 = error1 + error2;
            errsum = errsum + erro12 - errmax;
            area = area + area12 - rlist[maxerr];
            if ((defab1 != error1) && (defab2 != error2)) {
                if ((fabs(rlist[maxerr] - area12) <= 1.0e-5 * fabs(area12)) &&
                    (erro12 >= .99 * errmax))
                    iroff1++;
                if ((*last > 9) && (erro12 > errmax))
                    iroff2++;
            }
            rlist[maxerr] = area1;
            rlist[*last] = area2;
            errbnd = max(epsabs, epsrel * fabs(area));
            if (errsum > errbnd) {

                /* Test for roundoff error and eventually set error flag. */
                if ((iroff1 > 6) || (iroff2 > 20))
                    *ier = 2;

                /* Set error flag in the case that the number of subintervals
                    equals the limit. */
                if (*last == limit)
                    *ier = 1;

                /* Set error flag in the case of bad integrand behavior at a
                    point of the integration range. */
                if (max(fabs(a1), fabs(b2)) <= (1.0 + c * 1000.0 * epmach) *
                    (fabs(a2) + 1.0e4 * uflow))
                    *ier = 3;
            }
            /* Append the newly-created intervals to the list. */

            if (error2 <= error1) {
                alist[*last] = a2;
                blist[maxerr] = b1;
                blist[*last] = b2;
                elist[maxerr] = error1;
                elist[*last] = error2;
            }
            else {
                alist[maxerr] = a2;
                alist[*last] = a1;
                blist[*last] = b1;
                rlist[maxerr] = area2;
                rlist[*last] = area1;
                elist[maxerr] = error2;
                elist[*last] = error1;
            }

            /* Call DQSORT to maintain the descending ordering in the list of
                error estimates and select the subinterval with the
                largest error estimate (to be bisected next). */

            qsort(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);
            if ((*ier != 0) || (errsum <= errbnd)) break;
        }

        /* Compute final result. */

        result = 0.0;
        for (k = 0; k <= *last; k++) {
            result += rlist[k];
        }
        *abserr = errsum;
    _60:
        if (keyf != 1)
            *neval = (10 * keyf + 1) * (2 * (*neval) + 1);
        else
            *neval = 30 * (*neval) + 15;

        return result;
    }
}