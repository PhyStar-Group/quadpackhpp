#pragma once
#include"quadpack.hpp"
namespace quadpack {
    template<typename T_fun, realtype T_real>
    T_real Quadpack<T_fun,T_real>::gk15(T_fun f, T_real user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc)
    {
        static T_real XGK15[8] = {
            0.99145537112081263921,
            0.94910791234275852453,
            0.86486442335976907279,
            0.74153118559939443986,
            0.58608723546769113029,
            0.40584515137739716691,
            0.20778495500789846760,
            0.00000000000000000000 };
        static T_real WGK15[8] = {
            0.02293532201052922496,
            0.06309209262997855329,
            0.10479001032225018384,
            0.14065325971552591875,
            0.16900472663926790283,
            0.19035057806478540991,
            0.20443294007529889241,
            0.20948214108472782801 };
        static T_real WG7[4] = {
            0.12948496616886969327,
            0.27970539148927666790,
            0.38183005050511894495,
            0.41795918367346938776 };
        T_real fv1[7], fv2[7];
        T_real absc, centr, dhlgth;
        T_real fc, fsum, fval1, fval2, hlgth;
        T_real resg, resk, reskh, result;
        int j, jtw, jtwm1;

        centr = 0.5 * (a + b);
        hlgth = 0.5 * (b - a);
        dhlgth = fabs(hlgth);

        fc = (*f)(centr, user_data);
        resg = fc * WG7[3];
        resk = fc * WGK15[7];
        *resabs = fabs(resk);
        for (j = 0; j < 3; j++) {
            jtw = 2 * j + 1;
            absc = hlgth * XGK15[jtw];
            fval1 = (*f)(centr - absc, user_data);
            fval2 = (*f)(centr + absc, user_data);
            fv1[jtw] = fval1;
            fv2[jtw] = fval2;
            fsum = fval1 + fval2;
            resg += WG7[j] * fsum;
            resk += WGK15[jtw] * fsum;
            *resabs = *resabs + WGK15[jtw] * (fabs(fval1) + fabs(fval2));
        }
        for (j = 0; j < 4; j++) {
            jtwm1 = j * 2;
            absc = hlgth * XGK15[jtwm1];
            fval1 = (*f)(centr - absc, user_data);
            fval2 = (*f)(centr + absc, user_data);
            fv1[jtwm1] = fval1;
            fv2[jtwm1] = fval2;
            fsum = fval1 + fval2;
            resk = resk + WGK15[jtwm1] * fsum;
            *resabs = (*resabs) + WGK15[jtwm1] * (fabs(fval1) + fabs(fval2));
        }
        reskh = resk * 0.5;
        *resasc = WGK15[7] * fabs(fc - reskh);
        for (j = 0; j < 7; j++)
            *resasc = (*resasc) + WGK15[j] * (fabs(fv1[j] - reskh) +
                fabs(fv2[j] - reskh));
        result = resk * hlgth;
        *resabs = (*resabs) * dhlgth;
        *resasc = (*resasc) * dhlgth;
        *abserr = fabs((resk - resg) * hlgth);
        if ((*resasc != 0.0) && (*abserr != 0.0))
            *abserr = (*resasc) * min(1.0, pow((200.0 * (*abserr) / (*resasc)), 1.5));
        if (*resabs > uflow / (50.0 * epmach))
            *abserr = max(epmach * 50.0 * (*resabs), (*abserr));
        return result;
    }
}