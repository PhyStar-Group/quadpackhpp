// quadpackhpp.h: 标准系统包含文件的包含文件
// 或项目特定的包含文件。

#pragma once
#include <iostream>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <concepts>

namespace quadpack {
    constexpr auto uflow = DBL_MIN;
    constexpr auto oflow = DBL_MAX;
    constexpr auto epmach = DBL_EPSILON;
    constexpr auto LIMIT = 1000;
    constexpr auto MAXP1 = 21;
    constexpr auto Pi = 3.14159265358979323846;
    constexpr auto COSINE = 1;
    constexpr auto SINE = 2;
    constexpr auto FALSE = 0;
    constexpr auto TRUE = 1;

#ifndef min
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#endif
#ifndef max
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#endif


    template<typename T>
    concept realtype = std::is_floating_point<T>::value;//only allow float number

    //    typedef T_real(*dq_function_type)(T_real, void*);
    template<realtype T_real>
    using qwgt = T_real(*)(T_real, T_real, T_real, T_real, T_real, int);
    template<realtype T_real>
    static T_real qwgtc(T_real x, T_real c, T_real p2, T_real p3, T_real p4, int kp);
    template<realtype T_real>
    static T_real qwgto(T_real x, T_real omega, T_real p2, T_real p3, T_real p4, int integr);
    template<realtype T_real>
    static T_real qwgts(T_real x, T_real a, T_real b, T_real alpha, T_real beta, int integr);

    /* Integration routines */
    template<typename T_fun, typename T_param, realtype T_real = double>
    struct Quadpack {
        /* Gauss-Kronrod for integration over finite range. */
        static T_real gk15(T_fun f, T_param user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);
        static T_real gk21(T_fun f, T_param user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);
        static T_real gk31(T_fun f, T_param user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);
        static T_real gk41(T_fun f, T_param user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);
        static T_real gk51(T_fun f, T_param user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);
        static T_real gk61(T_fun f, T_param user_data[], T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);

        /* Gauss-Kronrod for integration over infinite range. */
        static T_real gk15i(T_fun f, T_param user_data[], T_real boun, int inf, T_real a, T_real b, T_real* abserr, T_real* resabs, T_real* resasc);

        /* Gauss-Kronrod for integration of weighted function. */
        static T_real gk15w(T_fun f, T_param user_data[], qwgt<T_real> w, T_real p1, T_real p2, T_real p3, T_real p4, int kp, T_real a, T_real b, T_real* abserr,
            T_real* resabs, T_real* resasc);
        static T_real qext(int* n, T_real epstab[], T_real* abserr, T_real res3la[], int* nres);
        static void qsort(int limit, int last, int* maxerr, T_real* ermax, T_real elist[], int iord[], int* nrmax);
        static T_real qagi(T_fun f, T_param user_data[], T_real bound, int inf, T_real epsabs, T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qags(T_fun f, T_param user_data[], T_real a, T_real b, T_real epsabs, T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qagp(T_fun f, T_param user_data[], T_real a, T_real b, int npts2, T_real* points, T_real epsabs, T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qng(T_fun f, T_param user_data[], T_real a, T_real b, T_real epsabs, T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qag(T_fun f, T_param user_data[], T_real a, T_real b, T_real epsabs, T_real epsrel, int irule, T_real* abserr, int* neval, int* ier);
        static T_real qage(T_fun f, T_param user_data[], T_real a, T_real b, T_real epsabs, T_real epsrel, int irule, T_real* abserr, int* neval, int* ier, int* last);
        static void qcheb(T_real* x, T_real* fval, T_real* cheb12, T_real* cheb24);
        static  T_real qc25o(T_fun f, T_param user_data[], T_real a, T_real b, T_real omega, int integr, int nrmom, int maxp1, int ksave, T_real* abserr, int* neval,
            T_real* resabs, T_real* resasc, int* momcom, T_real** chebmo);
        static T_real qfour(T_fun f, T_param user_data[], T_real a, T_real b, T_real omega, int integr, T_real epsabs, T_real epsrel, int icall, int maxp1,
            T_real* abserr, int* neval, int* ier, int* momcom, T_real** chebmo);
        static T_real qawfe(T_fun f, T_param user_data[], T_real a, T_real omega, int integr, T_real epsabs, int limlst, int maxp1, T_real* abserr, int* neval, int* ier,
            T_real* rslst, T_real* erlist, int* ierlst, T_real** chebmo);
        static T_real qawf(T_fun f, T_param user_data[], T_real a, T_real omega, int integr, T_real epsabs, T_real* abserr, int* neval, int* ier);
        static T_real qawo(T_fun f, T_param user_data[], T_real a, T_real b, T_real omega, int integr, T_real epsabs,
            T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qaws(T_fun f, T_param user_data[], T_real a, T_real b, T_real alfa, T_real beta, int wgtfunc,
            T_real epsabs, T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qawse(T_fun f, T_param user_data[], T_real a, T_real b, T_real alfa, T_real beta,
            int wgtfunc, T_real epsabs, T_real epsrel, T_real* abserr,
            int* neval, int* ier);
        static void qmomo(T_real alfa, T_real beta, T_real ri[], T_real rj[], T_real rg[], T_real rh[], int wgtfunc);
        static T_real qc25s(T_fun f, T_param user_data[], T_real a, T_real b, T_real bl, T_real br, T_real alfa,
            T_real beta, T_real ri[], T_real rj[], T_real rg[], T_real rh[],
            T_real* abserr, T_real* resasc, int wgtfunc, int* nev);
        static T_real qc25c(T_fun f, T_param user_data[], T_real a, T_real b, T_real c, T_real* abserr,
            int* krul, int* neval);
        static T_real qawc(T_fun f, T_param user_data[], T_real a, T_real b, T_real c, T_real epsabs,
            T_real epsrel, T_real* abserr, int* neval, int* ier);
        static T_real qawce(T_fun f, T_param user_data[], T_real a, T_real b, T_real c, T_real epsabs,
            T_real epsrel, T_real* abserr, int* neval, int* ier);
    };
}
#include"gk15.hpp"
#include"gk21.hpp"
#include"gk31.hpp"
#include"gk41.hpp"
#include"gk51.hpp"
#include"gk61.hpp"
#include"gk15i.hpp"
#include"gk15w.hpp"
#include"qext.hpp"
#include"qsort.hpp"
#include"qagi.hpp"
#include"qags.hpp"
#include"qagp.hpp"
#include"qng.hpp"
#include"qag.hpp"
#include"qage.hpp"
#include"qwgt.hpp"
#include"qcheb.hpp"
#include"qc25o.hpp"
#include"qfour.hpp"
#include"qawfe.hpp"
#include"qawf.hpp"
#include"qawo.hpp"
#include"qaws.hpp"
#include"qawse.hpp"
#include"qmomo.hpp"
#include"qc25s.hpp"
#include"qc25c.hpp"
#include"qawc.hpp"
#include"qawce.hpp"
namespace quadpack {
    template <realtype T_real, typename T_param>
    using pfun = T_real(*)(T_real, T_param*);
    template <typename T_param>
    using Quadpack_d = Quadpack<pfun<double, T_param>,T_param, double>;
    template <typename T_param>
    using Quadpack_f = Quadpack<pfun<float, T_param>,T_param, float>;
    template <typename T_param>
    using Quadpack_ld = Quadpack<pfun<long double, T_param>,T_param, long double>;
    template <realtype T_real, typename T_param>
    using Quadpack_ = Quadpack<pfun<T_real, T_param>,T_param, T_real>;
}

