#pragma once
#include"quadpack.hpp"
namespace quadpack {
    template<realtype T_real>
    T_real qwgtc(T_real x, T_real c, T_real p2, T_real p3, T_real p4, int kp)
    {
        return 1.0 / (x - c);
    }
    template<realtype T_real>
    T_real qwgto(T_real x, T_real omega, T_real p2, T_real p3, T_real p4, int wgtfunc)
    {
        T_real omx;

        omx = omega * x;
        if (wgtfunc == 1)
            return cos(omx);
        else
            return sin(omx);
    }
    template<realtype T_real>
    T_real qwgts(T_real x, T_real a, T_real b, T_real alpha, T_real beta, int wgtfunc)
    {
        T_real bmx, xma, result;

        xma = x - a;
        bmx = b - x;
        result = pow(xma, alpha) * pow(bmx, beta);
        switch (wgtfunc) {
        case 1:
            return result;
        case 2:
            return result * log(xma);
        case 3:
            return result * log(bmx);
        case 4:
            return result * log(xma) * log(bmx);
        default:
            return result;
        }
    }
}