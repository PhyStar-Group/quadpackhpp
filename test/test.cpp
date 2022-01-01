// quadpackhpp.cpp: 定义应用程序的入口点。
//
#include<iostream>
#include<math.h>
#include "../quadpack.hpp"

using namespace quadpack;
using namespace std;
using Real = double;
Real f1(Real x, Real user_data[])
{
	return exp(-x);
}
Real f2(Real x, Real user_data[])
{
	return 1.0 / sqrt(fabs(x * x + 2.0 * x - 2.0));
}
Real f3_1(Real x, Real user_data[])
{
	return log(x) / (1.0 + 100.0 * x * x);
}
Real f3_2(Real x, Real user_data[])
{
	return exp(-x * x);
}
Real f4(Real x, Real user_data[])
{
	return pow(x, 3.0) * log(fabs((x * x - 1.0) * (x * x - 2.0)));
}
Real f5(Real x, Real user_data[])
{
	return log(x) / sqrt(x);
}
Real f6_1(Real x, Real user_data[])
{
	return exp(-x * x) * log(1 - x);
}
Real f6_2(Real x, Real user_data[])
{
	return exp(-x);
}
Real f6_3(Real x, Real user_data[])
{
	return cos(100.0 * sin(x));
}
Real f7(Real x, Real user_data[])
{
	Real result;
	result = 1.0 / (5.0 * pow(x, 3.0) + 6.0);
	return result;
}
Real f8_1(Real x, Real user_data[])
{
	return ((x > 0.0) ? (1.0 / sqrt(x)) : 0.0);
}
Real f8_2(Real x, Real user_data[])
{
	return exp(-x);
}
Real f9_1(Real x, Real user_data[])
{
	if (x > 0.0) return log(x);
	return 0.0;
}
Real f9_2(Real x, Real user_data[])
{
	return log(x) * sin(10.0 * x);
}
Real f10(Real x, Real user_data[])
{
	Real result;
	result = 0.0;
	if (x > 0.0) result = (1.0 / pow(1.0 + log(x) * log(x), 2.0));
	return result;
}
Real f11(Real x, Real user_data[])
{
	return exp(-x);
}
Real f12_1(Real x, Real user_data[])
{
	if ((x == 0.0) || (x == 1.0)) return 0.0;
	return sqrt(x / (1.0 - x)) * log(x);
}
Real f12_2(Real x, Real user_data[])
{
	return sqrt(x) * log(x);
}
//using pfun = Real(*)(Real, Real*);
using qd = Quadpack_<Real>;
int main()
{
	{
		Real a, b, omega, result, abserr, epsabs, epsrel;
		Real** chebmo;
		int i, ier, icall, n, neval, momcom;
		n = 21;
		chebmo = (Real**)calloc(n, sizeof(Real*));
		for (i = 0; i < n; i++)
			chebmo[i] = (Real*)calloc(25, sizeof(Real));

		a = 0.0;
		b = 11.0;
		omega = 20.0;
		icall = 1;
		momcom = 0;
		epsabs = 1e-8;
		epsrel = 1e-12;
		result = qd::qfour(f1,nullptr, a, b, omega, COSINE, epsabs, epsrel,
			icall, MAXP1, &abserr, &neval, &ier, &momcom, chebmo);
		printf("\nresult = %.17lg\n", result);
		printf("abserr = %.17lg\n", abserr);
		printf("neval = %d\n", neval);
		printf("momcom = %d\n", momcom);
		for (i = 0; i < n; i++)
			free(chebmo[i]);
		free(chebmo);
	
	}
	{
		Real a, b, epsabs, epsrel, abserr;
		Real y;
		int irule, neval, ier, last;

		a = 0.0;
		b = 1.0;
		epsabs = 0.0;
		epsrel = 1e-8;

		for (irule = 1; irule <= 6; ++irule) {
			y = qd::qage(f2,nullptr, a, b, epsabs, epsrel, irule, &abserr, &neval, &ier, &last);

			printf("G/K rule = %i\n", irule);
			printf("dqage integral = %.17lg\n", y);
			printf("abserr = %.17lg, neval = %d, ier = %d\n",
				abserr, neval, ier);
			printf("\n");
		}
	}
	{
		Real bound, epsabs, epsrel, abserr;
		Real result;
		int inf, neval, ier;

		bound = 0.0;
		inf = 1;
		epsabs = 0.0;
		epsrel = 1e-8;

		result = qd::qagi(f3_1,nullptr, bound, inf, epsabs, epsrel, &abserr, &neval, &ier);

		printf("dqagi integral approximation = %.17lg\n", result);
		printf("estimate of absolute error = %.17lg\n", abserr);
		printf("number of function evaluations = %d\n", neval);
		printf("error code = %d\n", ier);
	}

	{
	Real a, b, epsabs, epsrel, abserr, points[4];
	Real  result;

	int neval, npts2, ier;

	a = 0.0;
	b = 3.0;
	npts2 = 4;
	points[0] = 1.0;	/* location of singularity #1 */
	points[1] = sqrt(2.0);	/* location of singularity #2 */
	epsabs = 0.0;
	epsrel = 1e-6;

	result = qd::qagp(f4,nullptr, a, b, npts2, points, epsabs, epsrel, &abserr, &neval, &ier);

	cout<<"qagp integral = "<<result<<endl;
	cout << "abserr = " << abserr <<", neval ="<< neval <<",ier="<<ier<<endl;
	}
	{
		Real a, b, epsabs, epsrel, abserr;
		Real y;
		int neval, ier;

		a = 0.0;
		b = 1.0;
		epsabs = 0.0;
		epsrel = 1e-6;
		y = qd::qags(f5,nullptr, a, b, epsabs, epsrel, &abserr, &neval, &ier);
		cout<<"qags integral = "<<y<<endl;
		cout << "abserr = " << abserr <<", neval ="<< neval <<",ier="<<ier<<endl;
	}
	{
		Real a, b, epsabs, epsrel, abserr;
		Real y;
		int irule, neval, ier;

		//    a = -1.0;
		//    b = 1.0;
		//    a = 0;
		//    b = 1;
		a = 0.0;
		b = 3.14159265358979;
		epsabs = 0.0;
		epsrel = 1e-3;

		for (irule = 1; irule <= 6; ++irule) {
			y = qd::qag(f6_3, nullptr, a, b, epsabs, epsrel, irule, &abserr, &neval, &ier);

			printf("G/K rule = %i\n", irule);
			printf("dqag integral = %.17lg\n", y);
			printf("abserr = %.17lg, neval = %d, ier = %d\n",
				abserr, neval, ier);
			printf("\n");
		}
	}
	{
		Real a, b, c, epsabs, epsrel, abserr, result;
		int neval, ier;

		/*  a and b are the integration limits */
		a = -1.0;
		b = 5.0;

		/*  c is the parameter of the weight function */
		c = 0.0;

		/*  epsabs and epsrel determine the accuracy requirement */
		epsabs = 0.0;
		epsrel = 1.0e-3;

		result = qd::qawc(f7, nullptr, a, b, c, epsabs, epsrel, &abserr, &neval, &ier);

		printf("Integral approximation = %.12lf\n", result);
		printf("Estimate of absolute error = %.12lf\n", abserr);
		printf("Number of function evaluations = %d\n", neval);
		printf("Error code = %d\n", ier);
	}
	{	
		Real a, omega, result, abserr, epsabs;
		int ier, neval;

		a = 0.0;
		omega = 0.5 * Pi;
		epsabs = 1.0e-8;
		result = qd::qawf(f8_2, nullptr, a, omega, COSINE, epsabs, &abserr, &neval, &ier);
		printf("\nresult = %.17lg\n", result);
		printf("abserr = %.17lg\n", abserr);
		printf("neval = %d\n", neval);
	}
	{
		Real a, b, omega, result, abserr, epsabs, epsrel;
		int ier, neval;

		a = 0.0;
		b = 1.0;
		omega = 10.0 * Pi;
		epsabs = 0.0;
		epsrel = 1e-6;
		result = qd::qawo(f9_1,nullptr, a, b, omega, SINE, epsabs, epsrel, &abserr, &neval, &ier);
		printf("\nresult = %.18lg\n", result);
		printf("abserr 0= %.18lg\n", abserr);
		printf("neval = %d\n", neval);
		printf("ier = %d\n", ier);
	}
	{
		Real a, b, alfa, beta, epsabs, epsrel, abserr, result;
		int wgtfunc, neval, ier;

		/*  a and b are the integration limits */
		a = 0.0;
		b = 1.0;

		/*  alfa, beta and wgtfunc determine the weight function */
		alfa = 0.0;
		beta = 0.0;
		wgtfunc = 2;

		/*  epsabs and epsrel determine the accuracy requirement */
		epsabs = 0.0;
		epsrel = 1.0e-5;

		result = qd::qaws(f10,nullptr, a, b, alfa, beta, wgtfunc, epsabs, epsrel, &abserr,
			&neval, &ier);

		printf("Integral approximation = %.12lf\n", result);
		printf("Estimate of absolute error = %.12lf\n", abserr);
		printf("Number of function evaluations = %d\n", neval);
		printf("Error code = %d\n", ier);
	}
	{
		Real a, b, omega, result, abserr, resabs, resasc;
		Real** chebmo;
		int i, nrmom, ksave, n, neval, momcom;

		n = 21;
		chebmo = (Real**)calloc(n, sizeof(Real*));
		for (i = 0; i < n; i++)
			chebmo[i] = (Real*)calloc(25, sizeof(Real));

		a = 0.0;
		b = 11.0;
		omega = 20.0;
		nrmom = 0;
		ksave = 0;
		momcom = 0;
		result = qd::qc25o(f11,nullptr, a, b, omega, COSINE, nrmom, MAXP1,
			ksave, &abserr, &neval, &resabs, &resasc,
			&momcom, chebmo);
		printf("\nresult = %.14lg\n", result);
		printf("abserr = %lg\n", abserr);
		printf("neval = %d\n", neval);
		printf("momcom = %d\n", momcom);
		for (i = 0; i < n; i++)
			free(chebmo[i]);
		free(chebmo);
	}
	{	
		Real a, b, epsabs, epsrel, abserr;
		Real y;
		int neval, ier;

		a = 0.0;
		b = 1.0;
		epsabs = 0.0e-8;
		epsrel = 1e-12;

		y = qd::qng(f12_2,nullptr, a, b, epsabs, epsrel, &abserr, &neval, &ier);

		printf("dqng integral = %.17lg\n", y);
		printf("abserr = %.17lg, neval = %d, ier = %d\n",
			abserr, neval, ier); 
	}
	cout << "Hello CMake." << endl;
	return 0;
}
