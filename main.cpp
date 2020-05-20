#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "winbgi2.h"
#include "rk4.h"

double y_0 = -5;
double t0 = 0;
double lambda = -2;

double fun(double t, double y);
double anl(double t);
double euler(double t0, double y_0, double t_k, double h, double(*fun)(double, double));
double rk4r(double t0, double y_0, double t_k, double h, double(*fun)(double, double));

void main()
{
	double h, ye,t_k,eps_e,eps_rk,ya,yrk,t,n;

	t_k = 5;
	h = 0.1;
	t = t0;

	printf("t\t\t ya\t\t ye\t\t eps_e\t\t yrk\t\t eps_rk\n");
	for (int i = 0; i < t_k / h; i++)
	{

		t += h;

		ya = anl(t);
		printf("%lf\t%lf\t", t, ya);

		ye = euler(t0, y_0, t, h, fun);
		eps_e = fabs(ye - ya) / fabs(ya);
		printf("%lf\t%lf\t", ye, eps_e);

		yrk = rk4r(t0, y_0, t, h, fun);
		eps_rk = fabs(yrk - ya) / fabs(ya);
		printf("%lf\t%lf\n", yrk, eps_rk);

	}

	graphics(800, 600);
	scale(0, -5, t_k, 7);

	FILE* f = fopen("RRZ_Dane.txt", "w");

	fprintf(f, "N\th\teps_e\teps_rk\n");
	for (int j = 0; j < 7; j++)
	{
		n = pow(2., j);
		h = t_k / n;

		ya = anl(t_k);

		ye = euler(t0, y_0, t_k, h, fun);
		eps_e = fabs(ye - ya) / fabs(ya);

		yrk = rk4r(t0, y_0, t_k, h, fun);
		eps_rk= fabs(yrk - ya) / fabs(ya);

		fprintf(f,"%lf\t%lf\t%lf\t%lf\n",n,h,eps_e,eps_rk);
		setcolor(0.1);
		circle(h, log10(eps_e),2);
		setcolor(0.9);
		circle(h, log10(eps_rk),2);
	}

	fclose(f);
	wait();
}

double euler(double t, double y, double h, double(*fun)(double, double))
{

	return y + h * fun(t, y);
}

double anl(double t)
{
	return y_0 * exp(lambda * (t - t0));
}

double fun(double t, double y)
{
	return lambda* y;
}

double euler(double t0, double y_0, double t_k, double h, double(*fun)(double, double))
{
	double y_1 = y_0 + h * fun(t0, y_0);
	if (t0 + h >= t_k)
	{
		return y_1;
	}
	return euler(t0 + h, y_1, t_k, h, fun);
}

double rk4r(double t0, double y_0, double t_k, double h, double(*fun)(double, double))
{
	double y_1 = rk4(t0, y_0, h, fun);
	if (t0 + h >= t_k)
	{
		return y_1;
	}
	return rk4r(t0 + h, y_1, t_k, h, fun);
}