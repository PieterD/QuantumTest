#ifndef COMPLEX_H_
#define COMPLEX_H_

typedef struct {
	double re,im;
} Complex;

Complex CNew(double re, double im);
Complex CAdd(Complex a, Complex b);
Complex CSca(Complex a, double b);
Complex CMul(Complex a, Complex b);
Complex CCon(Complex a);
int CNear(Complex a, Complex b);
Complex PhaseShift(double theta);
Complex RootOfUnity(long n, long k);
void PauliPhase(double theta, Complex *alpha, Complex *beta);

#endif

