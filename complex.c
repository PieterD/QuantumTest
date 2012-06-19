#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "integer.h"
#include "double.h"

#include "complex.h"

Complex CNew(double re, double im) {
	Complex c;
	c.re = re;
	c.im = im;
	return c;
}

Complex CAdd(Complex a, Complex b) {
	return CNew(a.re+b.re, a.im+b.im);
}

Complex CSca(Complex a, double b) {
	return CNew(a.re*b, a.im*b);
}

Complex CMul(Complex a, Complex b) {
	return CNew(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re);
}

Complex CCon(Complex a) {
	return CNew(a.re, -a.im);
}

int CNear(Complex a, Complex b) {
	if (!Near(a.re, b.re) || !Near(a.im, b.im)) {
		return 0;
	}
	return 1;
}

Complex PhaseShift(double theta) {
	return CNew(cos(theta), sin(theta));
}

Complex RootOfUnity(long n, long k) {
	long gcd = GCD(n, k);
	double theta;
	n /= gcd;
	k /= gcd;

	theta = 2*Pi*(double)(k)/(double)(n);
	return CNew(cos(theta), sin(theta));
}

void PauliPhase(double theta, Complex *alpha, Complex *delta) {
	Complex p = PhaseShift(theta);
	p = CAdd(p, CNew(-1,0));
	p = CSca(p, -0.5);
	*delta = p;
	p = CSca(p, -1);
	p = CAdd(p, CNew(1,0));
	*alpha = p;
}



