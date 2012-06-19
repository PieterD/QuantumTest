#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "integer.h"
#include "double.h"
#include "complex.h"
#include "matrix.h"
#include "gate.h"

#include "register.h"

Register NewRegister(unsigned long N, unsigned long base) {
	Register r;
	unsigned long i;
	assert(N>0);
	r.N = N;
	r.statenum = 1<<N;
	assert(base<r.statenum);
	r.state = malloc(sizeof*r.state * r.statenum);
	assert(r.state!=NULL);
	for (i=0; i<r.statenum; i++) {
		r.state[i] = CNew(0,0);
	}
	r.state[base] = CNew(1,0);
	return r;
}

SubReg NewSubReg(unsigned long N, ...) {
	va_list ap;
	unsigned long i;
	SubReg sr;
	sr.N = N;
	va_start(ap, N);
	for (i=0; i<sizeof sr.qb / sizeof*sr.qb; i++) {
		if (i<N) {
			sr.qb[i] = va_arg(ap, unsigned long);
		}
	}
	va_end(ap);
	return sr;
}

void FreeRegister(Register *r) {
	if (r != NULL) {
		free(r->state);
		r->state = NULL;
	}
}

void PrintRegister(Register *r, FILE *out) {
	unsigned long i;
	char buf[16];
	for (i=0; i<r->statenum; i++) {
		fprintf(out, "%7lu[%s]: %10f%+10fi\n", i, Binary(i, buf, r->N), r->state[i].re, r->state[i].im);
	}
	fprintf(out, "\n");
}

unsigned long Measure(Register *r) {
	unsigned long i;
	unsigned long m = MeasureCheat(r);

	for (i=0; i<r->statenum; i++) {
		r->state[i] = CNew(0,0);
	}
	r->state[m] = CNew(1,0);

	return m;
}
unsigned long MeasureCheat(Register *r) {
	unsigned long m = 0;
	unsigned long i;
	double rnd = Rnd();
	Complex sum = CNew(0,0);
	for (i=0; i<r->statenum; i++) {
		Complex a = r->state[i];
		Complex b = CCon(a);
		Complex c = CMul(a, b);
		Complex newsum = CAdd(sum, c);
		if (rnd >= sum.re && rnd < newsum.re) {
			m = i;
		}
		sum = newsum;
	}
	assert(CNear(sum, CNew(1,0)));

	return m;
}

unsigned long MeasureQubit(Register *r, unsigned long qb) {
	unsigned long i;
	unsigned long himask = ~0ul<<qb;
	unsigned long lomask = ~himask;
	unsigned long bit = 1<<qb;
	double norm;
	unsigned long m;
	Complex sum0 = CNew(0, 0);
	Complex sum1 = CNew(0, 0);
	Complex sum2;
	for (i=0; i<r->statenum/2; i++) {
		unsigned long z = ((i&himask)<<1) | (i&lomask);
		unsigned long o = z | bit;
		Complex a = r->state[z];
		Complex b = CCon(a);
		Complex c = CMul(a, b);
		sum0 = CAdd(sum0, c);

		a = r->state[o];
		b = CCon(a);
		c = CMul(a, b);
		sum1 = CAdd(sum1, c);
	}

	sum2 = CAdd(sum0, sum1);
	#if 0
	printf("sum0=%f, sum1=%f, sum2=%f\n", sum0.re, sum1.re, sum2.re);
	#endif
	assert(CNear(sum2,CNew(1,0)));
	
	if (Rnd()<sum0.re) {
		m = 0;
		norm = sum0.re;
	} else {
		m = 1;
		norm = sum1.re;
	}

	/* Normalize */
	norm = 1/sqrt(norm);

	for (i=0; i<r->statenum/2; i++) {
		unsigned long s = ((i&himask)<<1) | (i&lomask) | (bit*m);
		unsigned long n = s ^ bit;
		r->state[s] = CSca(r->state[s], norm);
		r->state[n] = CNew(0,0);
	}
	return m;
}

unsigned long MeasureSubReg(Register *r, SubReg sr) {
	unsigned long m = 0;
	unsigned long i;
	for (i=0; i<sr.N; i++) {
		m <<= 1;
		m |= MeasureQubit(r, sr.qb[i]);
	}
	return m;
}


