#include <assert.h>
#include <math.h>

#include "double.h"
#include "complex.h"
#include "matrix.h"
#include "register.h"

#include "gate.h"

/* Generic gates */
void Gate1(Register *r, Matrix m, unsigned long qb) {
	unsigned long i;
	unsigned long himask = ~0ul<<qb;
	unsigned long lomask = ~himask;
	unsigned long bit = himask & ~((~lomask)<<1);

	Matrix v = MNew(1, 2);
	Matrix mv = MNew(1, 2);
	Complex *v1 = MCel(&v, 0, 0);
	Complex *v2 = MCel(&v, 0, 1);
	Complex *mv1 = MCel(&mv, 0, 0);
	Complex *mv2 = MCel(&mv, 0, 1);

	assert(m.w==2);
	assert(m.h==2);
	assert(qb<r->statenum);

	#if 0
	printf("%08lX %08lX %08lX\n", himask, lomask, bit);
	#endif

	for (i=0; i<r->statenum/2; i++) {
		unsigned long z = ((i&himask)<<1) | (i&lomask);
		unsigned long o = z | bit;
		*v1 = r->state[z];
		*v2 = r->state[o];
		mv = MMul(m, v);
		r->state[z] = *mv1;
		r->state[o] = *mv2;
		#if 0
		printf("%lu %lu [%f%+fi %f%+fi\n", z, o, v1->re, v1->im, v2->re, v2->im);
		#endif
	}
}

void Gate2(Register *r, Matrix m, unsigned long qb1, unsigned long qb0) {
	unsigned long i;
	unsigned long qbhi = qb0>qb1?qb0:qb1;
	unsigned long qblo = qb0<qb1?qb0:qb1;
	unsigned long himask = ~0ul<<qbhi;
	unsigned long mdmask = ~0ul<<qblo;
	unsigned long bit0 = 1<<qb0;
	unsigned long bit1 = 1<<qb1;

	Matrix v = MNew(1, 4);
	Matrix mv = MNew(1, 4);
	Complex *v1 = MCel(&v, 0, 0);
	Complex *v2 = MCel(&v, 0, 1);
	Complex *v3 = MCel(&v, 0, 2);
	Complex *v4 = MCel(&v, 0, 3);
	Complex *mv1 = MCel(&mv, 0, 0);
	Complex *mv2 = MCel(&mv, 0, 1);
	Complex *mv3 = MCel(&mv, 0, 2);
	Complex *mv4 = MCel(&mv, 0, 3);

	assert(m.w==4);
	assert(m.h==4);
	assert(r->N>=2);
	assert(qbhi<r->statenum);
	assert(qblo<qbhi);

	#if 0
	printf("%08lX %08lX %08lX %08lX %lu %lu\n", himask, mdmask, bit0, bit1, qb0, qb1);
	#endif

	for (i=0; i<r->statenum/4; i++) {
		unsigned long s00,s01,s10,s11;
		s00 = i;
		s00 = (s00&~mdmask)|(s00&mdmask)<<1;
		s00 = (s00&~himask)|(s00&himask)<<1;
		s01 = s00 | bit0;
		s10 = s00 | bit1;
		s11 = s00 | bit0 | bit1;
		*v1 = r->state[s00];
		*v2 = r->state[s01];
		*v3 = r->state[s10];
		*v4 = r->state[s11];
		mv = MMul(m, v);
		r->state[s00] = *mv1;
		r->state[s01] = *mv2;
		r->state[s10] = *mv3;
		r->state[s11] = *mv4;
		#if 0
		printf("%ld:", i);
		printf(" %s", Binary(s00, buf, 3));
		printf(" %s", Binary(s01, buf, 3));
		printf(" %s", Binary(s10, buf, 3));
		printf(" %s", Binary(s11, buf, 3));
		printf("\n");
		#endif
	}
}

void GateC(Register *r, Matrix sm, unsigned long cb, unsigned long qb) {
	Matrix m = MNew(4, 4);
	Complex *sm00 = MCel(&sm, 0, 0);
	Complex *sm01 = MCel(&sm, 0, 1);
	Complex *sm10 = MCel(&sm, 1, 0);
	Complex *sm11 = MCel(&sm, 1, 1);

	assert(sm.w==2);
	assert(sm.h==2);

	m = MSet(m, CNew(1,0), CNew(0,0), CNew(0,0), CNew(0,0),
							CNew(0,0), CNew(1,0), CNew(0,0), CNew(0,0),
							CNew(0,0), CNew(0,0), *sm00, *sm01,
							CNew(0,0), CNew(0,0), *sm10, *sm11);

	Gate2(r, m, cb, qb);
}

/* Specific gates */
void GateHadamard(Register *r, unsigned long qb) {
	Matrix hadamard = MSca(MSet(MNew(2,2),
															CNew(1,0), CNew(1,0),
															CNew(1,0), CNew(-1,0)),
												 CNew(1.0/sqrt(2), 0));
	Gate1(r, hadamard, qb);
}

void GatePauliX(Register *r, unsigned long qb) {
	Matrix x = MSet(MNew(2,2),
									CNew(0,0), CNew(1,0),
									CNew(1,0), CNew(0,0));
	Gate1(r, x, qb);
}

void GatePauliY(Register *r, unsigned long qb) {
	Matrix y = MSet(MNew(2,2),
									CNew(0,0), CNew(0,-1),
									CNew(0,1), CNew(0,0));
	Gate1(r, y, qb);
}

void GatePauliZ(Register *r, unsigned long qb) {
	Matrix z = MSet(MNew(2,2),
									CNew(1,0), CNew(0,0),
									CNew(0,0), CNew(-1,0));
	Gate1(r, z, qb);
}

void GateSwap(Register *r, unsigned long qb1, unsigned long qb0) {
	GateCNot(r, qb1, qb0);
	GateCNot(r, qb0, qb1);
	GateCNot(r, qb1, qb0);
}


void GatePhase(Register *r, double theta, unsigned long qb) {
	Matrix phase = MSet(MNew(2,2),
											CNew(1,0), CNew(0,0),
											CNew(0,0), PhaseShift(theta));
	Gate1(r, phase, qb);
}

void GateCPhase(Register *r, double theta, unsigned long qb1, unsigned long qb0) {
	Matrix phase = MSet(MNew(2,2),
											CNew(1,0), CNew(0,0),
											CNew(0,0), PhaseShift(theta));
	GateC(r, phase, qb1, qb0);
}

void GateNot(Register *r, unsigned long qb) {
	GatePauliX(r, qb);
}

void GateCNot(Register *r, unsigned long qb1, unsigned long qb0) {
	Matrix not = MSet(MNew(2,2),
										CNew(0,0), CNew(1,0),
										CNew(1,0), CNew(0,0));
	GateC(r, not, qb1, qb0);
}

void GateCCNot(Register *r, unsigned long qb2, unsigned long qb1, unsigned long qb0) {
	Matrix not = MSet(MNew(2,2),
										CNew(0,0), CNew(1,0),
										CNew(1,0), CNew(0,0));
	Matrix hadamard = MSca(MSet(MNew(2,2),
															CNew(1,0), CNew(1,0),
															CNew(1,0), CNew(-1,0)),
												 CNew(1.0/sqrt(2), 0));
	Matrix phasehalf = MSet(MNew(2,2),
													CNew(1,0), CNew(0,0),
													CNew(0,0), PhaseShift(Pi/2));
	Matrix phasetriphalf = MSet(MNew(2,2),
															CNew(1,0), CNew(0,0),
															CNew(0,0), PhaseShift(3*Pi/2));

	Gate1(r, hadamard, qb0);
	GateC(r, phasehalf, qb1, qb0);
	GateC(r, not, qb2, qb1);
	GateC(r, phasetriphalf, qb1, qb0);
	GateC(r, not, qb2, qb1);
	GateC(r, phasehalf, qb2, qb0);
	Gate1(r, hadamard, qb0);
}


