#include <assert.h>

#include "register.h"
#include "gate.h"

#include "arith.h"

void ArithConstAdd1b(Register *r, unsigned long x, unsigned long s, unsigned long c, unsigned long z) {
	GateCNot(r, x, s);
	GateCNot(r, c, s);
	GateCCNot(r, x, c, z);
	GateSwap(r, c, z);
	GateCNot(r, s, z);
	GateCNot(r, x, z);
}

void ArithConstAdd1bFinish(Register *r, unsigned long s, unsigned long c) {
	GateCNot(r, s, c);
	GateNot(r, s);
}

void ArithAddConst(Register *r, SubReg x, unsigned long y, SubReg s, unsigned long z) {
	unsigned long c = s.qb[0];
	unsigned long i;
	assert(s.N == x.N+1);
	for (i=0; i<x.N; i++) {
		ArithConstAdd1b(r, x.qb[x.N-1-i], s.qb[s.N-1-i], c, z);
		if (y & (1<<i)) {
			ArithConstAdd1bFinish(r, s.qb[s.N-1-i], c);
		}
	}
}

/* A new quantum ripple-carry addition circuit (2008)
 * Steven A. Cuccaro, Thomas G. Draper, Samuel A. Kutin, David Petrie Moulton
 * http://arxiv.org/pdf/quant-ph/0410184v1.pdf */
static void MAJ(Register *r, unsigned long c, unsigned long b, unsigned long a) {
	GateCNot(r, a, b);
	GateCNot(r, a, c);
	GateCCNot(r, c, b, a);
}
static void UMS(Register *r, unsigned long c, unsigned long b, unsigned long a) {
	GateCCNot(r, c, b, a);
	GateCNot(r, a, c);
	GateCNot(r, c, b);
}
void ArithAddV(Register *r, SubReg A, SubReg B, unsigned long g, unsigned long c) {
	int i;
	r=r; g=g;
	assert(A.N == B.N);
	MAJ(r, g, B.qb[0], A.qb[0]);
	for (i=1; i<(int)A.N; i++) {
		MAJ(r, A.qb[i-1], B.qb[i], A.qb[i]);
	}
	GateCNot(r, A.qb[A.N-1], c);
	for (i=A.N-1; i>0; i--) {
		UMS(r, A.qb[i-1], B.qb[i], A.qb[i]);
	}
	UMS(r, g, B.qb[0], A.qb[0]);
}



