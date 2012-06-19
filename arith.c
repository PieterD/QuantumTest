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




