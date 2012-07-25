#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include "test.h"

#include "integer.h"
#include "double.h"
#include "complex.h"
#include "matrix.h"
#include "register.h"
#include "gate.h"
#include "arith.h"

STARTTEST(GCD) {
	long data[][3] = {
		{48, 18, 6},
		{19, 8, 1},
		{56, 132, 4},
		{62353, 723902, 23}
	};
	int i;
	for (i=0; i<(long)(sizeof data / sizeof*data); i++) {
		if (data[i][2] != GCD(data[i][0], data[i][1])) {
			FAIL("Invalid GCD!");
		}
	}
} ENDTEST(GCD)

STARTTEST(Near) {
	struct { double a,b; int c; } data[] = {
		{0.1, 0.0, 0},
		{0.2, 2*0.1, 1},
		{32, 2*2*2*2*2, 1}
	};
	int i;
	for (i=0; i<(long)(sizeof data/sizeof*data); i++) {
		if (Near(data[i].a, data[i].b) != data[i].c) {
			FAIL("Invalid Near()!");
		}
	}
} ENDTEST(Near)

STARTTEST(RootOfUnity) {
	Complex a = CMul(RootOfUnity(3, 2), RootOfUnity(3,1));
	if (!CNear(a, CNew(1,0))) {
		FAIL("Invalid root of unity multiplication!");
	}
	a = CMul(RootOfUnity(5,2), RootOfUnity(10,4));
	a = CMul(a, RootOfUnity(15, 3));
	if (!CNear(a, CNew(1,0))) {
		FAIL("Invalid root of unity multiplication 2!");
	}
	a = RootOfUnity(2,1);
	if (!CNear(a, CNew(-1,0))) {
		FAIL("Invalid root of unity 3!");
	}
	a = RootOfUnity(4,1);
	if (!CNear(a, CNew(0,1))) {
		FAIL("Invalid root of unity 4!");
	}
} ENDTEST(RootOfUnity)


STARTTEST(MSet) {
	Matrix m = MNew(4, 4);
	int i;
	m = MSet(m, CNew( 1, 2), CNew( 3, 4), CNew( 5, 6), CNew( 7, 8),
							CNew( 9,10), CNew(11,12), CNew(13,14), CNew(15,16),
							CNew(17,18), CNew(19,20), CNew(21,22), CNew(23,24),
							CNew(25,26), CNew(27,28), CNew(29,30), CNew(31,32));
	for (i=0; i<16; i++) {
		if (!CNear(m.d[i], CNew((double)(2*i+1), (double)(2*i+2)))) {
			FAIL("Invalid MSet()!");
		}
	}
} ENDTEST(MSet)

STARTTEST(MSca) {
	int x,y;
	Matrix m = MNew(4, 4);
	Complex c = CNew(3.5, 2);
	m = MSca(m, c);
	for (y=0; y<m.h; y++) {
		for (x=0; x<m.w; x++) {
			Complex *cel = MCel(&m, x, y);
			if (x==y) {
				if (!CNear(*cel, c)) {
					FAIL("Matrix scalar multiplication failed!");
				}
			} else {
				if (!CNear(*cel, CNew(0,0))) {
					FAIL("Matrix scalar multiplication failed!");
				}
			}
		}
	}
} ENDTEST(MSca)

STARTTEST(Pauli) {
	int runs=1000;
	while (runs-- > 0) {
		double theta = Rnd()*2*Pi;
		Matrix I = MNew(2, 2);
		Matrix Z = MSet(MNew(2, 2), CNew(1,0), CNew(0,0),
																CNew(0,0), CNew(-1,0));
		Matrix R;
		Complex direct = PhaseShift(theta);
		Complex alpha=CNew(0,0);
		Complex gamma=CNew(0,0);
		PauliPhase(theta, &alpha, &gamma);
		I = MSca(I, alpha);
		Z = MSca(Z, gamma);
		R = MAdd(I, Z);
		if (!CNear(*MCel(&R,0,0), CNew(1,0))) {
			printf("direct: %f+i%f\n", direct.re, direct.im);
			printf("alpha : %f+i%f\n", alpha.re, alpha.im);
			printf("gamma : %f+i%f\n", gamma.re, gamma.im);
			printf("matrix: %f+i%f\n", MCel(&R,0,0)->re, MCel(&R,1,1)->im);
			MPrint(I);
			MPrint(Z);
			MPrint(R);
			FAIL("Expected 1!");
		}
		if (!CNear(*MCel(&R,1,1), direct)) {
			printf("direct: %f+i%f\n", direct.re, direct.im);
			printf("alpha : %f+i%f\n", alpha.re, alpha.im);
			printf("gamma : %f+i%f\n", gamma.re, gamma.im);
			printf("matrix: %f+i%f\n", MCel(&R,0,0)->re, MCel(&R,1,1)->im);
			MPrint(I);
			MPrint(Z);
			MPrint(R);
			FAIL("Expected e^itheta!");
		}
	}
} ENDTEST(Pauli)

STARTTEST(MMul) {
	int i;
	Matrix m = MNew(4, 4);
	Matrix v = MNew(1, 4);
	Matrix v2 = MNew(1, 4);
	Matrix p;
	m = MSet(m, CNew(1,0), CNew(0,0), CNew(0,0), CNew(0,0),
							CNew(0,0), CNew(0,0), CNew(1,0), CNew(0,0),
							CNew(0,0), CNew(1,0), CNew(0,0), CNew(0,0),
							CNew(0,0), CNew(0,0), CNew(0,0), CNew(1,0));
	v = MSet(v, CNew(4.5,2.3), CNew(6.3,4.1), CNew(0,0), CNew(6.4,0));
	v2 = MSet(v2, CNew(4.5,2.3), CNew(0,0), CNew(6.3,4.1), CNew(6.4,0));
	p = MMul(m, v);
	if (p.w != 1 || p.h != 4) {
		printf("%dx%d\n", p.w, p.h);
		FAIL("Invalid matrix product size!");
	}
	for (i=0; i<4; i++) {
		if (!CNear(p.d[i], v2.d[i])) {
			FAIL("Invalid matrix product!");
		}
	}
} ENDTEST(MMul)

STARTTEST(Measure) {
	Register r = NewRegister(7, 99);
	unsigned long s = Measure(&r);
	if (s != 99) {
		FAIL("Measure failed!");
	}
	FreeRegister(&r);
} ENDTEST(Measure)

STARTTEST(Fredkin) {
	int i;
	for (i=0; i<100; i++) {
		unsigned long orig = rand()&7;
		unsigned long m;
		unsigned long new = orig;
		Register r = NewRegister(3, orig);

		if ((new&1) != 0) {
			new = (new&1) | ((new&2)<<1) | ((new&4)>>1);
		}
		GateFredkin(&r, 0,1,2);
		m = Measure(&r);
		if (m != new) {
			FAIL("Invalid result from Fredkin gate!");
		}
	}
} ENDTEST(Fredkin)

STARTTEST(Superposition) {
	int i;
	struct {int gate; int state[16];} data[] = {
		{-1, {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}},
		{ 0, {0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0}},
		{ 1, {0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0}},
		{ 2, {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}},
		{ 3, {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}},
		{ 1, {0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1}},
		{ 3, {0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1}},
	};
	Register r = NewRegister(4, 11);

	for (i=0; i<(int)(sizeof data/sizeof*data); i++) {
		int j;
		int run;
		unsigned long counts[16];
		if (data[i].gate >= 0) {
			GateHadamard(&r, data[i].gate);
		}
		for (j=0; j<16; j++) {
			counts[j] = 0;
		}
		for (run=0; run<1000; run++) {
			unsigned long state = MeasureCheat(&r);
			counts[state]++;
		}
		for (j=0; j<16; j++) {
			#if 0
			printf("%5d %d %5lu %10f%+10f\n", j, data[i].state[j], counts[j], r.state[j].re, r.state[j].im);
			#endif
			if (data[i].state[j]) {
				if (counts[j] == 0) {
					FAIL("Expected non-zero state count!");
				}
			} else {
				if (counts[j] != 0) {
					FAIL("Expected zero state count!");
				}
			}
		}
		/*puts("");*/
	}

	FreeRegister(&r);
} ENDTEST(Superposition)

STARTTEST(Entanglement) {
	int i;
	unsigned long counts[4] = {0,0,0,0};
	Register r = NewRegister(2, 0);

	GateHadamard(&r, 0);
	GateCNot(&r, 0, 1);

	for (i=0; i<1000; i++) {
		counts[MeasureCheat(&r)]++;
	}
	if (counts[0] == 0 || counts[3] == 0) {
		FAIL("Unexpected zero count!");
	}
	if (counts[1] != 0 || counts[2] != 0) {
		FAIL("Unexpected non-zero count!");
	}
} ENDTEST(Entanglement)

STARTTEST(SuperdenseCoding) {
	int i;
	int run;

	for (run=0; run<100; run++) {
		for (i=0; i<4; i++) {
			Register r = NewRegister(2, 0);
			unsigned long state;

			/* Prepare bell state */
			GateHadamard(&r, 0);
			GateCNot(&r, 0, 1);

			/* Encode qbit 0 */
			if (i == 1) {
				GatePauliZ(&r, 0);
			} else if (i == 2) {
				GatePauliX(&r, 0);
			} else if (i == 3) {
				GatePauliZ(&r, 0);
				GatePauliX(&r, 0);
			}

			/* Reverse bell state */
			GateCNot(&r, 0, 1);
			GateHadamard(&r, 0);

			state = Measure(&r);
			if (state != (unsigned long)i) {
				FAIL("Invalid superdense coding!");
			}
		}
	}
} ENDTEST(SuperdenseCoding)

STARTTEST(MeasureQubit) {
	Register r = NewRegister(4, 0);

	GateHadamard(&r, 0);
	GateHadamard(&r, 1);
	GateHadamard(&r, 2);
	GateHadamard(&r, 3);
#if 0
	{
		unsigned long m;
		PrintRegister(&r, stdout);

		m = MeasureQubit(&r, 0);
		printf("%lu\n", m);
		PrintRegister(&r, stdout);

		m = MeasureQubit(&r, 1);
		printf("%lu\n", m);
		PrintRegister(&r, stdout);

		m = MeasureQubit(&r, 2);
		printf("%lu\n", m);
		PrintRegister(&r, stdout);

		m = MeasureQubit(&r, 3);
		printf("%lu\n", m);
		PrintRegister(&r, stdout);
	}
#endif
} ENDTEST(MeasureQubit)

STARTTEST(Teleportation) {
	int run;

	for (run=0; run<100; run++) {
		unsigned long m, m0, m1;
		Register r = NewRegister(3, 0);

		/* Perturb qubit 0 */
		r.state[0].re = sqrt(0.25);
		r.state[1].re = sqrt(0.75);

		/* Prepare Bell state for shared entangled qubits 1,2 */
		GateHadamard(&r, 1);
		GateCNot(&r, 1, 2);

		/* Bell state measurement on qubit 0,1 */
		GateCNot(&r, 0, 1);
		GateHadamard(&r, 0);

		m0 = MeasureQubit(&r, 0);
		m1 = MeasureQubit(&r, 1);
		#if 0
		printf("Measured: %lu%lu\n", m1,m0);
		#endif
		m = (m1<<1) | m0;

		if (m == 1) {
			GatePhase(&r, Pi, 2);
		} else if (m == 2) {
			GateNot(&r, 2);
		} else if (m == 3) {
			GateNot(&r, 2);
			GatePhase(&r, Pi, 2);
		}

		/* Check if qubit 0 teleported to qubit 2 */
		if (!CNear(r.state[m], CNew(sqrt(0.25),0)) || !CNear(r.state[m|4], CNew(sqrt(0.75),0))) {
			FAIL("Quantum teleportation failed!");
		}
	}
} ENDTEST(Teleportation)

STARTTEST(CCNot) {
	int i;
	for (i=0; i<1000; i++) {
		unsigned int start = rand()&7;
		unsigned int end;
		Register r = NewRegister(3, start);
		end = MeasureCheat(&r);
		if (start != end) {
			FAIL("Invalid Pre-CCNot measurement!");
		}
		GateCCNot(&r, 2,1,0);
		end = MeasureCheat(&r);
		if (end == 7) {
			end = 6;
		} else if (end == 6) {
			end = 7;
		}
		if (start != end) {
			FAIL("Invalid CCNot measurement!");
		}
	}
} ENDTEST(CCNot)

STARTTEST(Increment) {
	Register r = NewRegister(10, 0);
	unsigned long m;

	m = MeasureCheat(&r);
	if (m != 0) {
		printf("m=%lu\n", m);
		FAIL("Expected 0!");
	}

	GateCNot(&r, 0, 1);
	GateNot(&r, 0);
	m = MeasureCheat(&r);
	if (m != 1) {
		FAIL("Expected 1!");
	}

	GateCNot(&r, 0, 1);
	GateNot(&r, 0);
	m = MeasureCheat(&r);
	if (m != 2) {
		FAIL("Expected 2!");
	}

	GateCNot(&r, 0, 1);
	GateNot(&r, 0);
	m = MeasureCheat(&r);
	if (m != 3) {
		FAIL("Expected 3!");
	}
} ENDTEST(Increment)

STARTTEST(Add) {
	int runs=100;
	for (; runs>0; runs--) {
		unsigned long a = rand()&((1<<5)-1);
		unsigned long b = rand()&((1<<5)-1);
		unsigned long m;
		Register r = NewRegister(12, a);
		SubReg X = NewSubReg(5, 4,3,2,1,0);
		SubReg S = NewSubReg(6, 10,9,8,7,6,5);

		ArithAddConst(&r, X, b, S, 11);

		m = MeasureSubReg(&r, X);
		if (m != a) {
			FAIL("Failed to re-measure original X");
		}
		m = MeasureSubReg(&r, S);
		if (m != a+b) {
			FAIL("Addition failed");
		}
		m = MeasureQubit(&r, 11);
		if (m != 0) {
			FAIL("Zero bit is not zero");
		}
	}
} ENDTEST(Add)

int main(void) {
	RUNTEST(GCD);
	RUNTEST(Near);
	RUNTEST(RootOfUnity);
	RUNTEST(MSet);
	RUNTEST(MSca);
	RUNTEST(Pauli);
	RUNTEST(MMul);
	RUNTEST(Measure);
	RUNTEST(Fredkin);
	RUNTEST(Superposition);
	RUNTEST(Entanglement);
	RUNTEST(SuperdenseCoding);
	RUNTEST(MeasureQubit);
	RUNTEST(Teleportation);
	RUNTEST(CCNot);
	RUNTEST(Increment);
	RUNTEST(Add);
	return 0;
}


