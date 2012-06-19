#ifndef REGISTER_H_
#define REGISTER_H_

#include <stdio.h>
#include "complex.h"
#include "matrix.h"

typedef struct {
	unsigned long N, statenum;
	Complex *state;
} Register;

typedef struct {
	unsigned long N;
	unsigned char qb[64];
} SubReg;

Register NewRegister(unsigned long N, unsigned long base);
SubReg NewSubReg(unsigned long N, ...);
void FreeRegister(Register *r);
void PrintRegister(Register *r, FILE *out);

unsigned long Measure(Register *r);
unsigned long MeasureCheat(Register *r);
unsigned long MeasureQubit(Register *r, unsigned long qb);
unsigned long MeasureSubReg(Register *r, SubReg sr);

#endif

