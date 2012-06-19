#ifndef GATE_H_
#define GATE_H_

#include "matrix.h"
#include "register.h"

void Gate1(Register *r, Matrix m, unsigned long qb);
void Gate2(Register *r, Matrix m, unsigned long qb1, unsigned long qb0);
void GateC(Register *r, Matrix sm, unsigned long cb, unsigned long qb);

void GateHadamard(Register *r, unsigned long qb);
void GatePauliX(Register *r, unsigned long qb);
void GatePauliY(Register *r, unsigned long qb);
void GatePauliZ(Register *r, unsigned long qb);

void GateSwap(Register *r, unsigned long qb1, unsigned long qb0);

void GatePhase(Register *r, double theta, unsigned long qb);
void GateCPhase(Register *r, double theta, unsigned long qb1, unsigned long qb0);

void GateNot(Register *r, unsigned long qb);
void GateCNot(Register *r, unsigned long qb1, unsigned long qb0);
void GateCCNot(Register *r, unsigned long qb2, unsigned long qb1, unsigned long qb0);

#endif

