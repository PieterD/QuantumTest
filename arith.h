#ifndef ARITH_H_
#define ARITH_H_

void ArithConstAdd1b(Register *r, unsigned long x, unsigned long s, unsigned long c, unsigned long z);
void ArithConstAdd1bFinish(Register *r, unsigned long s, unsigned long c);
void ArithAddConst(Register *r, SubReg x, unsigned long y, SubReg s, unsigned long z);
void ArithAddV(Register *r, SubReg A, SubReg B, unsigned long g, unsigned long c);

#endif

