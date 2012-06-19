#ifndef MATRIX_H_
#define MATRIX_H_

#include "complex.h"

typedef struct {
	int w,h;
	/*Complex *d;*/
	Complex d[16];
} Matrix;

Matrix MNew(int w, int h);
Matrix MSet(Matrix m, ...);
Complex* MCel(Matrix *m, int x, int y);
Matrix MAdd(Matrix a, Matrix b);
Matrix MSca(Matrix a, Complex b);
Matrix MMul(Matrix a, Matrix b);
void MPrint(Matrix m);

#endif

