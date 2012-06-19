#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>

#include "complex.h"

#include "matrix.h"

Matrix MNew(int w, int h) {
	Matrix m;
	int x,y;
	m.w = w;
	m.h = h;
	#if 0
	m.d = malloc(sizeof*m.d * w*h);
	assert(m.d!=NULL);
	#endif
	for (y=0; y<h; y++) {
		for (x=0; x<w; x++) {
			Complex *c = MCel(&m, x, y);
			if (x == y) {
				*c = CNew(1,0);
			} else {
				*c = CNew(0,0);
			}
		}
	}
	return m;
}

Complex* MCel(Matrix *m, int x, int y) {
	assert(x<m->w);
	assert(y<m->h);
	return &m->d[y*m->w+x];
}

Matrix MSet(Matrix m, ...) {
	va_list ap;
	int x,y;
	va_start(ap, m);
	for (y=0; y<m.h; y++) {
		for (x=0; x<m.w; x++) {
			Complex src = va_arg(ap, Complex);
			Complex *dst = MCel(&m, x, y);
			*dst = src;
		}
	}
	va_end(ap);
	return m;
}

Matrix MAdd(Matrix a, Matrix b) {
	Matrix m = MNew(a.w, a.h);
	int x,y;
	assert(a.w==b.w);
	assert(a.h==b.h);
	for (y=0; y<m.h; y++) {
		for (x=0; x<m.w; x++) {
			Complex *dst = MCel(&m, x, y);
			Complex *ca = MCel(&a, x, y);
			Complex *cb = MCel(&b, x, y);
			*dst = CAdd(*ca, *cb);
		}
	}
	return m;
}

Matrix MSca(Matrix a, Complex b) {
	Matrix m = MNew(a.w, a.h);
	int x,y;
	for (y=0; y<a.h; y++) {
		for(x=0; x<a.w; x++) {
			Complex *dst = MCel(&m, x, y);
			Complex *src = MCel(&a, x, y);
			*dst = CMul(*src, b);
		}
	}
	return m;
}

Matrix MMul(Matrix a, Matrix b) {
	Matrix m;
	int x,y;
	assert(a.w==b.h);
	m = MNew(b.w, a.h);
	for (y=0; y<a.h; y++) {
		for (x=0; x<b.w; x++) {
			int n;
			Complex *sum = MCel(&m, x, y);
			#if 0
			printf("[%d,%d] = ", x, y);
			#endif
			*sum = CNew(0,0);
			for (n=0; n<a.w; n++) {
				Complex *aa = MCel(&a, n, y);
				Complex *bb = MCel(&b, x, n);
				#if 0
				printf("+[%d,%d]*[%d,%d] ", n, y, x, n);
				#endif
				*sum = CAdd(*sum, CMul(*aa, *bb));
			}

			#if 0
			printf("\n");
			#endif
		}
	}
	return m;
}

void MPrint(Matrix m) {
	int x,y;
	for (y=0; y<m.h; y++) {
		for (x=0; x<m.w; x++) {
			Complex *c = MCel(&m, x,y);
			printf("[%.5f+i%.5f] ", c->re, c->im);
		}
		printf("\n");
	}
}

