#include "integer.h"

long GCD(long a, long b) {
	if (a < b) {
		a ^= b;
		b ^= a;
		a ^= b;
	}
	if (b == 0) {
		return a;
	}
	return GCD(b, a%b);
}

char *Binary(unsigned long a, char *buf, int len) {
	unsigned long mask = 1<<(len-1);
	int i = 0;
	for (; mask>0; mask>>=1) {
		if (a&mask) {
			buf[i++] = '1';
		} else {
			buf[i++] = '0';
		}
	}
	buf[len] = '\0';
	return buf;
}

