#include <stdlib.h>
#include <time.h>

int Near(double a, double b) {
	double neardiff = 0.00000001;
	if (a-neardiff < b && a+neardiff > b) {
		return 1;
	}
	return 0;
}

double Rnd(void) {
	static int doseed=1;
	if (doseed) {
		srand((unsigned int)time(NULL));
		doseed=0;
	}
	return (double)rand()/(double)RAND_MAX;
}

