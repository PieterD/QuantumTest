#ifndef TEST_H_
#define TEST_H_

#include <stdio.h>

#define STARTTEST(NAME) \
	static int Test ## NAME (void) {

#define ENDTEST(NAME) \
	return 0; \
	}

#define RUNTEST(NAME) \
	LOG("Testing " #NAME); \
	if (Test ## NAME ()) \
		LOG("FAIL"); \
	else \
		LOG("SUCCESS");

#define LOG(STRING) \
	fprintf(stderr, "log %s\n", STRING)

#define FAIL(STRING) \
	LOG(STRING); \
	return -1;


#endif
