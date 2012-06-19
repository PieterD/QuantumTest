all: quantum

FILES=integer.c double.c complex.c matrix.c register.c gate.c arith.c

quantum: main.c $(FILES)
	gcc $^ -o $@ -lm -W -Wall -Werror -ansi -pedantic -O2

test: test.c $(FILES)
	gcc $^ -o $@ -lm -W -Wall -Werror -ansi -pedantic -O2
	./test

clean:
	-rm test quantum

rebuild: clean all

.PHONY: clean rebuild
