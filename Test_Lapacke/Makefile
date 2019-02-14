CC=gcc
CFLAGS=-I/opt/lapack-3.8.0/include
LDFLAGS=-L/opt/lapack-3.8.0/lib -llapacke -lcblas -llapack -lblas -lgfortran -lm

all: cgeqrf dgetrs test_qr

cgeqrf:cgeqrf.c
	$(CC) -O3 -o $@ $< $(CFLAGS) $(LDFLAGS)

dgetrs:dgetrs.c
	$(CC) -O3 -o $@ $< $(CFLAGS) $(LDFLAGS)

test_qr:test_qr.c
	$(CC) -O3 -o $@ $< $(CFLAGS) $(LDFLAGS)

.PHONY:clean
clean:
	-rm cgeqrf dgetrs test_qr
