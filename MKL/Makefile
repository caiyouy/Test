CC=gcc
CFLAGS=  -m64 -I${MKLROOT}/include
LDFLAGS=  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

test_qr_mkl:test_qr_mkl.c
	${CC} -O3 -o $@ $< ${LDFLAGS} ${CFLAGS}