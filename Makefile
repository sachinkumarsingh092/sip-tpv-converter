.DEFAULT_GOAL := compile

CC := gcc
CC := ${CC}
CFLAGS := -Wall -O3 -g
INCLUDES := -I/usr/local/include
LIBS := /usr/local/lib/libgnuastro.so.10  /usr/local/lib/libwcs-7.1.a \
	   /usr/local/lib/libcfitsio.so.8 /usr/local/lib/libgsl.a /usr/local/lib/libgslcblas.a -lm

# See here for examples: https://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
# $< is patterned to match prerequisites and $@ matches the target.


compile: temp-wcs.c
	${CC} ${CFLAGS} temp-wcs.c ${INCLUDES} -o temp-wcs -pthread ${LIBS} && ./temp-wcs

clean_compile: temp-wcs.c
	rm ./test-fits/test-sip.fits && ${CC} ${CFLAGS} temp-wcs.c ${INCLUDES} -o temp-wcs -pthread ${LIBS} && ./temp-wcs

with_valgrind: temp-wcs.c
	${CC} ${CFLAGS} temp-wcs.c ${INCLUDES} -o temp-wcs -pthread ${LIBS} && valgrind ./temp-wcs

run_bscript: 
	./scripts/genpveq.sh
