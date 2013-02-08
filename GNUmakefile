#
# GNUmakefile for Galois field library
#
#

SRCS = gf_w4.c gf_w8.c gf_w16.c gf_w32.c gf_w64.c gf_w128.c gf_wgen.c gf.c gf_unit.c gf_time.c gf_mult.c gf_method.c gf_54.c gf_methods.c gf_div.c gf_rand.c gf_general.c gf_poly.c
HDRS = gf_complete.h gf_int.h
EXECUTABLES = gf_mult gf_div gf_unit gf_time gf_54 gf_methods gf_poly
CFLAGS = -O3 -msse4 -DINTEL_SSE4
LDFLAGS = -O3 -msse4
RM = /bin/rm -f

OBJS = $(addsuffix .o, $(basename $(SRCS)))

DEFAULT = $(EXECUTABLES)

default: $(DEFAULT)

all: $(OBJS)

gf_methods: gf_methods.o gf.o gf_method.o gf_wgen.o gf_w4.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o
gf_time: gf_time.o gf.o gf_method.o gf_wgen.o gf_w4.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o gf_rand.o gf_general.o
gf_unit: gf_unit.o gf.o gf_method.o gf_wgen.o gf_w4.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o gf_rand.o gf_general.o
gf_mult: gf_mult.o gf.o gf_wgen.o gf_w4.o gf_method.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o
gf_div: gf_div.o gf.o gf_wgen.o gf_w4.o gf_method.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o
gf_54: gf_54.o gf.o gf_wgen.o gf_w4.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o
gf_poly: gf_poly.o gf.o gf_method.o gf_w4.o gf_w8.o gf_w16.o gf_w32.o gf_w64.o gf_w128.o gf_rand.o gf_wgen.o

clean:
	$(RM) $(OBJS) gf_div.c

spotless: clean
	$(RM) *~ $(EXECUTABLES)

gf_div.o: gf_complete.h gf_method.h
gf_methods.o: gf_complete.h gf_method.h
gf_time.o: gf_complete.h gf_method.h gf_rand.h gf_general.h
gf_wgen.o: gf_int.h gf_complete.h
gf_w4.o: gf_int.h gf_complete.h
gf_w8.o: gf_int.h gf_complete.h
gf_w16.o: gf_int.h gf_complete.h
gf_w32.o: gf_int.h gf_complete.h
gf_w64.o: gf_int.h gf_complete.h
gf_54.o: gf_complete.h
gf_unit.o: gf_complete.h gf_method.h gf_rand.h gf_general.h
gf_general.o: gf_complete.h gf_int.h gf_general.h gf_rand.h
gf_mult.o: gf_complete.h gf_method.h
gf_method.o: gf_complete.h

gf_div.c: gf_mult.c
	sed 's/multiply/divide/g' gf_mult.c > gf_div.c
