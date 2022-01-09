OBJ=deps/obl/oblas_lite.o rs.o

TEST_UTILS=\
t/00util/test\
t/00util/bench

CFLAGS   = -O3 -g -std=c11 -Wall -I. -Ideps/obl -fno-inline -Wvla
CFLAGS  += -march=native -funroll-loops -ftree-vectorize 

all: $(TEST_UTILS)

t/00util/test: t/00util/test.o $(OBJ)

t/00util/bench: t/00util/bench.o $(OBJ)

test: CPPFLAGS+=-D_DEFAULT_SOURCE
test: clean $(TEST_UTILS)
	prove -I. -v t/*.t

clean:
	$(RM) *.o *.a $(TEST_UTILS) $(OBJ)

indent:
	find -name '*.[h,c]' | xargs clang-format -i

scan:
	scan-build $(MAKE) clean all

gperf: LDLIBS = -lprofiler -ltcmalloc
gperf: clean t/00util/bench
	CPUPROFILE_FREQUENCY=100000000 CPUPROFILE=gperf.prof ./t/00util/bench 200 20 512
	pprof ./t/00util/bench gperf.prof --callgrind > callgrind.gperf
	gprof2dot --format=callgrind callgrind.gperf -z main | dot -T svg > gperf.svg

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=undefined,implicit-conversion
ubsan: LDLIBS += -lubsan
ubsan: clean t/00util/test
	./t/00util/test
