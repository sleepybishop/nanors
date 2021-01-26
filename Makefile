OBJ=rs.o

CPPFLAGS = -D_DEFAULT_SOURCE
CFLAGS   = -O3 -g -std=c11 -Wall -I. -fno-inline -Wvla
CFLAGS  += -march=native -funroll-loops -ftree-vectorize 

all: test bench 

test: test.o $(OBJ)

bench: bench.o $(OBJ)

clean:
	$(RM) test *.o *.a

indent:
	clang-format -style=LLVM -i *.c *.h

scan:
	scan-build $(MAKE) clean all

ubsan: CC=clang
ubsan: CFLAGS += -fsanitize=undefined,implicit-conversion,integer
ubsan: LDFLAGS += -lubsan
ubsan: clean test
	./test

gperf: LDFLAGS = -lprofiler -ltcmalloc
gperf: clean bench
	CPUPROFILE_FREQUENCY=100000000 CPUPROFILE=gperf.prof ./bench 200 20 512
	pprof ./bench gperf.prof --callgrind > callgrind.gperf
	gprof2dot --format=callgrind callgrind.gperf -z main | dot -T svg > gperf.svg

