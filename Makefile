SHELL=/bin/bash
EDITOR=vim
CC=icc
export LD_LIBRARY_PATH=/usr/local/lib
CFLAGS_DEBUG=-ansi -pedantic -W -Wall -Wdouble-promotion -Wformat -Winline\
	   -Wmissing-prototypes -Wstrict-prototypes\
	   -Wshadow -Wpointer-arith\
	   -Wcast-qual -Wcast-align\
	   -Wwrite-strings \
	   -fshort-enums -fno-common -Dinline= -g
CFLAGS=-O3 
#CPPFLAGS=-I/home/encad/berthoumieux/OpenCL/ocl_4.1/include
#LDFLAGS=-L/home/encad/berthoumieux/OpenCL/ocl_4.1/lib -lOpenCL -loclUtil_x86_64 -lshrutil_x86_64 -lgmp

LDFLAGS=-L/usr/local/lib -L/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64/ -lm  cblas -lmkl_rt
EXECUTABLE=test

#
# Compilation options :
#
mklbenchmark:
	icc -D BENCHMARK $(CFLAGS) -o test main.c -lmkl_rt 

mkl:
	icc $(CFLAGS) -o test main.c -lmkl_rt  -lpthread

mkldebug:
	icc -DDEBUG $(CFLAGS_DEBUG) -o test main.c -lmkl_rt 

all: compile_production clean_production exec

debug: compile_debug clean_debug exec

benchmark: compile_bench clean_bench exec

testing: compile_production clean_production exec_testing

#
# Compilation :
#
main.prod: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@ 

main.debug: main.c
	$(CC) $(CFLAGS_DEBUG) $(CPPFLAGS) -D DEBUG -c $^ -o $@ 

main.bench: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -D BENCHMARK -c $^ -o $@ 

#
# Linking :
#
compile_production: main.prod
	$(CC) $^ -o $(EXECUTABLE)  $(LDFLAGS)

compile_debug: main.debug
	$(CC) $^ -o $(EXECUTABLE)  $(LDFLAGS)

compile_bench: main.bench
	$(CC) $^ -o $(EXECUTABLE)  $(LDFLAGS)

#
# Cleaning :
#
clean_production: 
	rm main.prod

clean_debug: 
	rm main.debug

clean_bench: 
	rm main.bench

#
# Execution :
#
exec:
	./$(EXECUTABLE) $(ARG)
exec_testing:
	./$(EXECUTABLE) 1 1 1

#
# Get preprocessor output (make pp) :
#
pp: compile.pp echo

compile.pp: 
	cp main.c main_pp.c;
	$(CC) $(CFLAGS) $(CPPFLAGS)  -E main_pp.c > preprocessor.c 

echo: 
	$(EDITOR) preprocessor.c
	rm preprocessor.c
	rm main_pp.c
