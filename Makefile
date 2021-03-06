HELL=/bin/bash
EDITOR=vim
CC=icc
export LD_LIBRARY_PATH=/usr/local/lib
CFLAGS_DEBUG=-ansi -pedantic -W -Wall -Wdouble-promotion -Wformat -Winline\
			 -Wmissing-prototypes -Wstrict-prototypes\
			 -Wshadow -Wpointer-arith\
			 -Wcast-qual -Wcast-align\
			 -Wwrite-strings \
			 -fshort-enums -fno-common -Dinline= -g
CFLAGS=-ipo -O3 -xHost
CPPFLAGS=-openmp -mkl
#LDFLAGS=-L/home/encad/berthoumieux/OpenCL/ocl_4.1/lib -lOpenCL -loclUtil_x86_64 -lshrutil_x86_64 -lgmp
INCLUDE=-Iinclude/
LDFLAGS=-static-intel -mkl -openmp
EXE_SZ=sz/sz.bin
EXE_ERROR=error/error.bin
EXE_VOTE=vote/vote.bin

#
## Compilation options :
#
all: compile_production clean_production

debug: compile_debug clean_debug

vote: compile_vote clean_vote

sz: compile_sz clean_sz

error: compile_error clean_error

benchmark: compile_bench clean_bench

testing: compile_production clean_production exec_testing

#
## Compilation :
#
main.prod: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -c $^ -o $@

main.debug: main.c
	$(CC) $(CFLAGS_DEBUG) $(CPPFLAGS) $(INCLUDE) -D DEBUG -c $^ -o $@ 

main.bench: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -D BENCHMARK -c $^ -o $@ 

main.error: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -D BENCHMARK -D RATIO_ -c $^ -o $@ 

main.sz: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -D BENCHMARK -D SZ_ -c $^ -o $@ 

main.vote: main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -D BENCHMARK -D MAJ_VOTE_ -c $^ -o $@ 

#
## Linking :
#
compile_production: main.prod
	$(CC) $^ -o $(EXECUTABLE)  $(LDFLAGS)

compile_debug: main.debug
	$(CC) $^ -o $(EXECUTABLE)  $(LDFLAGS)

compile_bench: main.bench
	$(CC) $^ -o $(EXECUTABLE)  $(LDFLAGS)

compile_sz: main.sz
	$(CC) $^ -o $(EXE_SZ)  $(LDFLAGS)

compile_error: main.error
	$(CC) $^ -o $(EXE_ERROR)  $(LDFLAGS)

compile_vote: main.vote
	$(CC) $^ -o $(EXE_VOTE)  $(LDFLAGS)

#
## Cleaning :
#
clean_production: 
	rm main.prod

clean_debug: 
	rm main.debug

clean_bench: 
	rm main.bench

clean_sz: 
	rm main.sz

clean_error: 
	rm main.error

clean_vote: 
	rm main.vote

#
## Execution :
#
exec:
	./$(EXECUTABLE) $(ARG)
exec_testing:
	./$(EXECUTABLE) 1 1 1 1000 --ans=u --graphs=y

#
## Get preprocessor output (make pp) :
#
pp: compile.pp echo

compile.pp: 
	cp main.c main_pp.c;
	$(CC) $(CFLAGS) $(CPPFLAGS) -D MAJ_VOTE_  -E main_pp.c > preprocessor.c 

echo: 
	$(EDITOR) preprocessor.c
	rm preprocessor.c
	rm main_pp.c

#
## Delete data files
#
mrproper: clean clean-data
clean:
	@rm *.dat
	@rm grover
	@rm graphs.eps
	@rm gnuplot.gp
	@rm *.eps

clean-data:
	rm data/*
