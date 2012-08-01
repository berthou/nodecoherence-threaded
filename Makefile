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
CFLAGS=-ipo -O3 -no-prec-div -xHost
CPPFLAGS=-openmp
#LDFLAGS=-L/home/encad/berthoumieux/OpenCL/ocl_4.1/lib -lOpenCL -loclUtil_x86_64 -lshrutil_x86_64 -lgmp
INCLUDE=-Iinclude/
LDFLAGS=-static-intel -mkl -openmp
EXECUTABLE=grover

#
# Compilation options :
#
all: compile_production clean_production

debug: compile_debug clean_debug

benchmark: compile_bench clean_bench

testing: compile_production clean_production exec_testing

#
# Compilation :
#
main.prod: src/main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -c $^ -o $@

main.debug: src/main.c
	$(CC) $(CFLAGS_DEBUG) $(CPPFLAGS) $(INCLUDE) -D DEBUG -c $^ -o $@ 

main.bench: src/main.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) -D BENCHMARK -c $^ -o $@ 

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
	./$(EXECUTABLE) 1 1 1 1000 --ans=u --graphs=y

#
# Get preprocessor output (make pp) :
#
pp: compile.pp echo

compile.pp: 
	cp src/main.c src/main_pp.c;
	$(CC) $(CFLAGS) $(CPPFLAGS)  -E src/main_pp.c > preprocessor.c 

echo: 
	$(EDITOR) preprocessor.c
	rm preprocessor.c
	rm src/main_pp.c

#
# Delete data files
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
