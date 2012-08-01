#! /bin/bash

icc -O3 -D BENCHMARK -mkl -o grover main.c -static-intel -openmp

if [ -f M$3-fast_osc.dat ]
then
	rm M$3-fast_osc.dat
fi

./grover $@ 

if [ -f fast-oscillations.eps ]
then
	rm fast-oscillations.eps
fi


if [ -f gnuplot.gp ]
then
	rm gnuplot.gp
fi

echo "set term postscript enhanced color" >> gnuplot.gp
echo "set output 'fast-oscillations.eps'" >> gnuplot.gp
#echo "set yrange [-1:1]" >> gnuplot.gp

if [ -f M$3-fast_osc.dat ]
then
	echo "plot 'M$3-fast_osc.dat' smooth unique title 'mean for N=$1..$2, M=$3' with lines" >> gnuplot.gp
fi
gnuplot gnuplot.gp

rm gnuplot.gp

evince fast-oscillations.eps
