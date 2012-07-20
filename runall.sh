#! /bin/bash

icc -O3 -mkl -o test main.c -static-intel -openmp

if [ -f graphs.eps ]
then
	rm graphs.eps
fi

if [ -f gnuplot.gp ]
then
	rm gnuplot.gp
fi

max=20;
max3=10;
max4=10;
max5=5;

echo "set term postscript enhanced color" >> gnuplot.gp
echo "set output 'graphs.eps'" >> gnuplot.gp
echo "set yrange [0:1]" >> gnuplot.gp

for (( m=1; m<=5; m++ ))
do
	if (( $m == 3 ));
	then
		max=$max3;
	elif (( $m == 4 ));
	then
		max=$max4;
	elif (( $m == 5 ));
	then
		max=$max5;
	fi
	echo "m= $m, max = $max";
	./test 1 $max $m
done

max=20;

for (( m=1; m<=5; m++ ))
do
	if (( $m == 3 || $m == 4));
	then
		max=$max4;
	elif (( $m == 5 ));
	then
		max=$max5;
	fi
	for (( n=1; n<=$max; n++ ))
	do
		echo "plot 'M"$m"N"$n".dat' using 1:2 with points, 'M"$m"N"$n".dat' using 1:2 smooth unique title 'Np= $n, M= $m' with lines" >> gnuplot.gp
	done
done

gnuplot gnuplot.gp
rm gnuplot.gp
