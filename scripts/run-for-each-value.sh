#! /bin/bash

make -Bs compile_production;
make -Bs clean_production;

if [ -f graphs.eps ]
then
	rm graphs.eps
fi

if [ -f gnuplot.gp ]
then
	rm gnuplot.gp
fi


echo "set term postscript enhanced color" >> gnuplot.gp
echo "set output 'graphs.eps'" >> gnuplot.gp
echo "set yrange [0:1]" >> gnuplot.gp

for (( i=$1; i<=$2; i++ ))
do
	./test $i $i $3 
	echo "plot 'M$3N$i.dat' using 1:2 with points, 'M$3N$i.dat' using 1:2 smooth unique title 'Np= $i, M= $3' with lines" >> gnuplot.gp
done

gnuplot gnuplot.gp
rm gnuplot.gp

#for (( i=$1; i<=$2; i++ ))
#do
	#rm M$3N$i.dat
#done

evince graphs.eps
