#! /bin/bash

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
	for (( j=0; j<=$3; j++ ))
	do
		if [ -f M$3N$i-$j.dat ]
		then
			echo "plot 'M$3N$i-$j.dat' using 1:2 smooth unique title 'N=$i, M=$3 Boson $j' with lines" >> gnuplot.gp
		fi

	done
done
gnuplot gnuplot.gp


#for (( i=$1; i<=$2; i++ ))
#do
#rm M$3N$i.dat
#done

evince graphs.eps
