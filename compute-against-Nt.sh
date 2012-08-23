#! /bin/bash

if [ -f sz/sz.bin ]
then
	rm sz/sz.bin
fi

if [ -f vote/vote.bin ]
then
	rm vote/vote.bin
fi

if [ -f error/error.bin ]
then
	rm error/error.bin
fi

if [ -f gnuplot.gp ]
then
	rm gnuplot.gp
fi

if [ -f $9 ]
then
	rm $9
fi

echo "set term postscript enhanced color" >> gnuplot.gp
echo "set output '"${10}"'" >> gnuplot.gp
echo "set yrange [$8:$9]" >> gnuplot.gp
echo "set xrange [0:$4]" >> gnuplot.gp
echo "set xlabel 'Nt'" >> gnuplot.gp

case "$7" in

	sz)
		make sz
		cd sz
		./sz.bin $1 $2 $3 $4 $5 $6 
		echo "set ylabel '<Sz/N>'" >> ../gnuplot.gp
		for (( i=$1; i<=$2; i=i+1 ))
		do
			for (( j=0; j<=$3; j++ ))
			do
				if [ -f M$3N$i-$j.dat ]
				then
					echo "plot 'sz/M$3N$i-$j.dat' using 1:2 smooth unique title 'N=$i, M=$3' with lines;" >> ../gnuplot.gp
				fi
			done
		done
		;;

	error)
		make error
		cd error
		./error.bin $1 $2 $3 $4 $5 $6 
		echo "set ylabel '<Sz>/sqrt(<Sz^2>-<Sz>^2)'" >> ../gnuplot.gp
		echo -n "plot " >> ../gnuplot.gp
		for (( i=$1; i<=$2-1; i=i+2 ))
		do
			echo -n " 'error/M$3N$i-1.dat' using 1:2 smooth unique title 'N=$i, M=$3' with lines," >> ../gnuplot.gp
		done
		i=$2
		echo -n " 'error/M$3N$i-1.dat' using 1:2 smooth unique title 'N=$i, M=$3' with lines;" >> ../gnuplot.gp
		;;

	vote)
		make vote
		cd vote
		./vote.bin $1 $2 $3 $4 $5 $6 
		echo "set ylabel '<Psuc>'" >> ../gnuplot.gp
		echo -n "plot " >> ../gnuplot.gp
		for (( i=$1; i<=$2-1; i=i+3 ))
		do
			echo -n " 'vote/M$3N$i-1.dat' using 1:2 smooth unique title 'N=$i, M=$3' with lines," >> ../gnuplot.gp
		done
		i=$2
		echo -n " 'vote/M$3N$i-1.dat' using 1:2 smooth unique title 'N=$i, M=$3' with lines;" >> ../gnuplot.gp
		;;
esac
cd ..
gnuplot gnuplot.gp
#rm gnuplot.gp
