#! /bin/bash

make mklbenchmark;

for (( i=$1; i<=$2; i++ ))
do
	./test $@; 
done
