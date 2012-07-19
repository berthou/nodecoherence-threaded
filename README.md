Simulating Grover's algorithm with Bose-Einstein condensates with no decoherence.
=============

Download
--------

Run 

	git clone https://github.com/berthou/nodecoherence.git

or download & extract 

	https://github.com/berthou/nodecoherence/zipball/master

then go to the main folder :

	cd nodecoherence/


Before compiling and running the programm you need to install dependencies.

Dependencies
------------

The following dependencies are required :

* [GNU Scientific Library](http://www.gnu.org/software/gsl/) -- run `make dependencies` to install the library (unix only) or install it by yourself
* [gnuplot](www.gnuplot.info/) -- (Optional)


Compilation & execution
-----------------------

To make sure that everything went fine run `make -s testing`. The output should be :

	Np = 1
	Agrees on at least 13 digits
	Done.

If you have gnuplot installed you can view the output graph by runnning 

	./plots.sh 1 1 1

You might change the arguments but the 2nd must be greater or equal than the 1st one and the last one should'nt be greater than 4 (too much computation needed).
