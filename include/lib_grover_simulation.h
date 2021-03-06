#ifndef LIB_GROVER_SIMULATION_H
#define LIB_GROVER_SIMULATION_H

#include "reference_data.h"
#include "inner20.h"

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>


void free2D(double **array,int row);
double **malloc2D(int row, int column);
double  *get_innerproduct_pointer(int Np);
int verification(double *points,int M,int Np);
void write_to_file(double **data,int Np,int M,int size, int *graphs);
double overlap(double **szmatelem,double *tempvector,double *inversevectors, double *energies, double t,unsigned long matrix_size,int M);
double **compute_overlap(double **szmatelem,double *tempvector,double *inversevectors, double *energies,unsigned long matrix_size,int Np,int notpoints,double tmax,int M, int *graphs, int *ANSlist);
double *compute_tempvector(double *combfactor, double *inversevector,unsigned long matrix_size);
void print_vector_double(double *v,int size);
void print_vector(double *v,int size);
void print_matrix(double *H,int size);
double binomial_coeff(int n, int k);
void compute_state_list(int number,int base, int M, int *data);
int compute_matrix_size(int n,int p);
double inner_product(double k,double kx,double np);

#include "lib_grover_simulation.cc"

#endif
