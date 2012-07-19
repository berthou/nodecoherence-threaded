#include <stdio.h>
#include <math.h>
#include "mkl_scalapack.h"

#include "lib_grover_simulation.h"

/* 
 *
 *		!!! check return values
 *
 */

/* Auxiliary routine: printing a matrix */
void print_matrix_( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %.10f", a[i+j*lda] );
		printf( "\n" );
	}
}

int main(int argc, char **argv)
{
	/* number of bosons : M
	 * number of sites  : Np
	 * computing from Np_start to Np_stop (arguments from the command line)
	 * The GNU Scientific Library is used for operations on matrices
	 * binomial_iterator and combfactor_iterator are respectively used in the BINOMIAL_COEFF and COMPUTE_COMBFACTOR macros as iterators
	 * combfactor_bin, combfactor_product are variables of the COMPUTE_COMBFACTOR macro
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		matrix_size,
		i,j,n,knx,
		notpoints=40,
		binomial_iterator,
		combfactor_iterator;

	int	   *statens    = NULL,
		   *statensp   = NULL;

	double *points     = NULL,
		   *combfactor = NULL,
		   *szmatelem  = NULL,
		   *tempvector = NULL,
		   diag_var,
		   matelem,
		   sum,
		   tmax=30,
		   combfactor_bin,
		   combfactor_product;

	double *innerprod  = NULL;


	int info,
		lwork;
	double wkopt,
		   *work = NULL,
		   *a    = NULL,
		   *w    = NULL;


#ifdef BENCHMARK
	printf("Starting at : ");
	fflush(stdout);
	system("date");
#endif

	/* Allocate statens and statensp arrays */
	statens  = (int*)malloc(M*sizeof(int));
	statensp = (int*)malloc(M*sizeof(int));


	for(Np=Np_start;Np<=Np_stop;Np++) {
#ifndef BENCHMARK
		printf("Np = %d\n",Np);
#endif
		/* Get the inner_product pointer */
		innerprod=get_innerproduct_pointer(Np);
		if (innerprod==NULL)
			return -1;
		/* Get matrix size */
		matrix_size=pow(Np+1,M);
		/* Allocate matrix */
		a=(double*)realloc(a,(matrix_size*matrix_size)*sizeof(double));
		w=(double*)realloc(w,matrix_size*sizeof(double));
		memset(a,0,matrix_size*matrix_size*sizeof(double));
		memset(w,0,matrix_size*sizeof(double));
		if(a == NULL || w == NULL) {
			printf("malloc error, exiting...\n");
			exit(-1);
		}
		/* Allocate combfactor, szmatelemarrays */
		combfactor = (double*)realloc(combfactor,matrix_size*sizeof(double));
		szmatelem  = (double*)realloc(szmatelem,matrix_size*sizeof(double));
		/* Iteration over rows of the matrix */
		for(i=0;i<matrix_size;i++)
		{
			compute_state_list(i,Np+1,M,statens);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[i]);
			szmatelem[i]=(2.0*statens[0]-(double)Np)/(double)Np;

			/* Iteration over columns of the matrix */
			for(j=0;j<matrix_size;j++)
			{
				compute_state_list(j,Np+1,M,statensp);

				/* If we are on the diagonal of H : */
				if(i==j)
				{
					diag_var=1.0;
					for(n=0;n<M;n++)
					{
						diag_var*=(statens[n]/(double)Np);
					}
				}
				/* For the upper part of the matrix : */
				if(j>=i) 
				{
					matelem=1.0;
					for(n=0;n<M;n++)
					{
						sum=0.0;
						for(knx=0;knx<=Np;knx++)
						{
							sum+=((double)knx/(double)Np)*innerprod[statens[n]*(Np+1)+knx]*innerprod[statensp[n]*(Np+1)+knx];
						}
						matelem*=sum;
					}
					if (i==j)
					{
						a[j*matrix_size+i] = pow(Np,2)*(matelem+diag_var); 
					}
					else
					{
						a[j*matrix_size+i] = pow(Np,2)*matelem; 
						a[i*matrix_size+j] = pow(Np,2)*matelem; 
					}
				}
			}
		}
#ifdef DEBUG
		print_matrix_( "Matrix", matrix_size, matrix_size, a, matrix_size);
#endif

		/* LAPACK computation */
		/* Query and allocate the optimal workspace */
		lwork = -1;
		dsyev("Vectors", "Upper", &matrix_size,a,&matrix_size, w,&wkopt,&lwork,&info);
		lwork = (int)wkopt;
		work = (double*)realloc(work,lwork*sizeof(double) );
		if (work == NULL) {
			printf("error malloc work\n");
				return -1;
		}
		memset(work,0,lwork*sizeof(double));
		/* Solve eigenproblem */
		dsyev("Vectors", "Upper", &matrix_size,a,&matrix_size, w,work,&lwork,&info);
		/* Check for convergence */
		if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues.\n" );
			exit( 1 );
		}

		/* End lapack */

		tempvector = compute_tempvector(combfactor,a,matrix_size);

#ifdef DEBUG
		print_matrix_( "Eigenvalues", 1, matrix_size, w, 1 );
		print_matrix_( "Eigenvectors (stored columnwise)", matrix_size, matrix_size, a, matrix_size);
		printf("\ncombfactor = \n\n");
		print_vector_double(combfactor,matrix_size);
		printf("\ntempvector = \n\n");
		print_vector_double(tempvector,matrix_size);
		printf("szmatelem = \n\n");
		print_vector_double(szmatelem,matrix_size);
#endif
		points = compute_overlap(szmatelem,tempvector,a,w,matrix_size,Np,notpoints,tmax);

		free(tempvector);

		write_to_file(points,Np,M,notpoints);

#ifdef BENCHMARK
		printf("Np=%d completed at : ",Np);
		fflush(stdout);
		system("date");
#endif
		verification(points,M,Np);
		free(points);
	}
	free(a);
	free(w);
	free(work);
	free(combfactor);
	free(szmatelem);
	free(statens);
	free(statensp);
	printf("Done.\n");
	return 0;
}

