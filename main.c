#include <stdio.h>
#include <math.h>
#include <omp.h>
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
		info,lwork,
		notpoints=300,
		binomial_iterator,
		combfactor_iterator;

	int	   *statens    = NULL,
		   *statensp   = NULL;

	double *points     = NULL,
		   *combfactor = NULL,
		   *szmatelem  = NULL,
		   *tempvector = NULL,
		   *work       = NULL,
		   *matrix     = NULL,
		   *energies   = NULL,
		   *innerprod  = NULL,
		   diag_var,
		   wkopt,
		   matelem,
		   sum,
		   tmax=150,
		   combfactor_bin,
		   combfactor_product;

	struct timeval start,
				   stop;

	/* openMP */
	int nb_procs =omp_get_num_procs();
	omp_set_num_threads(nb_procs);

#ifdef BENCHMARK
	printf("%d thread(s) available\n",nb_procs);
	printf("Starting at : ");
	fflush(stdout);
	system("date");
#endif

	/* Allocate statens and statensp arrays */
	statens  = (int*)malloc(M*sizeof(int));
	statensp = (int*)malloc(M*sizeof(int));

	for(Np=Np_start;Np<=Np_stop;Np++) {
#ifdef BENCHMARK
		gettimeofday(&start,NULL);
#endif
#ifdef DEBUG
		printf("Np = %d\n",Np);
#endif
		/* Get the inner_product pointer */
		innerprod=get_innerproduct_pointer(Np);
		if (innerprod==NULL) {
			printf("No pre-computed innerprod\n");
		}
		/* Get matrix size */
		matrix_size=pow(Np+1,M);
		/* Allocate matrix */
		matrix=(double*)realloc(matrix,(matrix_size*matrix_size)*sizeof(double));
		energies=(double*)realloc(energies,matrix_size*sizeof(double));
		/*memset(matrix,0,matrix_size*matrix_size*sizeof(double));
		  memset(energies,0,matrix_size*sizeof(double));
		  */
		if(matrix == NULL || energies == NULL) {
			printf("malloc error, exiting...\n");
		}
		/* Allocate combfactor, szmatelemarrays */
		combfactor = (double*)realloc(combfactor,matrix_size*sizeof(double));
		szmatelem  = (double*)realloc(szmatelem,matrix_size*sizeof(double));
		/* Iteration over rows of the matrix */
//#pragma omp parallel for shared(matrix_size,statens,statensp,Np,M,combfactor,szmatelem,innerprod,matrix) private(i,j,diag_var,n,matelem,sum,knx)
		for(i=0;i<matrix_size;i++)
		{
			compute_state_list(i,Np+1,M,statens);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[i]);
			szmatelem[i]=(2.0*statens[0]-(double)Np)/(double)Np;

			/* Iteration over columns of the matrix */
			for(j=i;j<matrix_size;j++)
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
					matrix[j*matrix_size+i] = pow(Np,2)*(matelem+diag_var); 
				}
				else
				{
					matrix[j*matrix_size+i] = pow(Np,2)*matelem; 
				}
			}
		}
#ifdef BENCHMARK
		gettimeofday(&stop,NULL);
		printf("Matrix computed in %f sec.\n",(stop.tv_sec+stop.tv_usec*1.0e-6)-(start.tv_sec+start.tv_usec*1.0e-6));
#endif
#ifdef DEBUG
		print_matrix_( "Matrix", matrix_size, matrix_size, matrix, matrix_size);
#endif

		/* LAPACK computation */
		/* Query and allocate the optimal workspace */
		lwork = -1;
		dsyev("Vectors", "Upper", &matrix_size,matrix,&matrix_size, energies,&wkopt,&lwork,&info);
		lwork = (int)wkopt;
		work = (double*)realloc(work,lwork*sizeof(double) );
		if (work == NULL) {
			printf("error malloc work\n");
			return -1;
		}
		memset(work,0,lwork*sizeof(double));
		/* Solve eigenproblem */
		gettimeofday(&start,NULL);
		dsyev("Vectors", "Upper", &matrix_size,matrix,&matrix_size, energies,work,&lwork,&info);
#ifdef BENCHMARK
		gettimeofday(&stop,NULL);
		printf("Eigensystem computed in %f sec.\n",(stop.tv_sec+stop.tv_usec*1.0e-6)-(start.tv_sec+start.tv_usec*1.0e-6));
#endif
		/* Check for convergence */
		if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues.\n" );
			exit( 1 );
		}

		/* End lapack */

		tempvector = compute_tempvector(combfactor,matrix,matrix_size);

#ifdef DEBUG
		print_matrix_( "Eigenvalues", 1, matrix_size, energies, 1 );
		print_matrix_( "Eigenvectors (stored columnwise)", matrix_size, matrix_size, matrix, matrix_size);
		printf("\ncombfactor = \n\n");
		print_vector_double(combfactor,matrix_size);
		printf("\ntempvector = \n\n");
		print_vector_double(tempvector,matrix_size);
		printf("szmatelem = \n\n");
		print_vector_double(szmatelem,matrix_size);
#endif
		points = compute_overlap(szmatelem,tempvector,matrix,energies,matrix_size,Np,notpoints,tmax);
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
	free(matrix);
	free(energies);
	free(work);
	free(combfactor);
	free(szmatelem);
	free(statens);
	free(statensp);
	printf("Done.\n");
	return 0;
}

