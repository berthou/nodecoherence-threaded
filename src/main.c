#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include "mkl_scalapack.h"

#define READ  0
#define WRITE 1
#define GRAPH 10
#define OSCILLATIONS 100
#define SLOW_OSCILLATION 101
#define FAST_OSCILLATION 102
#define NUMBER_OF_LOCAL_MAXIMUM 400
#define NOTPOINTS 1000

#include "lib_grover_simulation.h"


/* 
 *
 *		!!! check return values
 *
 */

int main(int argc, char **argv)
{
	/* number of bosons : M
	 * number of sites  : Np
	 * computing from Np_start to Np_stop (arguments from the command line)
	 * binomial_iterator and combfactor_iterator are respectively used in the BINOMIAL_COEFF and COMPUTE_COMBFACTOR macros as iterators
	 * combfactor_bin, combfactor_product are variables of the COMPUTE_COMBFACTOR macro
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		matrix_size,
		i,j,n,k,knx,
		info,lwork,
		number_of_extrema,
		notpoints = NOTPOINTS,
		binomial_iterator,
		combfactor_iterator;

	int	   *statens      = NULL,
		   *statensp     = NULL,
		   *answer	     = NULL,
		   *graphs       = NULL;

	double **points      = NULL,
		   *combfactor   = NULL,
		   **szmatelem   = NULL,
		   *tempvector   = NULL,
		   *work         = NULL,
		   *matrix       = NULL,
		   *energies     = NULL,
		   *innerprod    = NULL,
		   *oscillations = NULL,
		   tmax=atof(argv[4]),
		   diag_var,
		   wkopt,
		   matelem,
		   sum,
		   combfactor_bin,
		   combfactor_product;

	struct timeval start,
				   stop;

	/* Process command line arguments */
	graphs = (int*) malloc((M+1)*sizeof(int));
	answer = (int*) malloc(M*sizeof(int));
	if ((process_arguments(argc,argv,answer,graphs,M)) < 0)
	{
		free(graphs);
		free(answer);
		return -1;
	}

	if (graphs[M] == 0)
	{
		printf("No graphs queried, exiting.\n");
		free(graphs);
		free(answer);
		return -1;
	}

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
	if(statens == NULL || statensp == NULL) {
		printf("malloc error for statens, exiting...\n");
		exit(-1);
	}

	for(Np=Np_start;Np<=Np_stop;Np++) 
	{
		/* Get matrix size */
		matrix_size=pow(Np+1,M);
		/* Allocate matrix */
		matrix=(double*)realloc(matrix,(matrix_size*matrix_size)*sizeof(double));
		energies=(double*)realloc(energies,matrix_size*sizeof(double));
		if(matrix == NULL || energies == NULL) {
			printf("malloc error matrix & energies, exiting...\n");
			exit(-1);
		}
		/* Get the inner_product pointer */
		innerprod=get_innerproduct_pointer(Np);
		if (innerprod==NULL) {
			printf("No pre-computed innerprod\n");
			exit(-1);
		}
		/* Allocate combfactor, szmatelemarrays */
		combfactor = (double*)realloc(combfactor,matrix_size*sizeof(double));
		szmatelem  = alloc_szmatelem(szmatelem,matrix_size,M);
		tempvector = (double*)malloc(matrix_size*sizeof(double));
		if(combfactor == NULL || szmatelem == NULL || tempvector == NULL) {
			printf("malloc error for combfactor or szmatelem or tempvector, exiting...\n");
			exit(-1);
		}
		if(eigensystem_to_file(matrix,energies,tempvector,szmatelem,Np,M,READ) < 0)
		{
			notpoints = 1;

#ifdef BENCHMARK
			gettimeofday(&start,NULL);
#endif
#ifdef DEBUG
			printf("Np = %d\n",Np);
#endif
			/* Iteration over rows of the matrix */
			for(i=0;i<matrix_size;i++)
			{
				compute_state_list(i,Np+1,M,statens);
				COMPUTE_COMBFACTOR(Np,statens,M,combfactor[i]);
				for(k=0;k<M;k++)
					szmatelem[k][i]=(2.0*statens[k]-(double)Np)/(double)Np;

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
							if(answer[n] == 1)
								diag_var*=(statens[n]/(double)Np);
							else
								diag_var*=(1-(statens[n]/(double)Np));
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
				exit(-1);
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
			if(tempvector == NULL) {
				printf("malloc error tempvector, exiting...\n");
				return -1;
			}
			eigensystem_to_file(matrix,energies,tempvector,szmatelem,Np,M,WRITE);
			printf("Data saved, Np=%d completed at : ",Np);
		fflush(stdout);
		system("date");
			continue;
		}

#ifdef DEBUG
		print_matrix_( "Eigenvalues", 1, matrix_size, energies, 1 );
		print_matrix_( "Eigenvectors (stored columnwise)", matrix_size, matrix_size, matrix, matrix_size);
		printf("\ncombfactor = \n\n");
		print_vector_double(combfactor,matrix_size);
		printf("\ntempvector = \n\n");
		print_vector_double(tempvector,matrix_size);
		printf("szmatelem = \n\n");
		print2D(szmatelem,M,matrix_size);
#endif
		points = compute_overlap(szmatelem,tempvector,matrix,energies,matrix_size,Np,notpoints,tmax,M,graphs);
		free(tempvector);

		write_to_file(points,oscillations,Np,M,notpoints,graphs,GRAPH,number_of_extrema);
		/* Slow oscillations */
		if(M > 4)
		{
		oscillations = compute_oscillations(points,notpoints,&number_of_extrema,Np,SLOW_OSCILLATION);
		write_to_file(points,oscillations,Np,M,notpoints,graphs,SLOW_OSCILLATION,number_of_extrema);
		free(oscillations);
		}
		/* Fast oscillations */
		oscillations = compute_oscillations(points,notpoints,&number_of_extrema,Np,FAST_OSCILLATION);
		write_to_file(points,oscillations,Np,M,notpoints,graphs,FAST_OSCILLATION,number_of_extrema);
		free(oscillations); 

#ifdef BENCHMARK
		printf("Np=%d completed at : ",Np);
		fflush(stdout);
		system("date");
#endif
		free(points);
	}
	free(matrix);
	free(energies);
	free(work);
	free(combfactor);
	free2D(szmatelem,M);
	free(statens);
	free(statensp);
	free(graphs);
	free(answer);
	printf("Done.\n");
	return 0;
}

