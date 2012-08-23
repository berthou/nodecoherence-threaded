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
	/* computing from Np_start to Np_stop (arguments from the command line)
	 * binomial_iterator and combfactor_iterator are respectively used in the BINOMIAL_COEFF and COMPUTE_COMBFACTOR macros as iterators
	 * combfactor_bin, combfactor_product are variables of the COMPUTE_COMBFACTOR macro
	 */

	int M=atoi(argv[3]),
		Np,
		Np_start=atoi(argv[1]),
		Np_stop=atoi(argv[2]),
		i,j,n,k,knx,
		info,lwork,anslist_prod,ans_it=0,
		notpoints=200,
		binomial_iterator,
		combfactor_iterator;
	unsigned long	matrix_size;

	int	   *statens    = NULL,
		   *statensp   = NULL,
		   *answer     = NULL,
		   *ANSlist    = NULL,
		   *graphs     = NULL;

	double **points     = NULL,
		   *combfactor = NULL,
		   **szmatelem  = NULL,
		   *tempvector = NULL,
		   *work       = NULL,
		   *matrix     = NULL,
		   *energies   = NULL,
		   *innerprod  = NULL,
		   diag_var,
		   wkopt,
		   matelem,
		   sum,
		   tmax=atof(argv[4]),
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
			exit(-1);
		}
		ans_it=0;
		/* Get matrix size */
		matrix_size=pow(Np+1,M);
		/* Allocate matrix */
		matrix=(double*)realloc(matrix,(matrix_size*matrix_size)*sizeof(double));
		energies=(double*)realloc(energies,matrix_size*sizeof(double));
		ANSlist=(int*)calloc(matrix_size,sizeof(int));
		if(matrix == NULL || energies == NULL || ANSlist == NULL) {
			printf("malloc error matrix & energies, exiting...\n");
			exit(-1);
		}
		/* Allocate combfactor, szmatelemarrays */
		combfactor = (double*)realloc(combfactor,matrix_size*sizeof(double));
		szmatelem  = alloc_szmatelem(szmatelem,matrix_size,M);
		if(combfactor == NULL || szmatelem == NULL) {
			printf("malloc error for combfactor or szmatelem, exiting...\n");
			exit(-1);
		}
		/* Iteration over rows of the matrix */
		for(i=0;i<matrix_size;i++)
		{
			compute_state_list(i,Np+1,M,statens);
			COMPUTE_COMBFACTOR(Np,statens,M,combfactor[i]);
			for(k=0;k<M;k++)
				szmatelem[k][i]=(2.0*statens[k]-(double)Np)/(double)Np;

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
						if(answer[n] == 1)
							diag_var*=(statens[n]/(double)Np);
						else
							diag_var*=(1-(statens[n]/(double)Np));
					}
					/* process ANSlist */
					anslist_prod=1;
					for(n=0;n<M;n++)
					{
						anslist_prod *= ((double)statens[n]/(double)Np > 0.5) ? 1 : 0;
					}
					if (anslist_prod == 1)
						ANSlist[ans_it++]=i;
					else
						ANSlist[ans_it++]=-1;	
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
		printf("\nMatrix computed in %f sec.\n",(stop.tv_sec+stop.tv_usec*1.0e-6)-(start.tv_sec+start.tv_usec*1.0e-6));
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
		points = compute_overlap(szmatelem,tempvector,matrix,energies,matrix_size,Np,notpoints,tmax,M,graphs,ANSlist);
		free(tempvector);
		write_to_file(points,Np,M,notpoints,graphs);
		free(ANSlist);

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

