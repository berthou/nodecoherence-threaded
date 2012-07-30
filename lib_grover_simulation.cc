#ifndef LIB_GROVER_SIMULATION_CC
#define LIB_GROVER_SIMULATION_CC

#include <sys/time.h>

double *compute_fast_oscillations(double **points,int size,int *number_of_extrema) 
{
	int count = (NUMBER_OF_LOCAL_MAXIMUM > size+1) ? size+1 : NUMBER_OF_LOCAL_MAXIMUM ;
	double *result = (double*)malloc((count+1)*sizeof(double));
	double current_point=0,
		   last_local_extrema=0;
	int result_index=0,
		first_found=1,
		i=size+1+1;

	while (count > 0 && (i < 2*(size+1)))
	{
		while ((points[0][i] > current_point) && (i < 2*(size+1)))
		{
			current_point=points[0][i];
			i++;
		}
		/* Do not consider the last point as a local extrema */
		if(i-1 == 2*(size+1)-1)
			break;

		/* points[0][i-1] is a local maximum */
		//printf("Local maximum found : %f\n",points[0][i-1]);
		if(first_found)
		{
			first_found = 0;
			last_local_extrema = points[0][i-(size+1)-1];
		}
		else
		{
			result[result_index]=points[0][i-1-(size+1)]-last_local_extrema;
#ifndef BENCHMARK
			printf("Local maximum found at %f\n",result[result_index]);
#endif
			last_local_extrema = points[0][i-(size+1)-1];
			result_index++;
		}

		/* Go to the next point : */
		i++;

		while (points[0][i] < current_point && (i < 2*(size+1)))
		{
			current_point=points[0][i];
			i++;
		}
		/* Do not consider the last point as a local extrema */
		if(i-1 == 2*(size+1)-1)
			break;
		/* points[0][i-1] is a local minimum */ 
		//printf("Local minimum found : %f\n",points[0][i-1]);
		/* Go to the next point : */
		i++;
		/* Iteration is over */
		count--;
	};

	/* Compute the mean */
	result[result_index]=0;
	for(i=0;i<result_index;i++)
		result[result_index]+=result[i];

	result[result_index]/=result_index;

	*number_of_extrema = result_index;

#ifdef BENCHMARK
	printf("Mean = %f\n",result[result_index]);
#endif

	return result;
}

int process_arguments(int argc, char **argv,int *answer,int *graphs,int M)
{
	int i;
	graphs[M]=0;
	if(argc != 7)
	{
		printf("Usage : %s N_start N_stop M Tmax --ans=<u/d> --graphs=<y/n>\n",argv[0]);
		return -1;
	}
	if(strlen(argv[5]) != 6+M)
	{
		printf("Answer length is not correct\nanswer is a string of %dx u or d\n",M);
		return -2;
	}
	if(strlen(argv[6]) != 9+M)
	{
		printf("The list of graphs is not correct\nGraphs a string of %dx y or n\n",M);
		return -3;
	}
	for (i = 0; i < M; i++)
	{
		if(argv[5][6+i] == 'u')
			answer[i] = 1;
		else if (argv[5][6+i] == 'd')
			answer[i] = -1;
		else
		{
			printf("Answer is not correct, only 'u' or 'd'\n");
			return -4;
		}
		if(argv[6][9+i] == 'y') {
			graphs[M]++;
			graphs[i] = 1;
		}
		else if (argv[6][9+i] == 'n')
			graphs[i] = 0;
		else
		{
			printf("Graph is not correct, only 'y' or 'n'\n");
			return -5;
		}
	}
	return 0;
}



double **alloc_szmatelem(double **szmatelem, int matrix_size, int M)
{
	int i,j;
	if (szmatelem == NULL)
		szmatelem = (double**) malloc(M*sizeof(double*));
	for (i = 0; i < M; i++)
	{
		szmatelem[i] = (double*) realloc(szmatelem[i],matrix_size*sizeof(double));
	}
	return szmatelem;
}

double *get_innerproduct_pointer(int Np)
{
	double *pointer;

	if(Np==1)
		pointer=innerprodN1;
	else if(Np==2)
		pointer=innerprodN2;
	else if(Np==3)
		pointer=innerprodN3;
	else if(Np==4)
		pointer=innerprodN4;
	else if(Np==5)
		pointer=innerprodN5;
	else if(Np==6)
		pointer=innerprodN6;
	else if(Np==7)
		pointer=innerprodN7;
	else if(Np==8)
		pointer=innerprodN8;
	else if(Np==9)
		pointer=innerprodN9;
	else if(Np==10)
		pointer=innerprodN10;
	else if(Np==11)
		pointer=innerprodN11;
	else if(Np==12)
		pointer=innerprodN12;
	else if(Np==13)
		pointer=innerprodN13;
	else if(Np==14)
		pointer=innerprodN14;
	else if(Np==15)
		pointer=innerprodN15;
	else if(Np==16)
		pointer=innerprodN16;
	else if(Np==17)
		pointer=innerprodN17;
	else if(Np==18)
		pointer=innerprodN18;
	else if(Np==19)
		pointer=innerprodN19;
	else if(Np==20)
		pointer=innerprodN20;
	else
	{
		return NULL;
	}

	return pointer;
}

int verification(double *points,int M,int Np)
{
	int i,
		count_five=0,
		count_thirteen=0,
		count_ten=0;

	double diff;

	point *data;

	if (M==1)
		data = M1N_table;
	else if (M==2 && Np==1)
		data = M2N1_table;
	else if (M==2 && Np==10)
		data = M2N10_table;
	else if (M==3 && Np==1)
		data = M3N1_table;
	else if (M==3 && Np==5)
		data = M3N5_table;
	else if (M==4 && Np==1)
		data = M4N1_table;
	else if (M==4 && Np==2)
		data = M4N2_table;
	else if (M==4 && Np==3)
		data = M4N3_table;
	else {
#ifndef BENCHMARK
		printf("No data to compare.\n");
#endif
		return -1;
	}

	for(i=0;i<41;i++)
	{
		diff = fabs(points[i+41] - data[i].y);

		if( diff > 0.00001) {
			count_five++;
		}
		if( diff > 0.0000000001) {
			count_ten++;
		}
		if( diff > 0.0000000000001) {
			/*if (M>2)
			  printf("Differs between 10th and 13th digit : [%d]%.20f\n",i,diff);*/
			count_thirteen++;
		}
	}
	if (count_thirteen)
	{
		if (count_ten) 
		{
			/*printf("Errors on %d points on 10 digits\n",count_ten);*/
			if (count_five)
			{
				/*printf("Errors on %d points on 5 digits\n",count_five);*/
				/*return count_ten;*/
			}
			else 
			{
				printf("Agrees on at least 5 digits\n");
			}
		}
		else 
		{
			printf("Agrees on at least 10 digits\n");
		}
	}
	else 
	{
		printf("Agrees on at least 13 digits\n");
	}
	fflush(stdout);
	return count_thirteen;
}

void write_to_file(double **data,double *oscillations,int Np,int M,int size, int *graphs, int type,int number_of_extrema)
{
	int i,
		graph_id,
		index=0;

	int *file_desc;

	char x[20],
		 y[20];

	char filename[15],
		 m[2],
		 n[5],
		 char_index[2];

	file_desc = (int*)malloc(graphs[M]*sizeof(int));

	for(graph_id = 0 ; graph_id < M ; graph_id++) 
	{
		if(graphs[graph_id])
		{
			/* Build filename = "M(value of M)N(value of N)-graph_id.dat"
			*/
			memset(filename,0,15);
			memset(m,0,2);
			memset(char_index,0,2);
			memset(n,0,5);
			strcat(filename,"M");
			snprintf(m,sizeof(m),"%d",M);
			strcat(filename,m);
			if (type == GRAPH)
			{
				strcat(filename,"N");
				snprintf(n,sizeof(n),"%d",Np);
				strcat(filename,n);
				strcat(filename,"-");
				snprintf(char_index,sizeof(char_index),"%d",graph_id+1);
				strcat(filename,char_index);
				strcat(filename,".dat");
			}
			else
			{
				strcat(filename,"-osc");
				strcat(filename,".dat");
			}
			/* Open file */
			if(type == OSCILLATIONS)
			{
				if ((file_desc[graph_id]=open(filename,O_WRONLY | O_CREAT | O_APPEND, 00600))== -1)
				{
					printf("Unable to open the file %s\n",filename);
					return;
				}
			}
			else 
			{
				if((file_desc[graph_id]=open(filename,O_WRONLY | O_CREAT | O_TRUNC, 00600))== -1)
				{
					printf("Unable to open the file %s\n",filename);
					return;
				}
			}

			/* Write file */
			if(type == GRAPH)
			{
				for(i=0;i<size+1;i++)
				{
					snprintf(x,sizeof(x)+1,"%.15f\t",data[index][i]);
					write(file_desc[graph_id],x,strlen(x));
					snprintf(y,sizeof(y)+1,"%.15f\n",data[index][i+size+1]);
					write(file_desc[graph_id],y,strlen(y));
				}
			}
			if(type == OSCILLATIONS)
			{
				//for(i=0;i<number_of_extrema;i++)
				//{
				snprintf(x,sizeof(x)+1,"%.15f\n",oscillations[number_of_extrema]);
				write(file_desc[graph_id],x,strlen(x));
				//}

			}
			close(file_desc[graph_id]);
		}
		index++;
	}
}

	inline double 
overlap(double **szmatelem,double *tempvector,double *inversevectors, double *energies, double t,int matrix_size,int M)
{
	int nstilde,
		ms;
	double var   = 0.0,
		   sine  = 0.0,
		   cosine= 0.0,
		   temp;

	for(nstilde=0;nstilde<matrix_size;nstilde++)
	{
		sine=0.0;
		cosine=0.0;
		for(ms=0;ms<matrix_size;ms++)
		{
			temp   =  tempvector[ms]*inversevectors[ms*matrix_size+nstilde];
			cosine += temp*cos(t*energies[ms]);
			sine   += temp*sin(t*energies[ms]);
		}
		sine*=sine;
		cosine*=cosine;
		var+=szmatelem[M][nstilde]*(sine+cosine);
	}
	return var;

}

double **compute_overlap(double **szmatelem,double *tempvector,double *inversevectors, double *energies,int matrix_size,int Np,int notpoints,double tmax,int M, int *graphs)
{
	struct timeval start,
				   stop;

	gettimeofday(&start,NULL);

	int i,j,index=0;
	double **result = malloc2D(graphs[M],2*(notpoints+1));
	for(j=0;j<M;j++)
	{
		if(graphs[j]) 
		{
#pragma omp parallel for shared(szmatelem,tempvector,inversevectors,energies,matrix_size,Np,notpoints,tmax,result,index) private(i)
			for(i=0;i<notpoints+1;i++)
			{
				result[index][i]=i*(tmax/notpoints);
				result[index][i+notpoints+1]=overlap(szmatelem,tempvector,inversevectors,energies,((double)i*(tmax/notpoints))/Np,matrix_size,j);
				/* printf("t= %f\tSz/Nt= %.10f\n",result[i],result[i+notpoints+1]); */
			}
		}
		index++;
	}
	gettimeofday(&stop,NULL);

#ifdef BENCHMARK
	printf("overlap computed in %f sec.\n",(stop.tv_sec+stop.tv_usec*1.0e-6)-(start.tv_sec+start.tv_usec*1.0e-6));
#endif
	return result;
}



double *compute_tempvector(double *combfactor, double *inversevector,int matrix_size)
{
	double *result = (double*)malloc(matrix_size*sizeof(double));
	int i,j;
	for(i=0;i<matrix_size;i++)
	{
		result[i]=0.0;
		for(j=0;j<matrix_size;j++)
		{
			result[i]+=combfactor[j]*inversevector[i*matrix_size+j];
		}
	}
	return result;
}

double **malloc2D(int row, int column)
{
	int i;
	double **result = malloc(row*sizeof(double*));
	for (i = 0; i < row; i++)
	{
		result[i] = (double*)malloc(column*sizeof(double));
	}
	return result;
}

void free2D(double **array,int row)
{
	int i;
	for(i=0;i<row;i++)
		free(array[i]);
	free(array);
}
void print2D(double **array,int row, int column)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<column;j++)
			printf("%f\t",array[i][j]);
		printf("\n");
	}
	printf("\n");
}

void print_vector_double(double *v,int size)
{
	int i;
	for(i=0;i<size;i++)
	{
		printf("%.10f\t",v[i]);
	}
	printf("\n");
}

#define COMPUTE_COMBFACTOR(Np,statens,M,var)											\
	do {																				\
		var = 1.0/(sqrt(pow(2,Np*M)));													\
		combfactor_product=1.0;													\
		for(combfactor_iterator=0;combfactor_iterator<M;combfactor_iterator++) \
		{ 																		\
			BINOMIAL_COEFF(Np,statens[combfactor_iterator],combfactor_bin);			\
			combfactor_product*=combfactor_bin;										\
		}																		\
		combfactor_product=sqrt(combfactor_product);								\
		var*=combfactor_product;												\
	}while(0)



#define BINOMIAL_COEFF(n,k,var)				\
	do {										\
		var = 1.0;								\
		if(n-2*k>0)								\
		{										\
			for(binomial_iterator=n;binomial_iterator>=n-k+1;binomial_iterator--)				\
			var=var*binomial_iterator/(n-binomial_iterator+1);		\
		}										\
		else									\
		{										\
			for(binomial_iterator=n;binomial_iterator>=k+1;binomial_iterator--)						\
			var=var*binomial_iterator/(n-binomial_iterator+1);			\
		}											\
	}while(0)

void compute_state_list(int number,int base, int M, int *data)
{
	/* Convert number 'number' in base 'base'
	 * the result is stored in an array of size 'array_length'
	 * with LSB representation
	 */

	int iterator = M;

	/* Start to fill from the end of the array : */
	int index=M-1;

	while(iterator != 0)
	{
		data[index] = (number % base);
		number= (number / base);
		index--;

		iterator--;
	}
}

int compute_matrix_size(int n,int p)
{
	/*
	 * Computes n^p and return the result as a int
	 */
	int i;
	int result=n;
	for(i=1;i<p;i++) {
		result*=n;
	}
	return result;
}

#define FACTORIAL(x) 		\
	do {					\
		if(x<=1)			\
		x=1.0;			\
		else				\
		for(fact_iterator=x-1;fact_iterator>1;fact_iterator--){x=(x*fact_iterator);} \
	}while(0)


double inner_product(double k,double kx,double np)
{
	/* fact_iterator is used in the FACTORIAL macro */
	int i,
		start,
		stop,
		fact_iterator;

	double result=0,
		   k_fact    = k,
		   npk_fact  = np-k,
		   kx_fact   = kx,
		   npkx_fact = np-kx,
		   kn_fact,
		   npkxn_fact,
		   n_fact,
		   nkxk_fact;

	FACTORIAL(k_fact);
	FACTORIAL(npk_fact);
	FACTORIAL(kx_fact);
	FACTORIAL(npkx_fact);

	/* start = max(k-kx)
	 * stop = min(k,Np-kx)
	 */
	start = (k-kx > 0) ? k-kx : 0; 
	stop = (k <= np-kx) ? k : np-kx;

	/* Sum */
	for(i=start;i<=stop;i++){ 

		kn_fact=k-i;
		npkxn_fact=np-kx-i;
		n_fact=i;
		nkxk_fact=i+kx-k;

		FACTORIAL(kn_fact);
		FACTORIAL(npkxn_fact);
		FACTORIAL(n_fact);
		FACTORIAL(nkxk_fact);  

		result+=(pow(-1,i)/(kn_fact*npkxn_fact*n_fact*nkxk_fact));
	}
	/* Sum * sqrt */
	result*=sqrt((k_fact*npk_fact*kx_fact*npkx_fact)/pow(2,np));

	return result;
}

#endif
