#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define _MIN(a, b) (a<=b ? a : b)

// peut mettre des lignes et colonnes soit differentes soit meme
#define LIGNES					3
#define COLONNES				4

// array for the B column in every matrix
// lignes has to be less than the number of elements in b_array
// add as many elements as needed to b_array
double b_array[]=
{
	8.0f, 3.0f, 8.0f, 5.0f, 4.0f, 8.0f, 7.0f, 10.0f, 0.0f, 1.0f, 3.0f, 11.0f
};

double **allocation_dynamic(int lig, int col);
void affiche(double **matrix, int lig, int col);
void get_b(double **matrix, int lig);
double *gauss_elimination(double **matrix, int lig, int col);
void free_matrix(double ***matrix, int lig);
void free_matrix_solution(double **matrix);
void affiche_solution(double *solution, int lig);

// TEST MATRICES
void matrice_bord(double **matrix, int lig);
void matrice_ding_dong(double **matrix, int lig);
void matrice_franc(double **matrix, int lig);
void matrice_hilbert(double **matrix, int lig);
void matrice_kms(double **matrix, int lig, double p);
void matrice_lehmer(double **matrix, int lig);
void matrice_lotkin(double **matrix, int lig);
void matrice_moler(double **matrix, int lig);
double calculate_time_in_seconds(clock_t start, clock_t end);
void print_time_elapsed(const char *matrix_name, double elapsed);
double error_in_matrix_calculation(double **matrix, double *solution, int lig);
void print_error_in_calculation(const char *matrix, double error);

int main()
{
	double **mat=NULL;
	double *solutions=NULL;
	double time_elapsed=0.0f, error_calculation=0.0f;
	clock_t start, end;
	
	start=clock();
	printf("MATRICE BORD\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE BORD SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("bord", error_calculation);
	print_time_elapsed("bord", time_elapsed);
	
	start=clock();
	printf("MATRICE DING DONG\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_ding_dong(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE DING DONG SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("ding dong", error_calculation);
	print_time_elapsed("ding dong", time_elapsed);
	
	start=clock();
	printf("MATRICE FRANC\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_franc(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE FRANC SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("franc", error_calculation);
	print_time_elapsed("franc", time_elapsed);
	
	start=clock();
	printf("MATRICE HILBERT\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_hilbert(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE HILBERT SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("hilbert", error_calculation);
	print_time_elapsed("hilbert", time_elapsed);
	
	start=clock();
	printf("MATRICE KMS\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_kms(mat, LIGNES, 0.005f);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE KMS SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("kms", error_calculation);
	print_time_elapsed("kms", time_elapsed);
	
	start=clock();
	printf("MATRICE LEHMER\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_lehmer(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE LEHMER SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("lehmer", error_calculation);
	print_time_elapsed("lehmer", time_elapsed);
	
	start=clock();
	printf("MATRICE LOTKIN\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_lotkin(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE LOTKIN SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("lotkin", error_calculation);
	print_time_elapsed("lotkin", time_elapsed);
	
	start=clock();
	printf("MATRICE MOLER\n\n");
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_moler(mat, LIGNES);
	get_b(mat, LIGNES);
	affiche(mat, LIGNES, COLONNES);
	solutions=gauss_elimination(mat, LIGNES, COLONNES);
	printf("MATRICE MOLER SOLUTIONS\n");
	affiche_solution(solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	end=clock();
	time_elapsed=calculate_time_in_seconds(start, end);
	mat=allocation_dynamic(LIGNES, COLONNES);
	matrice_bord(mat, LIGNES);
	get_b(mat, LIGNES);
	error_calculation=error_in_matrix_calculation(mat, solutions, LIGNES);
	free_matrix(&mat, LIGNES);
	free_matrix_solution(&solutions);
	print_error_in_calculation("moler", error_calculation);
	print_time_elapsed("moler", time_elapsed);

	getchar();
	return 0;
}

double **allocation_dynamic(int lignes, int colonnes)
{
	// allocation de memoire pour les lignes
	double **a_retourner=(double**)malloc(sizeof(double*)*lignes);
	int i=0;

	for( i=0; i<lignes; ++i )
	{
		// allocation de memoire pour des colonnes
		a_retourner[i]=(double*)calloc(colonnes, sizeof(double));
	}

	return a_retourner;
}

void affiche(double **matrix, int lig, int col)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		printf("[");
		for( j=0; j<col; ++j )
		{
			printf(" %f ", matrix[i][j]);
		}
		printf("]\n");
	}
	printf("\n");

	return;
}

void get_b(double **matrix, int lig)
{
	int i=0, k=0;

	for( i=0; i<lig; ++i )
	{
		matrix[i][lig]=b_array[k];
		k++;
	}

	return;
}

double *gauss_elimination(double **matrix, int lig, int col)
{
	int i=0, j=0, k=0;
	double sum=0.0f;
	double *solution=(double*)malloc(sizeof(double)*lig);
	double pivot=0.0f;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<col; ++j )
		{
			if( i>j )
			{
				pivot=matrix[i][j]/matrix[j][j];
				for( k=0; k<col; ++k )
				{
					matrix[i][k]=matrix[i][k]-(pivot*matrix[j][k]);
				}
			}
		}
	}

	for( i=lig-1; i>=0; --i )
	{
		solution[i]=matrix[i][col-1];
		for( j=i; j<col; ++j )
		{
			if( j!=i )
			{
				solution[i]=solution[i]-(matrix[i][j]*solution[j]);
			}
		}
		solution[i]=solution[i]/matrix[i][i];
	}

	return solution;
}

void free_matrix(double ***matrix, int lig)
{
	int i=0;
	double **mat=*matrix;

	for( i=0; i<lig; ++i )
	{
		if( mat[i] )
		{
			free(mat[i]);
		}
	}

	if( mat )
	{
		free(mat);
		*matrix=NULL;
	}
}

void free_matrix_solution(double **matrix)
{
	double *mat=*matrix;

	if( mat )
	{
		free(mat);
		*matrix=NULL;
	}
}

void affiche_solution(double *matrix, int lig)
{
	int i=0;

	for( i=0; i<lig; ++i )
	{
		printf("[ x%d = %lf ]\n", (i+1), matrix[i]);
	}
	printf("\n");
}

/*MATRICE TEST*/
void matrice_bord(double **matrix, int lig)
{
	int i=0;

	// je recopie ce code a partir de ton mail
	for( i=0; i<lig; ++i )
	{
		matrix[i][i]=1;
	}

	for( i=0; i<lig; ++i )
	{
		matrix[0][i]=pow(2.0f, (1-(i+1)));
		matrix[i][0]=matrix[0][i];
	}

	return;
}

void matrice_ding_dong(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=1/(2.0f*lig-(i+1)-(j+1)+1.5f);
		}
	}

	return;
}

void matrice_franc(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			if( (i+1)>=((j+1)+2) )
			{
				matrix[i][j]=0;
			}
			else
			{
				matrix[i][j]=_MIN(i+1, j+1);
			}
		}
	}

	return;
}


void matrice_hilbert(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			matrix[i][j]=(1.0f/((i+1)+(j+1)-1));
		}
	}

	return;
}

void matrice_kms(double **matrix, int lig, double p)
{
	if( p<= 0 || p>= 1 )
	{
		printf("error in matrice_kms\n");
	}
	else
	{
		int i=0, j=0;

		for( i=0; i<lig; ++i )
		{
			for( j=0; j<lig+1; ++j )
			{
				matrix[i][j]=pow(p, fabs((double)(i+1)-(j+1)));
			}
		}
	}

	return;
}

void matrice_lehmer(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i)
	{
		for( j=0; j<lig; ++j)
		{
			if (i<=j)
			{
				matrix[i][j]=((double)(i+1)/(j+1));
			}
			else
			{
				matrix[i][j]=((double)(j+1)/(i+1));
			}
		}
	}

	return;
}

void matrice_lotkin(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i)
	{
		for( j=0; j<lig; ++j)
		{
			if( i==1 )
			{
				matrix[i][j]=1;
			}
			else
			{
				matrix[i][j]=(1.0f/((i+1)+(j+1)-1));
			}
		}
	}

	return;
}

void matrice_moler(double **matrix, int lig)
{
	int i=0, j=0;

	for( i=0; i<lig; ++i)
	{
		for( j=0; j<lig; ++j )
		{
			if( i==j )
			{
				matrix[i][j]=((double)i+1);
			}
			else
			{
				matrix[i][j]=(_MIN(i+1 ,j+1)-2.0f);
			}
		}
	}

	return;
}

double calculate_time_in_seconds(clock_t start, clock_t end)
{
	return ((double)(end-start)/CLOCKS_PER_SEC);
}

void print_time_elapsed(const char *matrix_name, double elapsed)
{
	printf("matrix %s time elapsed : %f seconds\n\n", matrix_name, elapsed);

	return;
}

double error_in_matrix_calculation(double **matrix, double *solution, int lig)
{
	double *ret_val=(double*)calloc(lig, sizeof(double));
	double *B=(double*)calloc(lig, sizeof(double));
	double error=0.0f;
	int i=0, j=0, k=0;

	for( i=0; i<lig; ++i )
	{
		B[i]=matrix[i][lig];
	}

	for( i=0; i<lig; ++i )
	{
		for( j=0; j<lig; ++j )
		{
			ret_val[i]=0.0f;
			for( k=0; k<lig; ++k )
			{
				ret_val[i]=ret_val[i]+(matrix[i][k]*solution[k]);
				//printf("B[%d]=%f\n", k, B[k]);
				//printf("ret_val[%d]=%f\n", i, ret_val[i]);
			}
		}
	}

	for( i=0; i<lig; i++ )
	{
		ret_val[i]=ret_val[i]-B[i];
		error=error+pow(ret_val[i], 2);
	}

	free_matrix_solution(&ret_val);
	error=sqrt(error);

	return error;
}

void print_error_in_calculation(const char *matrix, double error)
{
	printf("error in calculation of matrix %s : %f\n\n", matrix, error);
}
