#ifndef CS_H_
#define CS_H_

#include <string.h>              // required for memcpy()
#include <float.h>               // required for DBL_EPSILON
#include <math.h>                // required for fabs(), sqrt();
#include <stdlib.h>
#include <stdio.h>

#define MAX_ITERATION_COUNT 100

//Número de filas y columnas:
#define M  3
#define N  4

void principal(double A[M*N], int nrows, int ncols, double y[M], double x[N], int K);
										
int Singular_Value_Decomposition(double A[M*N], int nrows, int ncols,
								 double B[M], double x[N]);

 void Householders_Reduction_to_Bidiagonal_Form(double A[M*N], int nrows,
                                                      int ncols, double U[M*N], double V[N*N],
                                                      double diagonal[N], double superdiagonal[N]);

 int Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,
                                              double U[M*N], double V[N*N],
                                              double diagonal[N], double superdiagonal[N]);
                                              
 void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,
                                               double singular_values[N],
                                               double U[M*N], double V[N*N]);

void Singular_Value_Decomposition_Solve( double U[M*N],  double D[N],
		 	 	 	 	 	 	 	 	 double V[N*N], double tolerance,
										int nrows, int ncols,
										double B[M], double x[N]);

void matrixW(double A[M*N], double r[M], double W[N]);

void matrixR(double A[M*N], double x[N], double R[M]);

void cambioAtom(double recibe[M*N], double entrega[M*N], int j, int W, int tot);

void resta(double A[M], double B[M], double C[M]);

void atom0(double A[M*N], int j, int W, int tot);

int maximoW(double matrix[N], int n);

int buscar_maximo(double valores[N], int n);



#endif
