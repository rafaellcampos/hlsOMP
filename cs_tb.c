#include <stdio.h>
#include "cs.h"
//Banco de pruebas para algoritmo OMP en Compressed Sensing
//Para probar diferentes tama�os de matrices cambie el tama�o de A
//Tambi�n se puede cambiar la se�al y comprimida
//Para cambiar el tam�o de matrices h�galo desde el archivo cs.h



//Funci�n para imprimir columnas en el banco de pruebas
static void printRophi(int rophi, double *matrix, int cols)
{

	if (cols >= 1) printf(" %3.3f", matrix[rophi*cols+0]);
	if (cols >= 2) printf(" %3.3f", matrix[rophi*cols+1]);
	if (cols >= 3) printf(" %3.3f", matrix[rophi*cols+2]);
	if (cols >= 6) printf(" ...");
	if (cols >= 5) printf(" %3.3f", matrix[rophi*cols+(cols-2)]);
	if (cols >= 4) printf(" %3.3f", matrix[rophi*cols+(cols-1)]);
	printf(" \n");
}
//Funci�n para imprimir matrices (para visualizar resultados y debuggeo)
static void printMatrix(double *matrix, int rophis, int cols)
{
	if (rophis >= 1) printRophi(0, matrix, cols);
	if (rophis >= 2) printRophi(1, matrix, cols);
	if (rophis >= 3) printRophi(2, matrix, cols);
	if (rophis >= 6) printf("  ...\n");
	if (rophis >= 5) printRophi(rophis-2, matrix, cols);
	if (rophis >= 4) printRophi(rophis-1, matrix, cols);
}


//Principal funci�n del banco de pruebas
int main(int argc, char **argv)
{

	double A[M*N] = {
			-80, 30, 100, 40,
			-20, 40, -30, -40,
			 20, 100, -10, 80
	};

	double y[M] ={
			 270,
			  10,
			 450
	};

	double Rec[N];
	//Indique la densidad de la se�al
	int k = 3;


   principal(A, M, N, y,Rec,k);

   printf("Se�al recuperada: = \n"); printMatrix(Rec, 1, N);

   return 0;
}
