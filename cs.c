#include "cs.h"
#include <float.h>

// Módulo SVD basado en el trabajo de
// http://www.mymathlib.com/matrices/linearsystems/singular_value.html


void principal(double A[M*N], int nrows, int ncols, double y[M], double x[N], int K){

	static double W[N];
	static double An[M*N];
	static double r[M];
	static double y_l[M];
	static double x_l[N];
	static double aux[M];
	int test;

	copyYtoR: for(int cR = 0; cR<M; cR++) r[cR] = y[cR];
	Anto0:    for(int cA = 0; cA<M*N; cA++) An[cA] = 0;
	copyYtoY: for(int cY = 0; cY<M; cY++) y_l[cY] = y[cY];


	mainLoop: for(int it = 0; it<K; it++){

		int maximo;
		matrixW(A, r, W);
		maximo = buscar_maximo(W, N);
		cambioAtom(An, A, maximo, N, M*N);
		atom0(A, maximo, N, M*N);
		test = Singular_Value_Decomposition(An, M, N, y_l, x_l);
		matrixR(An, x_l, aux);
		resta(y_l, aux, r);
	}
	cpyXtoX: for(int cX = 0; cX<N; cX++) x[cX] = x_l[cX];


}

int Singular_Value_Decomposition(double A[M*N], int nrows, int ncols,
								 double B[M], double x[N]){

	static double U[M*N];
	static double singular_values[N];
	static double V[N*N];
	static double dummy_array[N];

   Householders_Reduction_to_Bidiagonal_Form( A, nrows, ncols, U, V,
                                                singular_values, dummy_array);

   if (Givens_Reduction_to_Diagonal_Form( nrows, ncols, U, V,
                                singular_values, dummy_array ) < 0) return -1;

   Sort_by_Decreasing_Singular_Values(nrows, ncols, singular_values, U, V);
   double tolerance = 0.01;

   Singular_Value_Decomposition_Solve(U,singular_values,V,tolerance,nrows,ncols,B,x);
  
   return 0;
}

void Householders_Reduction_to_Bidiagonal_Form(double A[M*N], int nrows,
                                                      int ncols, double U[M*N], double V[N*N],
                                                      double diagonal[N], double superdiagonal[N]){
   int i,j,k,ip1;
   double s, s2, si, scale;
   double dum;
   
   double *pu = 0;
   double *pv = 0;
   double *pui = 0;
   double *pvi = 0;
   
   double half_norm_squared;

// Copy A to U

   //memcpy(U,A, sizeof(double) * nrows * ncols);
   cpyA:for(int copy = 0; copy<M*N; copy ++)U[copy] = A[copy];
   
//
 
   diagonal[0] = 0.0;
   s = 0.0;
   scale = 0.0;
   primer:for ( i = 0, pui = U, ip1 = 1; i < ncols; pui += ncols, i++, ip1++ ) {
      superdiagonal[i] = scale * s;
//       
//                  Perform Householder transform on coluM*Ns.
//
//       Calculate the normed squared of the i-th coluM*N vector starting at
//       row i.
//
     segundo:for (j = i, pu = pui, scale = 0.0; j < nrows; j++, pu += ncols)
         scale += fabs( *(pu + i) );
       
      if (scale > 0.0) {
         for2:for (j = i, pu = pui, s2 = 0.0; j < nrows; j++, pu += ncols) {

            *(pu + i) /= scale;
            s2 += *(pu + i) * *(pu + i);
         }
//
//    
//       Chose sign of s which maximizes the norm
//  
         s = ( *(pui + i) < 0.0 ) ? sqrt(s2) : -sqrt(s2);
//
//       Calculate -2/u'u
//
         half_norm_squared = *(pui + i) * s - s2;
//
//       Transform remaining coluM*Ns by the Householder transform.
//
         *(pui + i) -= s;
         
          tercero:for (j = ip1; j < ncols; j++) {
           cuarto:for (k = i, si = 0.0, pu = pui; k < nrows; k++, pu += ncols)
               si += *(pu + i) * *(pu + j);
            si /= half_norm_squared;
            house1:for (k = i, pu = pui; k < nrows; k++, pu += ncols) {
               *(pu + j) += si * *(pu + i);
            }
         }
      }
      house2:for (j = i, pu = pui; j < nrows; j++, pu += ncols) *(pu + i) *= scale;
      diagonal[i] = s * scale;
//       
//                  Perform Householder transform on rows.
//
//       Calculate the normed squared of the i-th row vector starting at 
//       coluM*N i.
//
      s = 0.0;
      scale = 0.0;
      if (i >= nrows || i == (ncols - 1) ) continue;
      house3:for (j = ip1; j < ncols; j++) scale += fabs ( *(pui + j) );
      if ( scale > 0.0 ) {
         house4:for (j = ip1, s2 = 0.0; j < ncols; j++) {
            *(pui + j) /= scale;
            s2 += *(pui + j) * *(pui + j);
         }
         s = ( *(pui + ip1) < 0.0 ) ? sqrt(s2) : -sqrt(s2);
//
//       Calculate -2/u'u
//
         half_norm_squared = *(pui + ip1) * s - s2;
//
//       Transform the rows by the Householder transform.
//
         *(pui + ip1) -= s;
         house5:for (k = ip1; k < ncols; k++)
            superdiagonal[k] = *(pui + k) / half_norm_squared;
         if ( i < (nrows - 1) ) {
            house6:for (j = ip1, pu = pui + ncols; j < nrows; j++, pu += ncols) {
               house7:for (k = ip1, si = 0.0; k < ncols; k++)
                  si += *(pui + k) * *(pu + k);
               house8:for (k = ip1; k < ncols; k++) {
                  *(pu + k) += si * superdiagonal[k];
               }
            }
         }
         house9:for (k = ip1; k < ncols; k++) *(pui + k) *= scale;
      }
   }

// Update V
   pui = U + ncols * (ncols - 2);
   pvi = V + ncols * (ncols - 1);
   *(pvi + ncols - 1) = 1.0;
   s = superdiagonal[ncols - 1];
   pvi -= ncols;
   house10:for (i = ncols - 2, ip1 = ncols - 1; i >= 0; i--, pui -= ncols,
                                                      pvi -= ncols, ip1-- ) {
      if ( s != 0.0 ) {
         pv = pvi + ncols;
         house11:for (j = ip1; j < ncols; j++, pv += ncols)
            *(pv + i) = ( *(pui + j) / *(pui + ip1) ) / s;
         house12:for (j = ip1; j < ncols; j++) {
            si = 0.0;
            house13:for (k = ip1, pv = pvi + ncols; k < ncols; k++, pv += ncols)
               si += *(pui + k) * *(pv + j);
            house14:for (k = ip1, pv = pvi + ncols; k < ncols; k++, pv += ncols)
               *(pv + j) += si * *(pv + i);                  
         }
      }
      pv = pvi + ncols;
      house15:for ( j = ip1; j < ncols; j++, pv += ncols ) {
         *(pvi + j) = 0.0;
         *(pv + i) = 0.0;
      }
      *(pvi + i) = 1.0;
      s = superdiagonal[i];
   }

// Update U

   pui = U + ncols * (ncols - 1);
   house16:for (i = ncols - 1, ip1 = ncols; i >= 0; ip1 = i, i--, pui -= ncols ) {
      s = diagonal[i];
      house17:for ( j = ip1; j < ncols; j++) *(pui + j) = 0.0;
      if ( s != 0.0 ) {
         house18:for (j = ip1; j < ncols; j++) {
            si = 0.0;
            pu = pui + ncols;
            house19:for (k = ip1; k < nrows; k++, pu += ncols)
               si += *(pu + i) * *(pu + j);
            si = (si / *(pui + i) ) / s;
            house20:for (k = i, pu = pui; k < nrows; k++, pu += ncols)
               *(pu + j) += si * *(pu + i);                  
         }
         house21:for (j = i, pu = pui; j < nrows; j++, pu += ncols){
            *(pu + i) /= s;
         }
      }
      else 
         house22:for (j = i, pu = pui; j < nrows; j++, pu += ncols) *(pu + i) = 0.0;
      *(pui + i) += 1.0;
   }
}

int Givens_Reduction_to_Diagonal_Form( int nrows, int ncols,
                                              double U[M*N], double V[N*N],
                                              double diagonal[N], double superdiagonal[N]){

   double epsilon;
   double c, s;
   double f,g,h;
   double x,y,z;
   double *pu = 0;
   double *pv = 0;
   int i,j,k,m;
   int rotation_test;
   int iteration_count;
  
   red0:for (i = 0, x = 0.0; i < ncols; i++) {
      y = fabs(diagonal[i]) + fabs(superdiagonal[i]);
      if ( x < y ) x = y;
   }
   epsilon = x * DBL_EPSILON;
   red1:for (k = ncols - 1; k >= 0; k--) {
      iteration_count = 0;
      red2:while(1) {
         rotation_test = 1;
         red3:for (m = k; m >= 0; m--) {
            if (fabs(superdiagonal[m]) <= epsilon) {rotation_test = 0; break;}
            if (fabs(diagonal[m-1]) <= epsilon) break;
         }
         if (rotation_test) {
            c = 0.0;
            s = 1.0;
            red4:for (i = m; i <= k; i++) {
               f = s * superdiagonal[i];
               superdiagonal[i] *= c;
               if (fabs(f) <= epsilon) break;
               g = diagonal[i];
               h = sqrt(f*f + g*g);
               diagonal[i] = h;
               c = g / h;
               s = -f / h; 
               red5:for (j = 0, pu = U; j < nrows; j++, pu += ncols) {
                  y = *(pu + m - 1);
                  z = *(pu + i);
                  *(pu + m - 1 ) = y * c + z * s;
                  *(pu + i) = -y * s + z * c;
               }
            }
         }
         z = diagonal[k];
         if (m == k ) {
            if ( z < 0.0 ) {
               diagonal[k] = -z;
               red6:for ( j = 0, pv = V; j < ncols; j++, pv += ncols)
                  *(pv + k) = - *(pv + k);
            }
            break;
         }
         else {
            if ( iteration_count >= MAX_ITERATION_COUNT ) return -1;
            iteration_count++;
            x = diagonal[m];
            y = diagonal[k-1];
            g = superdiagonal[k-1];
            h = superdiagonal[k];
            f = ( (y - z) * ( y + z ) + (g - h) * (g + h) )/(2.0 * h * y);
            g = sqrt( f * f + 1.0 );
            if ( f < 0.0 ) g = -g;
            f = ( (x - z) * (x + z) + h * (y / (f + g) - h) ) / x;
// Next QR Transformtion
            c = 1.0;
            s = 1.0;
            red7:for (i = m + 1; i <= k; i++) {
               g = superdiagonal[i];
               y = diagonal[i];
               h = s * g;
               g *= c;
               z = sqrt( f * f + h * h );
               superdiagonal[i-1] = z;
               c = f / z;
               s = h / z;
               f =  x * c + g * s;
               g = -x * s + g * c;
               h = y * s;
               y *= c;
               red8:for (j = 0, pv = V; j < ncols; j++, pv += ncols) {
                  x = *(pv + i - 1);
                  z = *(pv + i);
                  *(pv + i - 1) = x * c + z * s;
                  *(pv + i) = -x * s + z * c;
               }
               z = sqrt( f * f + h * h );
               diagonal[i - 1] = z;
               if (z != 0.0) {
                  c = f / z;
                  s = h / z;
               } 
               f = c * g + s * y;
               x = -s * g + c * y;
               red9:for (j = 0, pu = U; j < nrows; j++, pu += ncols) {
                  y = *(pu + i - 1);
                  z = *(pu + i);
                  *(pu + i - 1) = c * y + s * z;
                  *(pu + i) = -s * y + c * z;
               }
            }
            superdiagonal[m] = 0.0;
            superdiagonal[k] = f;
            diagonal[k] = x;
         }
      } 
   }
   return 0;
}

void Sort_by_Decreasing_Singular_Values(int nrows, int ncols,
                                               double singular_values[N],
                                               double U[M*N], double V[N*N]){
   int i,j,max_index;
   double temp;
   double *p1 = 0;
   double *p2 = 0;
   //double *p1 = (double *)malloc(sizeof(double)*M*N);
   //double *p2 = (double *)malloc(sizeof(double)*M*N);

Sort0:for (i = 0; i < ncols - 1; i++) {
      max_index = i;
Sort1:for (j = i + 1; j < ncols; j++)
         if (singular_values[j] > singular_values[max_index] ) 
            max_index = j;
      if (max_index == i) continue;
      temp = singular_values[i];
      singular_values[i] = singular_values[max_index];
      singular_values[max_index] = temp;
      p1 = U + max_index;
      p2 = U + i;
Sort2:for (j = 0; j < nrows; j++, p1 += ncols, p2 += ncols) {
         temp = *p1;
         *p1 = *p2;
         *p2 = temp;
      } 
      p1 = V + max_index;
      p2 = V + i;
Sort3:for (j = 0; j < ncols; j++, p1 += ncols, p2 += ncols) {
         temp = *p1;
         *p1 = *p2;
         *p2 = temp;
      }
   } 
}

void Singular_Value_Decomposition_Solve( double U[M*N],  double D[N],
		 	 	 	 	 	 	 	 	 double V[N*N], double tolerance,
										int nrows, int ncols,
										double B[M], double x[N]){
	   int i,j,k;
	   double *pu = 0;
	   double *pv = 0;

	   double dum;
	   dum = DBL_EPSILON * D[0] * (double) ncols;
	   if (tolerance < dum) tolerance = dum;
	   pv = V;
	S1:  for ( i = 0; i < ncols; i++) {
	         x[i] = 0.0;
	S2:     for (j = 0; j < ncols; j++){
	           if (D[j] > tolerance ) {
	        	   	dum = 0.0;
	        	   	pu = U;
	S3:			for (k = 0; k < nrows; k++){
	            			dum += *(pu + j) * B[k];
	            			pu += ncols;
	            		}
	            		x[i] += dum * *(pv + j) / D[j];
	            }
			}
	         pv += ncols;
	      }
	}

void matrixW(double A[M*N], double r[M], double W[N]){
    int cols = N;
    int rows = M;
    double res;
rowsW:for(int j = 0; j < cols; j++){
	   res = 0;
colsW:   for(int i = 0; i < rows; i++){
	       res = res + A[i*N+j] * r[i];
	    }
	    W[j] = res;
	}
}

void matrixR(double A[M*N], double x[N], double R[M]){
	int cols = N;
	int rows = M;
	double res;
rows:	for(int i = 0; i < rows; i++){
	   res = 0;
	  cols: for(int j = 0; j < cols; j++){
	       res = res + A[i*N+j] * x[j];
	    }
	    R[i] = res;
		}
}

void cambioAtom(double recibe[M*N], double entrega[M*N], int j, int W, int tot){

	int res = 0;
 atomC: for(int i = j; i < tot; i = i + W){
    	recibe[i] = entrega[j + res];
    	res = res + W;
    }
}

void atom0(double A[M*N], int j, int W, int tot){

	atom0: for(int i = j; i < tot; i = i + W){
	   	A[i] = 0;
	   }
}

void resta(double A[M], double B[M], double C[M]){
resta: for(int i = 0; i<M; i++){
		C[i] = A[i] - B[i];
	}
}

int buscar_maximo(double valores[N], int n){

	int maximo_pos = 0;
max: for (int i = 0; i < n; i++) {
		if (valores[i] > valores[maximo_pos]) maximo_pos = i;
	}
	return maximo_pos;
}

