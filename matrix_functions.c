//
//  matrix_functions.c
//  qvz2
//  Code mostly taken from:
//  https://www.cs.rochester.edu/~brown/Crypto/assts/projects/adj.html
//

#include "matrix_functions.h"

/*
 Recursive definition of determinate using expansion by minors.
 */
double Determinant(double **a,int n)
{
    int i,j,j1,j2;
    double det = 0;
    double **m = NULL;
    
    if (n < 1) { /* Error */
        
    } else if (n == 1) { /* Shouldn't get used */
        det = a[0][0];
    } else if (n == 2) {
        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {
        det = 0;
        for (j1=0;j1<n;j1++) {
            m = malloc((n-1)*sizeof(double *));
            for (i=0;i<n-1;i++)
                m[i] = malloc((n-1)*sizeof(double));
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                    if (j == j1)
                        continue;
                    m[i-1][j2] = a[i][j];
                    j2++;
                }
            }
            det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
            for (i=0;i<n-1;i++)
                free(m[i]);
            free(m);
        }
    }
    return(det);
}

/*
 Find the cofactor matrix of a square matrix
 */
void CoFactor(double **a,int n,double **b)
{
    int i,j,ii,jj,i1,j1;
    double det;
    double **c;
    
    c = malloc((n-1)*sizeof(double *));
    for (i=0;i<n-1;i++)
        c[i] = malloc((n-1)*sizeof(double));
    
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) {
            
            /* Form the adjoint a_ij */
            i1 = 0;
            for (ii=0;ii<n;ii++) {
                if (ii == i)
                    continue;
                j1 = 0;
                for (jj=0;jj<n;jj++) {
                    if (jj == j)
                        continue;
                    c[i1][j1] = a[ii][jj];
                    j1++;
                }
                i1++;
            }
            
            /* Calculate the determinate */
            det = Determinant(c,n-1);
            
            /* Fill in the elements of the cofactor */
            b[i][j] = pow(-1.0,i+j+2.0) * det;
        }
    }
    for (i=0;i<n-1;i++)
        free(c[i]);
    free(c);
}

/*
 Transpose of a square matrix, do it in place
 */
void Transpose(double **a,int n)
{
    int i,j;
    double tmp;
    
    for (i=1;i<n;i++) {
        for (j=0;j<i;j++) {
            tmp = a[i][j];
            a[i][j] = a[j][i];
            a[j][i] = tmp;
        }
    }
}
/*
Inverse of a Matrix
*/
void Inverse(double** a, double **inv, int n){
    int i,j;
    
    double det;
    
    double **c;
    
    c = malloc((n)*sizeof(double *));
    for (i=0;i<n;i++)
        c[i] = malloc((n)*sizeof(double));
    
    det=Determinant(a,n);
    if (det == 0) {
        printf("ERROR: Matrix not invertible");
        return ;
    }
    else{
        CoFactor(a, n, c);
        Transpose(c, n);
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                inv[i][j] = c[i][j]/det;
            }
        }
    }
}

void MxM(double** a, double** c, double**results, int n, int m, int k){
    int i,j,ii;
    for (i = 0; i < m; i++) {
        for (ii=0; ii<k; ii++) {
            results[i][ii]=0;
            for (j=0; j<n; j++) {
                results[i][ii]+=a[i][j]*c[j][ii];
            }
        }
    }
}

void Mxv(double** a, double* v, double*results, int n, int m){
    int i,j;
    for (i = 0; i < m; i++) {
        results[i]=0;
        for (j=0; j<n; j++) {
            results[i]+=a[i][j]*v[j];
        }
    }
}