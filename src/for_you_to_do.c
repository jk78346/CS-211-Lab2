#include "../include/for_you_to_do.h"

#define MIN(x, y) (((x) < (y))?(x):(y))

int get_block_size(){
    //return the block size you'd like to use 
    /*add your code here */
    return 3;//128;
  
}

/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/


int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
    int i, j, k, t, max_idx, tmp;
    double max_val;
    double* tmp_row = (double*) malloc(n*sizeof(double));
    //ipiv            = (int*)    malloc(n*sizeof(int));
    //for(i = 0 ; i < n ; i++){ipiv[i] = i; printf("ipiv[%d]: %d\n", i, ipiv[i]);}
    double d_tmp;
    //printf("init\n");
    //print_matrix(A, n, n);
    for(i = 0 ; i < n-1 ; i++){
      max_idx = i;
      max_val = fabs(A[i*n+i]);
      for(t = i+1 ; t < n ; t++){
        if(fabs(A[t*n+i]) > max_val){
          max_idx = t;
          max_val = fabs(A[t*n+i]);
        }
      }
      if(max_val == 0){
        free(tmp_row);
        return -1;
      }else{
        if(max_idx != i){ // swap values of ipiv[i] and ipiv[max_idx]
          tmp = ipiv[i]; 
          ipiv[i] = ipiv[max_idx];
          ipiv[max_idx] = tmp;  
          // swap rows
 //         memcpy(tmp_row,       &A[i*n], n*sizeof(double));
 //         memcpy(&A[i*n], &A[max_idx*n], n*sizeof(double));
 //         memcpy(&A[max_idx*n], tmp_row, n*sizeof(double));
          for(t = 0 ; t < n ; t++){
            d_tmp = A[i*n+t];
            A[i*n+t] = A[max_idx*n+t];
            A[max_idx*n+t] = d_tmp;
          } 
        }
      }
      for(j = i+1;  j < n ; j++){
        A[j*n+i] = A[j*n+i] / A[i*n+i];
        for(k = i+1; k < n ; k++){
          A[j*n+k] -= A[j*n+i] * A[i*n+k];
        }
    //    A[j*n+i] = 0;
      }
        //printf("i= %d\n", i);
        //print_matrix(A, n, n);
        //printf("%d %d %d %d\n", ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
    }
    free(tmp_row);
    //printf("after dgetrf\n");
    //print_matrix(A, n, n);
    //printf("ipiv:\n");
    //printf("%d %d %d %d\n", ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
    return 0;
}

/**
 * 
8 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    int i, j;
    double* x = (double*) malloc(n*sizeof(double));
    if(UPLO == 'L'){  
      for(i = 0 ; i < n ; i++){
        x[i] = B[ipiv[i]];
        for(j = 0 ; j < i ; j++){
          x[i] -= A[i*n+j] * x[j];
        }
      }
      for(i = 0 ; i < n ; i++){
         B[i] = x[i];
      }
    }else if(UPLO == 'U'){
      for(i = n-1 ; i >= 0 ; i--){
        for(j = i+1; j < n ; j++){
          B[i] -= A[i*n+j] * B[j];
        }
        B[i] /= A[i*n+i]; 
     }
    }else{
      printf("undefined UPLO: %c\n", UPLO);
      exit(0);
    }
    free(x);
    return;
}

/**
 ( 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
/* add your code here */
	/* please just copy from your lab1 function optimal( ... ) */
  int i1, j1, k1;
  for(k = 0 ; k < n; k += b){
    for(i = 0 ; i < n ; i += b){
      for(j = 0 ; j < n ; j += b){
/*mini block*/
        for(k1 = k ; k1 < k+b ; k1++){
          for(i1 = i ; i1 < i+b ; i1++){
            for(j1 = j ; j1 < j+b ; j1++){
              register double r = C[i1*n+j1];
              r += A[i1*n+k1]*B[k1*n+j1];
              C[i1*n+j1] = r;
            }
          }
        }
/*end mini block*/
      }
    }
  }
  return;
}

/**
 * 
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      
 * 
 *      A     n by n     , square matrix
 * 
 *    
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *     
 *      b                , block size   
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
    int ib, end, i, j, k, t, max_idx, tmp;
    b = get_block_size();
    double max_val;
    double d_tmp;
    for(ib = 0 ; ib < n ; ib+=b){
      end = MIN(ib+b-1, n); // ending idx of this block inclusively
      for(i = ib ; i <= end ; i++){
        // find pivot row t
        max_idx = i;
        max_val = fabs(A[i*n+i]);
        for(t = i+1 ; t < n ; t++){
          if(fabs(A[t*n+i]) > max_val){
            max_idx = t;
            max_val = fabs(A[t*n+i]);
          }
        }
        if(max_val == 0){
          return -1;
        }else{
          if(max_idx != i){
            tmp = ipiv[i];
            ipiv[i] = ipiv[max_idx];
            ipiv[max_idx] = tmp;
            for(k = 0 ; k < n ; k++){
              d_tmp          = A[i*n+k];
              A[i*n+k]       = A[max_idx*n+k];
              A[max_idx*n+k] = d_tmp;
            }    
          }
        }
        //printf("after swap: ib: %d, i: %d\n", ib, i);
        //print_matrix(A, n, n);
        //getchar();
        for(j = i+1 ; j < n ;j++){
          A[j*n+i] /= A[i*n+i];
          for(k = i+1 ; k <= end ; k++){
            A[j*n+k] -= A[j*n+i]*A[i*n+k];
          }
        }
        //printf("after /= and -=, ib: %d, i: %d\n", ib, i);
        //print_matrix(A, n, n);
        //getchar();
      } // end i
      // update right block
      // solve LL*X = A, where A(ib:end, end+1:n) , X is the new A(ib:end, end+1:n)
      for(j = end+1 ; j < n ; j++){
        for(i = ib ; i<= end ;i++){
          for(k = ib ; k < i ; k++){
            A[i*n+j] -= A[i*n+k]/*LL*/*A[k*n+j]; 
          }    
        }
      }
      //printf("after -=, ib: %d, i: %d\n", ib, i);
      //print_matrix(A, n, n);
      //getchar();
      // the GEMM
      for(i = end+1 ; i < n ; i++){
        for(j = end+1 ; j < n ; j++){
          register double the_tmp = 0;
          for(k = ib ; k <= end ; k++){
            the_tmp += A[i*n+k]*A[k*n+j];
          }
          A[i*n+j] -= the_tmp;
          //printf("the_tmp: %f, i: %d, j: %d\n", the_tmp, i, j);
        }
      }
      //printf("after gemm, ib: %d, i: %d\n", ib, i);
      //print_matrix(A, n, n);
      //getchar();
    } // end ib
    //printf("after block\n");
    //print_matrix(A, n, n );
    return 0;
}

