/*
 * orig_cg.c: This file is part of the IEC project.
 *
 * IEC : Inspector/Executor Compiler
 *
 * Copyright (C) 2011 the Ohio State University
 *
 * This program can be redistributed and/or modified under the terms
 * of the license specified in the LICENSE.txt file at the root of the
 * project.
 *
 * Contact: P. Sadayappan <saday@cse.ohio-state.edu>
 *
 */
/**
 * @file: orig_cg.c
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>, John Eisenlohr <eisenloh@cse.ohio-state.edu>
 */
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "sys/time.h"
#include "ie.h"

#define MAX_ITER 500000

int readRBfile( FILE *fpIn, int *n, int *nz, int **rows, int **cols, double **vals )
{
  char header[128];
  fgets( header, sizeof(header), fpIn );
  fgets( header, sizeof(header), fpIn );
  fscanf(fpIn,"%s",header);

  // temporary values for compressed format
  int i,j,k;
  int nzC;
  int *rowsC, *colsC, *rowCounts;

  fscanf(fpIn,"%d",n);
  fscanf(fpIn,"%d",n);
  fscanf(fpIn,"%d",&nzC);
#ifdef DEBUG
  printf("n = %d\n", *n);
  printf("original nz = %d\n", nzC);
#endif

  fgets( header, sizeof(header), fpIn );
  fgets( header, sizeof(header), fpIn );

  rowsC = (int*)malloc( (*(n)+1)*sizeof(int) );
  colsC = (int*)malloc(nzC*sizeof(int));
  rowCounts = (int*)malloc( (*n)*sizeof(int) );
  

  // read in the rows array and initialize the
  // rowCounts rows will be corrected later
  for ( j = 0; j < (*n)+1; j++ ) {
    fscanf( fpIn, "%d", &(rowsC[j]) );
    rowsC[j]--;
    if ( j > 0 ) {
      rowCounts[j-1] = rowsC[j] - rowsC[j-1];
#ifdef DEBUG
      printf("rows[%d]=%d, rowCounts[%d]=%d\n", j, rowsC[j], j-1, rowCounts[j-1]);
#endif
    }
  }

  // read in the original compressed cols array
  // (no duplicates for symmetry)
  for (  j = 0; j < nzC; j++ ) {
    fscanf( fpIn, "%d", &(colsC[j]) );
    colsC[j]--;
    //    printf("colsC[%d]=%d\n", j, colsC[j]);
  }

#ifdef DEBUG
  printf("got original cols\n");
  fflush(stdout);
#endif

  // increment rowCounts for cols not recorded in compressed format
  for (  j = 0; j < (*n); j++ ) {
    fflush(stdout);
    for (  k = rowsC[j]; k < rowsC[j+1]; k++ ) {
      int col = colsC[k];
      if ( col > j )
	rowCounts[col]++;
    }
  }

#ifdef DEBUG
  printf("incremented rowCounts\n");
  fflush(stdout);
#endif

  (*rows) = (int*)malloc( ((*n)+1)*sizeof(int) );
  (*rows)[0] = 0;
  for ( j = 1; j < (*n)+1; j++ ) {
    (*rows)[j] = (*rows)[j-1] + rowCounts[j-1];
    rowCounts[j-1] = 0;
  }

  (*nz) = (*rows)[(*n)];
#ifdef DEBUG
  printf("uncompressed nz = %d\n", (*nz));
  fflush(stdout);
#endif

  (*cols) = (int*)malloc((*nz)*sizeof(int));
  (*vals) = (double*)malloc((*nz)*sizeof(double));  

  // in the following loop, we are going to be updating
  // the cols and vals arrays with duplicates for symmetry
  // We need two global indices into the cols and vals
  // arrays -- one for [i,j] and one for [j,i].
  // These are named thisRowGlobalIndex and laterRowGlobalIndex
  // These are computed as follows -- the rows array gives us
  // a base index into cols and vals for any row; we also
  // keep track of how many elements we have written into
  // cols and vals so far for each row so this gives an
  // offset index into the global array. The sum of the base
  // and offset is the global index.
  int thisRowGlobalIndex, laterRowGlobalIndex;
  double elt;
  for ( i = 0; i < (*n); i++ ) {
    for ( k = rowsC[i]; k < rowsC[i+1]; k++ ) {
      j = colsC[k];
      thisRowGlobalIndex = (*rows)[i]+rowCounts[i];
      rowCounts[i]++;
      // transfer the next column index in the row
      (*cols)[thisRowGlobalIndex] = j;
      // read the value and put in this row's data
      fscanf( fpIn, "%lf", &elt );
      (*vals)[thisRowGlobalIndex] = elt;
      // update later row for symmetry
      if ( j > i ) {
	laterRowGlobalIndex = (*rows)[j] + rowCounts[j];
	(*cols)[laterRowGlobalIndex] = i;
	(*vals)[laterRowGlobalIndex] = elt;
	rowCounts[j]++;
      }
    }
  }

#ifdef DEBUG
  // write out a text file
  FILE *fp = fopen("checkFile.txt", "w");
  fprintf(fp,"%d\n",(*n));
  fprintf(fp,"%d\n",(*nz));

  for ( i = 0; i < (*nz); i++ )
    fprintf(fp,"%lf\n",(*vals)[i]);

  for ( i = 0; i < (*nz); i++ )
    fprintf(fp,"%d\n",(*cols)[i]);

  for ( i = 0; i < (*n)+1; i++ )
    fprintf(fp,"%d\n",(*rows)[i]);
#endif

  return 1;
}


double rtclock() 
{ 
  struct timezone Tzp; 
  struct timeval Tp; 
  int stat; 
  stat = gettimeofday (&Tp, &Tzp); 
  if (stat != 0) printf("Error return from gettimeofday: %d",stat); 
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6); 
}
 
void csr_matrix_vector_multiply( int n,
				 double *vals, int *rows, int *cols,
				 double *x, double *y )
{
  int i, j, k;
  for ( i = 0; i < n; i++ ) {
    y[i] = 0;
    for ( k = rows[i]; k < rows[i+1]; k++ ) {
      j = cols[k];
      y[i] += vals[k]*x[j];
    }
  }
}

void generate_initial_guess( int n,
			     double *vals, int *rows, int *cols,
			     double *x, double *y )
{
  int j;
  for ( j = 0; j < n; j++ )
    x[j] = 0;
}


int main(int argc, char* argv[])
{
  if ( argc != 2 ) {
    printf("***************************\n***************************\n\n");
    printf("        NO INPUT FILE \n\n");
    printf("  Usage: ./cg matNNN.txt\n\n");
    printf("***************************\n");
    exit(1);
  }

  char* fname = argv[1];
  FILE *fpIn;
  int n, nz;
  double elt;
  int i,j, k, index;
  double tol = 0.0000000000000001;
  fpIn = fopen(fname,"r");

  int *rows, *cols;
  double *vals;

  size_t fnl = strlen(fname);
  if ( fname[fnl-2]=='r' && fname[fnl-1]=='b' ) {
    printf("File is in RB format\n");
    int ioStatus = readRBfile( fpIn, &n, &nz, &rows, &cols, &vals );
    fclose(fpIn);
    if ( ioStatus < 0 ) {
      printf("Error reading file\n");
      exit(1);
    }
  }
  else{
    fscanf(fpIn,"%d",&n);
    fscanf(fpIn,"%d",&nz);
    vals = (double*)malloc(nz*sizeof(double));
    rows = (int*)malloc((n+1)*sizeof(int));
    cols = (int*)malloc(nz*sizeof(int));
    for ( j = 0; j < nz; j++ )
      fscanf( fpIn, "%lf", &(vals[j]) );
    for ( j = 0; j < nz; j++ )
      fscanf( fpIn, "%d", &(cols[j]) );
    for ( j = 0; j < n+1; j++ )
      fscanf( fpIn, "%d", &(rows[j]) );
  }


#ifdef DEBUG
  printf("first row cols:\n");
  for ( j = rows[0]; j < rows[1]; j++ )
    printf("%d ",cols[j]);
  printf("\n");
  printf("last row cols:\n");
  for ( j = rows[n-1]; j < rows[n]; j++ )
    printf("%d ",cols[j]);
  printf("\n");

  printf("first non-zero elt %lf\n", vals[0]);
  printf("last non-zero elt %lf\n", vals[nz-1]);
#endif

  //  exit(1);

  double *x = (double*)malloc(n*sizeof(double));
  double *y = (double*)malloc(n*sizeof(double));

  time_t seconds;
  time(&seconds);
  //srand((unsigned int) seconds);
  for (  j = 0; j < n; j++ )
    y[j] = -1 + 2*( (double)rand() / (double)RAND_MAX );

  /* .................. */
  
  int result = 0;
  double rho, old_rho, beta, alpha, qdotp;
  double *r = (double*)malloc(n*sizeof(double));
  double *p_new = (double*)malloc(n*sizeof(double));
  double *p_old = (double*)malloc(n*sizeof(double));
  double *q = (double*)malloc(n*sizeof(double));
  double *Ax = (double*)malloc(n*sizeof(double));
  double temp;

  generate_initial_guess( n, vals, rows, cols, x, y );

  csr_matrix_vector_multiply( n, vals, rows, cols, x, Ax );
  
  for ( j = 0; j < n; j++ )
    r[j] = y[j] - Ax[j];
  int iter = 0;
  
  rho = 0.0;
  for(  j = 0 ; j < n ; j++ )
    {
      rho += r[j]*r[j];
      p_new[j] = r[j];
    }  

  int i_1,j_1,k_1,j_2,j_3,j_4;

  double start_t = rtclock();
  while ( rho >= tol && iter < MAX_ITER ) {
    iter++;
    /* printf("iter %d, rho = %lf\n", iter, rho); */

    for ( i_1 = 0; i_1 < n; i_1++ ) {
      temp =  0;
      for ( k_1 = rows[i_1]; k_1 < rows[i_1+1]; k_1++ ) {
	temp += vals[k_1]*p_new[cols[k_1]];
      }
      q[i_1] = temp;
    }

    // <q,p>
    qdotp = 0;
    for( j_1 = 0; j_1 < n; j_1++ )
      qdotp += q[j_1]*p_new[j_1];
    alpha = rho / qdotp;
    
    for( j_2 = 0; j_2 < n; j_2++ ) {
      x[j_2] += alpha*p_new[j_2];
      r[j_2] -= alpha*q[j_2];
      p_old[j_2] = p_new[j_2];
    }

    old_rho = rho;
    // <r,r>
    rho = 0;

    for( j_3 = 0; j_3 < n; j_3++ )
      rho += r[j_3]*r[j_3];
    
    beta = rho / old_rho;
    for( j_4 = 0; j_4 < n; j_4++ )
      p_new[j_4] = r[j_4] + beta * p_old[j_4] ;
  }
  double stop_t = rtclock();
  
  fprintf(stderr,"[IEC]:Original:%lf\n",stop_t - start_t);
  fflush(stderr);

  /* .................. */
  if ( iter > 0 ) {
    printf("CG converged in %d iterations\n",iter);
    double *checkY = (double*)malloc(n*sizeof(double));
    double error, maxError = 0;
    int maxErrorIndex = 0;
    csr_matrix_vector_multiply( n, vals, rows, cols, x, checkY );
    for ( j = 0; j < n; j++ ) {
      error = fabs( y[j] - checkY[j] );
      if ( error > maxError ) {
	maxError = error;
	maxErrorIndex = j;
      }
    }
    printf("Max error at index %d, error %lf\n", maxErrorIndex, maxError);
  }
  else
    printf("CG did not converge\n");

}
