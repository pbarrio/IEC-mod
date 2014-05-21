/*
 * orig_smvp.c: This file is part of the IEC project.
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
 * @file: orig_smvp.c
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "stdlib.h"
#include "sys/time.h"
#include "stdio.h"

double* y;
double* A;
double* x;
int N,alen;

/* double rtclock()  */
/* {  */
/*   struct timezone Tzp;  */
/*   struct timeval Tp;  */
/*   int stat;  */
/*   stat = gettimeofday (&Tp, &Tzp);  */
/*   if (stat != 0) printf("Error return from gettimeofday: %d",stat);  */
/*   return(Tp.tv_sec + Tp.tv_usec*1.0e-6);  */
/* } */

void print_output(FILE* outfile)
{
  fprintf(outfile,"N = %d\ny:\n",N);
  int i;
  for( i = 0 ; i  < N*N ; i++ )
    fprintf(outfile," %lf",y[i]);
  fprintf(outfile,"\nA:\n");
  for( i = 0 ; i < alen ; i++ )
    fprintf(outfile," %lf",A[i]);
  fprintf(outfile,"\nx:\n");
  for( i = 0 ; i  < N*N ; i++ )
    fprintf(outfile," %lf",x[i]);
  fprintf(outfile,"\n");
}


int main(int argc, char** argv)
{
  if( argc != 2 ){
    printf("Usage: ./smvp <N>\n");
    return 1;
  }
  
  int *ia,*col;
  N = atoi(argv[1]);
  alen = 4*N*(N-1);

  y = (double*) malloc(N*N*sizeof(double)); 
  x = (double*) malloc(N*N*sizeof(double)); 
  ia = (int*) malloc((N*N+1)*sizeof(int));
  A = (double*) malloc(alen*sizeof(double));
  col = (int*) malloc(alen*sizeof(int));
  
  int i,j,t,k;
  double r;
  int temp;
  for( i = 0 ;i < N*N ; i++ )
    {
      y[i] = (double)rand() / (double)(RAND_MAX-1);
      x[i] = (double)rand() / (double)(RAND_MAX-1);
    }
  for( i = 0 ; i < alen ; i++ )
    A[i] = (double)rand() / (double)(RAND_MAX - 1);

  ia[0] = 0;
  int count = 0;
  for (i = 0; i < N; i++) 
    for (j = 0; j < N; j++) {
      int neighs = 0;
      if (i > 0) {
        neighs++;
        col[count++] = (((i - 1) * N) + j);
      }
      if (j > 0) {
        neighs++;
        col[count++] = (((i * N) + j) - 1);
      }
      if (j < (N - 1)) {
        neighs++;
        col[count++] = (((i * N) + j) + 1);
      }
      if (i < (N - 1)) {
        neighs++;
        col[count++] = (((i + 1) * N) + j);
      }
      ia[((i * N) + j) + 1] = (ia[(i * N) + j] + neighs);
    }


  FILE* outfile = fopen("orig_output.dat","w");
  print_output(outfile);
  for( t = 0 ; t < 1 ; t++ )
    {
      for( i = 0 ; i < N*N ; i++ )
	{
	  r = 0;
	  for( j = ia[i] ; j <  ia[i+1] ; j++ ){
	    temp = col[j];
	    r += A[j] * x[temp];
	  }
	  y[i] = r;
	}
      for( k = 0 ; k < N*N ; k++ )
	x[k] = y[k];
      print_output(outfile);
    }
  fclose(outfile);
  free(y);
  free(x);
  free(ia);
  free(A);
  free(col);
   
  return 0;
}
