/*
 * alloc.h: This file is part of the IEC project.
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
 * @file: alloc.h
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 * @acknowledgements : Sandip Mazumder <mazumder.2@osu.edu>
 */
#ifndef __ALLOC_H__
#define __ALLOC_H__

double** malloc_2d(int dim0, int dim1)
{
  double* base = (double*)malloc(sizeof(double)*dim0*dim1);
  double** final = (double**) malloc(sizeof(double*)*dim0);
  int i;
  for( i = 0 ; i < dim0 ; i++ )
    final[i] = base + dim1 * i;
  return final;
}

void free_2d(double** final)
{
  free(final[0]);
  free(final);
}

int** malloc_2d_int(int dim0, int dim1)
{
  int* base = (int*)malloc(sizeof(int)*dim0*dim1);
  int** final = (int**) malloc(sizeof(int*)*dim0);
  int i;
  for( i = 0 ; i < dim0 ; i++ )
    final[i] = base + dim1 * i;
  return final;
}

void free_2d_int(int** final)
{
  free(final[0]);
  free(final);
}


#endif
