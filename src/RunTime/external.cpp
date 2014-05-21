/*
 * external.cpp: This file is part of the IEC project.
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
 * @file: external.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/external.hpp"
#include <cassert>

using namespace std;

petsc_solve::petsc_solve(int pid, int np, int nt, int s) : proc_id(pid), nprocs(np), nthreads(nt), size(s), orig_net(new net*[s]), new_index(new int[s]), nlocal(new int[nt])
{ 
  for( int i = 0 ; i < s ; i++ )
    orig_net[i] = NULL;
}
 

petsc_solve::~petsc_solve()
{
  delete[] orig_net;
  delete[] new_index;
  delete[] nlocal;
}


void petsc_solve::FindNewRowNumbers()
{
  int *nlocal_rows = new int[nprocs*nthreads];
  int *nproc_offset = new int[nprocs*nthreads];
  
  for( int i = 0; i < nprocs*nthreads ; i++ )
    nlocal_rows[i] = 0;

  for( int i = 0 ; i < size ; i++ ){
    assert(orig_net[i] != NULL && orig_net[i]->home != -1 );
    nlocal_rows[orig_net[i]->home]++;
  }
  
  int total = 0;
  for( int i = 0 ; i < nprocs*nthreads ; i++ ){
    if( i / nthreads == proc_id )
      nlocal[i%nthreads] = nlocal_rows[i];
    nproc_offset[i] = total;
    total = total + nlocal_rows[i];
    nlocal_rows[i] = 0;
  }
  
  for( int i = 0 ; i < size ;i++ ){
    int home = orig_net[i]->home;
    new_index[i] = nproc_offset[home] + nlocal_rows[home];
    nlocal_rows[home]++;
  }
 
  delete[] nlocal_rows;
  delete[] nproc_offset;
}

void petsc_solve::RenumberGlobalRows(int* orig_nums, int array_size) const
{
  for( int i = 0 ; i < array_size ; i++ ){
    assert(orig_nums[i] < size) ;
    orig_nums[i] = new_index[orig_nums[i]];
  }
}


void petsc_solve::print(FILE* outfile) const
{
  fprintf(outfile,"Matrix Size = %d, Nlocal = ",size);
  for( int i =0 ; i < nthreads ; i++ )
    fprintf(outfile," %d",nlocal[i]);
  fprintf(outfile,"\n");
#ifdef HIGH_DETAILS
  fprintf(outfile,"Old -> New :\n");
  for( int i = 0 ; i < size ; i++ )
    fprintf(outfile,"Net:%d, %d -> %d\n",orig_net[i]->my_num,i,new_index[i]);
#endif
}
