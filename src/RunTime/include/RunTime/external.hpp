/*
 * external.hpp: This file is part of the IEC project.
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
 * @file: external.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __EXTERNAL_HPP__
#define __EXTERNAL_HPP__

#include "RunTime/hypergraph.hpp"
#include <cstdio>
#include <set>

class inspector;

class petsc_solve{

private:

  net** const orig_net;
  
  int* new_index;

  const int size;

  const int nprocs;

  const int nthreads;

  const int proc_id;

  int* const nlocal;

public:

  petsc_solve(int,int,int,int);

  ~petsc_solve();  

  inline int GetNLocalRows( int tid) const { return nlocal[tid]; }

  void FindNewRowNumbers();

  void RenumberGlobalRows(int*,int) const;

  void print(FILE*) const;

  inline int GetLocalRows(int tid) const { return nlocal[tid]; }
  
  friend class inspector;
};


#endif
