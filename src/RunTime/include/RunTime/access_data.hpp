/*
 * access_data.hpp: This file is part of the IEC project.
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
 * @file: access_data.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef ACCESS_DATA_HPP
#define ACCESS_DATA_HPP

#include "stdio.h"
#include <map>
#include <set>
#ifdef USE_HSQPHASH
#include "DataStructures/qphash.hpp"
#endif

class inspector;

class access_data{
 private:
#ifdef USE_HSQPHASH
  qphash* have_set;
#else
  std::map<int,int> have_set;
#endif

  /// List of positions in this indirection array that we don't have
  std::set<int> dont_have_set;

  int array_size;
  const int nprocs;
  const int myid;
  const int my_num;
  const int myTeamID;
  const int myTeamSize;
  int local_start_index;
  int local_end_index;
  int stride;
  int* orig_array;

 public:

  access_data(int mn, int md, int np, int teamID, int teamSize);

  ~access_data();

  void SetParams(int,int,int*);
  
  bool HaveIndex(int);

  int GetIndex(int) const;

  ///Get the count of number of elements needed from other processes
  ///All indirection arrays assumed to be blocked partitioned initially
  ///Arguments : 
  /// 1 )Buffer, viewed as a 2D array [nprocs][n_ind_arrays]
  /// 2) n_ind_arrays
  void GetSendCounts(int*,int) const;

  ///Send indices whose values are to be retrieved
  /// Arguments :
  ///  1) Buffer used for communication 
  ///  2) the size of the buffer
  ///  3) the offset in the buffer for the message sent to each process
  void PopulateBuffer(int*,int,int*) const;

  ///Send the values of corresponding indices
  ///Arguments
  ///  1) The buffer to be used to read the indices, the same buffer is replaced with the values
  ///  2) The size of the buffer
  ///  3) The offsets in the buffer which pointing to the start of the requested elements of the current indirection array
  ///  4) The number of elements that are requested for this array [nprocs][n_ind_arrays]
  ///  5) n_ind_arrays
  void GetRequestedValue(int*,int,int*,int*,int) const;

  void AddToHaveSet(int*,int,int*);

#ifndef NDEBUG
  void print_access(FILE*);
#endif

  friend class inspector;
};  


#endif
