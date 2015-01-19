/*
 * global_data.cpp: This file is part of the IEC project.
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
 * @file: global_data.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/global_data.hpp"
#include <cassert>
#include <cstdio>

using namespace std;

global_data::global_data(int mn, int oas, bool iro):
  id(mn),
  orig_array_size(oas),
  is_read_only(iro)
{ 
  data_net_info = new net*[oas];
  is_constrained = false;
}


global_data::~global_data()
{
  if( data_net_info ){
    for( int i = 0 ; i < orig_array_size ; i++ )
      delete data_net_info[i];
    delete[] data_net_info;
  }
}


/**
 * \brief Constructor
 *
 * \param mn Identifier for this data
 * \param oas Size of the original array
 * \param of Offset: address of the first elem if we put all global data in a single buffer in ID order.
 * \param iro True if read-only array
 */
global_data_double::global_data_double(int mn, int oas, int of, bool iro):
  global_data(mn,oas,iro)
{
  stride_size = sizeof(double);
  for( int i = 0 ; i < oas ; i++ )
    data_net_info[i] = new net(mn,i,stride_size,of+i);
}

void global_data_double::SetStride(int st)
{
  stride_size = st*sizeof(double);
  for( int i =0 ; i < orig_array_size ; i++ )
    data_net_info[i]->weight = stride_size;
}

global_data_double::~global_data_double()
{}
