/*
 * data_scalars.cpp: This file is part of the IEC project.
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
 * @file: data_scalars.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/data_scalars.hpp"
#include "CompileTime/loops.hpp"

using namespace std;

void data_scalar::AddAccess(partitionable_loop* curr_loop, bool is_lhs, VariantT curr_update)
{
  if( is_lhs ){
    assert( read_loops.size() == 0 || read_loops.back() != curr_loop );
    if( write_loops.size() != 0 && write_loops.back().first == curr_loop ){
      if( write_loops.back().second != curr_update ){
	printf("Scalar %s is updated with multiple operators within loop %d\n",scalar_var->get_name().getString().c_str(),curr_loop->GetLoopNum());
	exit(1);
      }
    }
    else{
      write_loops.push_back(pair<partitionable_loop*,VariantT>(curr_loop,curr_update));
    }
  }
  else{
    assert(write_loops.size() == 0 || write_loops.back().first != curr_loop);
    if( read_loops.size() == 0 || read_loops.back() != curr_loop )
      read_loops.push_back(curr_loop);
  }
}
