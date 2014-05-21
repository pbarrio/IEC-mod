/*
 * global_loop.hpp: This file is part of the IEC project.
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
 * @file: global_loop.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __GLOBAL_LOOP_HPP__
#define __GLOBAL_LOOP_HPP__

#include <map>
#include "RunTime/hypergraph.hpp"

class inspector;

class global_loop{
 private:

  vertex** iter_vertex;

  const int my_num;

  const int num_iters;

  int nproc_local;
  //int* home_info;

 public:
  
  global_loop(int,int,int);

  ~global_loop();
  
  inline int GetVertexHome(int iter_value) const{
    return iter_vertex[iter_value]->home;
  }

  friend class inspector;
  
};

#endif
