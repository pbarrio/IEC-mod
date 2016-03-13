/*
 * global_loop.cpp: This file is part of the IEC project.
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
 * @file: global_loop.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/global_loop.hpp"
#include <cassert>

using namespace std;

global_loop::global_loop(int mn, int nit, int offset):
  my_num(mn),
  num_iters(nit),
  iter_vertex(new vertex * [nit]),
  nproc_local(0){

  for (int i = 0; i < num_iters; i++)
    iter_vertex[i] = new vertex(my_num, i, offset + i);
}

global_loop::~global_loop(){

  for (int i = 0; i < num_iters; i++)
    delete iter_vertex[i];
  delete[] iter_vertex;
}
