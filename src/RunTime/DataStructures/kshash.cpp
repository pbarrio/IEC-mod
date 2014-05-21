/*
 * kshash.cpp: This file is part of the IEC project.
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
 * @file: kshash.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 * @acknowledgements : Kevin Stock <stockk@cse.ohio-state.edu>
 */
#include "DataStructures/kshash.hpp"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cassert>


using namespace std;

kshash::kshash(int size):
  hash_size(197),
  hash_table(new int[197*2])
{
  for( int i = 0 ; i < hash_size ; i++ ){
    hash_table[2*i] = -1;
    hash_table[2*i+1] = 0;
  }
  n_elems = 0;
}

kshash::~kshash()
{
  delete[] hash_table;
}

void kshash::return_keys_vals(int* all_keys, int* all_vals)
{
  int counter =0;
  for( int i = 0; i < hash_size ; i++ )
    if( hash_table[2*i] != -1 ){
      all_keys[counter] = hash_table[2*i];
      all_vals[counter] = hash_table[2*i+1];
      counter++;
    }
}


void kshash::return_keys_vals(int* all_keys_vals)
{
  int counter =0;
  for( int i = 0; i < hash_size ; i++ )
    if( hash_table[2*i] != -1 ){
      all_keys_vals[counter++] = hash_table[2*i];
      all_keys_vals[counter++] = hash_table[2*i+1];
    }
}
