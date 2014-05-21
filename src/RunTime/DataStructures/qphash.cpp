/*
 * qphash.cpp: This file is part of the IEC project.
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
 * @file: qphash.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 * @acknowledgements : Qingpeng Niu <niuq@cse.ohio-state.edu>
 */
#include "DataStructures/qphash.hpp"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cassert>


using namespace std;

tt::tt(unsigned int k , int v, struct tt* nx)
{
  key = k;
  value = v;
  next = nx;
};


tt::~tt()
{
  if( next != NULL )
    delete next;
}

qphash::qphash(unsigned int size):
  hash_size(size),
  // log_hash_size(6),
  // A(0.5*sqrt(5.0)-1),
  hash_table(new struct tt*[size])
{
  for( int i = 0 ; i < hash_size ; i++ )
    hash_table[i] = NULL;
  n_elems = 0;
}

qphash::~qphash()
{
  for( int i = 0; i < hash_size ; i++ )
    if( hash_table[i] != NULL )
      delete hash_table[i];
  delete [] hash_table;
}


double qphash::metric() const
{
  int conflicts = 0;
  struct tt* iter;
  for( unsigned short i = 0 ; i < hash_size ; i++ )
    for( iter = hash_table[i] ; iter != NULL ; iter = iter->next )
      conflicts++;
  return (double)conflicts / (double)hash_size;
}


int qphash::return_key_val(unsigned int target_key)
{
  unsigned short pos = hash(target_key);
  struct tt* iter;
  for( iter = hash_table[pos] ; iter != NULL ; iter = iter->next)
    if( iter->key == target_key ){
      return iter->value;
    }
  return -1;
}

void qphash::hash_insert_only(unsigned int target_key, int target_value)
{
  assert( hash_check(target_key,target_value) == hash_size) ;
  unsigned int pos = hash(target_key);
  struct tt* new_elem = new struct tt(target_key,target_value,hash_table[pos]);
  hash_table[pos] = new_elem;
  n_elems++;
}


void qphash::return_keys_vals(int* all_keys, int* all_vals)
{
  int counter =0;
  for( unsigned short i = 0; i < hash_size ; i++ )
    if( hash_table[i] != NULL )
      for( struct tt* iter = hash_table[i] ; iter != NULL ; iter = iter->next ){
	all_keys[counter] = iter->key;
	all_vals[counter] = iter->value;
	counter++;
      }
}


void qphash::return_keys_vals(int* all_keys_vals)
{
  int counter =0;
  for( unsigned short i = 0; i < hash_size ; i++ )
    if( hash_table[i] != NULL )
      for( struct tt* iter = hash_table[i] ; iter != NULL ; iter = iter->next ){
	all_keys_vals[counter++] = iter->key;
	all_keys_vals[counter++] = iter->value;
      }
}
