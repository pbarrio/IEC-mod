/*
 * qphash.hpp: This file is part of the IEC project.
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
 * @file: qphash.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 * @acknowledgements : Qingpeng Niu <niuq@cse.ohio-state.edu>
 */
#ifndef __QP_HASH_HPP__
#define __QP_HASH_HPP__
#include <cstdio>

struct tt {
  struct tt* next;
  unsigned int key;
  int value;
  tt(unsigned int,int,tt*);
  ~tt();
};


class qphash{
 private:
  int n_elems;

  /* const float A; */

  const unsigned short hash_size;

  /* const unsigned int log_hash_size; */
  
  struct tt** const hash_table;
  
  inline unsigned short hash(unsigned int key) const{
    // int s = floor(A*32);
    // int x = key*s;
    // return x >> (32-log_hash_size);
    // key = (key+0x7ed55d16) + (key<<12);
    // key = (key^0xc761c23c) ^ (key>>19);
    // key = (key+0x165667b1) + (key<<5);
    // key = (key+0xd3a2646c) ^ (key<<9);
    // key = (key+0xfd7046c5) + (key<<3);
    // key = (key^0xb55a4f09) ^ (key>>16);
    //return a;
    return key % hash_size;
  }

  inline unsigned short hash_check(unsigned int target_key,int target_value){
    unsigned short pos = hash(target_key);
    struct tt* iter;
    for( iter = hash_table[pos] ; iter != NULL ; iter = iter->next)
      if( iter->key == target_key ){
	iter->value += target_value;
	return hash_size;
      }
  
    return pos;
  }
    
  
 public:
  
  qphash(unsigned int);
  
  ~qphash();

  inline int num_elems() const{ return n_elems; };
  
  int return_key_val(unsigned int);

  void hash_insert_only(unsigned int,int);

  inline void hash_insert(unsigned int target_key ,int target_value){
    unsigned short pos = hash_check(target_key,target_value);
    if( pos != hash_size ){
      //printf("Didnt Find\n");
      struct tt* new_elem = new struct tt(target_key,target_value,hash_table[pos]);
      hash_table[pos] = new_elem;
      n_elems++;
    }
  }

  void return_keys_vals(int*, int*);

  void return_keys_vals(int*);

  double metric() const;

};  
  
  



#endif
