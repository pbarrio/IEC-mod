/*
 * kshash.hpp: This file is part of the IEC project.
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
 * @file: kshash.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 * @acknowledgements : Kevin Stock <stockk@cse.ohio-state.edu>
 */
#ifndef __KS_HASH_HPP__
#define __KS_HASH_HPP__
#include <cstdio>
#include <cstdlib>


class kshash{
 private:
  int n_elems;

  const int hash_size;

  int* const hash_table;
  
  inline int hash(int key){
    return key % hash_size;
  }
  
 public:
  
  kshash(int);
  
  ~kshash();

  inline int num_elems() const{ return n_elems; };

  inline void hash_insert(int target_key, int target_value){
    int hash_pos = hash(target_key);
    if ( hash_table[2*hash_pos] == target_key )
      hash_table[2*hash_pos+1] += target_value;
    else{
      int counter = 0;
    
      while( hash_table[2*hash_pos] != -1 && hash_table[2*hash_pos] != target_key && counter < hash_size ){
	hash_pos = (hash_pos + 1)%hash_size;
	counter++;
      }
      
      if( counter == hash_size )
	exit(1);

      if( hash_table[2*hash_pos] == -1 ){
	n_elems++;
	hash_table[2*hash_pos] = target_key;
	hash_table[2*hash_pos+1] = target_value;
      }
      else{
	hash_table[2*hash_pos+1] += target_value;
      }
    }
  }

  void return_keys_vals(int*, int*);

  void return_keys_vals(int*);
};  
  
  



#endif
