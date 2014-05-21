/*
 * array_details.hpp: This file is part of the IEC project.
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
 * @file: array_details.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __ARRAY_DETAILS_HPP__
#define __ARRAY_DETAILS_HPP__

#include <iostream>
#include <string>
#include <rose.h>
#include <deque>
#include "CompileTime/loops.hpp"

class access_details;

/* enum OP_TYPE{ */
/*   PLUS_EQ, */
/*   ASSIGN_EQ, */
/*   STAR_EQ */
/* } */

struct array_string{
  std::string array_name;
  std::string array_size;
  std::string array_stride;
  array_string();
  array_string(std::string);
};

class array_details: public AstAttribute{

 private:
  
  const SgVariableSymbol* array_symbol;

  const int my_num;

  const std::string array_name;

  const std::string stride_string;

  const std::string size_string;

  static int num_arrays;

  SgVariableDeclaration* local_array;

  SgVariableDeclaration* local_array_size;

  std::deque<access_details*> all_accesses;

  std::deque<partitionable_loop*> read_loops;

  std::deque< std::pair<partitionable_loop*,VariantT> > write_loops;

  std::deque<partitionable_loop*> local_writes;

  void CheckConflicts(access_details*);

  access_details* last_access;

 public:
  
  array_details(const SgVariableSymbol*,array_string*);

  ~array_details();

  inline const std::string GetArraySize() const { 
    const std::string new_string = size_string; 
    return new_string; 
  }

  inline const SgVariableSymbol* GetArraySymbol() { return array_symbol; }

  inline int GetArrayNum() const { return my_num; }

  inline int GetROMask() const{
    if(write_loops.size() == 0 && local_writes.size() == 0 )
      return 1;
    else
      return 0;
  }

  void GenerateStrideInfo(SgStatement*);

  access_details* AddAccess(SgExpression*, branch*, descriptor_type, bool, partitionable_loop*, VariantT);

  void GetCommunicationInfo(SgStatement*);

  void InitAccessArrays(SgStatement*);

  void DeleteArrays(SgStatement*);

  void InitLocalArrays(SgStatement*);

  inline SgVarRefExp* GetLocalArrayRef() const { return SageBuilder::buildVarRefExp(local_array) ; }

  static SgFunctionSymbol* malloc_fn;
  static SgFunctionSymbol* malloc_2d_double_fn;
  static SgFunctionSymbol* malloc_2d_float_fn;
  static SgFunctionSymbol* free_fn;
  static SgFunctionSymbol* free_2d_double_fn;
  static SgFunctionSymbol* free_2d_float_fn;
  static SgFunctionSymbol* get_size_fn;
  static SgFunctionSymbol* populate_fn;
  static SgFunctionSymbol* set_array_stride_fn;
  static SgFunctionSymbol* communicate_reads_for_fn;
  static SgFunctionSymbol* communicate_writes_for_fn;
};



#endif
