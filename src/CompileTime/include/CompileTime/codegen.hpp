/*
 * codegen.hpp: This file is part of the IEC project.
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
 * @file: codegen.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __CODEGEN_HPP__
#define __CODEGEN_HPP__

#include "CompileTime/opts.hpp"
#include "CompileTime/array_details.hpp"
#include "CompileTime/loops.hpp"
#include "CompileTime/access_details.hpp"
#include "CompileTime/indarray_details.hpp"
#include "CompileTime/safe_scalars.hpp"
#include "CompileTime/data_scalars.hpp"
#include <rose.h>
#include <Cxx_Grammar.h>
#include <deque>
#include <string>
#include <map>

class driver;

class codegen{

 private:
  
  SgProject* project;

  std::deque<array_details*> data_arrays;

  std::deque< safe_scalar_details* > safe_scalars;

  std::deque< partitionable_loop* > all_loops;

  int num_partitionable_loops;

  int num_partitionable_groups;

  std::deque< indarray_details* > indirection_arrays;
  
  std::deque<SgVariableSymbol*> iterators;
  
  std::deque<SgVariableSymbol*> parameters;

  std::deque<data_scalar*> data_scalars;

  std::deque<SgVariableSymbol*> outer_vars;

  std::deque<SgVariableSymbol*> privatizable;

  SgTreeCopy copy_helper;

  void GenerateCCode(SgGlobal*);

  void GenerateDataSizeArray(SgStatement*);

  void GetCommunicationInfo(SgStatement*);
  
  static SgFunctionSymbol* is_known_fn;

  static SgFunctionSymbol* get_elem_fn;

  SgVariableDeclaration* data_num_count;

  int num_data_arrays;

  SgVariableDeclaration* ro_mask;

  SgName proc_id_name;

  SgName nprocs_name;

  SgName nthreads_name;

  SgName myid_name;

  SgName thread_id_name;

  void InsertOMPParallel(SgBasicBlock*);

  SgStatement* root_node;

 public:
  
  codegen(SgStatement*,std::deque<SgForStatement*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,const std::deque<array_string*>& );

  void InitDataArrays(std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,const std::deque<array_string*>&);

  void InitIndirectionArrays(std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,const std::deque<array_string*>&);

  void AnalyseDataAccessExpressions();

  ~codegen();

  friend class driver;
};

#endif
