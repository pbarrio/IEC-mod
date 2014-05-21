/*
 * indarray_details.hpp: This file is part of the IEC project.
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
 * @file: indarray_details.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __INDARRAY_DETAILS_HPP__
#define __INDARRAY_DETAILS_HPP__

#include <string>
#include <rose.h>
#include "CompileTime/array_details.hpp"

class indarray_details: public AstAttribute{
  
 private:
  
  const std::string indarray_name;
  
  const std::string indarray_size;

  const std::string indarray_stride;

  const SgVariableSymbol* indarray_symbol;

  const int my_num;

  static int num_indarrays;

 public:
  
  indarray_details(SgVariableSymbol*,array_string*);

  ~indarray_details() { };

  const SgVariableSymbol* GetArraySymbol() const { return indarray_symbol; }

  int GetIndarrayNum() const { return my_num; }

  void GenerateSetAccessParam(SgStatement*) const;

  static SgFunctionSymbol* set_access_param_fn;

};


#endif
