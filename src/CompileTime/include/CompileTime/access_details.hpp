/*
 * access_details.hpp: This file is part of the IEC project.
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
 * @file: access_details.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __ACCESS_DETAILS_HPP__
#define __ACCESS_DETAILS_HPP__

#include "CompileTime/loops.hpp"
#include "CompileTime/conditional.hpp"
#include <rose.h>
#include <string>
#include <cassert>

class array_details;


class access_details : public AstAttribute{

 private:
  
  SgExpression* access_exp;

  descriptor_type access_type;

  const branch* my_branch;

  SgTreeCopy copy_helper;

  SgVariableDeclaration* access_array;

  SgVariableDeclaration* offset_var;
  
  static int total_num;

  const int my_num;

  bool inspector_gen;

  SgExpression* executor_expr;

 public:

  static SgFunctionSymbol* add_pin_fn;

  static SgFunctionSymbol* renumber_offset_fn;
  
  static SgFunctionSymbol* renumber_fn;

  static SgFunctionSymbol* renumber_const_offset_fn;
  
  access_details(SgExpression*,const branch*,descriptor_type);

  inline const SgExpression* GetExpression() const { return access_exp; }

  void InitAccessArray(SgStatement*,std::string,SgFunctionSymbol*);

  void DeleteAccessArray(SgStatement*,SgFunctionSymbol*) const ;
  
  friend class array_details;

  std::pair<SgStatement*,SgStatement*> GenerateCArrayStore(const int) ;

  void RenumberAccess(SgStatement*,const int) const ;

  SgExpression* GenerateExecutorExpression();

};


#endif
