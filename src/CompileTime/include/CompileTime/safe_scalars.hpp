/*
 * safe_scalars.hpp: This file is part of the IEC project.
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
 * @file: safe_scalars.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __SAFE_SCALARS_HPP__
#define __SAFE_SCALARS_HPP__

#include <rose.h>
#include <string.h>
#include <cassert>
#include "CompileTime/loops.hpp"
#include "CompileTime/access_details.hpp"

namespace OmpSupport{
  class OmpAttribute;
};


struct last_defined_attribute : public AstAttribute{
  SgStatement* last_defined;
  last_defined_attribute(SgStatement* a) { last_defined = a; }
  last_defined_attribute() { last_defined = NULL; }
};
  

class safe_scalar_details : public AstAttribute{

 private:
  static int n_safe_scalars;

  const int my_num;

  SgVariableDeclaration* shadow_scalar;
  
  branch* curr_counter;

  descriptor_type curr_type;

  SgStatement* last_defined;

  SgVariableSymbol* orig_var;

 public:
  
  safe_scalar_details(SgVariableSymbol*);

  ~safe_scalar_details() { }

  void InitShadowScalars(SgStatement*);

  void SetSafeScalarValues(SgIfStmt*) const;

  SgVariableDeclaration* GetShadowScalar() const { assert(shadow_scalar); return shadow_scalar; }

  void SetScalarUnknown(SgStatement*) const;

  void AddPrivateVariables(OmpSupport::OmpAttribute*) const; 
  
  inline SgVariableSymbol* getOriginalVar() const { return orig_var; }

  friend class partitionable_loop;
};


#endif
