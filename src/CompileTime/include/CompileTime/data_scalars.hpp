/*
 * data_scalars.hpp: This file is part of the IEC project.
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
 * @file: data_scalars.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __DATA_SCALARS_HPP__
#define __DATA_SCALARS_HPP__

#include <rose.h>
#include <deque>

class partitionable_loop;

class data_scalar : public AstAttribute{

private:

  SgVariableSymbol* scalar_var;
  
  std::deque<partitionable_loop*> read_loops;
  
  std::deque< std::pair<partitionable_loop*,VariantT> > write_loops;

public:

  data_scalar(SgVariableSymbol* cv) : scalar_var(cv) { }
  
  ~data_scalar() {}

  void AddAccess(partitionable_loop*, bool, VariantT);

  inline SgVariableSymbol* getOriginalVar() const { return scalar_var; }

};




#endif
