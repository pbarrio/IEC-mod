/*
 * conditional.hpp: This file is part of the IEC project.
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
 * @file: conditional.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __CONDITIONAL_HPP__
#define __CONDITIONAL_HPP__

#include "CompileTime/loops.hpp"

#include <rose.h>
#include <assert.h>

class conditional : public branch{

protected:
  
  static int nconditional;

  const int cond_num;

  SgVariableDeclaration* if_counter;

  SgVariableDeclaration* if_array;

  SgExpression* cond_exp;

  SgIfStmt* orig_if;

public:

  virtual ~conditional() { }
  
  conditional(SgIfStmt*);

  inline bool IsLoop() const { return false; }

  inline SgExpression* GetCondnExp() const { return cond_exp; }

  inline SgVariableDeclaration* GetCondnArray() const { return if_array; }

  virtual bool IsThen() const=0;

  virtual SgVariableDeclaration* GetCondnCounter() const=0;

  inline SgVariableDeclaration* GetIfCounter() const { return if_counter; }
  
  void InitCondnArray(SgStatement*,SgFunctionSymbol*);

  void DeleteCondnArray(SgStatement*,SgFunctionSymbol*) const;

  void GenerateIfCounterIncrement(SgStatement*) const;

  void GenerateIfCounterIncrementAfter(SgStatement*) const;

  virtual void GenerateBranchCounterIncrement(SgStatement*) const=0;

  void ResetIfCounter(SgStatement*) const;

  virtual void ResetBranchCounter(SgStatement*) const=0;

  friend class partitionable_loop;
};


class then_branch : public conditional{

private: 
  
  SgVariableDeclaration* then_counter;

public:

  ~then_branch() { }
  
  inline bool IsThen() const { return true; }

  then_branch(SgIfStmt* c) : conditional(c), then_counter(NULL) {} 

  SgVariableDeclaration* GetCondnCounter() const { assert(then_counter); return then_counter; }

  void InitIfThenCounter(SgStatement*);

  void GenerateBranchCounterIncrement(SgStatement*) const;

  void ResetBranchCounter(SgStatement*) const;
};


class else_branch : public conditional{

private: 

  SgVariableDeclaration* else_counter;

public:

  ~else_branch() { }

  inline bool IsThen() const { return false; }
  
  else_branch(SgIfStmt* c) : conditional(c), else_counter(NULL) {}

  SgVariableDeclaration* GetCondnCounter() const { assert(else_counter) ; return else_counter; }

  void InitElseCounter(SgStatement*);

  void GenerateBranchCounterIncrement(SgStatement*) const;
  
  void ResetBranchCounter(SgStatement*) const;
};




#endif
