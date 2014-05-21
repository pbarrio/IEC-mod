/*
 * conditional.cpp: This file is part of the IEC project.
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
 * @file: conditional.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/conditional.hpp"
#include <sstream>

using namespace std;

int conditional::nconditional = 0;

conditional::conditional(SgIfStmt* c) : orig_if(c), cond_num(nconditional++), if_counter(NULL), if_array(NULL), branch(true)
{
  cond_exp = isSgExprStatement(orig_if->get_conditional())->get_expression();
}


void conditional::InitCondnArray(SgStatement* insert_before, SgFunctionSymbol* malloc_fn)
{
  ostringstream if_array_oss;
  if_array_oss << "__if_cond_" << cond_num << "__" ;
  SgName if_array_name(if_array_oss.str());
  assert(if_counter);
  SgMultiplyOp* if_array_size = SageBuilder::buildBinaryExpression<SgMultiplyOp>(SageBuilder::buildSizeOfOp(SageBuilder::buildIntType()),SageBuilder::buildVarRefExp(if_counter));
  SgFunctionCallExp* if_array_malloc = SageBuilder::buildFunctionCallExp(malloc_fn,SageBuilder::buildExprListExp(if_array_size));
  SgAssignInitializer* if_array_init = SageBuilder::buildAssignInitializer(SageBuilder::buildCastExp(if_array_malloc,SageBuilder::buildPointerType(SageBuilder::buildIntType())));
  if_array = SageBuilder::buildVariableDeclaration(if_array_name,SageBuilder::buildPointerType(SageBuilder::buildIntType()),if_array_init);
  SageInterface::insertStatementBefore(insert_before,if_array);
}

void conditional::DeleteCondnArray(SgStatement* insert_before, SgFunctionSymbol* free_fn) const
{
  assert(if_array);
  SgFunctionCallExp* free_exp = SageBuilder::buildFunctionCallExp(free_fn,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(if_array)));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(free_exp));
}

void conditional::GenerateIfCounterIncrement(SgStatement* insert_before) const
{
  assert(if_counter);
  SgAssignOp* inc_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(if_counter),SageBuilder::buildBinaryExpression<SgAddOp>(SageBuilder::buildVarRefExp(if_counter),SageBuilder::buildIntVal(1)));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(inc_exp));
}

void conditional::GenerateIfCounterIncrementAfter(SgStatement* insert_after) const
{
  assert(if_counter);
  SgAssignOp* inc_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(if_counter),SageBuilder::buildBinaryExpression<SgAddOp>(SageBuilder::buildVarRefExp(if_counter),SageBuilder::buildIntVal(1)));
  SageInterface::insertStatementAfter(insert_after,SageBuilder::buildExprStatement(inc_exp));
}


void conditional::ResetIfCounter(SgStatement* insert_before) const
{
  assert(if_counter);
  SgAssignOp* reset_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(if_counter),SageBuilder::buildIntVal(0));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(reset_exp));
}


void then_branch::InitIfThenCounter(SgStatement* insert_before)
{
  assert(if_counter == NULL);
  ostringstream if_oss;
  if_oss << "__if_" << cond_num << "__" ;
  SgName if_counter_name(if_oss.str());
  if_counter = SageBuilder::buildVariableDeclaration(if_counter_name,SageBuilder::buildIntType(),SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0)));
  SageInterface::insertStatementBefore(insert_before,if_counter);

  assert(then_counter == NULL);
  ostringstream then_oss;
  then_oss << "__then_" << cond_num << "__";
  SgName then_counter_name(then_oss.str());
  then_counter = SageBuilder::buildVariableDeclaration(then_counter_name,SageBuilder::buildIntType(),SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0)));
  SageInterface::insertStatementBefore(insert_before,then_counter);
}


void then_branch::GenerateBranchCounterIncrement(SgStatement* insert_before) const
{
  assert(then_counter);
  SgAssignOp* inc_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(then_counter),SageBuilder::buildBinaryExpression<SgAddOp>(SageBuilder::buildVarRefExp(then_counter),SageBuilder::buildIntVal(1)));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(inc_exp));
}

void then_branch::ResetBranchCounter(SgStatement* insert_before) const
{
  assert(then_counter);
  SgAssignOp* reset_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(then_counter),SageBuilder::buildIntVal(0));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(reset_exp));
}


void else_branch::InitElseCounter(SgStatement* insert_before)
{
  assert(else_counter == NULL);
  ostringstream else_oss;
  else_oss << "__else_" << cond_num << "__";
  SgName else_counter_name(else_oss.str());
  else_counter = SageBuilder::buildVariableDeclaration(else_counter_name,SageBuilder::buildIntType(),SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0)));
  SageInterface::insertStatementBefore(insert_before,else_counter);
}

void else_branch::GenerateBranchCounterIncrement(SgStatement* insert_before) const
{
  assert(else_counter);
  SgAssignOp* inc_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(else_counter),SageBuilder::buildBinaryExpression<SgAddOp>(SageBuilder::buildVarRefExp(else_counter),SageBuilder::buildIntVal(1)));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(inc_exp));
}

void else_branch::ResetBranchCounter(SgStatement* insert_before) const
{
  assert(else_counter);
  SgAssignOp* reset_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(else_counter),SageBuilder::buildIntVal(0));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(reset_exp));
}

