/*
 * indarray_details.cpp: This file is part of the IEC project.
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
 * @file: indarray_details.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/indarray_details.hpp"
#include <boost/regex.hpp>

using namespace std;

int indarray_details::num_indarrays = 0;
SgFunctionSymbol* indarray_details::set_access_param_fn = NULL;

indarray_details::indarray_details(SgVariableSymbol* sym, array_string* st):
  indarray_symbol(sym),
  indarray_name(st->array_name),
  indarray_size(st->array_size),
  indarray_stride(st->array_stride),
  my_num(num_indarrays)
{
  num_indarrays++;
}


void indarray_details::GenerateSetAccessParam(SgStatement* insert_before) const
{
  assert(set_access_param_fn);
  SgIntVal* arg1 = SageBuilder::buildIntVal(my_num);
  SgVarRefExp* arg2 = SageBuilder::buildVarRefExp(indarray_size);
  SgName proc_id_var_name("__proc_id__");
  SgName nprocs_var_name("__nprocs__");
  SgVarRefExp* arg3 = SageBuilder::buildVarRefExp(indarray_stride);
  string indarray_size_var_string = indarray_name;
  indarray_size_var_string.insert(0,"__size_").append("__");
  SgName indarray_size_var_name(indarray_size_var_string);
  SgVariableDeclaration* indarray_size_var = SageBuilder::buildVariableDeclaration(indarray_size_var_name,SageBuilder::buildIntType(),SageBuilder::buildAssignInitializer(SageBuilder::buildVarRefExp(indarray_size)));
  SageInterface::insertStatementBefore(insert_before,indarray_size_var);
  SgDivideOp* arg4_mult_rhs = SageBuilder::buildBinaryExpression<SgDivideOp>(SageBuilder::buildVarRefExp(indarray_size_var),SageBuilder::buildVarRefExp(nprocs_var_name));
  SgMultiplyOp* arg4_add_rhs = SageBuilder::buildBinaryExpression<SgMultiplyOp>(arg4_mult_rhs,SageBuilder::buildVarRefExp(proc_id_var_name));
  SgAddOp* arg4 = SageBuilder::buildBinaryExpression<SgAddOp>(arg4_add_rhs,SageBuilder::buildVarRefExp(indarray_name));
  SgExprListExp* args = SageBuilder::buildExprListExp(arg1,arg2,arg3,arg4);
  SgFunctionCallExp* init_fn_call = SageBuilder::buildFunctionCallExp(set_access_param_fn,args);
  SgExprStatement* init_fn_stmt = SageBuilder::buildExprStatement(init_fn_call);
  SageInterface::insertStatementBefore(insert_before,init_fn_stmt);
}
