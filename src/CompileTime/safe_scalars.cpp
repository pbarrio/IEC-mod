/*
 * safe_scalars.cpp: This file is part of the IEC project.
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
 * @file: safe_scalars.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/safe_scalars.hpp"
#include "OmpAttribute.h"
#include <sstream>

using namespace std;

int safe_scalar_details::n_safe_scalars = 0;

safe_scalar_details::safe_scalar_details(SgVariableSymbol* ov):
  my_num(n_safe_scalars++),
  shadow_scalar(NULL),
  orig_var(ov)
{ 
  curr_counter = NULL;
  curr_type = OA_DESCRIPTOR;
  last_defined = NULL;
}


void safe_scalar_details::InitShadowScalars(SgStatement* insert_before)
{
  assert(shadow_scalar == NULL );
  ostringstream shadow_scalar_name_stream;
  //string shadow_scalar_string = scalar_name;
  shadow_scalar_name_stream << "__shadow_" << my_num << "__" ;
  string shadow_scalar_string = shadow_scalar_name_stream.str();
  SgAssignInitializer* init_rhs = SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0));
  shadow_scalar = SageBuilder::buildVariableDeclaration(shadow_scalar_string,SageBuilder::buildIntType(),init_rhs);
  SageInterface::insertStatementBefore(insert_before,shadow_scalar);
}


void safe_scalar_details::SetSafeScalarValues(SgIfStmt* curr_if_stmt) const
{
  assert(shadow_scalar);
  SgStatement* then_stmt = curr_if_stmt->get_true_body();
  SgAssignOp* set_known_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(shadow_scalar),SageBuilder::buildIntVal(1));
  SgExprStatement* set_known_stmt = SageBuilder::buildExprStatement(set_known_op);
  SageInterface::insertStatementAfter(then_stmt,set_known_stmt);

  SgStatement* else_stmt = curr_if_stmt->get_false_body();
  SgAssignOp* set_unknown_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(shadow_scalar),SageBuilder::buildIntVal(0));
  SgExprStatement* set_unknown_stmt = SageBuilder::buildExprStatement(set_unknown_op);
  SageInterface::insertStatementAfter(else_stmt,set_unknown_stmt);
}

void safe_scalar_details::SetScalarUnknown(SgStatement* insert_before) const
{
  assert(shadow_scalar);
  SgAssignOp* set_unknown_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(shadow_scalar),SageBuilder::buildIntVal(0));
  SgExprStatement* set_unknown_stmt = SageBuilder::buildExprStatement(set_unknown_op);
  SageInterface::insertStatementBefore(insert_before,set_unknown_stmt);
}

void safe_scalar_details::AddPrivateVariables(OmpSupport::OmpAttribute* curr_attribute) const {
  // string new_name = scalar_name;
  // curr_attribute->addVariable(OmpSupport::e_firstprivate,scalar_name);
}
