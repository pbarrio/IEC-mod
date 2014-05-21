/*
 * aux.cpp: This file is part of the IEC project.
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
 * @file: aux.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include <rose.h>
#include <string.h>
#include "CompileTime/loops.hpp"
#include "CompileTime/access_details.hpp"
#include "CompileTime/safe_scalars.hpp"
#include <map>

using namespace std;

bool CompareExpression(const SgExpression* lhs, const SgExpression* rhs)
{
  bool isEqual = false;
  //cout << "Comparing " << lhs->unparseToString() << " and " << rhs->unparseToString() << endl;
  
  if( ( isSgBinaryOp(lhs) ) && ( isSgBinaryOp(rhs) ) ) {
    const SgBinaryOp* lhs_bin_op = isSgBinaryOp(lhs);
    const SgBinaryOp* rhs_bin_op = isSgBinaryOp(rhs);
    if( lhs_bin_op->variantT() == rhs_bin_op->variantT() )
      isEqual = CompareExpression(lhs_bin_op->get_lhs_operand(),rhs_bin_op->get_lhs_operand() ) && CompareExpression(lhs_bin_op->get_rhs_operand(),rhs_bin_op->get_rhs_operand() );
  }
  else if( ( isSgUnaryOp(lhs) ) && ( isSgUnaryOp(rhs) ) ) {
    const SgUnaryOp* lhs_unary_op = isSgUnaryOp(lhs);
    const SgUnaryOp* rhs_unary_op = isSgUnaryOp(rhs);
    if( lhs_unary_op->variantT() == rhs_unary_op->variantT() )
      isEqual = CompareExpression(lhs_unary_op->get_operand(),rhs_unary_op->get_operand());
  }
  else if( ( isSgValueExp(lhs) ) && ( isSgValueExp(rhs) ) ) {
    const SgValueExp* lhs_value = isSgValueExp(lhs);
    const SgValueExp* rhs_value = isSgValueExp(rhs);
    VariantT lhs_v_type = lhs_value->variantT();
    VariantT rhs_v_type = rhs_value->variantT();
    if( lhs_v_type == rhs_v_type )
      switch(lhs_v_type){
      case V_SgIntVal:
	if( isSgIntVal(lhs_value)->get_value() == isSgIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgLongIntVal:
	if( isSgLongIntVal(lhs_value)->get_value() == isSgLongIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgLongLongIntVal:
	if( isSgLongLongIntVal(lhs_value)->get_value() == isSgLongLongIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgShortVal:
	if( isSgShortVal(lhs_value)->get_value() == isSgShortVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedIntVal:
	if( isSgUnsignedIntVal(lhs_value)->get_value() == isSgUnsignedIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedLongLongIntVal:
	if( isSgUnsignedLongLongIntVal(lhs_value)->get_value() == isSgUnsignedLongLongIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedLongVal:
	if( isSgUnsignedLongVal(lhs_value)->get_value() == isSgUnsignedLongVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedShortVal:
	if( isSgUnsignedShortVal(lhs_value)->get_value() == isSgUnsignedShortVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      default:
	assert(0);
      }
  }
  else if( ( isSgVarRefExp(lhs) ) && ( isSgVarRefExp(rhs) ) ) {
    const SgVarRefExp* lhs_var = isSgVarRefExp(lhs);
    const SgVarRefExp* rhs_var = isSgVarRefExp(rhs);
    string lhs_string = lhs_var->get_symbol()->get_name().getString(), rhs_string = rhs_var->get_symbol()->get_name().getString();
    if( lhs_string.compare(rhs_string) == 0 )
      isEqual = true;
  }
  return isEqual;
}

  

const SgVariableSymbol* GetLoopIterator(const SgForInitStatement* for_init)
{
  const SgStatementPtrList& temp2 = for_init->get_init_stmt();
  const SgExprStatement* temp3 = isSgExprStatement(*(temp2.begin()));
  assert(temp3);
  const SgExpression* temp4 = temp3->get_expression();
  const SgExpression* temp5 = isSgBinaryOp(temp4)->get_lhs_operand();
  assert(temp5);
  const SgVarRefExp* iter_ref = isSgVarRefExp(temp5);
  assert(iter_ref);
  const SgVariableSymbol* loop_iterator = iter_ref->get_symbol();
  return loop_iterator;
}

SgVariableSymbol* GetLoopIterator( SgForInitStatement* for_init)
{
  SgStatementPtrList& temp2 = for_init->get_init_stmt();
  SgExprStatement* temp3 = isSgExprStatement(*(temp2.begin()));
  assert(temp3);
  SgExpression* temp4 = temp3->get_expression();
  SgExpression* temp5 = isSgBinaryOp(temp4)->get_lhs_operand();
  assert(temp5);
  SgVarRefExp* iter_ref = isSgVarRefExp(temp5);
  assert(iter_ref);
  SgVariableSymbol* loop_iterator = iter_ref->get_symbol();
  return loop_iterator;
}


bool isNumber( const SgExpression* exp )
{
  const string lb_string = exp->unparseToString();
  string::const_iterator it = lb_string.begin() ;
  for(  ; it != lb_string.end() ; it++ )
    if( !isdigit(*it) && !isspace(*it) )
      break;
  return it == lb_string.end() ;
}




bool CompareAccessExpression(const SgExpression* lhs, const SgExpression* rhs)
{
  bool isEqual = false;

  if( ( isSgBinaryOp(lhs) ) && ( isSgBinaryOp(rhs) ) ) {
    const SgBinaryOp* lhs_bin_op = isSgBinaryOp(lhs);
    const SgBinaryOp* rhs_bin_op = isSgBinaryOp(rhs);
    if( lhs_bin_op->variantT() == rhs_bin_op->variantT() )
      isEqual = CompareAccessExpression(lhs_bin_op->get_lhs_operand(),rhs_bin_op->get_lhs_operand() ) && CompareAccessExpression(lhs_bin_op->get_rhs_operand(),rhs_bin_op->get_rhs_operand() );
  }
  else if( ( isSgUnaryOp(lhs) ) && ( isSgUnaryOp(rhs) ) ) {
    const SgUnaryOp* lhs_unary_op = isSgUnaryOp(lhs);
    const SgUnaryOp* rhs_unary_op = isSgUnaryOp(rhs);
    if( lhs_unary_op->variantT() == rhs_unary_op->variantT() )
      isEqual = CompareAccessExpression(lhs_unary_op->get_operand(),rhs_unary_op->get_operand());
  }
  else if( ( isSgValueExp(lhs) ) && ( isSgValueExp(rhs) ) ) {
    const SgValueExp* lhs_value = isSgValueExp(lhs);
    const SgValueExp* rhs_value = isSgValueExp(rhs);
    VariantT lhs_v_type = lhs_value->variantT();
    VariantT rhs_v_type = rhs_value->variantT();
    if( lhs_v_type == rhs_v_type )
      switch(lhs_v_type){
      case V_SgIntVal:
	if( isSgIntVal(lhs_value)->get_value() == isSgIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgLongIntVal:
	if( isSgLongIntVal(lhs_value)->get_value() == isSgLongIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgLongLongIntVal:
	if( isSgLongLongIntVal(lhs_value)->get_value() == isSgLongLongIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgShortVal:
	if( isSgShortVal(lhs_value)->get_value() == isSgShortVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedIntVal:
	if( isSgUnsignedIntVal(lhs_value)->get_value() == isSgUnsignedIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedLongLongIntVal:
	if( isSgUnsignedLongLongIntVal(lhs_value)->get_value() == isSgUnsignedLongLongIntVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedLongVal:
	if( isSgUnsignedLongVal(lhs_value)->get_value() == isSgUnsignedLongVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      case V_SgUnsignedShortVal:
	if( isSgUnsignedShortVal(lhs_value)->get_value() == isSgUnsignedShortVal(rhs_value)->get_value() )
	  isEqual = true;
	break;
      default:
	assert(0);
      }
  }
  else if( ( isSgVarRefExp(lhs) ) && ( isSgVarRefExp(rhs) ) ) {
    const SgVarRefExp* lhs_var = isSgVarRefExp(lhs);
    const SgVarRefExp* rhs_var = isSgVarRefExp(rhs);

    const SgVariableSymbol* lhs_symbol = lhs_var->get_symbol();
    const SgVariableSymbol* rhs_symbol = rhs_var->get_symbol();

    if( lhs_symbol->get_name().getString().compare(rhs_symbol->get_name().getString()) == 0 ){
      const last_defined_attribute* lhs_last = static_cast<const last_defined_attribute*>(lhs_var->getAttribute("last_defined"));
      const last_defined_attribute* rhs_last = static_cast<const last_defined_attribute*>(rhs_var->getAttribute("last_defined"));

      if( lhs_last->last_defined == rhs_last->last_defined )
	isEqual = true;
    }
  }
  return isEqual;
}
