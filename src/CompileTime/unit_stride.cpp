/*
 * unit_stride.cpp: This file is part of the IEC project.
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
 * @file: unit_stride.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/safe_scalars.hpp"
#include "CompileTime/array_details.hpp"
#include "CompileTime/data_scalars.hpp"
#include "CompileTime/driver.hpp"

using namespace std;

void partitionable_loop::AnalyseTargetLoop()
{
  for( set<SgVariableSymbol*>::iterator it = indirection_scalars.begin() ; it != indirection_scalars.end() ; it++){
    //assert((*it)->attributeExists("indirection_scalar"));
    safe_scalar_details* curr_scalar_details = static_cast<safe_scalar_details*>((*it)->getAttribute("indirection_scalar"));
    curr_scalar_details->curr_counter = NULL;
    curr_scalar_details->curr_type = CONSTANT;
    curr_scalar_details->last_defined = NULL;
  }
  AnalyseTargetStatement(orig_for->get_loop_body());
}




void partitionable_loop::AnalyseTargetStatement(SgStatement* curr_stmt)
{
  VariantT v_body = curr_stmt->variantT();
  
  switch(v_body){
  case V_SgBasicBlock:
    AnalyseTargetBasicBlock(isSgBasicBlock(curr_stmt));
    break;
  case V_SgForStatement:
    AnalyseTargetFor(isSgForStatement(curr_stmt));
    break;
  case V_SgExprStatement:
    AnalyseTargetExprStatement(isSgExprStatement(curr_stmt));
    break;
  case V_SgIfStmt:
    AnalyseTargetIfStatement(isSgIfStmt(curr_stmt));
    break;
  default:
    ;
  }
}

void partitionable_loop::AnalyseTargetBasicBlock( SgBasicBlock* curr_basic_block)
{
  SgStatementPtrList& stmt_list = curr_basic_block->get_statements();
  SgStatementPtrList::const_iterator stmt_iter;

  for( stmt_iter = stmt_list.begin() ; stmt_iter != stmt_list.end() ; stmt_iter++ ){
    AnalyseTargetStatement(*stmt_iter);
  }
}


void partitionable_loop::AnalyseTargetFor( SgForStatement* curr_for)
{
  AnalyseTargetStatement(curr_for->get_loop_body());
  scope* curr_scope = static_cast<scope*>(curr_for->getAttribute("scope"));
  inner_loop* curr_inner = static_cast<inner_loop*>(curr_for->getAttribute("inner_loop"));

  for( set<SgVariableSymbol*>::iterator it = curr_inner->written_iscalars.begin() ; it != curr_inner->written_iscalars.end() ; it++ ){
    safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>((*it)->getAttribute("indirection_scalar"));
    curr_safe_scalar->curr_counter = curr_scope->getScope().back();
    curr_safe_scalar->curr_type = IA_DESCRIPTOR;
    curr_safe_scalar->last_defined = curr_for;
  }
}


void partitionable_loop::AnalyseTargetIfStatement( SgIfStmt* curr_if )
{
  AnalyseTargetStatement(curr_if->get_true_body());
  if( curr_if->get_false_body() )
    AnalyseTargetStatement(curr_if->get_false_body());
  scope* curr_scope = static_cast<scope*>(curr_if->getAttribute("scope"));
  then_branch* curr_then = static_cast<then_branch*>(curr_if->getAttribute("then_branch"));
  for( set<SgVariableSymbol*>::iterator it = curr_then->written_iscalars.begin() ; it != curr_then->written_iscalars.end() ; it++ ){
    safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>((*it)->getAttribute("indirection_scalar"));
    curr_safe_scalar->curr_counter = curr_scope->getScope().back();
    curr_safe_scalar->curr_type = IA_DESCRIPTOR;
    curr_safe_scalar->last_defined = curr_if;
  }
  if( curr_if->get_false_body() ){
    else_branch* curr_else = static_cast<else_branch*>(curr_if->getAttribute("else_branch"));
    for( set<SgVariableSymbol*>::iterator it = curr_then->written_iscalars.begin() ; it != curr_then->written_iscalars.end() ; it++ ){
      safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>((*it)->getAttribute("indirection_scalar"));
      curr_safe_scalar->curr_counter = curr_scope->getScope().back();
      curr_safe_scalar->curr_type = IA_DESCRIPTOR;
      curr_safe_scalar->last_defined = curr_if;
    }
  }
}

void partitionable_loop::AnalyseTargetExprStatement(SgExprStatement* curr_stmt)
{
  IAssignmentAttr* curr_type = static_cast<IAssignmentAttr*>(curr_stmt->getAttribute("assignment_type"));
  scope* curr_scope = static_cast<scope*>(curr_stmt->getAttribute("scope"));
  if( isSgFunctionCallExp(curr_stmt->get_expression() ) ){
    SgExprListExp* curr_args = isSgFunctionCallExp(curr_stmt->get_expression())->get_args();
    SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
    for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ ){
      AnalyseTargetExpression(*curr_arg_iter,curr_scope,false,curr_stmt->get_expression()->variantT());
    }    
  }
  else{
    if( curr_type->IsIAssignment() ){
      SgAssignOp* curr_exp = isSgAssignOp(curr_stmt->get_expression());
      assert(curr_exp);
      SgVarRefExp* lhs_exp = isSgVarRefExp(curr_exp->get_lhs_operand());
      assert(lhs_exp);
      safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>(lhs_exp->get_symbol()->getAttribute("indirection_scalar"));
      pair<branch*,descriptor_type> ret_val = AnalyseSafeExpression(curr_exp->get_rhs_operand(),curr_scope);
      curr_safe_scalar->curr_counter = ret_val.first;
      curr_safe_scalar->curr_type = ret_val.second;
      curr_safe_scalar->last_defined = curr_stmt;
      printf("IAssignment Statement %s : type = %d\n",curr_stmt->unparseToString().c_str(),ret_val.second);
    }
    else {
      SgBinaryOp* curr_assign = isSgBinaryOp(curr_stmt->get_expression());
      VariantT curr_assign_type = curr_assign->variantT();
      AnalyseTargetExpression(curr_assign->get_lhs_operand(),curr_scope,true,curr_assign_type);
      AnalyseTargetExpression(curr_assign->get_rhs_operand(),curr_scope,false,curr_assign_type);
    }
  }
}




void partitionable_loop::AnalyseTargetExpression(SgExpression* curr_exp, scope* curr_scope, bool is_lhs, VariantT assign_type)
{
  if( isSgBinaryOp(curr_exp) ){
    SgBinaryOp* curr_binaryop = isSgBinaryOp(curr_exp);
    
    if( isSgPntrArrRefExp(curr_binaryop) ){
      SgExpression* lhs_exp = curr_binaryop->get_lhs_operand();
      SgExpression* rhs_exp = curr_binaryop->get_rhs_operand();

      while( isSgPntrArrRefExp(lhs_exp) ){
	rhs_exp = isSgPntrArrRefExp(lhs_exp)->get_rhs_operand();
	lhs_exp = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
      }

      SgVarRefExp* curr_array_var = isSgVarRefExp(lhs_exp);
      assert(curr_array_var);
      if( driver::CheckArray(curr_array_var) ){
	pair<branch*,descriptor_type> rhs_ret = AnalyseSafeExpression(rhs_exp,curr_scope);
	array_details* curr_array = static_cast<array_details*>(curr_array_var->get_symbol()->getAttribute("data_array"));
	access_details* curr_access = curr_array->AddAccess(rhs_exp,rhs_ret.first,rhs_ret.second,is_lhs,this,assign_type);
	printf(" to array %s\n",curr_array_var->get_symbol()->get_name().getString().c_str());
	rhs_exp->setAttribute("access_details",curr_access);
      }
    }
    else{
      AnalyseTargetExpression(curr_binaryop->get_lhs_operand(),curr_scope,is_lhs,assign_type);
      AnalyseTargetExpression(curr_binaryop->get_rhs_operand(),curr_scope,is_lhs,assign_type);
    }
  }
  else if( isSgUnaryOp(curr_exp) ){
    AnalyseTargetExpression(isSgUnaryOp(curr_exp)->get_operand(),curr_scope,is_lhs,assign_type);
  }
  else if( isSgFunctionCallExp(curr_exp) ){
    SgExprListExp* curr_args = isSgFunctionCallExp(curr_exp)->get_args();
    SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
    
    for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ ){
      AnalyseTargetExpression(*curr_arg_iter,curr_scope,false,assign_type);
    }
  }
  else if( isSgVarRefExp(curr_exp) ){
    SgVariableSymbol* curr_var = isSgVarRefExp(curr_exp)->get_symbol();
    if(curr_var->attributeExists("data_scalar")){
      data_scalar* curr_data_scalar = static_cast<data_scalar*>(curr_var->getAttribute("data_scalar"));
      curr_data_scalar->AddAccess(this,is_lhs,assign_type);
    }
  }
}



pair<branch*,descriptor_type> partitionable_loop::AnalyseSafeExpression(SgExpression* curr_exp, scope* curr_scope)
{
  pair<branch*,descriptor_type> ret_val;
  ret_val.first = NULL;
  ret_val.second = CONSTANT;
  if( isSgVarRefExp(curr_exp) ){
    SgVariableSymbol* curr_symbol = isSgVarRefExp(curr_exp)->get_symbol();
    deque<branch*> scope_stack = curr_scope->getScope();
    if( curr_symbol->attributeExists("indirection_scalar") ){
      safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>(curr_symbol->getAttribute("indirection_scalar"));
      SgStatement* last_defined = curr_safe_scalar->last_defined;
      if( last_defined ){
	scope* last_defined_scope = static_cast<scope*>(last_defined->getAttribute("scope"));
	assert( last_defined_scope->getScope().size() <= curr_scope->getScope().size() );
	deque<branch*> scope_stack = curr_scope->getScope();
	deque<branch*> last_defined_scope_stack = last_defined_scope->getScope();
	assert ( last_defined_scope_stack.back() == scope_stack[last_defined_scope_stack.size()-1] );
	ret_val.first = curr_safe_scalar->curr_counter;
	ret_val.second = curr_safe_scalar->curr_type;
	int i = scope_stack.size() - 1;
	for( ; i > last_defined_scope_stack.size() - 1 ; i-- )
	  if( !scope_stack[i]->IsLoop() ){
	    ret_val.first = scope_stack[i];
	    ret_val.second = IA_DESCRIPTOR;
	    curr_exp->setAttribute("last_defined",new last_defined_attribute((static_cast<conditional*>(scope_stack[i]))->orig_if));
	    break;
	  }
	if( i == last_defined_scope_stack.size() - 1 )
	  curr_exp->setAttribute("last_defined",new last_defined_attribute(last_defined));
	}
      else
	curr_exp->setAttribute("last_defined",new last_defined_attribute());	
    }
    else{
      deque<branch*>::reverse_iterator it = scope_stack.rbegin();
      //Check if it is an iterator
      for( ; it != scope_stack.rend() ; it++ )
	if( (*it)->IsLoop()) {
	  loop* curr_loop = static_cast<loop*>(*it);
	  if( curr_loop->GetIterator() == curr_symbol ){
	    assert(curr_loop->IsUsed());
	    ret_val.first = curr_loop;
	    ret_val.second = OA_DESCRIPTOR;
	    curr_exp->setAttribute("last_defined",new last_defined_attribute(curr_loop->orig_for));
	    break;
	  }
	}
	else{
	  ret_val.first = (*it);
	  ret_val.second = IA_DESCRIPTOR;
	  curr_exp->setAttribute("last_defined",new last_defined_attribute((static_cast<conditional*>(*it))->orig_if));
	  break;
	}    
    }
  }
  else if( isSgBinaryOp(curr_exp) ){
    SgExpression* curr_lhs = isSgBinaryOp(curr_exp)->get_lhs_operand();
    SgExpression* curr_rhs = isSgBinaryOp(curr_exp)->get_rhs_operand();
    if( isSgPntrArrRefExp(curr_exp) ){
      assert( !isSgPntrArrRefExp(curr_lhs) );
      ret_val = AnalyseSafeExpression(curr_rhs,curr_scope);
      if( ret_val.second == OA_DESCRIPTOR )
	ret_val.second = IA_DESCRIPTOR;
    }
    else{
      pair<branch*,descriptor_type> lhs_ret,rhs_ret;
      lhs_ret = AnalyseSafeExpression(curr_lhs,curr_scope);
      rhs_ret = AnalyseSafeExpression(curr_rhs,curr_scope);
      deque<branch*> scope_stack = curr_scope->getScope();
      deque<branch*>::reverse_iterator it = scope_stack.rbegin();
      for( ; it != scope_stack.rend() ; it++ )
	if( *it == lhs_ret.first || *it == rhs_ret.first )
	  break;
      if( it != scope_stack.rend() )
	if( *it == lhs_ret.first ){
	  ret_val.first = lhs_ret.first;
	  if( lhs_ret.first != rhs_ret.first && lhs_ret.second == OA_DESCRIPTOR && isSgAddOp(curr_exp) )
	    ret_val.second = OA_DESCRIPTOR;
	  else
	    ret_val.second = IA_DESCRIPTOR;
	}
	else{
	  ret_val.first = rhs_ret.first;
	  if( lhs_ret.first != rhs_ret.first && rhs_ret.second == OA_DESCRIPTOR && isSgAddOp(curr_exp) )
	    ret_val.second = OA_DESCRIPTOR;
	  else
	    ret_val.second = IA_DESCRIPTOR;
	}
    }
  }
  else if( isSgUnaryOp(curr_exp) ){
    pair<branch*,descriptor_type> unary_ret = AnalyseSafeExpression(isSgUnaryOp(curr_exp)->get_operand(),curr_scope);
    ret_val.first = unary_ret.first;
    if( ret_val.second == OA_DESCRIPTOR)
      ret_val.second = IA_DESCRIPTOR;
  }
  return ret_val;
}
