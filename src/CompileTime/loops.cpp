/*
 * loops.cpp: This file is part of the IEC project.
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
 * @file: loops.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/loops.hpp"
#include "CompileTime/array_details.hpp"
#include "CompileTime/safe_scalars.hpp"
#include "OmpAttribute.h"
#include <sstream>
#include <vector>

extern bool CompareExpression(const SgExpression* , const SgExpression*);
extern bool isNumber(const SgExpression*);
using namespace std;

SgVariableDeclaration* partitionable_loop::iter_num_count = NULL;
deque<partitionable_loop*> partitionable_loop::all_loops;
SgFunctionSymbol* partitionable_loop::communicate_reads_fn = NULL;
SgFunctionSymbol* partitionable_loop::init_write_ghosts_fn = NULL;
SgFunctionSymbol* partitionable_loop::communicate_writes_fn = NULL;
SgFunctionSymbol* partitionable_loop::reduce_scalar_fn = NULL;
SgFunctionSymbol* inner_loop::is_not_known_scalar_fn = NULL;
SgFunctionSymbol* partitionable_loop::memset_fn = NULL;
SgFunctionSymbol* partitionable_loop::get_proc_iter_size_fn = NULL;
SgFunctionSymbol* partitionable_loop::add_vertex_fn = NULL;
SgFunctionSymbol* partitionable_loop::get_vertex_home_fn = NULL;
//SgFunctionSymbol* partitionable_loop::add_pin_to_net_fn = NULL;
SgFunctionSymbol* partitionable_loop::is_known_fn = NULL;
SgFunctionSymbol* partitionable_loop::get_elem_fn = NULL;

int inner_loop::total_inner = 0;

branch::branch( bool isu ) : is_used(isu) {}

void branch::GenerateSetGaurds(SgStatement* insert_before, SgExpression* status_exp) const
{
  if( is_used ){ 
    for( set<SgVariableSymbol*>::iterator it = written_iscalars.begin() ; it != written_iscalars.end() ; it++ ){
      assert((*it)->attributeExists("indirection_scalar"));
      safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>((*it)->getAttribute("indirection_scalar"));
      curr_safe_scalar->SetScalarUnknown(insert_before);
    }
    if( status_exp){
      SgStatement* set_iter_status = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(status_exp),SageBuilder::buildIntVal(0)));
      SageInterface::insertStatementBefore(insert_before,set_iter_status);
    }
  }
}

loop::loop(SgForStatement* of, const int dp):
  orig_for(of),
  executor_for(NULL),
  depth(dp),
  branch(false)
{
  SgForInitStatement*  orig_init = orig_for->get_for_init_stmt();
  SgStatementPtrList& orig_init_list = orig_init->get_init_stmt();
  SgExprStatement* orig_init_stmt = isSgExprStatement(*(orig_init_list.begin()));
  assert(orig_init_stmt);
  SgAssignOp* orig_expr = isSgAssignOp(orig_init_stmt->get_expression());
  assert(orig_expr);
    //Set the iterator
  SgVarRefExp* orig_iter_var = isSgVarRefExp(orig_expr->get_lhs_operand());
  assert(orig_iter_var);
  orig_iterator = orig_iter_var->get_symbol();

  //Set the lower bounds exp
  lb_exp = orig_expr->get_rhs_operand();

  //Get the upper bounds exp
  SgExprStatement* orig_ub_stmt = isSgExprStatement(orig_for->get_test());
  assert(orig_ub_stmt);
  SgLessThanOp* orig_ub_exp = isSgLessThanOp(orig_ub_stmt->get_expression());
  assert(orig_ub_exp);
  ub_exp = orig_ub_exp->get_rhs_operand();
}

void loop::AddPrivateVariables(OmpSupport::OmpAttribute* curr_attribute) const
{
  curr_attribute->addVariable(OmpSupport::e_firstprivate,orig_iterator->get_name().getString());
}


partitionable_loop::~partitionable_loop()
{
  DeleteAttributes();
  read_data_arrays.clear();
  write_data_arrays.clear();
  read_data_scalars.clear();
  write_data_scalars.clear();
  indirection_scalars.clear();
  indirection_arrays.clear();
  i_assignment_stmts.clear();
}

void partitionable_loop::DeleteAttributes()
{
  class deleteAttributes: public AstPrePostProcessing {
  public:
    void preOrderVisit( SgNode* node ){
      if( isSgForStatement(node) ){
	if( node->attributeExists("inner_loop") ){
	  inner_loop* curr_inner = static_cast<inner_loop*>(node->getAttribute("inner_loop"));
	  delete curr_inner;
	  node->removeAttribute("inner_loop");
	}
	if( node->attributeExists("scope")) {
	  scope* curr_scope = static_cast<scope*>(node->getAttribute("scope"));
	  delete curr_scope;
	  node->removeAttribute("scope");
	}
      }
      if( isSgIfStmt(node) ){
	if( node->attributeExists("then_branch") ){
	  then_branch* curr_then = static_cast<then_branch*>(node->getAttribute("then_branch"));
	  delete curr_then;
	  node->removeAttribute("then_branch");
	  if( node->attributeExists("else_branch") ){
	    else_branch* curr_else = static_cast<else_branch*>(node->getAttribute("else_branch"));
	    delete curr_else;
	    node->removeAttribute("else_branch");
	  }
	}
	if( node->attributeExists("scope")) {
	  scope* curr_scope = static_cast<scope*>(node->getAttribute("scope"));
	  delete curr_scope;
	  node->removeAttribute("scope");
	}
      }
      if( isSgExprStatement(node) ){
	if( node->attributeExists("assignment_type") ){
	  IAssignmentAttr* curr_type = static_cast<IAssignmentAttr*>(node->getAttribute("assignment_type"));
	  delete curr_type;
	  node->removeAttribute("assignment_type");
	}
	if( node->attributeExists("scope")) {
	  scope* curr_scope = static_cast<scope*>(node->getAttribute("scope"));
	  delete curr_scope;
	  node->removeAttribute("scope");
	}
      }
    }
    void postOrderVisit(SgNode* node) { }
  };

  deleteAttributes delete_visitor;
  delete_visitor.traverse(orig_for->get_loop_body());
}


void partitionable_loop::GenerateLocalBlock(SgStatement* insert_before, SgFunctionSymbol* malloc_fn)
{
  //Start Values
  SgSubtractOp* diff_exp = SageBuilder::buildBinaryExpression<SgSubtractOp>(SageInterface::copyExpression(ub_exp),SageInterface::copyExpression(lb_exp));
  SgName nprocs_var_name("__nprocs__");
  SgName proc_id_var("__proc_id__");
  SgMultiplyOp* start_init_split = SageBuilder::buildBinaryExpression<SgMultiplyOp>(SageBuilder::buildBinaryExpression<SgDivideOp>(diff_exp,SageBuilder::buildVarRefExp(nprocs_var_name)),SageBuilder::buildVarRefExp(proc_id_var));
  SgAddOp* start_init = SageBuilder::buildBinaryExpression<SgAddOp>(start_init_split,SageInterface::copyExpression(lb_exp));
  SgAssignInitializer* start_init_assign = SageBuilder::buildAssignInitializer(start_init);
  ostringstream local_block_start_stream;
  local_block_start_stream << "__local_block_" << my_num << "_start__" ;
  SgName local_block_start_name(local_block_start_stream.str());
  local_block_start = SageBuilder::buildVariableDeclaration(local_block_start_name,SageBuilder::buildIntType(),start_init_assign);
  SageInterface::insertStatementBefore(insert_before,local_block_start);

  ostringstream local_block_end_stream;
  local_block_end_stream << "__local_block_" << my_num << "_end__" ;
  SgName local_block_end_name(local_block_end_stream.str());
  local_block_end = SageBuilder::buildVariableDeclaration(local_block_end_name,SageBuilder::buildIntType());
  SageInterface::insertStatementBefore(insert_before,local_block_end);
  SgSubtractOp* condn_rhs = SageBuilder::buildBinaryExpression<SgSubtractOp>(SageBuilder::buildVarRefExp(nprocs_var_name),SageBuilder::buildIntVal(1));
  SgNotEqualOp* condn = SageBuilder::buildBinaryExpression<SgNotEqualOp>(SageBuilder::buildVarRefExp(proc_id_var),condn_rhs);
  SgSubtractOp* diff_exp_ub = SageBuilder::buildBinaryExpression<SgSubtractOp>(SageInterface::copyExpression(ub_exp),SageInterface::copyExpression(lb_exp));
  SgDivideOp* div_exp_ub = SageBuilder::buildBinaryExpression<SgDivideOp>(diff_exp_ub,SageBuilder::buildVarRefExp(nprocs_var_name));
  SgAddOp* mult_exp_ub_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(SageBuilder::buildVarRefExp(proc_id_var),SageBuilder::buildIntVal(1));
  SgMultiplyOp* ub_init_exp = SageBuilder::buildBinaryExpression<SgMultiplyOp>(div_exp_ub,mult_exp_ub_rhs);
  SgAssignOp* ub_init_then_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(local_block_end),SageBuilder::buildBinaryExpression<SgAddOp>(ub_init_exp,SageInterface::copyExpression(lb_exp)));
  SgExprStatement* ub_init_then = SageBuilder::buildExprStatement(ub_init_then_exp);
  SgAssignOp* ub_init_else_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(local_block_end),SageInterface::copyExpression(ub_exp));
  SgExprStatement* ub_init_else = SageBuilder::buildExprStatement(ub_init_else_exp);
  SgIfStmt* ub_stmt = SageBuilder::buildIfStmt(condn,ub_init_then,ub_init_else);
  SageInterface::insertStatementBefore(insert_before,ub_stmt);

  ostringstream status_array_stream;
  status_array_stream << "__" << my_num << "_status__" ;
  SgName status_array_name(status_array_stream.str());
  SgMultiplyOp* status_array_size = SageBuilder::buildBinaryExpression<SgMultiplyOp>(SageBuilder::buildSizeOfOp(SageBuilder::buildIntType()),SageBuilder::buildBinaryExpression<SgSubtractOp>(SageBuilder::buildVarRefExp(local_block_end),SageBuilder::buildVarRefExp(local_block_start)));
  SgAssignInitializer* status_array_init = SageBuilder::buildAssignInitializer(SageBuilder::buildCastExp(SageBuilder::buildFunctionCallExp(malloc_fn,SageBuilder::buildExprListExp(status_array_size)),SageBuilder::buildPointerType(SageBuilder::buildIntType())));
  status_array = SageBuilder::buildVariableDeclaration(status_array_name,SageBuilder::buildPointerType(SageBuilder::buildIntType()),status_array_init);
  SageInterface::insertStatementBefore(insert_before,status_array);

  SgFunctionCallExp* status_array_set0 = SageBuilder::buildFunctionCallExp(memset_fn,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(status_array),SageBuilder::buildIntVal(0),SageInterface::copyExpression(status_array_size)));
  SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(status_array_set0));
}


void partitionable_loop::ResetStatusArray(SgStatement* insert_before, SgFunctionSymbol* malloc_fn, SgFunctionSymbol* free_fn) const
{
  SgExprStatement* free_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(free_fn,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(status_array))));
  SageInterface::insertStatementBefore(insert_before,free_stmt);
  
  SgFunctionCallExp* size_exp = SageBuilder::buildFunctionCallExp(get_proc_iter_size_fn,SageBuilder::buildExprListExp(SageBuilder::buildIntVal(my_group)));
  SgMultiplyOp* array_size_exp = SageBuilder::buildBinaryExpression<SgMultiplyOp>(SageBuilder::buildSizeOfOp(SageBuilder::buildIntType()),size_exp);
  SgExprStatement* init_array = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(status_array),SageBuilder::buildCastExp(SageBuilder::buildFunctionCallExp(malloc_fn,SageBuilder::buildExprListExp(array_size_exp)),SageBuilder::buildPointerType(SageBuilder::buildIntType()))));
  SageInterface::insertStatementBefore(insert_before,init_array);
  
  SgExprStatement* set_array0 = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(memset_fn,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(status_array),SageBuilder::buildIntVal(0),SageInterface::copyExpression(array_size_exp))));
  SageInterface::insertStatementBefore(insert_before,set_array0);
}


void partitionable_loop::FreeStatusArray(SgStatement* insert_before, SgFunctionSymbol* free_fn) const
{
  SgExprStatement* free_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(free_fn,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(status_array))));
  SageInterface::insertStatementBefore(insert_before,free_stmt);
}

void partitionable_loop::ResetShadowScalars(SgStatement* insert_before) const
{
  for( set<SgVariableSymbol*>::iterator it = written_iscalars.begin() ; it != written_iscalars.end() ; it++ ){
    assert((*it)->attributeExists("indirection_scalar"));
    safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>((*it)->getAttribute("indirection_scalar"));
    curr_safe_scalar->SetScalarUnknown(insert_before);
  }
}

void partitionable_loop::InitCounters(SgStatement* insert_before)
{
  assert(body_counter == NULL);
  ostringstream iterator_name_stream ;
  iterator_name_stream << "__body_" << my_num << "__" ;
  string new_counter_name = iterator_name_stream.str();
  SgAssignInitializer* new_init = SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0));
  body_counter = SageBuilder::buildVariableDeclaration(new_counter_name,SageBuilder::buildIntType(),new_init);
  SageInterface::insertStatementBefore(insert_before,body_counter);

  inner_loop_nodes = NodeQuery::querySubTree(orig_for->get_loop_body(),V_SgForStatement);
  for( SgNodePtrList::const_iterator it = inner_loop_nodes.begin() ; it != inner_loop_nodes.end() ; it++ ){
    assert((*it)->attributeExists("inner_loop"));
    inner_loop* next_inner = static_cast<inner_loop*>((*it)->getAttribute("inner_loop"));
    next_inner->InitCounters(insert_before);
  }

  inner_condn_nodes = NodeQuery::querySubTree(orig_for->get_loop_body(),V_SgIfStmt);
  for( SgNodePtrList::const_iterator it = inner_condn_nodes.begin() ; it != inner_condn_nodes.end() ; it++ ){
    assert((*it)->attributeExists("then_branch"));
    then_branch* next_then = static_cast<then_branch*>((*it)->getAttribute("then_branch"));
    next_then->InitIfThenCounter(insert_before);
    if( (*it)->attributeExists("else_branch")){
      else_branch* next_else = static_cast<else_branch*>((*it)->getAttribute("else_branch"));
      next_else->InitElseCounter(insert_before);
    }
  }
}


void partitionable_loop::GenerateResetCounters(SgStatement* insert_before) const 
{
  assert(body_counter);
  SgAssignOp* body_reset_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(body_counter),SageBuilder::buildIntVal(0));
  SgExprStatement* body_reset_stmt = SageBuilder::buildExprStatement(body_reset_op);
  SageInterface::insertStatementBefore(insert_before,body_reset_stmt);
}

void partitionable_loop::InitBoundsArrays(SgStatement* insert_before, SgFunctionSymbol* a)
{
  assert(body_counter);
  ostringstream local_iter_stream;
  local_iter_stream << "__nlocal_" << my_num << "__" ;
  string nlocal_iter_string = local_iter_stream.str();
  SgName nlocal_iter_name(nlocal_iter_string);
  SgAssignInitializer* nlocal_iter_rhs = SageBuilder::buildAssignInitializer(SageBuilder::buildVarRefExp(body_counter));
  nlocal_iters = SageBuilder::buildVariableDeclaration(nlocal_iter_name,SageBuilder::buildIntType(),nlocal_iter_rhs);
  SageInterface::insertStatementBefore(insert_before,nlocal_iters);
  for( SgNodePtrList::const_iterator it = inner_loop_nodes.begin() ; it != inner_loop_nodes.end() ; it++ ){
    assert((*it)->attributeExists("inner_loop"));
    inner_loop* next_inner = static_cast<inner_loop*>((*it)->getAttribute("inner_loop"));
    next_inner->InitBoundsArrays(insert_before,a);
  }
  for( SgNodePtrList::const_iterator it = inner_condn_nodes.begin() ; it != inner_condn_nodes.end() ; it++ ){
    assert((*it)->attributeExists("then_branch"));
    then_branch* next_then = static_cast<then_branch*>((*it)->getAttribute("then_branch"));
    next_then->InitCondnArray(insert_before,a);
  }
}


void partitionable_loop::ResetAllCounters(SgStatement* insert_before)  const
{
  GenerateResetCounters(insert_before);
  for( SgNodePtrList::const_iterator it = inner_loop_nodes.begin() ; it != inner_loop_nodes.end() ; it++ ){
    assert((*it)->attributeExists("inner_loop"));
    inner_loop* next_inner = static_cast<inner_loop*>((*it)->getAttribute("inner_loop"));
    next_inner->GenerateResetCounters(insert_before);
  }
  for( SgNodePtrList::const_iterator it = inner_condn_nodes.begin() ; it != inner_condn_nodes.end() ; it++ ){
    assert((*it)->attributeExists("then_branch"));
    then_branch* next_then = static_cast<then_branch*>((*it)->getAttribute("then_branch"));
    next_then->ResetIfCounter(insert_before);
    next_then->ResetBranchCounter(insert_before);
    if( (*it)->attributeExists("else_branch")){
      else_branch* next_else = static_cast<else_branch*>((*it)->getAttribute("else_branch"));
      next_else->ResetBranchCounter(insert_before);
    }    
  }  
}

void partitionable_loop::GenerateBodyCounterIncrement(SgStatement* insert_before) const
{
  assert(body_counter);
  SgVarRefExp* body_counter_ref = SageBuilder::buildVarRefExp(body_counter);
  SgIntVal* body_counter_inc_val = SageBuilder::buildIntVal(1);
  SgAddOp* body_counter_inc_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(body_counter_ref,body_counter_inc_val);
  SgAssignOp* body_counter_inc = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(body_counter),body_counter_inc_rhs);
  SgExprStatement* body_counter_inc_stmt = SageBuilder::buildExprStatement(body_counter_inc);
  SageInterface::insertStatementBefore(insert_before,body_counter_inc_stmt);
}


void partitionable_loop::DeleteArrays(SgStatement* insert_before, SgFunctionSymbol* free_fn) const
{
  for( SgNodePtrList::const_iterator it = inner_loop_nodes.begin() ; it != inner_loop_nodes.end() ; it++ ){
    inner_loop* curr_inner = static_cast<inner_loop*>((*it)->getAttribute("inner_loop"));
    curr_inner->DeleteBoundsArrays(insert_before,free_fn);
  }
  for( SgNodePtrList::const_iterator it = inner_condn_nodes.begin() ; it != inner_condn_nodes.end() ; it++ ){
    then_branch* curr_condn = static_cast<then_branch*>((*it)->getAttribute("then_branch"));
    curr_condn->DeleteCondnArray(insert_before,free_fn);
  }
}


inner_loop::inner_loop(SgForStatement* of, const int dp):
  loop(of,dp),
  my_num(total_inner++)
{
  //printf("New Inner Loop: depth:%d, lb_exp:%s, ub_exp:%s\n",depth,lb_exp->unparseToString().c_str(),ub_exp->unparseToString().c_str());
  body_counter = NULL;
  loop_counter = NULL;
  lb_array = NULL;
  ub_array = NULL;
}

void inner_loop::InitCounters(SgStatement* insert_before)
{
  if( is_used ){
    assert(loop_counter == NULL);
    ostringstream loop_iterator_stream;
    loop_iterator_stream << "__loop_i_" << my_num << "__" ;
    string loop_counter_name = loop_iterator_stream.str();
    SgAssignInitializer* loop_init = SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0));
    loop_counter = SageBuilder::buildVariableDeclaration(loop_counter_name,SageBuilder::buildIntType(),loop_init);
    SageInterface::insertStatementBefore(insert_before,loop_counter);

    assert(body_counter == NULL);
    ostringstream body_iterator_stream;
    body_iterator_stream << "__body_i_" << my_num << "__" ;
    string body_counter_name = body_iterator_stream.str();
    SgAssignInitializer* body_init = SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0));
    body_counter = SageBuilder::buildVariableDeclaration(body_counter_name,SageBuilder::buildIntType(),body_init);
    SageInterface::insertStatementBefore(insert_before,body_counter);
  }
}


void inner_loop::InitBoundsArrays(SgStatement* insert_before, SgFunctionSymbol* malloc_fn)
{
  if( is_used ){
    SgPointerType* array_type = SageBuilder::buildPointerType(SageBuilder::buildIntType());
  
    //Lower bounds malloc fn args
    SgVarRefExp* lb_size = SageBuilder::buildVarRefExp(loop_counter);
    SgSizeOfOp* lb_type_size = SageBuilder::buildSizeOfOp(SageBuilder::buildIntType());
    SgMultiplyOp* lb_array_size = SageBuilder::buildBinaryExpression<SgMultiplyOp>(lb_size,lb_type_size);
    SgExprListExp* lb_malloc_fn_args = SageBuilder::buildExprListExp(lb_array_size);

    //Lower bound array init exp
    SgFunctionCallExp* lb_malloc_fn_exp = SageBuilder::buildFunctionCallExp(malloc_fn,lb_malloc_fn_args);
    SgCastExp* lb_init_exp = SageBuilder::buildCastExp(lb_malloc_fn_exp,array_type);
    SgAssignInitializer* lb_init = SageBuilder::buildAssignInitializer(lb_init_exp);
  
    //Lower Bound array declaration
    ostringstream lb_name_stream;
    lb_name_stream << "__lb_" << orig_iterator->get_name().getString() << "_i_" << my_num << "__" ;
    SgName lb_name(lb_name_stream.str());
    lb_array = SageBuilder::buildVariableDeclaration(lb_name,array_type,lb_init);
    SageInterface::insertStatementBefore(insert_before,lb_array);

    //Upper bounds malloc fn args
    SgVarRefExp* ub_size = SageBuilder::buildVarRefExp(loop_counter);
    SgSizeOfOp* ub_type_size = SageBuilder::buildSizeOfOp(SageBuilder::buildIntType());
    SgMultiplyOp* ub_array_size = SageBuilder::buildBinaryExpression<SgMultiplyOp>(ub_size,ub_type_size);
    SgExprListExp* ub_malloc_fn_args = SageBuilder::buildExprListExp(ub_array_size);

    //Upper bound array init exp
    SgFunctionCallExp* ub_malloc_fn_exp = SageBuilder::buildFunctionCallExp(malloc_fn,ub_malloc_fn_args);
    SgCastExp* ub_init_exp = SageBuilder::buildCastExp(ub_malloc_fn_exp,array_type);
    SgAssignInitializer* ub_init = SageBuilder::buildAssignInitializer(ub_init_exp);
  
    //Upper Bound array declaration
    ostringstream ub_name_stream;
    ub_name_stream << "__ub_" << orig_iterator->get_name().getString() << "_i_" << my_num << "__" ;
    SgName ub_name(ub_name_stream.str());
    ub_array = SageBuilder::buildVariableDeclaration(ub_name,array_type,ub_init);
    SageInterface::insertStatementBefore(insert_before,ub_array);  
  }
}


void inner_loop::DeleteBoundsArrays(SgStatement* insert_before, SgFunctionSymbol* free_fn) const 
{
  if( is_used ){
    if( lb_array){
      SgExprListExp* fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(lb_array));
      SgFunctionCallExp* fn_call_exp = SageBuilder::buildFunctionCallExp(free_fn,fn_args);
      SgExprStatement* fn_call_stmt = SageBuilder::buildExprStatement(fn_call_exp);
      SageInterface::insertStatementBefore(insert_before,fn_call_stmt);
    }
    if( ub_array ){
      SgExprListExp* fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(ub_array));
      SgFunctionCallExp* fn_call_exp = SageBuilder::buildFunctionCallExp(free_fn,fn_args);
      SgExprStatement* fn_call_stmt = SageBuilder::buildExprStatement(fn_call_exp);
      SageInterface::insertStatementBefore(insert_before,fn_call_stmt);
    }
  }
}


pair<SgStatement*,SgStatement*> inner_loop::GenerateCArrayStore() const 
{
  pair<SgStatement*,SgStatement*> ret_stmts;
  ret_stmts.first = NULL; ret_stmts.second = NULL;
  if( is_used ){
    //Update the lower bound
    SgPntrArrRefExp* lb_lhs = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(lb_array),SageBuilder::buildVarRefExp(loop_counter));
    SgAssignOp* lb_assign_op = SageBuilder::buildBinaryExpression<SgAssignOp>(lb_lhs,SageInterface::copyExpression(lb_exp));
    SgExprStatement* lb_stmt = SageBuilder::buildExprStatement(lb_assign_op);
    ret_stmts.first = lb_stmt;
    //SageInterface::insertStatementBefore(insert_before,lb_stmt);

    //Update the upper bound
    SgPntrArrRefExp* ub_lhs = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(ub_array),SageBuilder::buildVarRefExp(loop_counter));
    SgAssignOp* ub_assign_op = SageBuilder::buildBinaryExpression<SgAssignOp>(ub_lhs,SageInterface::copyExpression(ub_exp));
    SgExprStatement* ub_stmt = SageBuilder::buildExprStatement(ub_assign_op);
    ret_stmts.second = ub_stmt;
    //SageInterface::insertStatementBefore(insert_before,ub_stmt);
  }
  return ret_stmts;
}

SgIfStmt* inner_loop::GenerateIfFirstIterStmt(SgStatement* then_stmt) const
{
  if( is_used ){
    SgEqualityOp* condn = SageBuilder::buildBinaryExpression<SgEqualityOp>(SageBuilder::buildVarRefExp(orig_iterator),SageInterface::copyExpression(lb_exp));
    SgIfStmt* if_stmt = SageBuilder::buildIfStmt(condn,then_stmt,SageBuilder::buildNullStatement());
    return if_stmt;
  }
  else
    return NULL;
}


void inner_loop::GenerateLoopCounterIncrement(SgStatement* insert_before) const
{
  if( is_used ){
    assert(loop_counter);
    SgVarRefExp* loop_counter_ref = SageBuilder::buildVarRefExp(loop_counter);
    SgIntVal* loop_counter_inc_val = SageBuilder::buildIntVal(1);
    SgAddOp* loop_counter_inc_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(loop_counter_ref,loop_counter_inc_val);
    SgAssignOp* loop_counter_inc = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(loop_counter),loop_counter_inc_rhs);
    SgExprStatement* loop_counter_inc_stmt = SageBuilder::buildExprStatement(loop_counter_inc);
    SageInterface::insertStatementBefore(insert_before,loop_counter_inc_stmt);
  }
}

void inner_loop::GenerateBodyCounterIncrement(SgStatement* insert_before) const
{
  if( is_used ){
    assert(body_counter);
    SgVarRefExp* body_counter_ref = SageBuilder::buildVarRefExp(body_counter);
    SgIntVal* body_counter_inc_val = SageBuilder::buildIntVal(1);
    SgAddOp* body_counter_inc_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(body_counter_ref,body_counter_inc_val);
    SgAssignOp* body_counter_inc = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(body_counter),body_counter_inc_rhs);
    SgExprStatement* body_counter_inc_stmt = SageBuilder::buildExprStatement(body_counter_inc);
    SageInterface::insertStatementBefore(insert_before,body_counter_inc_stmt);
  }
}

void inner_loop::GenerateLoopCounterIncrementAfter(SgStatement* insert_after) const
{
  if( is_used ){
    assert(loop_counter);
    SgVarRefExp* loop_counter_ref = SageBuilder::buildVarRefExp(loop_counter);
    SgIntVal* loop_counter_inc_val = SageBuilder::buildIntVal(1);
    SgAddOp* loop_counter_inc_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(loop_counter_ref,loop_counter_inc_val);
    SgAssignOp* loop_counter_inc = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(loop_counter),loop_counter_inc_rhs);
    SgExprStatement* loop_counter_inc_stmt = SageBuilder::buildExprStatement(loop_counter_inc);
    SageInterface::insertStatementAfter(insert_after,loop_counter_inc_stmt);
  }
}

void inner_loop::GenerateBodyCounterIncrementAfter(SgStatement* insert_after) const
{
  if( is_used ){
    assert(body_counter);
    SgVarRefExp* body_counter_ref = SageBuilder::buildVarRefExp(body_counter);
    SgIntVal* body_counter_inc_val = SageBuilder::buildIntVal(1);
    SgAddOp* body_counter_inc_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(body_counter_ref,body_counter_inc_val);
    SgAssignOp* body_counter_inc = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(body_counter),body_counter_inc_rhs);
    SgExprStatement* body_counter_inc_stmt = SageBuilder::buildExprStatement(body_counter_inc);
    SageInterface::insertStatementAfter(insert_after,body_counter_inc_stmt);
  }
}


void inner_loop::GenerateResetCounters(SgStatement* insert_before) const 
{
  if( is_used ){
    assert(loop_counter && body_counter);
    SgAssignOp* loop_reset_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(loop_counter),SageBuilder::buildIntVal(0));
    SgExprStatement* loop_reset_stmt = SageBuilder::buildExprStatement(loop_reset_op);
    SageInterface::insertStatementBefore(insert_before,loop_reset_stmt);

    SgAssignOp* body_reset_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(body_counter),SageBuilder::buildIntVal(0));
    SgExprStatement* body_reset_stmt = SageBuilder::buildExprStatement(body_reset_op);
    SageInterface::insertStatementBefore(insert_before,body_reset_stmt);
  }
}

void inner_loop::GenerateResetCountersAfter(SgStatement*& insert_after) const 
{
  if( is_used ){
    assert(loop_counter && body_counter);
    SgAssignOp* loop_reset_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(loop_counter),SageBuilder::buildIntVal(0));
    SgExprStatement* loop_reset_stmt = SageBuilder::buildExprStatement(loop_reset_op);
    SageInterface::insertStatementAfter(insert_after,loop_reset_stmt);
    insert_after = loop_reset_stmt;
    
    SgAssignOp* body_reset_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(body_counter),SageBuilder::buildIntVal(0));
    SgExprStatement* body_reset_stmt = SageBuilder::buildExprStatement(body_reset_op);
    SageInterface::insertStatementAfter(insert_after,body_reset_stmt);
    insert_after = body_reset_stmt;
  }
}


