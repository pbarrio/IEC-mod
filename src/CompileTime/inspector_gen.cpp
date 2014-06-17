/*
 * inspector_gen.cpp: This file is part of the IEC project.
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
 * @file: inspector_gen.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/loops.hpp"
#include "CompileTime/conditional.hpp"
#include "CompileTime/safe_scalars.hpp"
#include "CompileTime/array_details.hpp"
#include "CompileTime/indarray_details.hpp"
#include "CompileTime/driver.hpp"

using namespace std;

void partitionable_loop::GenerateCInspector(SgStatement* graph_gen, SgStatement* init_counter, SgStatement* init_thread_counter, SgStatement* init_access_array)
{
  printf("Inspector for :%s\n",orig_for->unparseToString().c_str());

  //Graph loop header
  ResetShadowScalars(graph_gen);
  SgAssignOp* graph_loop_init_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(orig_iterator),SageBuilder::buildVarRefExp(local_block_start));
  SgExprStatement* graph_loop_init_stmt = SageBuilder::buildExprStatement(graph_loop_init_exp);
  SgBasicBlock* graph_loop_init_bb = SageBuilder::buildBasicBlock(graph_loop_init_stmt);
  SgForInitStatement* graph_loop_init = SageBuilder::buildForInitStatement(graph_loop_init_bb->get_statements());
  SgLessThanOp* graph_loop_test_exp = SageBuilder::buildBinaryExpression<SgLessThanOp>(SageBuilder::buildVarRefExp(orig_iterator),SageBuilder::buildVarRefExp(local_block_end));
  SgExprStatement* graph_loop_test = SageBuilder::buildExprStatement(graph_loop_test_exp);
  SgExpression* graph_loop_increment = (SageInterface::copyExpression(orig_for->get_increment()));

  SgPntrArrRefExp* curr_iter_ref = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(status_array),SageBuilder::buildBinaryExpression<SgSubtractOp>(SageBuilder::buildVarRefExp(orig_iterator),SageBuilder::buildVarRefExp(local_block_start)));
  SgEqualityOp* curr_iter_status_exp = SageBuilder::buildBinaryExpression<SgEqualityOp>(curr_iter_ref,SageBuilder::buildIntVal(0));
  SgNullStatement* graph_loop_body = SageBuilder::buildNullStatement();
  SgIfStmt* curr_iter_status = SageBuilder::buildIfStmt(curr_iter_status_exp,graph_loop_body,SageBuilder::buildNullStatement());
  SgForStatement* graph_loop = SageBuilder::buildForStatement(graph_loop_init,graph_loop_test,graph_loop_increment,curr_iter_status);
  SageInterface::insertStatementBefore(graph_gen,graph_loop);

  //Set status of curr iterations
  SgAssignOp* curr_iter_set_status_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(curr_iter_ref),SageBuilder::buildIntVal(1));
  SageInterface::insertStatementBefore(graph_loop_body,SageBuilder::buildExprStatement(curr_iter_set_status_exp));

  //Graph Gen body
  assert(add_vertex_fn);
  SgIntVal* add_vertex_arg1 = SageBuilder::buildIntVal(my_group);
  SgVarRefExp* add_vertex_arg2 = SageBuilder::buildVarRefExp(orig_iterator);
  SgExprListExp* add_vertex_fn_args = SageBuilder::buildExprListExp(add_vertex_arg1,add_vertex_arg2);
  SgFunctionCallExp* add_vertex_fn_call = SageBuilder::buildFunctionCallExp(add_vertex_fn,add_vertex_fn_args);
  SgExprStatement* add_vertex_fn_stmt = SageBuilder::buildExprStatement(add_vertex_fn_call);
  SageInterface::insertStatementBefore(graph_loop_body,add_vertex_fn_stmt);

  //Declare Status counter
  ostringstream status_counter_stream;
  status_counter_stream << "__" << my_num << "_status_counter__" ;
  SgName status_counter_name(status_counter_stream.str());
  SgVariableDeclaration* status_counter_var = SageBuilder::buildVariableDeclaration(status_counter_name,SageBuilder::buildIntType(),SageBuilder::buildAssignInitializer(SageBuilder::buildIntVal(0)));
  SageInterface::insertStatementBefore(init_counter,status_counter_var);

  //Init counter Loop header
  ResetShadowScalars(init_counter);
  SgForInitStatement* init_counter_loop_init = isSgForInitStatement(SageInterface::copyStatement(orig_for->get_for_init_stmt()));
  SgStatement* init_counter_loop_test = (SageInterface::copyStatement(orig_for->get_test()));
  SgExpression* init_counter_loop_increment = (SageInterface::copyExpression(orig_for->get_increment()));
  SgNullStatement* init_counter_loop_body = SageBuilder::buildNullStatement();
  SgForStatement* init_counter_loop = SageBuilder::buildForStatement(init_counter_loop_init,init_counter_loop_test,init_counter_loop_increment,init_counter_loop_body);
  SageInterface::insertStatementBefore(init_counter,init_counter_loop);
  
  assert(get_vertex_home_fn);

  //Init Counters body
  //Condition to check if iteration maps to current process
  SgIntVal* init_counter_arg1 = SageBuilder::buildIntVal(my_group);
  SgVarRefExp* init_counter_arg2 = SageBuilder::buildVarRefExp(orig_iterator);
  SgExprListExp* init_counter_fn_args = SageBuilder::buildExprListExp(init_counter_arg1,init_counter_arg2);
  SgFunctionCallExp* init_counter_fn_call = SageBuilder::buildFunctionCallExp(get_vertex_home_fn,init_counter_fn_args);
//   SgName init_counter_nthreads_name("__nthreads__");
  SgName init_counter_procid_name("__proc_id__");
//   SgVarRefExp* init_counter_nthreads_var = SageBuilder::buildVarRefExp(init_counter_nthreads_name);
//   SgDivideOp* init_counter_lhs = SageBuilder::buildBinaryExpression<SgDivideOp>(init_counter_fn_call,init_counter_nthreads_var);
  SgVarRefExp* init_counter_procid_var = SageBuilder::buildVarRefExp(init_counter_procid_name);
//   SgEqualityOp* init_counter_condn_exp = SageBuilder::buildBinaryExpression<SgEqualityOp>(init_counter_lhs,init_counter_procid_var);
  SgEqualityOp* init_counter_condn_exp = SageBuilder::buildBinaryExpression<SgEqualityOp>(init_counter_fn_call,init_counter_procid_var);
  SgNullStatement* init_counter_local_then = SageBuilder::buildNullStatement();
  SgIfStmt* init_counter_local_iter_stmt = SageBuilder::buildIfStmt(init_counter_condn_exp,init_counter_local_then,SageBuilder::buildNullStatement());
  SageInterface::insertStatementBefore(init_counter_loop_body,init_counter_local_iter_stmt);

  //Condition to check the status of the iteration
  SgPntrArrRefExp* init_counter_iter_status_exp = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(GetStatusArray()),SageBuilder::buildVarRefExp(status_counter_var));
  SgEqualityOp* init_counter_iter_check = SageBuilder::buildBinaryExpression<SgEqualityOp>(init_counter_iter_status_exp,SageBuilder::buildIntVal(0));
  SgNullStatement* init_counter_then_body = SageBuilder::buildNullStatement();
  SgIfStmt* init_counter_iter_status_stmt = SageBuilder::buildIfStmt(init_counter_iter_check,init_counter_then_body,SageBuilder::buildNullStatement());
  SageInterface::insertStatementBefore(init_counter_local_then,init_counter_iter_status_stmt);

  //Statement to increment the status counter
  SgExprStatement* init_counter_iter_status_increment = SageBuilder::buildExprStatement(SageBuilder::buildUnaryExpression<SgPlusPlusOp>(SageBuilder::buildVarRefExp(status_counter_var)));
  SageInterface::insertStatementBefore(init_counter_local_then,init_counter_iter_status_increment);

  //Statement to set the status of the iterations
  SgExprStatement* init_counter_iter_status_set_stmt = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(init_counter_iter_status_exp),SageBuilder::buildIntVal(1)));
  SageInterface::insertStatementBefore(init_counter_then_body,init_counter_iter_status_set_stmt);


  //Init thread counters header and body
  SgNullStatement* init_thread_counter_then_body = NULL;
  if( init_thread_counter ){
    SgForInitStatement* init_thread_counter_loop_init = isSgForInitStatement(SageInterface::copyStatement(orig_for->get_for_init_stmt()));
    SgStatement* init_thread_counter_loop_test = (SageInterface::copyStatement(orig_for->get_test()));
    SgExpression* init_thread_counter_loop_increment = (SageInterface::copyExpression(orig_for->get_increment()));
    SgNullStatement* init_thread_counter_loop_body = SageBuilder::buildNullStatement();
    SgForStatement* init_thread_counter_loop = SageBuilder::buildForStatement(init_thread_counter_loop_init,init_thread_counter_loop_test,init_thread_counter_loop_increment,init_thread_counter_loop_body);
    SageInterface::insertStatementBefore(init_thread_counter,init_thread_counter_loop);
    
    SgIntVal* init_thread_counter_arg1 = SageBuilder::buildIntVal(my_group);
    SgVarRefExp* init_thread_counter_arg2 = SageBuilder::buildVarRefExp(orig_iterator);
    SgExprListExp* init_thread_counter_fn_args = SageBuilder::buildExprListExp(init_thread_counter_arg1,init_thread_counter_arg2);
    SgFunctionCallExp* init_thread_counter_fn_call = SageBuilder::buildFunctionCallExp(get_vertex_home_fn,init_thread_counter_fn_args);
    SgName init_thread_counter_myid_name("__myid__");
    SgVarRefExp* init_thread_counter_myid_var = SageBuilder::buildVarRefExp(init_thread_counter_myid_name);
    SgEqualityOp* init_thread_counter_condn_exp = SageBuilder::buildBinaryExpression<SgEqualityOp>(init_thread_counter_fn_call,init_thread_counter_myid_var);
    init_thread_counter_then_body = SageBuilder::buildNullStatement();
    SgNullStatement* init_thread_counter_else_body = SageBuilder::buildNullStatement();
    SgIfStmt* init_thread_counter_first_stmt = SageBuilder::buildIfStmt(init_thread_counter_condn_exp,init_thread_counter_then_body,init_thread_counter_else_body);
    SageInterface::insertStatementBefore(init_thread_counter_loop_body,init_thread_counter_first_stmt);
  }


  //Init Access arrays header
  SgForInitStatement* init_access_loop_init = isSgForInitStatement(SageInterface::copyStatement(orig_for->get_for_init_stmt()));
  SgStatement* init_access_loop_test = (SageInterface::copyStatement(orig_for->get_test()));
  SgExpression* init_access_loop_increment = (SageInterface::copyExpression(orig_for->get_increment()));
  SgNullStatement* init_access_loop_body = SageBuilder::buildNullStatement();
  SgForStatement* init_access_loop = SageBuilder::buildForStatement(init_access_loop_init,init_access_loop_test,init_access_loop_increment,init_access_loop_body);
  SageInterface::insertStatementBefore(init_access_array,init_access_loop);

  //Init Access arrays body
  SgIntVal* init_access_array_arg1 = SageBuilder::buildIntVal(my_group);
  SgVarRefExp* init_access_array_arg2 = SageBuilder::buildVarRefExp(orig_iterator);
  SgExprListExp* init_access_array_fn_args = SageBuilder::buildExprListExp(init_access_array_arg1,init_access_array_arg2);
  SgFunctionCallExp* init_access_array_fn_call = SageBuilder::buildFunctionCallExp(get_vertex_home_fn,init_access_array_fn_args);
  SgName init_access_array_myid_name("__myid__");
  SgVarRefExp* init_access_array_myid_var = SageBuilder::buildVarRefExp(init_access_array_myid_name);
  SgEqualityOp* init_access_array_condn_exp = SageBuilder::buildBinaryExpression<SgEqualityOp>(init_access_array_fn_call,init_access_array_myid_var);
  SgNullStatement* init_access_array_then_body = SageBuilder::buildNullStatement();
  SgNullStatement* init_access_array_else_body = SageBuilder::buildNullStatement();
  SgIfStmt* init_access_array_first_stmt = SageBuilder::buildIfStmt(init_access_array_condn_exp,init_access_array_then_body,init_access_array_else_body);
  SageInterface::insertStatementBefore(init_access_loop_body,init_access_array_first_stmt);


  GenerateCStatement(orig_for->get_loop_body(),graph_loop_body,init_counter_then_body,init_thread_counter_then_body,init_access_array_then_body,curr_iter_ref,init_counter_iter_status_exp);
  if( init_thread_counter_then_body )
    GenerateBodyCounterIncrement(init_thread_counter_then_body);
  else
    GenerateBodyCounterIncrement(init_counter_then_body);
  GenerateBodyCounterIncrement(init_access_array_then_body);
}



void partitionable_loop::GenerateCStatement(SgStatement* target_stmt, SgStatement* graph_stmt, SgStatement* init_counter_stmt, SgStatement* init_thread_counter_stmt, SgStatement* init_access_stmt, SgExpression* status_exp, SgExpression* counter_status_exp)
{
  VariantT v_stmt = target_stmt->variantT();

  switch(v_stmt){
  case V_SgBasicBlock:
    GenerateCBasicBlock(isSgBasicBlock(target_stmt),graph_stmt,init_counter_stmt,init_thread_counter_stmt,init_access_stmt,status_exp,counter_status_exp);
    break;
  case V_SgForStatement:
    GenerateCForStatement(isSgForStatement(target_stmt),graph_stmt,init_counter_stmt,init_thread_counter_stmt,init_access_stmt,status_exp,counter_status_exp);
    break;
  case V_SgExprStatement:
    GenerateCExprStatement(isSgExprStatement(target_stmt),graph_stmt,init_counter_stmt,init_thread_counter_stmt,init_access_stmt,status_exp,counter_status_exp);
    break;
  case V_SgIfStmt:
    GenerateCIfStatement(isSgIfStmt(target_stmt),graph_stmt,init_counter_stmt,init_thread_counter_stmt,init_access_stmt,status_exp,counter_status_exp);
    break;
  default:
    ;
    break;
  }
}



void partitionable_loop::GenerateCBasicBlock(SgBasicBlock* target_bb, SgStatement* graph_gen, SgStatement* init_counter, SgStatement* init_thread_counter, SgStatement* init_access, SgExpression* status_exp, SgExpression* counter_status_exp)
{
  SgStatementPtrList& target_list = target_bb->get_statements();
  SgStatementPtrList::const_iterator target_iter;
  
  //Basic block in graph generator
  SgNullStatement* graph_gen_stmt = SageBuilder::buildNullStatement();
  SgBasicBlock* graph_gen_bb = SageBuilder::buildBasicBlock(graph_gen_stmt);
  SageInterface::insertStatementBefore(graph_gen,graph_gen_bb);
  
  //Basic block in init_counter
  SgNullStatement* init_counter_stmt = SageBuilder::buildNullStatement();
  SgBasicBlock* init_counter_bb = SageBuilder::buildBasicBlock(init_counter_stmt);
  SageInterface::insertStatementBefore(init_counter,init_counter_bb);

  //Basic block in init_thread_counter
  SgNullStatement* init_thread_counter_stmt = NULL;
  if( init_thread_counter){
    init_thread_counter_stmt = SageBuilder::buildNullStatement();
    SgBasicBlock* init_thread_counter_bb = SageBuilder::buildBasicBlock(init_thread_counter_stmt);
    SageInterface::insertStatementBefore(init_thread_counter,init_thread_counter_bb);
  }

  //Basic block in init_access
  SgNullStatement* init_access_stmt = SageBuilder::buildNullStatement();
  SgBasicBlock* init_access_bb = SageBuilder::buildBasicBlock(init_access_stmt);
  SageInterface::insertStatementBefore(init_access,init_access_bb);

  for( target_iter = target_list.begin();target_iter != target_list.end() ; target_iter++ )
    GenerateCStatement(*target_iter,graph_gen_stmt,init_counter_stmt,init_thread_counter_stmt,init_access_stmt,status_exp,counter_status_exp);
}


void partitionable_loop::GenerateCIfStatement(SgIfStmt* target_if, SgStatement* graph_gen, SgStatement* init_counter, SgStatement* init_thread_counter, SgStatement* init_access, SgExpression* status_exp, SgExpression* counter_status_exp)
{
  SgExprStatement* target_condn_stmt = isSgExprStatement(target_if->get_conditional());

  then_branch* curr_then_branch = static_cast<then_branch*>(target_if->getAttribute("then_branch"));
  assert(curr_then_branch);

  //For statement in graph generator
  SgExpression* graph_check_conditional = GenerateCSafeStatement(target_condn_stmt->get_expression());
  SgExpression* graph_gen_conditional = SageInterface::copyExpression(target_condn_stmt->get_expression());
  ReplaceIndirectionArrays(graph_gen_conditional);
  SgStatement* graph_gen_then = SageBuilder::buildNullStatement();
  SgStatement* graph_gen_else = SageBuilder::buildNullStatement();
  SgIfStmt* graph_check_then = SageBuilder::buildIfStmt(graph_gen_conditional,graph_gen_then,graph_gen_else);
  SgStatement* graph_check_else = SageBuilder::buildNullStatement();
  SgIfStmt* graph_check_stmt = SageBuilder::buildIfStmt(graph_check_conditional,graph_check_then,graph_check_else);
  SageInterface::insertStatementBefore(graph_gen,graph_check_stmt);
  curr_then_branch->GenerateSetGaurds(graph_check_else,status_exp);

  //For statement in Init Counter
  SgExpression* init_counter_check_conditional = GenerateCSafeStatement(target_condn_stmt->get_expression());
  SgExpression* init_counter_conditional = SageInterface::copyExpression(target_condn_stmt->get_expression());
  ReplaceIndirectionArrays(init_counter_conditional);
  SgStatement* init_counter_then = SageBuilder::buildNullStatement();
  SgStatement* init_counter_else = SageBuilder::buildNullStatement();
  SgIfStmt* init_counter_check_then = SageBuilder::buildIfStmt(init_counter_conditional,init_counter_then,init_counter_else);
  SgStatement* init_counter_check_else = SageBuilder::buildNullStatement();
  SgIfStmt* init_counter_check_stmt = SageBuilder::buildIfStmt(init_counter_check_conditional,init_counter_check_then,init_counter_check_else);
  SageInterface::insertStatementBefore(init_counter,init_counter_check_stmt);
  curr_then_branch->GenerateSetGaurds(init_counter_check_else,counter_status_exp);

  //For statement in Init thread Counter
  SgStatement* init_thread_counter_then = NULL;
  SgStatement* init_thread_counter_else = NULL;
  if( init_thread_counter ){
    SgExpression* init_thread_counter_conditional = SageInterface::copyExpression(target_condn_stmt->get_expression());
    ReplaceIndirectionArrays(init_thread_counter_conditional);
    init_thread_counter_then = SageBuilder::buildNullStatement();
    init_thread_counter_else = SageBuilder::buildNullStatement();
    SgIfStmt* init_thread_counter_stmt = SageBuilder::buildIfStmt(init_thread_counter_conditional,init_thread_counter_then,init_thread_counter_else);
    SageInterface::insertStatementBefore(init_thread_counter,init_thread_counter_stmt);
  }

  //For statement in Populate array
  SgExpression* init_access_conditional = SageInterface::copyExpression(target_condn_stmt->get_expression());
  ReplaceIndirectionArrays(init_access_conditional);
  SgPntrArrRefExp* cond_array_exp = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(curr_then_branch->GetCondnArray()),SageBuilder::buildVarRefExp(curr_then_branch->GetIfCounter()));
  SgAssignOp* store_condn_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(cond_array_exp,init_access_conditional);
  SageInterface::insertStatementBefore(init_access,SageBuilder::buildExprStatement(store_condn_exp));
  SgStatement* init_access_then = SageBuilder::buildNullStatement();
  SgStatement* init_access_else = SageBuilder::buildNullStatement();
  SgIfStmt* init_access_stmt = SageBuilder::buildIfStmt(SageInterface::copyExpression(cond_array_exp),init_access_then,init_access_else);
  SageInterface::insertStatementBefore(init_access,init_access_stmt);

  GenerateCStatement(target_if->get_true_body(),graph_gen_then,init_counter_then,init_thread_counter_then,init_access_then,status_exp,counter_status_exp);
  if( init_thread_counter ){
    curr_then_branch->GenerateIfCounterIncrement(init_thread_counter);
    curr_then_branch->GenerateBranchCounterIncrement(init_thread_counter_then);
  }
  else{
    curr_then_branch->GenerateIfCounterIncrementAfter(init_counter_check_then);
    curr_then_branch->GenerateBranchCounterIncrement(init_counter_then);
  }
  curr_then_branch->GenerateIfCounterIncrement(init_access);
  curr_then_branch->GenerateBranchCounterIncrement(init_access_then);
  
  if( target_if->get_false_body() ){
    GenerateCStatement(target_if->get_false_body(),graph_gen_else,init_counter_else,init_thread_counter_else,init_access_else,status_exp,counter_status_exp);
    else_branch* curr_else_branch = static_cast<else_branch*>(target_if->getAttribute("else_branch"));
    assert(curr_else_branch);
    if( init_thread_counter )
      curr_else_branch->GenerateBranchCounterIncrement(init_thread_counter_else);
    else
      curr_else_branch->GenerateBranchCounterIncrement(init_counter_else);
    curr_else_branch->GenerateBranchCounterIncrement(init_access_else);
  }
}



void partitionable_loop::GenerateCForStatement(SgForStatement* target_for, SgStatement* graph_gen, SgStatement* init_counter, SgStatement* init_thread_counter, SgStatement* init_access, SgExpression* status_exp, SgExpression* counter_status_exp)
{
  loop* curr_loop = static_cast<loop*>(target_for->getAttribute("inner_loop"));
  assert(curr_loop);
  
  //For statement in graph generator
  if( curr_loop->IsUsed() ){
    SgForInitStatement* graph_gen_for_init = isSgForInitStatement(SageInterface::copyStatement(target_for->get_for_init_stmt()));
    SgExprStatement* graph_gen_test = isSgExprStatement(SageInterface::copyStatement(target_for->get_test()));
    assert(graph_gen_test);
    SgExpression* graph_gen_inc = SageInterface::copyExpression(target_for->get_increment());
    SgNullStatement* graph_gen_body = SageBuilder::buildNullStatement();
    SgForStatement* graph_gen_for = SageBuilder::buildForStatement(graph_gen_for_init,graph_gen_test,graph_gen_inc,graph_gen_body);
    SgExpression* graph_check_lower_bound = GenerateCSafeStatement(curr_loop->GetLBExp());
    SgExpression* graph_check_upper_bound = GenerateCSafeStatement(curr_loop->GetUBExp());
    SgExpression* graph_check_bounds_cond = NULL;
    if( graph_check_lower_bound ){
      if( graph_check_upper_bound){
  	graph_check_bounds_cond = SageBuilder::buildBinaryExpression<SgAndOp>(graph_check_lower_bound,graph_check_upper_bound);
      }
      else
  	graph_check_bounds_cond = graph_check_lower_bound;
    }
    else{
      if( graph_check_upper_bound )
  	graph_check_bounds_cond = graph_check_upper_bound;
    }
    SgStatement* graph_check_bounds_stmt;
    if( graph_check_bounds_cond ){
      graph_check_bounds_stmt = SageBuilder::buildIfStmt(graph_check_bounds_cond,graph_gen_for,SageBuilder::buildNullStatement());
      curr_loop->GenerateSetGaurds(isSgIfStmt(graph_check_bounds_stmt)->get_false_body(),status_exp);
    }      
    else
      graph_check_bounds_stmt = graph_gen_for;
    SageInterface::insertStatementBefore(graph_gen,graph_check_bounds_stmt);
    
    SgStatementPtrList& graph_gen_for_init_list = graph_gen_for_init->get_init_stmt();
    SgExprStatement* graph_gen_for_init_stmt = isSgExprStatement(*(graph_gen_for_init_list.begin()));
    assert(graph_gen_for_init_stmt);
    SgAssignOp* graph_gen_for_init_expr = isSgAssignOp(graph_gen_for_init_stmt->get_expression());
    assert(graph_gen_for_init_expr);
    ReplaceIndirectionArrays(graph_gen_for_init_expr->get_rhs_operand());
    SgLessThanOp* graph_gen_test_exp = isSgLessThanOp(graph_gen_test->get_expression());
    assert(graph_gen_test_exp);
    ReplaceIndirectionArrays(graph_gen_test_exp->get_rhs_operand());

    //For statement in Init counter
    SgForInitStatement* init_counter_for_init = isSgForInitStatement(SageInterface::copyStatement(target_for->get_for_init_stmt()));
    SgExprStatement* init_counter_test = isSgExprStatement(SageInterface::copyStatement(target_for->get_test()));
    assert(init_counter_test);
    SgExpression* init_counter_inc = SageInterface::copyExpression(target_for->get_increment());
    SgNullStatement* init_counter_body = SageBuilder::buildNullStatement();
    SgForStatement* init_counter_for = SageBuilder::buildForStatement(init_counter_for_init,init_counter_test,init_counter_inc,init_counter_body);
    SgExpression* init_counter_check_lower_bound = GenerateCSafeStatement(curr_loop->GetLBExp());
    SgExpression* init_counter_check_upper_bound = GenerateCSafeStatement(curr_loop->GetUBExp());
    SgExpression* init_counter_check_bounds_cond = NULL;
    if( init_counter_check_lower_bound ){
      if( init_counter_check_upper_bound){
  	init_counter_check_bounds_cond = SageBuilder::buildBinaryExpression<SgAndOp>(init_counter_check_lower_bound,init_counter_check_upper_bound);
      }
      else
  	init_counter_check_bounds_cond = init_counter_check_lower_bound;
    }
    else{
      if( init_counter_check_upper_bound )
  	init_counter_check_bounds_cond = init_counter_check_upper_bound;
    }
    SgStatement* init_counter_check_bounds_stmt;
    if( init_counter_check_bounds_cond ){
      init_counter_check_bounds_stmt = SageBuilder::buildIfStmt(init_counter_check_bounds_cond,init_counter_for,SageBuilder::buildNullStatement());
      curr_loop->GenerateSetGaurds(isSgIfStmt(init_counter_check_bounds_stmt)->get_false_body(),counter_status_exp);
    }
    else
      init_counter_check_bounds_stmt = init_counter_for;
    SageInterface::insertStatementBefore(init_counter,init_counter_check_bounds_stmt);

    SgStatementPtrList& init_counter_for_init_list = init_counter_for_init->get_init_stmt();
    SgExprStatement* init_counter_for_init_stmt = isSgExprStatement(*(init_counter_for_init_list.begin()));
    assert(init_counter_for_init_stmt);
    SgAssignOp* init_counter_for_init_expr = isSgAssignOp(init_counter_for_init_stmt->get_expression());
    assert(init_counter_for_init_expr);
    ReplaceIndirectionArrays(init_counter_for_init_expr->get_rhs_operand());
    SgLessThanOp* init_counter_test_exp = isSgLessThanOp(init_counter_test->get_expression());
    assert(init_counter_test_exp);
    ReplaceIndirectionArrays(init_counter_test_exp->get_rhs_operand());

    //For statement in Init thread counter
    SgNullStatement* init_thread_counter_body = NULL;
    if( init_thread_counter ){
      SgForInitStatement* init_thread_counter_for_init = isSgForInitStatement(SageInterface::copyStatement(target_for->get_for_init_stmt()));
      SgExprStatement* init_thread_counter_test = isSgExprStatement(SageInterface::copyStatement(target_for->get_test()));
      assert(init_thread_counter_test);
      SgExpression* init_thread_counter_inc = SageInterface::copyExpression(target_for->get_increment());
      init_thread_counter_body = SageBuilder::buildNullStatement();
      SgForStatement* init_thread_counter_for = SageBuilder::buildForStatement(init_thread_counter_for_init,init_thread_counter_test,init_thread_counter_inc,init_thread_counter_body);
      SageInterface::insertStatementBefore(init_thread_counter,init_thread_counter_for );
      SgStatementPtrList& init_thread_counter_for_init_list = init_thread_counter_for_init->get_init_stmt();
      SgExprStatement* init_thread_counter_for_init_stmt = isSgExprStatement(*(init_thread_counter_for_init_list.begin()));
      assert(init_thread_counter_for_init_stmt);
      SgAssignOp* init_thread_counter_for_init_expr = isSgAssignOp(init_thread_counter_for_init_stmt->get_expression());
      assert(init_thread_counter_for_init_expr);
      ReplaceIndirectionArrays(init_thread_counter_for_init_expr->get_rhs_operand());
      SgLessThanOp* init_thread_counter_test_exp = isSgLessThanOp(init_thread_counter_test->get_expression());
      assert(init_thread_counter_test_exp);
      ReplaceIndirectionArrays(init_thread_counter_test_exp->get_rhs_operand());
    }

    //Statement to Populate bound arrays
    pair<SgStatement*,SgStatement*> bounds_store_stmts = curr_loop->GenerateCArrayStore();
    if( bounds_store_stmts.first ){
      SageInterface::insertStatementBefore(init_access,bounds_store_stmts.first);
      assert(isSgExprStatement(bounds_store_stmts.first));
      assert(isSgBinaryOp(isSgExprStatement(bounds_store_stmts.first)->get_expression()));
      SgExprStatement* lb_store_stmt = isSgExprStatement(bounds_store_stmts.first);
      SgBinaryOp* lb_store_exp = isSgBinaryOp(lb_store_stmt->get_expression());
      ReplaceIndirectionArrays(lb_store_exp->get_rhs_operand());
    }
    if( bounds_store_stmts.second ){
      SageInterface::insertStatementBefore(init_access,bounds_store_stmts.second);
      assert(isSgExprStatement(bounds_store_stmts.second));
      assert(isSgBinaryOp(isSgExprStatement(bounds_store_stmts.second)->get_expression()));
      ReplaceIndirectionArrays(isSgBinaryOp(isSgExprStatement(bounds_store_stmts.second)->get_expression())->get_rhs_operand());
    }
    SgForInitStatement* init_access_for_init = isSgForInitStatement(SageInterface::copyStatement(target_for->get_for_init_stmt()));
    SgExprStatement* init_access_test = isSgExprStatement(SageInterface::copyStatement(target_for->get_test()));
    assert(init_access_test);
    SgExpression* init_access_inc = SageInterface::copyExpression(target_for->get_increment());
    SgNullStatement* init_access_body = SageBuilder::buildNullStatement();
    SgForStatement* init_access_for = SageBuilder::buildForStatement(init_access_for_init,init_access_test,init_access_inc,init_access_body);
    SageInterface::insertStatementBefore(init_access,init_access_for);
    SgStatementPtrList& init_access_for_init_list = init_access_for_init->get_init_stmt();
    SgExprStatement* init_access_for_init_stmt = isSgExprStatement(*(init_access_for_init_list.begin()));
    assert(init_access_for_init_stmt);
    SgAssignOp* init_access_for_init_expr = isSgAssignOp(init_access_for_init_stmt->get_expression());
    assert(init_access_for_init_expr);
    ReplaceIndirectionArrays(init_access_for_init_expr->get_rhs_operand());
    SgLessThanOp* init_access_test_exp = isSgLessThanOp(init_access_test->get_expression());
    assert(init_access_test_exp);
    ReplaceIndirectionArrays(init_access_test_exp->get_rhs_operand());


    // curr_loop->SetInspectorLoops(graph_gen_for,init_counter_for,init_access_for);
    GenerateCStatement(target_for->get_loop_body(),graph_gen_body,init_counter_body,init_thread_counter_body,init_access_body,status_exp,counter_status_exp);
    
    if( init_thread_counter ){
      curr_loop->GenerateBodyCounterIncrement(init_thread_counter_body);
      curr_loop->GenerateLoopCounterIncrement(init_thread_counter);
    }
    else{
      curr_loop->GenerateBodyCounterIncrement(init_counter_body);
      curr_loop->GenerateLoopCounterIncrementAfter(init_counter_for);
    }
    curr_loop->GenerateBodyCounterIncrement(init_access_body);
    curr_loop->GenerateLoopCounterIncrement(init_access);
  }
  else
    GenerateCStatement(target_for->get_loop_body(),graph_gen,init_counter,init_thread_counter,init_access,status_exp,counter_status_exp);
}



void partitionable_loop::GenerateCExprStatement(SgExprStatement* target_stmt, SgStatement* graph_gen, SgStatement* init_counter, SgStatement* init_thread_counter,  SgStatement* init_access, SgExpression* status_exp, SgExpression* counter_status_exp)
{
  SgExpression* target_exp = target_stmt->get_expression();
  IAssignmentAttr* target_type = static_cast<IAssignmentAttr*>(target_stmt->getAttribute("assignment_type"));
  if( target_type->IsIAssignment()) {
    SgBinaryOp* target_binaryop = isSgBinaryOp(target_exp);
    SgVarRefExp* lhs_var = isSgVarRefExp(target_binaryop->get_lhs_operand());
    assert(lhs_var);
    assert(lhs_var->get_symbol()->getAttribute("indirection_scalar"));
    safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>(lhs_var->get_symbol()->getAttribute("indirection_scalar"));
    //Graph gen statement
    SgExpression* new_expr = SageInterface::copyExpression(target_stmt->get_expression());
    SgExprStatement* new_stmt = SageBuilder::buildExprStatement(new_expr);
    // SageInterface::insertStatementBefore(graph_gen,graph_gen_stmt);
    SgExpression* graph_gen_cond = GenerateCSafeStatement(target_binaryop->get_rhs_operand());
    if( graph_gen_cond ){
      SgStatement* graph_gen_else = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(status_exp),SageBuilder::buildIntVal(0)));
      SgIfStmt* graph_gen_stmt = SageBuilder::buildIfStmt(graph_gen_cond,new_stmt,graph_gen_else);
      SageInterface::insertStatementBefore(graph_gen,graph_gen_stmt);
      curr_safe_scalar->SetSafeScalarValues(isSgIfStmt(graph_gen_stmt));
    }
    else
      SageInterface::insertStatementBefore(graph_gen,new_stmt);
    ReplaceIndirectionArrays(new_expr);


    //init counter statement
    SgExpression* init_counter_expr = SageInterface::copyExpression(target_stmt->get_expression());
    SgExprStatement* new_init_counter_stmt = SageBuilder::buildExprStatement(init_counter_expr);
    SgExpression* init_counter_cond = GenerateCSafeStatement(target_binaryop->get_rhs_operand());
    if( init_counter_cond ){
      SgStatement* init_counter_else = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(counter_status_exp),SageBuilder::buildIntVal(0)));
      SgIfStmt* init_counter_stmt = SageBuilder::buildIfStmt(init_counter_cond,new_init_counter_stmt,init_counter_else);
      SageInterface::insertStatementBefore(init_counter,init_counter_stmt);
      curr_safe_scalar->SetSafeScalarValues(isSgIfStmt(init_counter_stmt));
    }
    else
      SageInterface::insertStatementBefore(init_counter,new_init_counter_stmt);
    ReplaceIndirectionArrays(init_counter_expr);

    //init thread counter statement
    if( init_thread_counter ){
      SgExpression* init_thread_counter_expr = SageInterface::copyExpression(target_stmt->get_expression());
      SgExprStatement* new_init_thread_counter_stmt = SageBuilder::buildExprStatement(init_thread_counter_expr);
      SageInterface::insertStatementBefore(init_thread_counter,new_init_thread_counter_stmt);
      ReplaceIndirectionArrays(init_thread_counter_expr);
    }

    //init access statement
    SgExpression* init_access_expr = SageInterface::copyExpression(target_stmt->get_expression());
    SgExprStatement* init_access_stmt = SageBuilder::buildExprStatement(init_access_expr);
    SageInterface::insertStatementBefore(init_access,init_access_stmt);
    ReplaceIndirectionArrays(init_access_expr);

  }
  else{
    if( isSgFunctionCallExp(target_exp) ){
      SgExprListExp* target_args = isSgFunctionCallExp(target_exp)->get_args();
      SgExpressionPtrList& target_arg_list = target_args->get_expressions();
    
      for( SgExpressionPtrList::const_iterator target_arg_iter = target_arg_list.begin(); target_arg_iter != target_arg_list.end() ; target_arg_iter++ ){
  	GenerateCExpression(*target_arg_iter,graph_gen,init_counter,init_thread_counter,init_access,status_exp,counter_status_exp);
      }    
    }
    else{
      GenerateCExpression(target_exp,graph_gen,init_counter,init_thread_counter,init_access,status_exp,counter_status_exp);
    }
  }
}



void partitionable_loop::GenerateCExpression(SgExpression* target_exp, SgStatement* graph_gen, SgStatement* init_counter, SgStatement* init_thread_counter, SgStatement* init_access, SgExpression* status_exp, SgExpression* counter_status_exp)
{
  if( isSgUnaryOp(target_exp) )
    GenerateCExpression(isSgUnaryOp(target_exp)->get_operand(),graph_gen,init_counter,init_thread_counter,init_access,status_exp,counter_status_exp);
  else if( isSgBinaryOp(target_exp) ){
    SgBinaryOp* target_binaryop = isSgBinaryOp(target_exp);
    VariantT curr_v = target_binaryop->variantT();
    if( curr_v == V_SgPntrArrRefExp ){
      SgExpression* target_lhs_exp = target_binaryop->get_lhs_operand();
      SgExpression* target_rhs_exp = target_binaryop->get_rhs_operand();
      while( isSgPntrArrRefExp(target_lhs_exp) ){
	target_rhs_exp = isSgBinaryOp(target_lhs_exp)->get_rhs_operand();
	target_lhs_exp = isSgBinaryOp(target_lhs_exp)->get_lhs_operand();
      }
      assert(isSgVarRefExp(target_lhs_exp));
      //graph generation
      SgVariableSymbol* target_array_var = isSgVarRefExp(target_lhs_exp)->get_symbol();
      if( target_array_var->attributeExists("data_array")) {
	array_details* target_array = static_cast<array_details*>(target_array_var->getAttribute("data_array"));
	assert(target_rhs_exp->attributeExists("access_details"));
	access_details* target_exp = static_cast<access_details*>(target_rhs_exp->getAttribute("access_details"));
	pair<SgStatement*,SgStatement*> new_stmts = target_exp->GenerateCArrayStore(target_array->GetArrayNum());
	if( new_stmts.first ){
	  SgExprStatement* add_pin_stmt = isSgExprStatement(new_stmts.first);
	  assert(add_pin_stmt);
	  SgExpression* graph_gen_cond = GenerateCSafeStatement(target_rhs_exp);
	  if( graph_gen_cond ){
	    SgStatement* graph_gen_else = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(status_exp),SageBuilder::buildIntVal(0)));
	    SgStatement* graph_gen_stmt = SageBuilder::buildIfStmt(graph_gen_cond,add_pin_stmt,graph_gen_else);
	    SageInterface::insertStatementBefore(graph_gen,graph_gen_stmt);
	  }
	  else
	    SageInterface::insertStatementBefore(graph_gen,add_pin_stmt);
	  SgFunctionCallExp* add_pin_fn_call = isSgFunctionCallExp(add_pin_stmt->get_expression());
	  assert(add_pin_fn_call);
	  SgExpressionPtrList&  add_pin_fn_args = add_pin_fn_call->get_args()->get_expressions();
	  ReplaceIndirectionArrays(*(++add_pin_fn_args.begin()));
	    
	  //init counter
	  SgExpression* init_counter_cond = GenerateCSafeStatement(target_rhs_exp);
	  if( init_counter_cond ){
	    SgStatement* init_counter_else = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageInterface::copyExpression(counter_status_exp),SageBuilder::buildIntVal(0)));
	    SgStatement* init_counter_stmt = SageBuilder::buildIfStmt(init_counter_cond,SageBuilder::buildNullStatement(),init_counter_else);
	    SageInterface::insertStatementBefore(init_counter,init_counter_stmt);
	  }
	}
	if( new_stmts.second ){
	  SageInterface::insertStatementBefore(init_access,new_stmts.second);
	  SgExprStatement* array_store_stmt = isSgExprStatement(new_stmts.second);
	  if (array_store_stmt){
	    ReplaceIndirectionArrays(isSgBinaryOp(array_store_stmt->get_expression())->get_rhs_operand());
	  }
	  else{
	    SgIfStmt* array_store_stmt_if = isSgIfStmt(new_stmts.second);
	    assert(array_store_stmt_if);
	    array_store_stmt = isSgExprStatement(array_store_stmt_if->get_true_body());
	    assert(array_store_stmt);
	    ReplaceIndirectionArrays(isSgBinaryOp(array_store_stmt->get_expression())->get_rhs_operand());
	    SgExprStatement* array_store_if_condn = isSgExprStatement(array_store_stmt_if->get_conditional());
	    assert(array_store_if_condn);
	    SgBinaryOp* array_store_if_condn_op = isSgBinaryOp(array_store_if_condn->get_expression());
	    ReplaceIndirectionArrays(array_store_if_condn_op->get_rhs_operand());
	  }
	}
      }
    }
    else{
      GenerateCExpression(target_binaryop->get_lhs_operand(),graph_gen,init_counter,init_thread_counter,init_access,status_exp,counter_status_exp);
      GenerateCExpression(target_binaryop->get_rhs_operand(),graph_gen,init_counter,init_thread_counter,init_access,status_exp,counter_status_exp);
    }
  }
}


SgExpression* partitionable_loop::GenerateCSafeStatement(SgExpression* target_exp)
{
  if( isSgUnaryOp(target_exp) ){
    return GenerateCSafeStatement(isSgUnaryOp(target_exp)->get_operand());
  }
  if( isSgBinaryOp(target_exp) ){
    SgBinaryOp* target_binaryop = isSgBinaryOp(target_exp);
    VariantT v_exp = target_binaryop->variantT();
    
    if( v_exp == V_SgPntrArrRefExp ){
      assert(is_known_fn);
      SgVarRefExp* target_indarray = isSgVarRefExp(target_binaryop->get_lhs_operand());
      assert(target_indarray);
      SgVariableSymbol* target_indarray_var = target_indarray->get_symbol();
      assert(target_indarray_var->attributeExists("indirection_array"));
      indarray_details* curr_indarray = static_cast<indarray_details*>(target_indarray_var->getAttribute("indirection_array"));
      SgIntVal* arg1 = SageBuilder::buildIntVal(curr_indarray->GetIndarrayNum());
      SgExpression* arg2 = SageInterface::copyExpression(target_binaryop->get_rhs_operand());
      ReplaceIndirectionArrays(arg2);
      SgExprListExp* args = SageBuilder::buildExprListExp(arg1,arg2);
      SgFunctionCallExp* condition = SageBuilder::buildFunctionCallExp(is_known_fn,args);
      //SgIfStmt* new_stmt = SageBuilder::buildIfStmt(condition,curr_stmt,SageBuilder::buildNullStatement());
      SgExpression* rhs_condn = GenerateCSafeStatement(target_binaryop->get_rhs_operand());
      if( rhs_condn )
  	return SageBuilder::buildBinaryExpression<SgAndOp>(rhs_condn,condition);
      else
  	return condition;
    }
    else{
      //SgExpression* lhs_stmt = GenerateCSafeStatement(target_binaryop->get_lhs_operand(),curr_stmt);
      SgExpression* lhs_condn = GenerateCSafeStatement(target_binaryop->get_lhs_operand());
      SgExpression* rhs_condn = GenerateCSafeStatement(target_binaryop->get_rhs_operand());
      if( lhs_condn ){
  	if( rhs_condn )
  	  return SageBuilder::buildBinaryExpression<SgAndOp>(lhs_condn,rhs_condn);
  	else
  	  return lhs_condn;
      }
      if( rhs_condn )
  	return rhs_condn;
      else
  	return NULL;
    }
  }
  else if( isSgVarRefExp(target_exp) ){
    SgVariableSymbol* safe_scalar_var = isSgVarRefExp(target_exp)->get_symbol();
    if (safe_scalar_var->attributeExists("indirection_scalar")){
      safe_scalar_details* curr_safe_scalar = static_cast<safe_scalar_details*>(safe_scalar_var->getAttribute("indirection_scalar"));
      SgNotEqualOp* condition = SageBuilder::buildBinaryExpression<SgNotEqualOp>(SageBuilder::buildVarRefExp(curr_safe_scalar->GetShadowScalar()),SageBuilder::buildIntVal(0));
      //SgIfStmt* new_stmt = SageBuilder::buildIfStmt(condition,curr_stmt,SageBuilder::buildNullStatement());
      return condition;
    }
    else 
      return NULL;
  }
  else
    return NULL;
}



void partitionable_loop::ReplaceIndirectionArrays(SgExpression* target_exp)
{
  if( target_exp->attributeExists("last_defined") )
    target_exp->removeAttribute("last_defined");
  if( isSgUnaryOp(target_exp) )
    ReplaceIndirectionArrays(isSgUnaryOp(target_exp)->get_operand());
  if( isSgBinaryOp(target_exp) ){
    SgBinaryOp* target_binaryop = isSgBinaryOp(target_exp);
    if( isSgPntrArrRefExp(target_exp) ){
      assert(target_binaryop);
      SgExpression* access_exp = target_binaryop->get_rhs_operand();
      ReplaceIndirectionArrays(target_binaryop->get_rhs_operand());
      SgVarRefExp* target_array = isSgVarRefExp(target_binaryop->get_lhs_operand());
      assert(target_array);
      assert(get_elem_fn);
      SgVariableSymbol* curr_indarray_var = target_array->get_symbol();
      assert(curr_indarray_var->attributeExists("indirection_array"));
      indarray_details* curr_indarray = static_cast<indarray_details*>(curr_indarray_var->getAttribute("indirection_array"));
      assert(curr_indarray);
      SgIntVal* arg1 = SageBuilder::buildIntVal(curr_indarray->GetIndarrayNum());
      SgExpression* arg2 = SageInterface::copyExpression(access_exp);
      SgFunctionCallExp* new_exp = SageBuilder::buildFunctionCallExp(get_elem_fn,SageBuilder::buildExprListExp(arg1,arg2));
      SageInterface::replaceExpression(target_exp,new_exp);    
    }
    else{
      ReplaceIndirectionArrays(target_binaryop->get_lhs_operand());
      ReplaceIndirectionArrays(target_binaryop->get_rhs_operand());
    }
  }
}



