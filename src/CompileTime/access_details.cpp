/*
 * access_details.cpp: This file is part of the IEC project.
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
 * @file: access_details.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include <cstdio>
#include <sstream>
#include "CompileTime/access_details.hpp"

using namespace std;

SgFunctionSymbol* access_details::add_pin_fn = NULL;
SgFunctionSymbol* access_details::renumber_offset_fn = NULL;
SgFunctionSymbol* access_details::renumber_const_offset_fn = NULL;
SgFunctionSymbol* access_details::renumber_fn = NULL;
int access_details::total_num = 0;


access_details::access_details(SgExpression* ae, const branch* ml, descriptor_type mt):
  access_exp(ae),
  my_branch(ml),
  access_type(mt),
  my_num(total_num++),
  inspector_gen(false),
  executor_expr(NULL)
{
  if( my_branch->IsLoop() ){
    const loop* my_loop = static_cast<const loop*>(ml);
    printf(" %s of type %d from loop %s ",access_exp->unparseToString().c_str(),mt,my_loop->GetIterator()->get_name().getString().c_str());
  }
  else{
    const conditional* my_branch = static_cast<const conditional*>(ml);
    printf("Found Access %s of type 1 from conditional\n",access_exp->unparseToString().c_str());
  }
  access_array = NULL;
  offset_var = NULL;
}




void access_details::InitAccessArray(SgStatement* insert_before, string array_name_string, SgFunctionSymbol* malloc_fn)
{
  assert(access_array == NULL );
  if( my_branch->IsLoop() ){
    const loop* my_loop = static_cast<const loop*>(my_branch);
    if( my_loop->IsPartitionableLoop() && access_type == OA_DESCRIPTOR )
      return;
  }
  //For now assume all arrays are ints
  SgPointerType* array_type = SageBuilder::buildPointerType(SageBuilder::buildIntType());
  ostringstream access_oss;

  //Malloc size exp
  SgVarRefExp* mult1;
  if( my_branch->IsLoop() ){
    const loop* my_loop = static_cast<const loop*>(my_branch);
    if( access_type == IA_DESCRIPTOR ){
      mult1 = SageBuilder::buildVarRefExp(my_loop->GetBodyCounter());
      access_oss << "__i_" ;
    }
    else{
      mult1 = SageBuilder::buildVarRefExp(my_loop->GetLoopCounter());
      access_oss << "__o_";
    }
  }
  else{
    const conditional* my_condn = static_cast<const conditional*>(my_branch);
    mult1 = SageBuilder::buildVarRefExp(my_condn->GetCondnCounter());
    if( my_condn->IsThen() )
      access_oss << "__t_" ;
    else
      access_oss << "__e_" ;
  }

  SgSizeOfOp* mult2 = SageBuilder::buildSizeOfOp(SageBuilder::buildIntType());
  SgMultiplyOp* sizeof_exp = SageBuilder::buildBinaryExpression<SgMultiplyOp>(isSgExpression(mult1),isSgExpression(mult2));
  
  //malloc fn call
  SgExprListExp* malloc_args = SageBuilder::buildExprListExp(isSgExpression(sizeof_exp));
  SgFunctionCallExp* malloc_fn_call = SageBuilder::buildFunctionCallExp(malloc_fn,malloc_args);

  access_oss << array_name_string << "_" ;
  if( my_branch->IsLoop() ){
    const loop* my_loop = static_cast<const loop*>(my_branch);
    access_oss <<  my_loop->GetIterator()->get_name().getString() << "_" << my_num << "_";
  }
  else
    access_oss << my_num << "_" ;

  SgName array_name(access_oss.str());

  SgCastExp* rhs_cast_exp = SageBuilder::buildCastExp(isSgExpression(malloc_fn_call),array_type);
  SgAssignInitializer* array_init_exp = SageBuilder::buildAssignInitializer(isSgExpression(rhs_cast_exp));
  SgVariableDeclaration* array_init_stmt = SageBuilder::buildVariableDeclaration(array_name,array_type,array_init_exp);
  
  SageInterface::insertStatementBefore(insert_before,array_init_stmt);
  access_array = array_init_stmt;
  
  if( access_type == OA_DESCRIPTOR) {
    access_oss << "_offset" ;
    SgName offset_var_name(access_oss.str());
    offset_var = SageBuilder::buildVariableDeclaration(offset_var_name,SageBuilder::buildIntType());
    SageInterface::insertStatementBefore(insert_before,offset_var);
  }
  
}


void access_details::DeleteAccessArray(SgStatement* insert_before, SgFunctionSymbol* free_fn) const 
{
  if(access_array != NULL ){
    SgVarRefExp* free_arg1 = SageBuilder::buildVarRefExp(access_array);
    SgExprListExp* free_args = SageBuilder::buildExprListExp(free_arg1);
    SgFunctionCallExp* free_call = SageBuilder::buildFunctionCallExp(free_fn,free_args);
    SgExprStatement* free_stmt = SageBuilder::buildExprStatement(free_call);
    SageInterface::insertStatementBefore(insert_before,free_stmt);
  }
}

pair<SgStatement*,SgStatement*> access_details::GenerateCArrayStore( const int array_num) 
{
  pair<SgStatement*,SgStatement*> ret_stmts;
  ret_stmts.first = NULL; ret_stmts.second = NULL;
  if( !inspector_gen ){
  SgStatement* graph_gen_stmt = NULL;
    assert(add_pin_fn);
    //To add a pin 
    SgIntVal* graph_gen_arg1 = SageBuilder::buildIntVal(array_num);
    SgExpression* graph_gen_arg2 = SageInterface::copyExpression(access_exp);
    //SgExpression* graph_gen_arg2 = static_cast<SgExpression*>(copy_helper.copyAst(access_exp));
    SgIntVal* graph_gen_arg3;
    SgIntVal* graph_gen_arg4;
   
    if( access_type == OA_DESCRIPTOR ){
      assert( my_branch->IsLoop() );
      const loop* my_loop = static_cast<const loop*>(my_branch);
      graph_gen_arg3 = SageBuilder::buildIntVal(1);
      if( my_loop->IsPartitionableLoop() )
	graph_gen_arg4 = SageBuilder::buildIntVal(1);
      else
	graph_gen_arg4 = SageBuilder::buildIntVal(0);
    }
    else{
      graph_gen_arg3 = SageBuilder::buildIntVal(0);
      graph_gen_arg4 = SageBuilder::buildIntVal(0);
    }
    SgExprListExp* graph_gen_args = SageBuilder::buildExprListExp(graph_gen_arg1,graph_gen_arg2,graph_gen_arg3,graph_gen_arg4);
    SgFunctionCallExp* graph_gen_fn_call = SageBuilder::buildFunctionCallExp(add_pin_fn,graph_gen_args);
    graph_gen_stmt = SageBuilder::buildExprStatement(graph_gen_fn_call);
    ret_stmts.first = graph_gen_stmt;
    //SageInterface::insertStatementBefore(graph_gen,graph_gen_stmt);

    //For Init access
    if( access_array ){
      SgVarRefExp* init_access_index;
      if( my_branch->IsLoop() ){
	const loop* my_loop = static_cast<const loop*>(my_branch);
	if( access_type == OA_DESCRIPTOR )
	  init_access_index = SageBuilder::buildVarRefExp(my_loop->GetLoopCounter());
	else
	  init_access_index = SageBuilder::buildVarRefExp(my_loop->GetBodyCounter());
      }
      else{
	const conditional* my_condn = static_cast<const conditional*>(my_branch);
	init_access_index = SageBuilder::buildVarRefExp(my_condn->GetCondnCounter());
      }

      SgVarRefExp* init_access_array = SageBuilder::buildVarRefExp(access_array);
      SgPntrArrRefExp* init_access_array_ref = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(init_access_array,init_access_index);
      SgExpression* init_access_value = SageInterface::copyExpression(access_exp); //static_cast<SgExpression*>(copy_helper.copyAst(access_exp)); 
      SgAssignOp* init_access_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(init_access_array_ref,init_access_value);
      SgExprStatement* init_access_stmt = SageBuilder::buildExprStatement(init_access_exp);

      if( access_type == OA_DESCRIPTOR ){
	const loop* my_loop = static_cast<const loop*>(my_branch);
	SgIfStmt* init_access_firstiter = my_loop->GenerateIfFirstIterStmt(init_access_stmt);
	assert(init_access_firstiter);
	ret_stmts.second = init_access_firstiter;
	//SageInterface::insertStatementBefore(init_access,init_access_firstiter);
      }
      else
	ret_stmts.second = init_access_stmt;
      //SageInterface::insertStatementBefore(init_access,init_access_stmt);
    }
    inspector_gen = true;
  }
  return ret_stmts;
}


void access_details::RenumberAccess(SgStatement* insert_before, const int array_num) const 
{
  if( access_array ){
    SgName thread_id_name("__thread_id__");
    SgVarRefExp* thread_id_var = SageBuilder::buildVarRefExp(thread_id_name);
    SgIntVal* array_num_ref = SageBuilder::buildIntVal(array_num);
    SgVarRefExp* array_ref = SageBuilder::buildVarRefExp(access_array);
    SgExprStatement* renumber_stmt;
    if( access_type == OA_DESCRIPTOR ){
      const loop* my_loop = static_cast<const loop*>(my_branch);
      SgVarRefExp* access_size = SageBuilder::buildVarRefExp(my_loop->GetLoopCounter());
      SgExpression* lb_array = SageBuilder::buildVarRefExp(my_loop->GetLBArray());
      assert( lb_array != NULL );
      SgExprListExp* all_args = SageBuilder::buildExprListExp(thread_id_var,array_num_ref,access_size,array_ref,lb_array);
      assert(renumber_offset_fn);
      SgFunctionCallExp* renumber_fn_call = SageBuilder::buildFunctionCallExp(renumber_offset_fn,all_args);
      renumber_stmt = SageBuilder::buildExprStatement(renumber_fn_call);
    }
    else{
      SgVarRefExp* access_size ;
      if( my_branch->IsLoop() ){
	const loop* my_loop = static_cast<const loop*>(my_branch);
	access_size = SageBuilder::buildVarRefExp(my_loop->GetBodyCounter());
      }
      else{
	const conditional* my_condn = static_cast<const conditional*>(my_branch);
	access_size = SageBuilder::buildVarRefExp(my_condn->GetCondnCounter());
      }
      SgExprListExp* all_args = SageBuilder::buildExprListExp(thread_id_var,array_num_ref,access_size,array_ref);
      SgFunctionCallExp* renumber_fn_call = SageBuilder::buildFunctionCallExp(renumber_fn,all_args);
      renumber_stmt = SageBuilder::buildExprStatement(renumber_fn_call);
    }
    SageInterface::insertStatementBefore(insert_before,renumber_stmt);
  }
}


SgExpression* access_details::GenerateExecutorExpression() 
{
  SgExpression* new_access;
  new_access = NULL;
  if( executor_expr == NULL ){
    if( access_array ){
      if( access_type == OA_DESCRIPTOR){
	assert(offset_var);
	assert(my_branch->IsLoop());
	const loop* my_loop = static_cast<const loop*>(my_branch);
	SgVarRefExp* loop_iterator = SageBuilder::buildVarRefExp(my_loop->GetIterator());
	executor_expr =  SageBuilder::buildBinaryExpression<SgAddOp>(SageBuilder::buildVarRefExp(offset_var),loop_iterator);
	SgVarRefExp* access_counter = SageBuilder::buildVarRefExp(my_loop->GetLoopCounter());
	SgPntrArrRefExp* access_var_init_rhs = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(access_array),access_counter);
	SgAssignOp* access_var_init = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(offset_var),access_var_init_rhs);
	SgExprStatement* access_var_init_stmt = SageBuilder::buildExprStatement(access_var_init);
	my_loop->InsertBeforeExecutor(access_var_init_stmt);
      }
      else{
	SgVarRefExp* access_counter;
	if( my_branch->IsLoop() ){
	  const loop* my_loop = static_cast<const loop*>(my_branch);
	  if( my_loop->GetDepth() == 0 )
	    access_counter = SageBuilder::buildVarRefExp(my_loop->GetIterator());
	  else
	    access_counter = SageBuilder::buildVarRefExp(my_loop->GetBodyCounter());
	}
	else{
	  const conditional* my_condn = static_cast<const conditional*>(my_branch);
	  access_counter = SageBuilder::buildVarRefExp(my_condn->GetCondnCounter());
	}
	executor_expr = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(access_array),access_counter);
      }
      new_access = executor_expr;
    }
  }
  else
    if( access_array )
      new_access = SageInterface::copyExpression(executor_expr);
  return new_access;
}

