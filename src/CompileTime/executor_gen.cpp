/*
 * executor_gen.cpp: This file is part of the IEC project.
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
 * @file: executor_gen.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/loops.hpp"
#include "CompileTime/conditional.hpp"
#include "CompileTime/array_details.hpp"
#include "CompileTime/access_details.hpp"
#include "CompileTime/driver.hpp"

using namespace std;


void partitionable_loop::GenerateCExecutor(SgStatement* insert_before)
{
	ResetAllCounters(insert_before);
  
// 	SgName thread_id_name("__thread_id__");
	if( read_communication ){
		SgStatement* comm_read_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_reads_fn,SageBuilder::buildExprListExp(/*SageBuilder::buildVarRefExp(thread_id_name),*/SageBuilder::buildIntVal(my_num))));
		SageInterface::insertStatementBefore(insert_before,comm_read_stmt);
	}
  
	SgIntVal* new_lb = SageBuilder::buildIntVal(0);
	SgAssignOp*  new_lb_exp = SageBuilder::buildBinaryExpression<SgAssignOp>( SageBuilder::buildVarRefExp(orig_iterator), new_lb);
	SgStatement* new_lb_stmt = SageBuilder::buildExprStatement(new_lb_exp);

	assert(nlocal_iters);
	SgVarRefExp* new_ub = SageBuilder::buildVarRefExp(nlocal_iters);
	SgLessThanOp* new_ub_exp = SageBuilder::buildBinaryExpression<SgLessThanOp>( SageBuilder::buildVarRefExp(orig_iterator), new_ub);
	SgStatement* new_ub_stmt = SageBuilder::buildExprStatement(new_ub_exp);

	SgExpression* new_inc_exp = SageBuilder::buildUnaryExpression<SgPlusPlusOp>( SageBuilder::buildVarRefExp(orig_iterator) );
  
	SgNullStatement* new_body = SageBuilder::buildNullStatement();

	executor_for = SageBuilder::buildForStatement(new_lb_stmt,new_ub_stmt,new_inc_exp,new_body);
  
	SageInterface::insertStatementBefore(insert_before,executor_for);

	GenerateCExecutorStatement(orig_for->get_loop_body(),new_body);

	if( write_communication ){
		SgStatement* init_ghosts_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(init_write_ghosts_fn,SageBuilder::buildExprListExp(/*SageBuilder::buildVarRefExp(thread_id_name),*/SageBuilder::buildIntVal(my_num))));
		SageInterface::insertStatementBefore(executor_for,init_ghosts_stmt);
		SgStatement* comm_write_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_writes_fn,SageBuilder::buildExprListExp(/*SageBuilder::buildVarRefExp(thread_id_name),*/SageBuilder::buildIntVal(my_num))));
		SageInterface::insertStatementAfter(executor_for,comm_write_stmt);
	}    
  
	for( set<SgVariableSymbol*>::iterator it = write_data_scalars.begin() ; it != write_data_scalars.end() ; it++ ){
		SgAddressOfOp* reduce_fn_arg = SageBuilder::buildUnaryExpression<SgAddressOfOp>(SageBuilder::buildVarRefExp(*it));
		SgFunctionCallExp* reduce_fn_exp = SageBuilder::buildFunctionCallExp(reduce_scalar_fn,SageBuilder::buildExprListExp(/*SageBuilder::buildVarRefExp(thread_id_name),*/reduce_fn_arg));
		SageInterface::insertStatementBefore(insert_before,SageBuilder::buildExprStatement(reduce_fn_exp));
	}
}


void partitionable_loop::GenerateCExecutorStatement(SgStatement* curr_stmt, SgStatement* insert_before)
{
	VariantT v_stmt = curr_stmt->variantT();

	switch(v_stmt){
	case V_SgBasicBlock:
		GenerateCExecutorBasicBlock(isSgBasicBlock(curr_stmt),insert_before);
		break;
	case V_SgForStatement:
		GenerateCExecutorForStatement(isSgForStatement(curr_stmt),insert_before);
		break;
	case V_SgExprStatement:
		GenerateCExecutorExprStatement(isSgExprStatement(curr_stmt),insert_before);
		break;
	case V_SgIfStmt:
		GenerateCExecutorIfStatement(isSgIfStmt(curr_stmt),insert_before);
		break;
	default:
		;
		break;
	}
}


void partitionable_loop::GenerateCExecutorBasicBlock(SgBasicBlock* curr_bb, SgStatement* insert_before)
{
	SgStatementPtrList& curr_list = curr_bb->get_statements();
	SgStatementPtrList::const_iterator curr_iter;

	for( curr_iter = curr_list.begin() ; curr_iter != curr_list.end() ; curr_iter++ ){
		GenerateCExecutorStatement(*curr_iter,insert_before);
	}
}


void partitionable_loop::GenerateCExecutorForStatement(SgForStatement* curr_for, SgStatement* insert_before)
{
	inner_loop* curr_inner = static_cast<inner_loop*>(curr_for->getAttribute("inner_loop"));
	if( curr_inner->IsUsed() ){
		//lower bound
		SgVarRefExp* init_lhs = SageBuilder::buildVarRefExp(curr_inner->GetIterator());
		SgPntrArrRefExp* init_rhs = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(curr_inner->GetLBArray()),SageBuilder::buildVarRefExp(curr_inner->GetLoopCounter()) );
		SgAssignOp* init_op = SageBuilder::buildBinaryExpression<SgAssignOp>(init_lhs,init_rhs);
		SgExprStatement* init_stmt = SageBuilder::buildExprStatement(init_op);
  
		//upper bound
		SgVarRefExp* test_lhs = SageBuilder::buildVarRefExp(curr_inner->GetIterator());
		SgPntrArrRefExp* test_rhs = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(curr_inner->GetUBArray()),SageBuilder::buildVarRefExp(curr_inner->GetLoopCounter()));
		SgLessThanOp* test_exp = SageBuilder::buildBinaryExpression<SgLessThanOp>(test_lhs,test_rhs);
		SgExprStatement* test_stmt = SageBuilder::buildExprStatement(test_exp);
  
		//test expressions
		SgExpression* inc_exp = SageBuilder::buildUnaryExpression<SgPlusPlusOp>(SageBuilder::buildVarRefExp(curr_inner->GetIterator()));
    
		SgNullStatement* new_body = SageBuilder::buildNullStatement();
		SgForStatement* new_for = SageBuilder::buildForStatement(init_stmt,test_stmt,inc_exp,new_body);
		curr_inner->SetExecutorFor(new_for);
		SageInterface::insertStatementBefore(insert_before,new_for);
		GenerateCExecutorStatement(curr_for->get_loop_body(),new_body);

		curr_inner->GenerateLoopCounterIncrement(insert_before);
		curr_inner->GenerateBodyCounterIncrement(new_body);
	}
	else{
		//lower bound
		SgForInitStatement* new_lb = static_cast<SgForInitStatement*>(SageInterface::copyStatement(curr_for->get_for_init_stmt()));
		SgStatement* new_ub = SageInterface::copyStatement(curr_for->get_test());
		SgExpression* inc_exp = SageInterface::copyExpression(curr_for->get_increment());
		SgNullStatement* new_body = SageBuilder::buildNullStatement();
		SgForStatement* new_for = SageBuilder::buildForStatement(new_lb,new_ub,inc_exp,new_body);
		SageInterface::insertStatementBefore(insert_before,new_for);
		GenerateCExecutorStatement(curr_for->get_loop_body(),new_body);
	}
}


void partitionable_loop::GenerateCExecutorIfStatement( SgIfStmt* curr_if, SgStatement* insert_before)
{
	then_branch* curr_then = static_cast<then_branch*>(curr_if->getAttribute("then_branch"));
	SgExpression* new_condn = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(curr_then->GetCondnArray()),SageBuilder::buildVarRefExp(curr_then->GetIfCounter()));
	SgNullStatement* new_then = SageBuilder::buildNullStatement();
	SgNullStatement* new_else = SageBuilder::buildNullStatement();

	SgIfStmt* new_if = SageBuilder::buildIfStmt(SageBuilder::buildExprStatement(new_condn),new_then,new_else);
	SageInterface::insertStatementBefore(insert_before,new_if);

	GenerateCExecutorStatement(curr_if->get_true_body(),new_then);
  
	curr_then->GenerateIfCounterIncrement(insert_before);
	curr_then->GenerateBranchCounterIncrement(new_then);

	if( curr_if->attributeExists("else_branch")){
		else_branch* curr_else = static_cast<else_branch*>(curr_if->getAttribute("else_branch"));
		GenerateCExecutorStatement(curr_if->get_false_body(),new_else);
		curr_else->GenerateBranchCounterIncrement(new_else);
	}  
}


void partitionable_loop::GenerateCExecutorExprStatement(SgExprStatement* curr_stmt, SgStatement* insert_before)
{
	IAssignmentAttr* curr_type = static_cast<IAssignmentAttr*>(curr_stmt->getAttribute("assignment_type"));
	if( !curr_type->IsIAssignment() ){
		SgExprStatement* new_stmt = SageBuilder::buildExprStatement(SageInterface::copyExpression(curr_stmt->get_expression()));
		SageInterface::insertStatementBefore(insert_before,new_stmt);
		GenerateCExecutorExpression(new_stmt->get_expression(),curr_stmt->get_expression());
	}
}


void partitionable_loop::GenerateCExecutorExpression(SgExpression* new_exp , SgExpression* curr_exp)
{
	if( isSgBinaryOp(curr_exp) ){
		if( isSgPntrArrRefExp(curr_exp) ){
			SgPntrArrRefExp* curr_arrayref = isSgPntrArrRefExp(curr_exp);
			SgExpression* new_array_ref = isSgPntrArrRefExp(new_exp)->get_lhs_operand();
			SgExpression* new_access_exp = isSgPntrArrRefExp(new_exp)->get_rhs_operand();
			SgExpression* array_ref = curr_arrayref->get_lhs_operand();
			SgExpression* access_exp = curr_arrayref->get_rhs_operand();
			while( isSgPntrArrRefExp(array_ref) ){
				access_exp = isSgPntrArrRefExp(array_ref)->get_rhs_operand();
				new_access_exp = isSgPntrArrRefExp(new_array_ref)->get_rhs_operand();
				array_ref = isSgPntrArrRefExp(array_ref)->get_lhs_operand();
				new_array_ref = isSgPntrArrRefExp(new_array_ref)->get_lhs_operand();
			}
			SgVariableSymbol* array_var = isSgVarRefExp(array_ref)->get_symbol();
			if(array_var->attributeExists("data_array") ){
				array_details* curr_data_array = static_cast<array_details*>(array_var->getAttribute("data_array"));
				assert(access_exp->attributeExists("access_details"));
				access_details* curr_access = static_cast<access_details*>(access_exp->getAttribute("access_details"));

				SageInterface::replaceExpression(new_array_ref,curr_data_array->GetLocalArrayRef());

				SgExpression* new_access = curr_access->GenerateExecutorExpression();
				if( new_access )
					SageInterface::replaceExpression(new_access_exp,new_access);
			}
		}
		else{
			SgBinaryOp* curr_binaryop = isSgBinaryOp(curr_exp);
			SgBinaryOp* new_binaryop = isSgBinaryOp(new_exp);
			GenerateCExecutorExpression(new_binaryop->get_lhs_operand(),curr_binaryop->get_lhs_operand());
			GenerateCExecutorExpression(new_binaryop->get_rhs_operand(),curr_binaryop->get_rhs_operand());
		}
	}
	else if( isSgUnaryOp(curr_exp) ){
		GenerateCExecutorExpression(isSgUnaryOp(new_exp)->get_operand(),isSgUnaryOp(curr_exp)->get_operand());
	}
	else if( isSgFunctionCallExp(curr_exp)  ){
		SgExprListExp* curr_args = isSgFunctionCallExp(curr_exp)->get_args();
		SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
		SgExprListExp* new_args = isSgFunctionCallExp(new_exp)->get_args();
		SgExpressionPtrList& new_arg_list = new_args->get_expressions();
		SgExpressionPtrList::const_iterator new_arg_iter = new_arg_list.begin();

		for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ , new_arg_iter++){
			GenerateCExecutorExpression(*new_arg_iter,*curr_arg_iter);
		}    
	}
}
