/*
 * detector.cpp: This file is part of the IEC project.
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
 * @file: detector.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/loops.hpp"
#include "CompileTime/conditional.hpp"
#include "CompileTime/driver.hpp"

using namespace std;

extern bool CompareExpression(const SgExpression* , const SgExpression*);

partitionable_loop::partitionable_loop(SgForStatement* of):
  loop(of,0),
  write_communication(false),
  read_communication(false)
{
  local_block_start = NULL;
  local_block_end = NULL;
  nlocal_iters = NULL;
  body_counter = NULL;
  SetUsed();
}


bool partitionable_loop::CheckProperties()
{
  deque<branch*> curr_scope;
  curr_scope.push_back(this);
  iterators.insert(orig_iterator);
  bool check_property = true;
  vector<SgVarRefExp*> lb_params = SageInterface::querySubTree<SgVarRefExp>(lb_exp);
  for( vector<SgVarRefExp*>::iterator it = lb_params.begin(); it != lb_params.end()  ; it++ )
    if( iterators.find((*it)->get_symbol()) != iterators.end() ){
      check_property = false;
      break;
    }
    else
      parameters.insert((*it)->get_symbol());
  if( !check_property )
    return false;
  vector<SgVarRefExp*> ub_params = SageInterface::querySubTree<SgVarRefExp>(ub_exp);
  for( vector<SgVarRefExp*>::iterator it = ub_params.begin(); it != ub_params.end()  ; it++ )
    if( iterators.find((*it)->get_symbol()) != iterators.end() ){
      check_property = false;
      break;
    }
    else
      parameters.insert((*it)->get_symbol());
  if( !check_property )
    return false;

  if( !GetIScalarArrays(orig_for->get_loop_body(),curr_scope) )
    return false;

  //The worklist to complete the list
  int num_iscalars;
  int num_iarrays;
  do{
    num_iscalars = indirection_scalars.size();
    num_iarrays = indirection_arrays.size();
    
    for( deque<SgExprStatement*>::iterator it = assignment_stmts.begin() ; it != assignment_stmts.end() ; it++ ){
      SgBinaryOp* stmt_exp = isSgBinaryOp((*it)->get_expression());
      if( stmt_exp){
	scope* stmt_scope = static_cast<scope*>((*it)->getAttribute("scope"));
	SgExpression* lhs_exp = stmt_exp->get_lhs_operand();
	if( isSgPntrArrRefExp(lhs_exp) ){
	  SgExpression* lhs_array = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
	  while(isSgPntrArrRefExp(lhs_array) )
	    lhs_array = isSgPntrArrRefExp(lhs_array)->get_lhs_operand();
	  SgVarRefExp* lhs_array_var = isSgVarRefExp(lhs_array);
	  assert(lhs_array_var);
	  if( indirection_arrays.find(lhs_array_var->get_symbol()) != indirection_arrays.end() ){
	    printf("Changing Value of indirection array in %s\n",(*it)->unparseToString().c_str());
	    //exit(1);
	    check_property = false;
	    break;
	  }
	}
	else if( isSgVarRefExp(lhs_exp) ){
	  if( indirection_scalars.find(isSgVarRefExp(lhs_exp)->get_symbol()) != indirection_scalars.end() ){
	    deque<SgExprStatement*>::const_iterator kt;
	    for( kt = i_assignment_stmts.begin() ; kt != i_assignment_stmts.end() ; kt++ )
	      if( *kt == *it )
		break;
	    if( kt == i_assignment_stmts.end() ){
	      deque<branch*> const& statement_stack = stmt_scope->getScope();
	      for( deque<branch*>::const_iterator jt = statement_stack.begin() ; jt != statement_stack.end() ; jt++) 
		(*jt)->AddWriteIScalar(isSgVarRefExp(lhs_exp)->get_symbol());
	      (*it)->setAttribute("assignment_type",new IAssignmentAttr(true));
	      i_assignment_stmts.push_back((*it));
	      check_property = AddIArraysScalars(stmt_exp->get_rhs_operand(),stmt_scope->getScope());
	    }
	  }
	}
      }
    }
  }while( check_property && ( indirection_arrays.size() != num_iarrays || indirection_scalars.size() != num_iscalars ) );

  if( !check_property )
    return false;
  
  // printf("Indirection Arrays :");
  // for( set<SgVariableSymbol*>::const_iterator it = indirection_arrays.begin() ; it != indirection_arrays.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nIndirection Scalars :");
  // for( set<SgVariableSymbol*>::const_iterator it = indirection_scalars.begin() ; it != indirection_scalars.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\n");

  for( deque<SgExprStatement*>::iterator it = i_assignment_stmts.begin() ; it != i_assignment_stmts.end() ; it++ ){
    deque<SgExprStatement*>::iterator jt;
    for( jt = assignment_stmts.begin() ; jt != assignment_stmts.end() ; jt++ )
      if( *it == *jt ){
	assignment_stmts.erase(jt);
	break;
      }
    SgAssignOp* curr_assign = isSgAssignOp((*it)->get_expression());
    assert(curr_assign);
    SgVarRefExp* curr_lhs = isSgVarRefExp(curr_assign->get_lhs_operand());
    assert(curr_lhs);
    written_iscalars.insert(curr_lhs->get_symbol());
  }
  for( deque<SgExprStatement*>::iterator it = assignment_stmts.begin() ; it != assignment_stmts.end() ; it++ )
    (*it)->setAttribute("assignment_type",new IAssignmentAttr(false));


  //Check for "privatizability" property of the IScalars
  set<SgVariableSymbol*> defined_iscalars;
  set<SgVariableSymbol*> defined_dscalars;
  defined_iscalars.insert(indirection_scalars.begin(),indirection_scalars.end());
  for( set<SgVariableSymbol*>::iterator it = written_iscalars.begin() ; it != written_iscalars.end() ; it++ ){
    defined_iscalars.erase(*it);
  }
  //defined_iscalars.erase(written_iscalars.begin(),written_iscalars.end());

  if( !CheckPrivatizable(orig_for->get_loop_body(),defined_iscalars,defined_dscalars) )
    return false;

  privatizable.insert(defined_dscalars.begin(),defined_dscalars.end());

  for( set<SgVariableSymbol*>::iterator it = indirection_scalars.begin() ; it != indirection_scalars.end() ; it++ )
    if( written_iscalars.find(*it) == written_iscalars.end() )
      parameters.insert(*it);
  for( set<SgVariableSymbol*>::iterator it = parameters.begin() ; it != parameters.end() ; it++ )
    indirection_scalars.erase(*it);
  
  //Finally check that the inner dimension of "data arrays" use only those loops-iterators that are replicated or parametes
  if( !CheckInnerDimensions() )
    return false;

  //At this stage indirection_scalars and written_iscalars must have the same thing

  // printf("Read Data Arrays :");
  // for( set<SgVariableSymbol*>::const_iterator it = read_data_arrays.begin() ; it != read_data_arrays.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nWrite Data Arrays :");
  // for( set<SgVariableSymbol*>::const_iterator it = write_data_arrays.begin() ; it != write_data_arrays.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nLocal Read Arrays :");
  // for( set<SgVariableSymbol*>::const_iterator it = local_read_arrays.begin() ; it != local_read_arrays.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nLocal Write Arrays :");
  // for( set<SgVariableSymbol*>::const_iterator it = local_write_arrays.begin() ; it != local_write_arrays.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nRead Data Scalars :");
  // for( set<SgVariableSymbol*>::const_iterator it = read_data_scalars.begin() ; it != read_data_scalars.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nWrite Data Scalars :");
  // for( set<SgVariableSymbol*>::const_iterator it = write_data_scalars.begin() ; it != write_data_scalars.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\nPrivatizable Data Scalars :");
  // for( set<SgVariableSymbol*>::const_iterator it = defined_dscalars.begin() ; it != defined_dscalars.end() ; it++ )
  //   printf(" %s",(*it)->get_name().getString().c_str());
  // printf("\nIndirection Scalars :");
  // for( set<SgVariableSymbol*>::const_iterator it = written_iscalars.begin() ; it != written_iscalars.end() ; it++ )
  //   printf(" %s",(*it)->get_name().getString().c_str());
  // printf("\nParameters :");
  // for( set<SgVariableSymbol*>::const_iterator it = parameters.begin() ; it != parameters.end() ; it++ )
  //   printf(" %s",(*it)->get_name().getString().c_str());
  // printf("\nIterators :");
  // for( set<SgVariableSymbol*>::const_iterator it = iterators.begin() ; it != iterators.end() ; it++ )
  //   printf(" %s",(*it)->get_name().getString().c_str());
  // printf("\n");
  
  return true;
}


bool partitionable_loop::GetIScalarArrays(SgStatement* curr_stmt, deque<branch*>& curr_scope)
{
  VariantT curr_variant  =  curr_stmt->variantT();
  bool check_property = true;

  switch(curr_variant){
  case V_SgBasicBlock:
    check_property = GetIScalarArraysBasicBlock(isSgBasicBlock(curr_stmt),curr_scope);
    break;
  case V_SgForStatement:
    check_property = GetIScalarArraysFor(isSgForStatement(curr_stmt),curr_scope);
    break;
  case V_SgExprStatement:
    check_property = GetIScalarArraysExprStatement(isSgExprStatement(curr_stmt),curr_scope);
    break;
  case V_SgIfStmt:
    check_property = GetIScalarArraysIf(isSgIfStmt(curr_stmt),curr_scope);
    break;
  case V_SgVariableDeclaration:
    printf("Error: Cannot Handle Variable Declerations\n");
    check_property = false;
    break;
  default:
    printf("Unhandled statement type in Loop\n");
    check_property = false;
    break;
  }

  return check_property;
}


bool partitionable_loop::GetIScalarArraysBasicBlock(SgBasicBlock* curr_bb, deque<branch*>& curr_scope)
{
  SgStatementPtrList& curr_list = curr_bb->get_statements();
  SgStatementPtrList::iterator curr_iter;

  bool check_property = true;
  for( curr_iter = curr_list.begin() ; curr_iter != curr_list.end() ; curr_iter++ ){
    if( !GetIScalarArrays(*curr_iter,curr_scope) ){
      check_property = false;
      break;
    }
  }
  return check_property;
}


bool partitionable_loop::GetIScalarArraysFor(SgForStatement* curr_for, deque<branch*>& curr_scope)
{
  inner_loop* new_inner_loop = new inner_loop(curr_for,curr_scope.size());
  curr_for->setAttribute("inner_loop",new_inner_loop);
  curr_for->setAttribute("scope",new scope(curr_scope));
  bool lb_used = AddIArraysScalars(new_inner_loop->GetLBExp(),curr_scope);
  bool ub_used = AddIArraysScalars(new_inner_loop->GetUBExp(),curr_scope);
  if( !lb_used || !ub_used )
    return false;
  iterators.insert(new_inner_loop->GetIterator());

  curr_scope.push_back(new_inner_loop);
  bool check_property = GetIScalarArrays(curr_for->get_loop_body(),curr_scope);
  curr_scope.pop_back();
  return check_property;
}


bool partitionable_loop::GetIScalarArraysIf(SgIfStmt* curr_if, deque<branch*>& curr_scope)
{
  curr_if->setAttribute("scope",new scope(curr_scope));
  SgExprStatement* curr_condn = isSgExprStatement(curr_if->get_conditional());
  assert(curr_condn);
  if( !AddIArraysScalars(curr_condn->get_expression(),curr_scope) )
    return false;
  then_branch* new_then_stmt = new then_branch(curr_if);
  curr_if->setAttribute("then_branch",new_then_stmt);
  curr_scope.push_back(new_then_stmt);
  if( !GetIScalarArrays(curr_if->get_true_body(),curr_scope) )
    return false;
  curr_scope.pop_back();
  if( curr_if->get_false_body() != NULL ){
    else_branch* new_else_stmt = new else_branch(curr_if);
    curr_if->setAttribute("else_branch",new_else_stmt);
    curr_scope.push_back(new_else_stmt);
    bool check_property = GetIScalarArrays(curr_if->get_false_body(),curr_scope);
    curr_scope.pop_back();
    return check_property;
  }
  else
    return true;
}
		
bool partitionable_loop::GetIScalarArraysExprStatement(SgExprStatement* curr_stmt, deque<branch*>& curr_scope)
{
  scope* scope_attr = new scope(curr_scope);
  curr_stmt->setAttribute("scope",scope_attr);
  bool check_property = true;
  
  SgExpression* curr_exp = curr_stmt->get_expression();
  if( isSgFunctionCallExp(curr_exp) ){
    SgExprListExp* curr_args = isSgFunctionCallExp(curr_exp)->get_args();
    SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
    assignment_stmts.push_back(curr_stmt);    
    for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ ){
      if( !GetIScalarArraysExpression(*curr_arg_iter,curr_scope) ){
	check_property = false;
	break;
      }
    }    
  }
  else if( isSgBinaryOp(curr_exp) ){
    VariantT curr_variant = curr_exp->variantT();
    
    SgExpression* lhs_exp;
    deque<branch*>::iterator it;    
    loop* outer_loop;
    SgVariableSymbol* lhs_symbol;

    switch(curr_variant){
    case V_SgAssignOp:
    case V_SgAndAssignOp:
    case V_SgDivAssignOp:
    case V_SgLshiftAssignOp:
    case V_SgMinusAssignOp:
    case V_SgMultAssignOp:
    case V_SgPlusAssignOp:
    case V_SgRshiftAssignOp:
      lhs_exp = isSgBinaryOp(curr_exp)->get_lhs_operand();
      if( isSgVarRefExp(lhs_exp) ){
	lhs_symbol = isSgVarRefExp(lhs_exp)->get_symbol();
	for( it = curr_scope.begin() ; it != curr_scope.end() ; it++ )
	  if( (*it)->IsLoop() ){
	    outer_loop = static_cast<loop*>(*it);
	    if( outer_loop->GetIterator() == lhs_symbol ){
	      printf("Iterator assigned within the loop at %s\n",curr_stmt->unparseToString().c_str());
	      check_property = false;
	      break;
	    }
	  }
      }
      else if( !isSgPntrArrRefExp(lhs_exp) ){
	printf("Unhandled expression on LHS in %s\n",curr_stmt->unparseToString().c_str());
	check_property = false;
      }
      if( check_property ){
	assignment_stmts.push_back(curr_stmt);
	check_property = GetIScalarArraysExpression(curr_exp,curr_scope);
      }
      break;
    default:
      printf("Uhandled Expression Statement : %s\n",curr_stmt->unparseToString().c_str());
      check_property = false;
    };
  }
  else{
    printf("Unhandled Statement Type : %s\n",curr_stmt->unparseToString().c_str());
    check_property = false;
  }
  return check_property;
}												 


bool partitionable_loop::GetIScalarArraysExpression(SgExpression* curr_exp, deque<branch*>& curr_scope)
{
  bool check_property = true;
  if( isSgBinaryOp(curr_exp) ){
    SgExpression* lhs_exp = isSgBinaryOp(curr_exp)->get_lhs_operand();
    SgExpression* rhs_exp = isSgBinaryOp(curr_exp)->get_rhs_operand();
    
    if( isSgDotExp(curr_exp) ){
      printf("Use of Structs in %s not handled presently\n",curr_exp->unparseToString().c_str());
      check_property = false;
    }
    else if( isSgPntrArrRefExp(curr_exp) ){
      //Indirection used only for the outermost dimension (C-style). 
      //Inner dimensions are affine=, verified by scop_extractor.
      while( isSgPntrArrRefExp(lhs_exp) ){
	rhs_exp = isSgPntrArrRefExp(lhs_exp)->get_rhs_operand();
	lhs_exp = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
      }
      assert(isSgVarRefExp(lhs_exp));
      //Ideally the else of this conditional must break the property of partitionable loops
      // Essentially, this is ignoring use of arrays not specified in the arrays pragma.
      // This is a way to handle privitazable arrays. NOT RECOMMENDED!!! But is a hack for now
      if( driver::CheckArray(isSgVarRefExp(lhs_exp))){
	check_property = AddIArraysScalars(rhs_exp,curr_scope);
      }
    }
    else{
      check_property = GetIScalarArraysExpression(lhs_exp,curr_scope) && GetIScalarArraysExpression(rhs_exp,curr_scope);
    }
  }
  else if( isSgUnaryOp(curr_exp) ){
    check_property = GetIScalarArraysExpression(isSgUnaryOp(curr_exp)->get_operand(),curr_scope);
  }
  else if( isSgFunctionCallExp(curr_exp) ){
    //The Scop Extractor ensures this is a pure math function. All expressions are read-only
    SgExprListExp* curr_args = isSgFunctionCallExp(curr_exp)->get_args();
    SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
    
    for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ ){
      if( !GetIScalarArraysExpression(*curr_arg_iter,curr_scope) ){
	check_property = false;
	break;
      }
    }    
  }
  else{
    if( !isSgVarRefExp(curr_exp) && !isSgValueExp(curr_exp) )
      check_property = false;
  }
  return check_property;
}


bool partitionable_loop::AddIArraysScalars(SgExpression* curr_exp, const deque<branch*>& curr_scope)
{
  bool check_property = true;
  if( isSgVarRefExp(curr_exp) ){
    SgVariableSymbol* curr_symbol = isSgVarRefExp(curr_exp)->get_symbol();
    bool is_iterator = false;
    for( deque<branch*>::const_iterator it = curr_scope.begin() ; it != curr_scope.end() ; it++ )
      if( (*it)->IsLoop() ){
	const loop* outer_loop = static_cast<const loop*>(*it);
	if( outer_loop->GetIterator() == curr_symbol ){
	  is_iterator = true;
	  break;
	}
      }
    if( !is_iterator )
      indirection_scalars.insert(curr_symbol);
  }
  else if( isSgBinaryOp(curr_exp) ){
    if( isSgPntrArrRefExp(curr_exp) ){
      SgExpression* array_var = isSgPntrArrRefExp(curr_exp)->get_lhs_operand();
      if( !isSgVarRefExp(array_var) ){
	printf("Cant handle multi-dimensional indirection arrays : %s\n",curr_exp->unparseToString().c_str());
	check_property = false;
      }
      if( check_property ){
	indirection_arrays.insert(isSgVarRefExp(array_var)->get_symbol());
	check_property = AddIArraysScalars(isSgBinaryOp(curr_exp)->get_rhs_operand(),curr_scope);
      }
    }
    else if( isSgCompoundAssignOp(curr_exp) || isSgAssignOp(curr_exp) || isSgDotExp(curr_exp) ){
      printf("Cant handle this <IExpr> : %s\n",curr_exp->unparseToString().c_str());
      check_property = false;
    }
    else{
      bool is_lhs = AddIArraysScalars(isSgBinaryOp(curr_exp)->get_lhs_operand(),curr_scope);
      bool is_rhs = AddIArraysScalars(isSgBinaryOp(curr_exp)->get_rhs_operand(),curr_scope);
      check_property = is_lhs && is_rhs;
    }
  }
  else if( isSgUnaryOp(curr_exp) ){
    check_property = AddIArraysScalars(isSgUnaryOp(curr_exp)->get_operand(),curr_scope);
  }
  else{
    if( !isSgValueExp(curr_exp) )
      check_property = false;
  }
  return check_property;
}



bool partitionable_loop::CheckPrivatizable(SgStatement* curr_stmt, set<SgVariableSymbol*>& defined_iscalars, set<SgVariableSymbol*>& defined_dscalars)
{
  VariantT curr_variant  =  curr_stmt->variantT();
  bool check_property = true;

  switch(curr_variant){
  case V_SgBasicBlock:
    check_property = CheckPrivatizableBasicBlock(isSgBasicBlock(curr_stmt),defined_iscalars,defined_dscalars);
    break;
  case V_SgForStatement:
    check_property = CheckPrivatizableFor(isSgForStatement(curr_stmt),defined_iscalars,defined_dscalars);
    break;
  case V_SgExprStatement:
    check_property = CheckPrivatizableExprStatement(isSgExprStatement(curr_stmt),defined_iscalars,defined_dscalars);
    break;
  case V_SgIfStmt:
    check_property = CheckPrivatizableIf(isSgIfStmt(curr_stmt),defined_iscalars,defined_dscalars);
    break;
  default:
    printf("Unhandled statement type in Loop\n");
    check_property = false;
    break;
  }
  return check_property;
}

bool partitionable_loop::CheckPrivatizableBasicBlock(SgBasicBlock* curr_bb, set<SgVariableSymbol*>& defined_iscalars, set<SgVariableSymbol*>& defined_dscalars)
{
  SgStatementPtrList& curr_list = curr_bb->get_statements();
  SgStatementPtrList::iterator curr_iter;
  bool check_property = true;

  for( curr_iter = curr_list.begin() ; curr_iter != curr_list.end() ; curr_iter++ ){
    if( !CheckPrivatizable(*curr_iter,defined_iscalars,defined_dscalars) ){
      check_property = false;
      break;
    }
  }
  return check_property;
}


bool partitionable_loop::CheckPrivatizableFor(SgForStatement* curr_for, set<SgVariableSymbol*>& defined_iscalars, set<SgVariableSymbol*>& defined_dscalars)
{
  set<SgVariableSymbol*> body_iscalars;
  set<SgVariableSymbol*> body_dscalars;

  inner_loop* curr_inner_loop = static_cast<inner_loop*>(curr_for->getAttribute("inner_loop"));
  scope* curr_scope = static_cast<scope*>(curr_for->getAttribute("scope"));
  if( !CheckVariables(curr_inner_loop->GetLBExp(),defined_iscalars,curr_scope) || !CheckVariables(curr_inner_loop->GetUBExp(),defined_iscalars,curr_scope) )
    return false;
  
  body_iscalars.insert(defined_iscalars.begin(),defined_iscalars.end());
  body_dscalars.insert(defined_dscalars.begin(),defined_dscalars.end());
  bool check_property = CheckPrivatizable(curr_for->get_loop_body(),body_iscalars,body_dscalars);
  // if( !curr_inner_loop->IsUsed() )
  //   printf("[IE-Debug] : Unused loop : %s\n",curr_for->unparseToString().c_str());
  return check_property;
}

bool partitionable_loop::CheckPrivatizableIf(SgIfStmt* curr_if, set<SgVariableSymbol*>& defined_iscalars, set<SgVariableSymbol*>& defined_dscalars)
{
  scope* curr_scope = static_cast<scope*>(curr_if->getAttribute("scope"));
  if( !CheckVariables(isSgExprStatement(curr_if->get_conditional())->get_expression(),defined_iscalars,curr_scope) )
    return false;

  set<SgVariableSymbol*> then_iscalars;
  set<SgVariableSymbol*> else_iscalars;
  set<SgVariableSymbol*> then_dscalars;
  set<SgVariableSymbol*> else_dscalars;

  then_iscalars.insert(defined_iscalars.begin(),defined_iscalars.end());
  then_dscalars.insert(defined_dscalars.begin(),defined_dscalars.end());
  if( !CheckPrivatizable(curr_if->get_true_body(),then_iscalars,then_dscalars) )
    return false;

  for( set<SgVariableSymbol*>::iterator it = defined_iscalars.begin() ; it != defined_iscalars.end() ; it++ )
    then_iscalars.erase(*it);
  for( set<SgVariableSymbol*>::iterator it = defined_dscalars.begin() ; it != defined_dscalars.end() ; it++ )
    then_dscalars.erase(*it);

  
  if( curr_if->get_false_body() ){
    //else_iscalars.insert(defined_iscalars.begin(),defined_iscalars.end());
    for( set<SgVariableSymbol*>::iterator it = defined_iscalars.begin() ; it != defined_iscalars.end() ; it++ ){
      SgVariableSymbol* curr_symbol = (*it);
      else_iscalars.insert(curr_symbol);
    }
    else_dscalars.insert(defined_dscalars.begin(),defined_dscalars.end());
    if( !CheckPrivatizable(curr_if->get_false_body(),else_iscalars,else_dscalars) )
      return false;
    for( set<SgVariableSymbol*>::iterator it = defined_iscalars.begin() ; it != defined_iscalars.end() ; it++ )
      else_iscalars.erase(*it);
    for( set<SgVariableSymbol*>::iterator it = defined_dscalars.begin() ; it != defined_dscalars.end() ; it++ )
      else_dscalars.erase(*it);
  }

  for( set<SgVariableSymbol*>::iterator it = then_iscalars.begin() ; it != then_iscalars.end() ; it++ ){
    set<SgVariableSymbol*>::iterator jt = else_iscalars.find(*it);
    if( jt != else_iscalars.end() ){
      defined_iscalars.insert(*it);
    }
  }
  for( set<SgVariableSymbol*>::iterator it = then_dscalars.begin() ; it != then_dscalars.end() ; it++ )
    if( else_dscalars.find(*it) != else_dscalars.end() )
      defined_dscalars.insert(*it);
  return true;
}


bool partitionable_loop::CheckPrivatizableExprStatement(SgExprStatement* curr_stmt, set<SgVariableSymbol*>& defined_iscalars, set<SgVariableSymbol*>& defined_dscalars)
{
  scope* curr_scope = static_cast<scope*>(curr_stmt->getAttribute("scope"));
  SgExpression* curr_exp = curr_stmt->get_expression();
  bool check_property = true;
  
  if( isSgFunctionCallExp(curr_exp) ){
    SgExprListExp* curr_args = isSgFunctionCallExp(curr_exp)->get_args();
    SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
    
    for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ ){
      SgExpression* curr_arg_exp = *curr_arg_iter;
      if( isSgVarRefExp(curr_arg_exp ) ){
	//All function arguements are considered to be read only;
	check_property = CheckPrivatizableExpression(curr_arg_exp,defined_iscalars,defined_dscalars,curr_scope);
      }
      else if( !isSgValueExp(curr_arg_exp) )
	check_property = false;
      if( !check_property )
	break;
    }    
  }
  else if( isSgBinaryOp(curr_exp) ){
    SgBinaryOp* curr_assign = isSgBinaryOp(curr_exp);
    assert(curr_assign);
    IAssignmentAttr* curr_type = static_cast<IAssignmentAttr*>(curr_stmt->getAttribute("assignment_type"));
    if( curr_type->IsIAssignment() ){
      if( CheckVariables(curr_assign->get_rhs_operand(),defined_iscalars,curr_scope) ){
	SgVarRefExp* curr_lhs = isSgVarRefExp(curr_assign->get_lhs_operand());
	assert(curr_lhs);
	SgVariableSymbol* new_defined_iscalar = curr_lhs->get_symbol();
	defined_iscalars.insert(new_defined_iscalar);
	deque<branch*>::const_iterator it = curr_scope->getScope().begin()++;
	for( ; it != curr_scope->getScope().end() ; it++ ){
	  written_iscalars.insert(curr_lhs->get_symbol());
	}
      }
      else
	check_property = false;
    }
    else{
      SgExpression* lhs_exp = curr_assign->get_lhs_operand();
      if ( isSgVarRefExp(lhs_exp) ){
	SgVariableSymbol* curr_symbol = isSgVarRefExp(lhs_exp)->get_symbol();
	if( isSgAssignOp(curr_assign) ){
	  if( CheckPrivatizableExpression(curr_assign->get_rhs_operand(),defined_iscalars,defined_dscalars,curr_scope) )
	    defined_dscalars.insert(curr_symbol);
	  else
	    check_property = false;
	}
	else{
	  if ( CheckPrivatizableExpression(curr_assign->get_rhs_operand(),defined_iscalars,defined_dscalars,curr_scope) ){
	    if( defined_dscalars.find(curr_symbol) == defined_dscalars.end() ){
	      if( read_data_scalars.find(curr_symbol) != read_data_scalars.end() ){
		printf("Unhandled dependence due to update of previously read data scalar %s\n",curr_symbol->get_name().getString().c_str());
		check_property = false;
	      }
	      else
		write_data_scalars.insert(curr_symbol);
	    }
	  }
	  else
	    check_property = false;
	}
      }
      else if( isSgPntrArrRefExp(lhs_exp) ){
	SgExpression* rhs_exp = NULL;
	while( isSgPntrArrRefExp(lhs_exp) ){
	  rhs_exp = isSgPntrArrRefExp(lhs_exp)->get_rhs_operand();
	  lhs_exp = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
	}
	assert(rhs_exp && isSgVarRefExp(lhs_exp));
	//This is again a hack to handle privatizable arrays. THis is NOT RECOMMEDED!!
	if( CheckPrivatizableExpression(curr_assign->get_rhs_operand(),defined_iscalars,defined_dscalars,curr_scope) ){
	  if( driver::CheckArray(isSgVarRefExp(lhs_exp)) ){
	    if( CheckVariables(rhs_exp,defined_iscalars,curr_scope) ){
	      SgVariableSymbol* curr_array = isSgVarRefExp(lhs_exp)->get_symbol();
	      if( isSgVarRefExp(rhs_exp) && isSgVarRefExp(rhs_exp)->get_symbol() == orig_iterator ){
		//This is special case for allowing read and write to direct access from partitionable loop
		// For now no offset allowed. 
		if( read_data_arrays.find(curr_array) != read_data_arrays.end() ){
		  printf("Unhandled dependence due to local write to a previously read array %s\n",curr_array->get_name().getString().c_str());
		  check_property = false;
		}
		else
		  local_write_arrays.insert(curr_array);
	      }
	      else{
		if( read_data_arrays.find(curr_array) != read_data_arrays.end() || local_read_arrays.find(curr_array) != local_read_arrays.end() ) {
		  printf("Unhandled dependence due to write to an array %s that previously read from\n",curr_array->get_name().getString().c_str());
		  check_property = false;
		}
		else
		  write_data_arrays.insert(curr_array);
	      }
	    }
	    else
	      check_property = false;
	  }
	}
	else
	  check_property = false;
      }
      else{
	printf("Unhandled expression on LHS of %s\n",curr_stmt->unparseToString().c_str());
	check_property = false;
      }
    }
  }
  else{
    printf("Unhandled expression statement %s\n",curr_stmt->unparseToString().c_str());
    check_property = false;
  }
  return check_property;
}


//This function evaluates only right-hand side expressions or function call arguements
bool partitionable_loop::CheckPrivatizableExpression(SgExpression* curr_exp, set<SgVariableSymbol*>& defined_iscalars, set<SgVariableSymbol*>& defined_dscalars, scope* curr_scope)
{
  bool check_property = true;
  if( isSgVarRefExp(curr_exp) ){
    SgVariableSymbol* curr_symbol = isSgVarRefExp(curr_exp)->get_symbol();
    if( !isSgPointerType(curr_symbol->get_type()) ){
      //Hack to allow privatizable arrays being passed into math functions
      if( isSgArrayType(curr_symbol->get_type()) ){
	if( driver::CheckArray(isSgVarRefExp(curr_exp)) )
	  check_property = false;
      }
      else{
	if( defined_dscalars.find(curr_symbol) == defined_dscalars.end()) {
	  if( write_data_scalars.find(curr_symbol) != write_data_scalars.end() ){
	    printf("Unhandled dependendence due to read from previously updated data scalar %s\n",curr_symbol->get_name().getString().c_str());
	    check_property = false;
	  }
	  else
	    read_data_scalars.insert(curr_symbol);
	}
      }
    }
    else
      check_property = false;
  }
  else if( isSgBinaryOp(curr_exp) ){
    SgExpression* lhs_exp = isSgBinaryOp(curr_exp)->get_lhs_operand();
    SgExpression* rhs_exp = isSgBinaryOp(curr_exp)->get_rhs_operand();
    
    if( isSgPntrArrRefExp(curr_exp) ){
      while( isSgPntrArrRefExp(lhs_exp) ){
	rhs_exp = isSgPntrArrRefExp(lhs_exp)->get_rhs_operand();
	lhs_exp = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
      }
      assert(isSgVarRefExp(lhs_exp) );
      //Hack to handle privatizable arrays.
      if( driver::CheckArray(isSgVarRefExp(lhs_exp)) ){
	if( CheckVariables(rhs_exp,defined_iscalars,curr_scope) ){
	  SgVariableSymbol* curr_array = isSgVarRefExp(lhs_exp)->get_symbol();
	  if( write_data_arrays.find(curr_array) != write_data_arrays.end() ){
	    printf("Unhandled dependence due to array read from data array %s that was previously written to\n",curr_array->get_name().getString().c_str());
	    check_property = false;
	  }
	  else{
	    if( isSgVarRefExp(rhs_exp) && isSgVarRefExp(rhs_exp)->get_symbol() == orig_iterator ){
	      local_read_arrays.insert(curr_array);
	    }
	    else if( local_write_arrays.find(curr_array) != local_write_arrays.end() ){
	      printf("Unhandled dependence due to array read from data array %s that was previously written to locally\n",curr_array->get_name().getString().c_str());
	      check_property = false;
	    }
	    else{
	      read_data_arrays.insert(curr_array);
	    }
	  }
	}
	else
	  check_property = false;
      }
    }
    else{
      check_property = CheckPrivatizableExpression(lhs_exp,defined_iscalars,defined_dscalars,curr_scope) &&
	CheckPrivatizableExpression(rhs_exp,defined_iscalars,defined_dscalars,curr_scope);
    }
  }
  else if( isSgUnaryOp(curr_exp) ){
    check_property = CheckPrivatizableExpression(isSgUnaryOp(curr_exp)->get_operand(),defined_iscalars,defined_dscalars,curr_scope);
  }
  else if( isSgFunctionCallExp(curr_exp) ){
    SgExprListExp* curr_args = isSgFunctionCallExp(curr_exp)->get_args();
    SgExpressionPtrList& curr_arg_list = curr_args->get_expressions();
    bool check_property = true;

    for( SgExpressionPtrList::const_iterator curr_arg_iter = curr_arg_list.begin(); curr_arg_iter != curr_arg_list.end() ; curr_arg_iter++ ){
      if( !CheckPrivatizableExpression(*curr_arg_iter,defined_iscalars,defined_dscalars,curr_scope) ){
	check_property = false;
	break;
      }
    }
  }
  else{
    if(!isSgValueExp(curr_exp) )
      check_property = false;
  }
  return check_property;
}


//This function checks that iscalars used are read only or privatizable
bool partitionable_loop::CheckVariables(SgExpression* curr_exp, set<SgVariableSymbol*>& defined_iscalars, scope* curr_scope)
{
  bool check_property = true;
  if( isSgVarRefExp(curr_exp) ){
    SgVariableSymbol* curr_symbol = isSgVarRefExp(curr_exp)->get_symbol();
    bool is_iterator = false;
    deque<branch*> scope_stack = curr_scope->getScope();
    for( deque<branch*>::iterator it = scope_stack.begin() ; it != scope_stack.end() ; it++ )
      if( (*it)->IsLoop() ){
	loop* curr_loop = static_cast<loop*>(*it);
	if( curr_loop->GetIterator() == curr_symbol ){
	  is_iterator = true;
	  curr_loop->SetUsed();
	  break;
	}
      }
    if( !is_iterator ){
      set<SgVariableSymbol*>::iterator posn = defined_iscalars.find(curr_symbol);
      if( posn == defined_iscalars.end() ){
	printf("Use of IScalar %s doesnt come from within the loop\n",curr_exp->unparseToString().c_str());
	check_property = false;
      }
    }
  }  
  else if( isSgBinaryOp(curr_exp) ){
    SgExpression* lhs_exp = isSgBinaryOp(curr_exp)->get_lhs_operand();
    SgExpression* rhs_exp = isSgBinaryOp(curr_exp)->get_rhs_operand();
    if( isSgPntrArrRefExp(curr_exp) ){
      check_property = CheckVariables(rhs_exp,defined_iscalars,curr_scope);
    }
    else{
      check_property = CheckVariables(lhs_exp,defined_iscalars,curr_scope) && CheckVariables(rhs_exp,defined_iscalars,curr_scope);
    }
  }
  else if( isSgUnaryOp(curr_exp) ){
    check_property = CheckVariables(isSgUnaryOp(curr_exp)->get_operand(),defined_iscalars,curr_scope);
  }
  else if( !isSgValueExp(curr_exp) )
    check_property = false;
  return check_property;
}

bool partitionable_loop::CheckInnerDimensions() const
{
  //A loop is replicated or is "not used" when it is not referenced within the leading dimension expression of 
  //an array reference and loop bounds are only parameters.
  vector<SgForStatement*> all_inner = SageInterface::querySubTree<SgForStatement>(orig_for->get_loop_body());
  set<SgVariableSymbol*> non_used_iterators;
  for( vector<SgForStatement*>::iterator it = all_inner.begin() ; it != all_inner.end() ; it++ ){
    SgForStatement* curr_inner_loop = *it;
    assert(curr_inner_loop->attributeExists("inner_loop"));
    inner_loop* curr_inner = static_cast<inner_loop*>(curr_inner_loop->getAttribute("inner_loop"));
    if( !curr_inner->IsUsed() ){
      bool is_inner_used = false;
      vector<SgVarRefExp*> lb_var_refs = SageInterface::querySubTree<SgVarRefExp>(curr_inner->GetLBExp());
      for( vector<SgVarRefExp*>::iterator it = lb_var_refs.begin() ; it != lb_var_refs.end() ; it++ )
	if( parameters.find((*it)->get_symbol()) == parameters.end() ){
	  is_inner_used = true;
	  break;
	}
      if( is_inner_used )
	curr_inner->SetUsed();
      else{
	vector<SgVarRefExp*> ub_var_refs = SageInterface::querySubTree<SgVarRefExp>(curr_inner->GetUBExp());
	for( vector<SgVarRefExp*>::iterator it = ub_var_refs.begin() ; it != ub_var_refs.end() ; it++ )
	  if( parameters.find((*it)->get_symbol()) == parameters.end() ){
	    is_inner_used = true;
	    break;
	  }
      }
      if( is_inner_used )
	curr_inner->SetUsed();
    }
  }

  // printf("Non-Used iterators: ");
  // for( set<SgVariableSymbol*>::iterator it = non_used_iterators.begin() ; it != non_used_iterators.end() ; it++ )
  //   printf(" %s,",(*it)->get_name().getString().c_str());
  // printf("\n");
  //Now ensure that all inner-dimensions of array expression refers only parameters and "unused" loop-iterators
  vector<SgPntrArrRefExp*> all_array_refs = SageInterface::querySubTree<SgPntrArrRefExp>(orig_for->get_loop_body());
  bool check_property = true;
  for( vector<SgPntrArrRefExp*>::iterator it = all_array_refs.begin() ; check_property && it != all_array_refs.end() ; it++ ){
    SgPntrArrRefExp* curr_exp = (*it);
    SgExprStatement* curr_stmt = SageInterface::getEnclosingNode<SgExprStatement>(curr_exp);
    if( !curr_stmt->attributeExists("scope") ){
      //printf("[IE-Debug] Statement %s has no scope, shouldnt happen\n",curr_stmt->unparseToString().c_str());
      //assert(0);
      break;
    }
    scope* curr_scope = static_cast<scope*>(curr_stmt->getAttribute("scope"));
    SgExpression* lhs_exp = curr_exp->get_lhs_operand();
    while( check_property &&  isSgPntrArrRefExp(lhs_exp) ){
      check_property = CheckInnerExpression(curr_exp->get_rhs_operand(),curr_scope);
      curr_exp = isSgPntrArrRefExp(lhs_exp);
      lhs_exp = curr_exp->get_lhs_operand();
    }
    assert(isSgVarRefExp(lhs_exp));
    if( check_property && !driver::CheckArray(isSgVarRefExp(lhs_exp) ) ){
      check_property = CheckInnerExpression(curr_exp->get_rhs_operand(),curr_scope);
    }
  }
  return check_property;
}


bool partitionable_loop::CheckInnerExpression(SgExpression* curr_exp, scope* curr_scope) const
{
  bool check_property = true;
  vector<SgVarRefExp*> all_vars = SageInterface::querySubTree<SgVarRefExp>(curr_exp);
  for( vector<SgVarRefExp*>::iterator jt = all_vars.begin() ; check_property && jt != all_vars.end() ; jt++ ){
    SgVariableSymbol* curr_symbol = (*jt)->get_symbol();
    if( parameters.find(curr_symbol) == parameters.end() ){
      deque<branch*>::const_iterator kt = curr_scope->getScope().begin() ;
      for( ; kt != curr_scope->getScope().end() ; kt++ ){
	if( (*kt)->IsLoop()) {
	  const loop* curr_loop = static_cast<const loop*>(*kt);
	  SgVariableSymbol* compare_iterator = curr_loop->GetIterator();
	  if( curr_loop->GetIterator() == curr_symbol) {
	    if( curr_loop->IsUsed() ){
	      check_property = false;
	      printf("Inner Dimension of %s expression uses a non-replicated inner-loop. Currently unhandled\n",curr_exp->unparseToString().c_str());
	    }
	    break;
	  }
	}
      }
      if( kt == curr_scope->getScope().end() ){
	printf("Inner Dimension of %s expression uses neither a replicated loop iterator or a parameter. Currently unhandled\n",curr_exp->unparseToString().c_str());
	check_property = false;
      }
    }
  }
  return check_property;
}
