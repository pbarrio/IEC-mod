/*
 * array_details.cpp: This file is part of the IEC project.
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
 * @file: array_details.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/array_details.hpp"
#include "CompileTime/loops.hpp"
#include <boost/regex.hpp>
#include <cassert>
#include <map>
#include "CompileTime/access_details.hpp"

using namespace std;

extern bool CompareExpression(const SgExpression*, const SgExpression* );
extern bool CompareAccessExpression(const SgExpression*, const SgExpression*);

int array_details::num_arrays = 0;
SgFunctionSymbol* array_details::malloc_fn = NULL;
SgFunctionSymbol* array_details::malloc_2d_double_fn = NULL;
SgFunctionSymbol* array_details::malloc_2d_float_fn = NULL;
SgFunctionSymbol* array_details::free_fn = NULL;
SgFunctionSymbol* array_details::free_2d_double_fn = NULL;
SgFunctionSymbol* array_details::free_2d_float_fn = NULL;
SgFunctionSymbol* array_details::get_size_fn = NULL;
SgFunctionSymbol* array_details::populate_fn = NULL;
SgFunctionSymbol* array_details::set_array_stride_fn = NULL;
SgFunctionSymbol* array_details::communicate_reads_for_fn = NULL;
SgFunctionSymbol* array_details::communicate_writes_for_fn = NULL;


array_string::array_string(string pragma_details)
{
  boost::regex array_details_regex("(\\w*)[^\\[]*\\[([^\\]]*)\\][^\\[]*\\[([^\\]]*)\\](.*)");// \\s*\\[([^\\]]*)\\]\\[([^\\]]*)\\]\\s*
  boost::cmatch array_details_match;
  // printf("String: %s\n",pragma_details.c_str());
  if( boost::regex_match(pragma_details.c_str(),array_details_match,array_details_regex)){
      array_name.assign(array_details_match[1].first,array_details_match[1].second);
      array_size.assign(array_details_match[2].first,array_details_match[2].second);  
      array_stride.assign(array_details_match[3].first,array_details_match[3].second);  
      printf("Array Details : %s %s %s\n",array_name.c_str(),array_size.c_str(),array_stride.c_str());
  }
  else
    assert(0);
}

array_string::array_string()
{
}

array_details::array_details(const SgVariableSymbol* sym, array_string* st):
  array_symbol(sym),
  array_name(st->array_name),
  size_string(st->array_size),
  stride_string(st->array_stride),
  my_num(num_arrays)
{
  num_arrays++;
  last_access = NULL;
}


array_details::~array_details()
{
  for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ )
    delete (*it);
  all_accesses.clear();
  write_loops.clear();
  read_loops.clear();
  local_writes.clear();
}


void array_details::GenerateStrideInfo(SgStatement* insert_before)
{
  if( stride_string.compare(" 1 ") != 0 ){
    assert(set_array_stride_fn);
    SgExpression* arg2;
    if( isdigit(stride_string.c_str()[0]) ){
      arg2 = SageBuilder::buildIntVal(atoi(stride_string.c_str()));
    }
    else{
      SgName stride_name(stride_string);
      arg2 = SageBuilder::buildVarRefExp(stride_name);
    }
    SgExprListExp* args = SageBuilder::buildExprListExp(SageBuilder::buildIntVal(my_num),arg2);
    SgExprStatement* stride_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(set_array_stride_fn,args));
    SageInterface::insertStatementBefore(insert_before,stride_stmt);
  }
}


access_details* array_details::AddAccess( SgExpression* access_exp , branch* curr_counter, descriptor_type curr_type, bool is_lhs, partitionable_loop* curr_outer, VariantT assign_type )
{
  access_details* ret_access = NULL;
  for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ )
    if( CompareAccessExpression((*it)->GetExpression(),access_exp) ){
      printf("Found duplicate access ");
      ret_access = (*it);
      break;
    }
  if( ret_access == NULL ){
    printf("Found new access ");
    access_details* new_access = new access_details(access_exp,curr_counter,curr_type);
    all_accesses.push_back(new_access);
    ret_access = new_access;
    CheckConflicts(new_access);

    //No communication needed for OA_DESCRIPTOR access from partitionable loops
    if( new_access->my_branch->IsLoop() &&  new_access->my_branch == curr_outer && new_access->access_type == OA_DESCRIPTOR ){
      if( is_lhs ){
	if( local_writes.size() != 0 ){
	  partitionable_loop* last_local = *(local_writes.rbegin() );
	  if( last_local != static_cast<partitionable_loop*>(curr_outer) )
	    local_writes.push_back(static_cast<partitionable_loop*>(curr_outer));
	}
	else
	  local_writes.push_back(static_cast<partitionable_loop*>(curr_outer));
      }
    }
    else{
      if( is_lhs ){
	if( read_loops.size() != 0 ){
	  partitionable_loop* last_read = *(read_loops.rbegin());
	  assert(last_read != static_cast<partitionable_loop*>(curr_outer));
	}
	if( write_loops.size() != 0 ){
	  pair<partitionable_loop*,VariantT> last_write = *(write_loops.rbegin());
	  if( last_write.first != static_cast<partitionable_loop*>(curr_outer))
	    write_loops.push_back(pair<partitionable_loop*,VariantT>(static_cast<partitionable_loop*>(curr_outer),assign_type));
	  else{
	    if( last_write.second != assign_type ){
	      printf("Array %s is updated with different operators in loop %s : Serializing dependence\n",array_name.c_str(),curr_outer->unparseToString().c_str() );
	      exit(1);
	    }
	  }
	}
	else
	  write_loops.push_back(pair<partitionable_loop*,VariantT>(static_cast<partitionable_loop*>(curr_outer),assign_type));
      }
      else{
	if( write_loops.size() != 0 ){
	  partitionable_loop* last_write = (*(write_loops.rbegin())).first;
	  assert(last_write != static_cast<partitionable_loop*>(curr_outer));
	}
	if( read_loops.size() != 0 ){
	  partitionable_loop* last_read = *(read_loops.rbegin());
	  if( last_read != static_cast<partitionable_loop*>(curr_outer))
	    read_loops.push_back(static_cast<partitionable_loop*>(curr_outer));
	}
	else
	  read_loops.push_back(static_cast<partitionable_loop*>(curr_outer));
      }
    }


  }
  return ret_access;
}


void array_details::CheckConflicts(access_details* new_access)
{
  //If new access is OA_DESCRIPTOR from partitionable loop, check if there are any other accesses of type OA_DESCRIPTOR
  if( new_access->my_branch->IsLoop() ){
    const loop* my_loop = static_cast<const loop*>(new_access->my_branch);
    if( my_loop->IsPartitionableLoop() && new_access->access_type == OA_DESCRIPTOR )
      for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ ){
	if( (*it)->my_branch->IsLoop()) {
	  access_details* curr_access = *it;
	  const loop* curr_loop = static_cast<const loop*>(curr_access->my_branch);
	  if( curr_access->access_type == OA_DESCRIPTOR ){
	    //If the previous access is an inner-loop change type to IA
	    if( !curr_loop->IsPartitionableLoop() ){
	      printf("Changed access %s from type 0 to type 1\n",new_access->access_exp->unparseToString().c_str());
	      new_access->access_type = IA_DESCRIPTOR;
	      break;
	    }
	    //If the previous access is from another partitionable loop check that the loop bounds are the same.
	    else{
	      if( curr_loop->GetGroupNum() != curr_loop->GetGroupNum() ){
		printf("Changed access %s from type 0 to type 1\n",new_access->access_exp->unparseToString().c_str());
		new_access->access_type = IA_DESCRIPTOR;
		break;
	      }
	    }
	  }
	}
      }
    //If the new access is from inner loop and is OA_DESCRIPTOR then change all prev. OA_DESCRIPTOR and partitionable loop to IA_DESCRIPTOR 
    else if( new_access->access_type == OA_DESCRIPTOR ){
      for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ )
	if( (*it)->my_branch->IsLoop() ) {
	  access_details* curr_access = *it;
	  const loop* curr_loop = static_cast<const loop*>(curr_access->my_branch);
	  if( curr_loop->IsPartitionableLoop() && curr_access->access_type == OA_DESCRIPTOR ){
	    printf("Changed access %s from type 0 to type 1\n",curr_access->access_exp->unparseToString().c_str());
	    curr_access->access_type = OA_DESCRIPTOR;
	  }
	}
    }
  }
}

void array_details::GetCommunicationInfo(SgStatement* insert_before)
{
  printf("Array %s, communication req., size = %d,%d,%d :",array_name.c_str(),(int)write_loops.size(),(int)local_writes.size(),(int)read_loops.size());
  SgName thread_id_name("__thread_id__");
  assert(communicate_reads_for_fn && communicate_writes_for_fn );
  if( write_loops.size() != 0 || local_writes.size() != 0 ){
    int* read_loop_num = new int[read_loops.size()];
    int count = 0;
    for( deque<partitionable_loop*>::iterator it = read_loops.begin() ; it != read_loops.end() ; it++, count++ )
      read_loop_num[count] = (*it)->GetLoopNum();

    for( deque< pair<partitionable_loop*,VariantT> >::iterator it = write_loops.begin() ; it != write_loops.end() ; it++ ){
      (*it).first->CommunicateWrites();
      SgExprListExp* write_fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildIntVal((*it).first->GetLoopNum()),SageBuilder::buildIntVal(my_num));
      SgExprStatement* write_fn_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_writes_for_fn,write_fn_args));
      SageInterface::insertStatementBefore(insert_before,write_fn_stmt);
      int write_loop_num = (*it).first->GetLoopNum();
      printf(" %d(write)",write_loop_num);
      if( read_loops.size() != 0 ){
	int pos;
	for( pos =0 ; pos < read_loops.size() ; pos++ )
	  if( write_loop_num < read_loop_num[pos] )
	    break;
	if( pos == read_loops.size() ){
	  read_loops[0]->CommunicateReads();
	  SgExprListExp* read_fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildIntVal(read_loops[0]->GetLoopNum()),SageBuilder::buildIntVal(my_num));
	  SgExprStatement* read_fn_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_reads_for_fn,read_fn_args));
	  SageInterface::insertStatementBefore(insert_before,read_fn_stmt);
	  printf(" %d(read)",read_loop_num[0]);
	}
	else{
	  read_loops[pos]->CommunicateReads();
	  SgExprListExp* read_fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildIntVal(read_loops[pos]->GetLoopNum()),SageBuilder::buildIntVal(my_num));
	  SgExprStatement* read_fn_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_reads_for_fn,read_fn_args));
	  SageInterface::insertStatementBefore(insert_before,read_fn_stmt);

	  printf(" %d(read)",read_loop_num[pos]);
	}
      }
    }
    
    for( deque<partitionable_loop*>::iterator it = local_writes.begin() ; it != local_writes.end() ; it++ ){
      int write_loop_num = (*it)->GetLoopNum();
      if( read_loops.size() != 0 ){
    	int pos;
    	for( pos =0 ; pos < read_loops.size() ; pos++ )
    	  if( write_loop_num < read_loop_num[pos] )
    	    break;
    	if( pos == read_loops.size() ){
    	  read_loops[0]->CommunicateReads();
	  SgExprListExp* read_fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildIntVal(read_loops[0]->GetLoopNum()),SageBuilder::buildIntVal(my_num));
	  SgExprStatement* read_fn_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_reads_for_fn,read_fn_args));
	  SageInterface::insertStatementBefore(insert_before,read_fn_stmt);
    	  printf(" %d(read)",read_loop_num[0]);
    	}
    	else{
    	  read_loops[pos]->CommunicateReads();
	  SgExprListExp* read_fn_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildIntVal(read_loops[pos]->GetLoopNum()),SageBuilder::buildIntVal(my_num));
	  SgExprStatement* read_fn_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(communicate_reads_for_fn,read_fn_args));
	  SageInterface::insertStatementBefore(insert_before,read_fn_stmt);

    	  printf(" %d(read)",read_loop_num[pos]);
    	}
      }
    }
    delete[] read_loop_num;
  }
  printf("\n");
}


void array_details::InitAccessArrays(SgStatement* insert_before)
{
  for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ ){
    (*it)->InitAccessArray(insert_before,array_name,malloc_fn);
  }
}

void array_details::DeleteArrays(SgStatement* insert_before)
{
  for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ ){
    (*it)->DeleteAccessArray(insert_before,free_fn);
  }
  SgType* base_type = isSgPointerType(array_symbol->get_type())->get_base_type();
  SgFunctionSymbol* free_to_use;
  if( isSgPointerType(base_type) ){
    SgType* base_type_2d = isSgPointerType(base_type)->get_base_type();
    if( isSgTypeFloat(base_type_2d) )
      free_to_use = free_2d_float_fn;
    else if( isSgTypeDouble(base_type_2d) )
      free_to_use = free_2d_double_fn;
    else
      assert(0);
  }
  else
    free_to_use = free_fn;
  SgFunctionCallExp* free_call = SageBuilder::buildFunctionCallExp(free_to_use,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(local_array)));
  SgExprStatement* free_stmt = SageBuilder::buildExprStatement(free_call);
  SageInterface::insertStatementBefore(insert_before,free_stmt);
}




void array_details::InitLocalArrays(SgStatement* insert_before)
{
  assert(get_size_fn);
  string local_array_size_string = array_name;
  local_array_size_string.insert(0,"__localsize_").append("__");
  SgName local_array_size_name(local_array_size_string);
  SgName thread_id_name("__thread_id__");
  SgVarRefExp* thread_id_var = SageBuilder::buildVarRefExp(thread_id_name);
  SgIntVal* local_array_num = SageBuilder::buildIntVal(my_num);
  SgFunctionCallExp* local_array_size_exp = SageBuilder::buildFunctionCallExp(get_size_fn,SageBuilder::buildExprListExp(thread_id_var,local_array_num));
  SgAssignInitializer* local_array_size_rhs = SageBuilder::buildAssignInitializer(local_array_size_exp);
  local_array_size = SageBuilder::buildVariableDeclaration(local_array_size_string,SageBuilder::buildIntType(),local_array_size_rhs);
  SageInterface::insertStatementBefore(insert_before,local_array_size);
  
  SgPointerType* local_array_type = isSgPointerType(array_symbol->get_type());
  assert(local_array_type);
  string local_array_string = array_name;
  local_array_string.insert(0,"__local_").append("__");
  SgType* base_type = local_array_type->get_base_type();

  SgName thread_id_name2("__thread_id__");
  SgVarRefExp* populate_arg1 = SageBuilder::buildVarRefExp(thread_id_name2);
  SgIntVal* populate_arg2 = SageBuilder::buildIntVal(my_num);
  SgExpression* populate_arg5;
  SgExpression* populate_arg3;
  SgExpression* populate_arg4;
  assert(populate_fn);
  
  if( stride_string.compare(" 1 ") == 0 ){
    assert( isSgPointerType(base_type) == NULL );
    assert(malloc_fn);
    SgSizeOfOp* sizeof_exp = SageBuilder::buildSizeOfOp(base_type);
    SgVarRefExp* local_size_exp = SageBuilder::buildVarRefExp(local_array_size);
    SgMultiplyOp* malloc_arg = SageBuilder::buildBinaryExpression<SgMultiplyOp>(local_size_exp,sizeof_exp);
    SgFunctionCallExp* malloc_exp = SageBuilder::buildFunctionCallExp(malloc_fn,SageBuilder::buildExprListExp(malloc_arg));
    SgCastExp* local_array_cast = SageBuilder::buildCastExp(malloc_exp,local_array_type);
    SgAssignInitializer* local_array_rhs = SageBuilder::buildAssignInitializer(local_array_cast);
    local_array = SageBuilder::buildVariableDeclaration(local_array_string,local_array_type,local_array_rhs);
    SageInterface::insertStatementBefore(insert_before,local_array);
    populate_arg5 = SageBuilder::buildIntVal(1);
    populate_arg3 = SageBuilder::buildVarRefExp(local_array);
    populate_arg4 = SageBuilder::buildVarRefExp(array_symbol->get_name());
  }
  else{
    assert(isSgPointerType(base_type));
    SgType* base_type_2d = isSgPointerType(base_type)->get_base_type();
    SgFunctionSymbol* malloc_to_use;
    if( isSgTypeFloat(base_type_2d) )
      malloc_to_use = malloc_2d_float_fn;
    else if( isSgTypeDouble(base_type_2d) )
      malloc_to_use = malloc_2d_double_fn;
    else
      assert(0);
    
    int stride_size = stride_string.size();
    SgVarRefExp* malloc_arg1 = SageBuilder::buildVarRefExp(local_array_size);
    SgExpression* malloc_arg2;
    if( isdigit(stride_string.c_str()[0]) ){
      malloc_arg2 = SageBuilder::buildIntVal(atoi(stride_string.c_str()));
      populate_arg5 = SageBuilder::buildIntVal(atoi(stride_string.c_str()));
    }
    else{
      SgName stride_name(stride_string);
      malloc_arg2 = SageBuilder::buildVarRefExp(stride_name);
      populate_arg5 = SageBuilder::buildVarRefExp(stride_name);
    }
    SgFunctionCallExp* malloc_exp = SageBuilder::buildFunctionCallExp(malloc_to_use,SageBuilder::buildExprListExp(malloc_arg1,malloc_arg2));
    SgAssignInitializer* local_array_rhs = SageBuilder::buildAssignInitializer(malloc_exp);
    local_array = SageBuilder::buildVariableDeclaration(local_array_string,local_array_type,local_array_rhs);
    SageInterface::insertStatementBefore(insert_before,local_array);
    populate_arg3 = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(local_array),SageBuilder::buildIntVal(0));
    populate_arg4 = SageBuilder::buildBinaryExpression<SgPntrArrRefExp>(SageBuilder::buildVarRefExp(array_symbol->get_name()),SageBuilder::buildIntVal(0));
  }
  
  SgExprListExp* populate_args = SageBuilder::buildExprListExp(populate_arg1,populate_arg2,populate_arg3,populate_arg4,populate_arg5);
  SgFunctionCallExp* populate_fn_call = SageBuilder::buildFunctionCallExp(populate_fn,populate_args);
  SgExprStatement* populate_fn_stmt = SageBuilder::buildExprStatement(populate_fn_call);
  SageInterface::insertStatementBefore(insert_before,populate_fn_stmt);

  for( deque<access_details*>::iterator it = all_accesses.begin() ; it != all_accesses.end() ; it++ )
    (*it)->RenumberAccess(insert_before,my_num);
}


