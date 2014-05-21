/*
 * analyseAST.cpp: This file is part of the IEC project.
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
 * @file: analyseAST.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/codegen.hpp"
#include <boost/regex.hpp>
#include "CompileTime/data_scalars.hpp"

using namespace std;

extern bool CompareExpression(const SgExpression* , const SgExpression*);

// codegen::codegen(SgPragmaDeclaration* origloop_pragma, SgProject* pr):
  // project(pr),
codegen::codegen(SgStatement* rn, deque<SgForStatement*>& target_loop_nodes, set<SgVariableSymbol*>& data_array_symbols, set<SgVariableSymbol*>& data_scalar_symbols, set<SgVariableSymbol*>& indirection_array_symbols, set<SgVariableSymbol*>& indirection_scalar_symbols, set<SgVariableSymbol*>& iterator_symbols, set<SgVariableSymbol*>& parameter_symbols, set<SgVariableSymbol*>& outer_var_symbols, set<SgVariableSymbol*>& privatizable_symbols, const deque<array_string*>& array_info ):
  proc_id_name("__proc_id__"),
  nprocs_name("__nprocs__"),
  myid_name("__myid__"),
  thread_id_name("__thread_id__"),
  nthreads_name("__nthreads__"),
  data_num_count(NULL),
  ro_mask(NULL),
  num_data_arrays(0),
  root_node(isSgStatement(rn))
{
  num_partitionable_loops = 0;
  num_partitionable_groups = 0;
  for( deque<SgForStatement*>::iterator it = target_loop_nodes.begin() ; it != target_loop_nodes.end() ; it++ ){
    SgForStatement* curr_loop = *it;
    partitionable_loop* curr_part_loop = static_cast<partitionable_loop*>((*it)->getAttribute("partitionable_loop"));
    bool set_group = false;
    for( deque<partitionable_loop*>::iterator jt = all_loops.begin() ; jt != all_loops.end() ; jt++ ){
      partitionable_loop* old_part_loop = *jt;
      if( CompareExpression(curr_part_loop->GetLBExp(),old_part_loop->GetLBExp()) && CompareExpression(curr_part_loop->GetUBExp(),old_part_loop->GetUBExp())){
	curr_part_loop->SetGroupLoopNum(num_partitionable_loops++,old_part_loop->GetGroupNum());
	set_group = true;
	break;
      }
    }
    if( !set_group )
      curr_part_loop->SetGroupLoopNum(num_partitionable_loops++,num_partitionable_groups++);
    printf("[IE] Partitionable Loop : LB: %s, UB: %s, Num : %d, Group: %d\n",curr_part_loop->GetLBExp()->unparseToString().c_str(),curr_part_loop->GetUBExp()->unparseToString().c_str(),curr_part_loop->GetLoopNum(),curr_part_loop->GetGroupNum());
    all_loops.push_back(curr_part_loop);
  }
  InitIndirectionArrays(indirection_array_symbols,indirection_scalar_symbols,array_info);
  InitDataArrays(data_array_symbols,data_scalar_symbols,array_info);
  iterators.insert(iterators.begin(),iterator_symbols.begin(),iterator_symbols.end());
  parameters.insert(parameters.begin(),parameter_symbols.begin(),parameter_symbols.end());
  outer_vars.insert(outer_vars.begin(),outer_var_symbols.begin(),outer_var_symbols.end());
  privatizable.insert(privatizable.begin(),privatizable_symbols.begin(),privatizable_symbols.end());
  printf("[IE] Parameters:");
  for( deque<SgVariableSymbol*>::const_iterator it = parameters.begin() ; it != parameters.end() ; it++ )
    printf(" %s",(*it)->get_name().getString().c_str());
  printf("\n[IE] Iterators:");
  for( deque<SgVariableSymbol*>::const_iterator it = iterators.begin() ; it != iterators.end() ; it++ )
    printf(" %s,",(*it)->get_name().getString().c_str());
  printf("\n[IE] Outer Variables:");
  for( deque<SgVariableSymbol*>::const_iterator it = outer_vars.begin() ; it != outer_vars.end() ; it++ )
    printf(" %s,",(*it)->get_name().getString().c_str());
  printf("\n[IE] Privatizable:");
  for( deque<SgVariableSymbol*>::const_iterator it = privatizable.begin() ; it != privatizable.end() ; it++ )
    printf(" %s,",(*it)->get_name().getString().c_str());
  printf("\n");

  for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ )
    (*it)->AnalyseTargetLoop();
}

codegen::~codegen()
{
  for( deque<array_details*>::iterator it = data_arrays.begin() ; it != data_arrays.end() ; it++ )
    delete (*it);
  data_arrays.clear();

  for( deque<safe_scalar_details*>::iterator it = safe_scalars.begin(); it != safe_scalars.end() ; it++ )
    delete (*it);
  safe_scalars.clear();
  
  for( deque<indarray_details*>::iterator it = indirection_arrays.begin() ; it != indirection_arrays.end() ; it++ )
    delete (*it);
  indirection_arrays.clear();

  for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ )
    delete (*it);
  all_loops.clear();
  
  parameters.clear();
  iterators.clear();
  outer_vars.clear();
  privatizable.clear();
}



void codegen::InitIndirectionArrays(set<SgVariableSymbol*>& indirection_array_symbols, set<SgVariableSymbol*>& indirection_scalar_symbols, const deque<array_string*>& array_info)
{
  printf("[IE] Indirection Arrays:");
  for( set<SgVariableSymbol*>::iterator it = indirection_array_symbols.begin() ; it != indirection_array_symbols.end() ; it++) {
    SgVariableSymbol* curr_array = *it;
    deque<array_string*>::const_iterator jt;
    string array_name = curr_array->get_name().getString();
    for( jt = array_info.begin() ; jt != array_info.end() ; jt++ )
      if( (*jt)->array_name.compare(array_name) == 0 )
	break;
    assert(jt != array_info.end() );
    indarray_details* new_indarray = new indarray_details(curr_array,*jt);
    (*it)->setAttribute("indirection_array",new_indarray);
    indirection_arrays.push_back(new_indarray);
    printf(" %s,",(*it)->get_name().getString().c_str());
  }
  printf("\n[IE] Indirection Scalars: ");
  for( set<SgVariableSymbol*>::iterator it = indirection_scalar_symbols.begin() ; it != indirection_scalar_symbols.end() ; it++) {
    SgVariableSymbol* curr_var = *it;
    safe_scalar_details* new_safe_scalar = new safe_scalar_details(curr_var);
    (*it)->setAttribute("indirection_scalar",new_safe_scalar);
    safe_scalars.push_back(new_safe_scalar);
    printf(" %s,",(*it)->get_name().getString().c_str());
  }
  printf("\n");
}

void codegen::InitDataArrays(set<SgVariableSymbol*>& data_array_symbols, set<SgVariableSymbol*>& data_scalar_symbols, const deque<array_string*>& array_info)
{
  printf("[IE] Data Arrays:");
  for( set<SgVariableSymbol*>::iterator it = data_array_symbols.begin() ; it != data_array_symbols.end() ; it++ ){
    SgVariableSymbol* curr_array = *it;
    deque<array_string*>::const_iterator jt;
    string array_name = curr_array->get_name().getString();
    for( jt = array_info.begin() ; jt != array_info.end() ; jt++ )
      if( (*jt)->array_name.compare(array_name) == 0 )
	break;
    assert(jt != array_info.end() );
    array_details* new_array_details = new array_details(curr_array,*jt);
    (*it)->setAttribute("data_array",new_array_details);
    data_arrays.push_back(new_array_details);
    printf(" %s,",(*it)->get_name().getString().c_str());
  }
  printf("\n");

  printf("[IE] Data Scalars:");
  for( set<SgVariableSymbol*>::iterator it = data_scalar_symbols.begin() ; it != data_scalar_symbols.end() ; it++ ){
    data_scalar* new_data_scalar = new data_scalar(*it);
    (*it)->setAttribute("data_scalar",new_data_scalar);
    printf(" %s,",(*it)->get_name().getString().c_str());
    data_scalars.push_back(new_data_scalar);
  }
  printf("\n");
}


