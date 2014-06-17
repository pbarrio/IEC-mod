/*
 * codegen.cpp: This file is part of the IEC project.
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
 * @file: codegen.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/codegen.hpp"
#include "OmpAttribute.h"

using namespace std;

extern const SgVariableSymbol* GetLoopIterator(const SgForInitStatement*);
extern SgVariableSymbol* GetLoopIterator(SgForInitStatement*);

SgFunctionSymbol* codegen::is_known_fn = NULL;
SgFunctionSymbol* codegen::get_elem_fn = NULL;

void codegen::GenerateCCode(SgGlobal* global_scope)
{
	//NULL statements to represent the three parts of the inspector.
	SgNullStatement* graph_gen = SageBuilder::buildNullStatement();
	SageInterface::insertStatementBefore(root_node,graph_gen);
  
	SgNullStatement* init_counter = SageBuilder::buildNullStatement();
	SageInterface::insertStatementBefore(root_node,init_counter);
 
	SgNullStatement* init_thread_counter = SageBuilder::buildNullStatement();
	SgNullStatement* init_access_array = SageBuilder::buildNullStatement();  
	SgNullStatement* executor_stmt = SageBuilder::buildNullStatement();
	SgBasicBlock* omp_block = NULL;
	if( CompilerOpts::IsOpenMPMode() ){
		omp_block = SageBuilder::buildBasicBlock(init_thread_counter,init_access_array,executor_stmt);
		SageInterface::insertStatementBefore(root_node,omp_block);
	}
	else{
		SageInterface::insertStatementBefore(root_node,init_thread_counter);
		SageInterface::insertStatementBefore(root_node,init_access_array);
		SageInterface::insertStatementBefore(root_node,executor_stmt);
	}

	SgExprStatement *init_thread_id,*init_myid;
  
	if( CompilerOpts::IsOpenMPMode() ){ 
		//Get number of threads
		SgName get_num_threads_fn_name("omp_get_max_threads");
		SgFunctionSymbol* get_num_threads_fn = global_scope->lookup_function_symbol("omp_get_max_threads");
		assert(get_num_threads_fn);
		SgFunctionCallExp* get_num_threads_exp = SageBuilder::buildFunctionCallExp(get_num_threads_fn);
		SgName nthreads_decl("__nthreads__");
		SgAssignOp* nthreads_assign_op = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(nthreads_decl),get_num_threads_exp);
		SgExprStatement* nthreads_stmt = SageBuilder::buildExprStatement(nthreads_assign_op);
		SageInterface::insertStatementBefore(graph_gen ,nthreads_stmt);
    
		SgName omp_thread_id_name("omp_get_thread_num");
		SgFunctionSymbol* omp_thread_id_fn = global_scope->lookup_function_symbol(omp_thread_id_name);
		init_thread_id = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildFunctionCallExp(omp_thread_id_fn)));
		SgMultiplyOp* init_myid_lhs = SageBuilder::buildBinaryExpression<SgMultiplyOp>(SageBuilder::buildVarRefExp(proc_id_name),SageBuilder::buildVarRefExp(nthreads_name));
		SgAddOp* init_myid_rhs = SageBuilder::buildBinaryExpression<SgAddOp>(init_myid_lhs,SageBuilder::buildVarRefExp(thread_id_name));
		SgAssignOp* init_myid_exp = SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(myid_name),init_myid_rhs);
		init_myid = SageBuilder::buildExprStatement(init_myid_exp);
	}
	else{
		SgExprStatement* nthreads_stmt = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(nthreads_name),SageBuilder::buildIntVal(1)));
		SageInterface::insertStatementBefore(graph_gen,nthreads_stmt);
		init_thread_id = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(thread_id_name),SageBuilder::buildIntVal(0)));
		init_myid = SageBuilder::buildExprStatement(SageBuilder::buildBinaryExpression<SgAssignOp>(SageBuilder::buildVarRefExp(myid_name),SageBuilder::buildVarRefExp(proc_id_name)));
	}
	SageInterface::insertStatementBefore(init_thread_counter,init_thread_id);
	SageInterface::insertStatementBefore(init_thread_counter,init_myid);

	//Create array with loop sizes
	SgArrayType* array_type = SageBuilder::buildArrayType(SageBuilder::buildIntType(),SageBuilder::buildIntVal(num_partitionable_groups));
	SgExprListExp* iter_sizes = NULL; 

	int next_group = 0;
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		if( curr_part_loop->GetGroupNum() == next_group ){
			if( iter_sizes )
				iter_sizes = SageBuilder::buildExprListExp(iter_sizes,SageBuilder::buildBinaryExpression<SgSubtractOp>(SageInterface::copyExpression(curr_part_loop->GetUBExp()),SageInterface::copyExpression(curr_part_loop->GetLBExp())));
			else
				iter_sizes =  SageBuilder::buildExprListExp(SageBuilder::buildBinaryExpression<SgSubtractOp>(SageInterface::copyExpression(curr_part_loop->GetUBExp()),SageInterface::copyExpression(curr_part_loop->GetLBExp())));
			next_group++;
		}
	}
	string var_string = "iter_num_count";
	SgVariableDeclaration* iter_num_count = SageBuilder::buildVariableDeclaration(var_string,array_type,SageBuilder::buildAggregateInitializer(iter_sizes));
	SageInterface::insertStatementBefore(graph_gen,iter_num_count);


	//Create array with data sizes
	GenerateDataSizeArray(graph_gen);
  
	//Create Inspector
	SgName create_inspector_fn_name("create_inspector");
	SgFunctionSymbol* create_inspector_fn = global_scope->lookup_function_symbol(create_inspector_fn_name);
	assert(create_inspector_fn);
	SgVarRefExp* arg1 = SageBuilder::buildVarRefExp(proc_id_name);
	SgVarRefExp* arg2 = SageBuilder::buildVarRefExp(nprocs_name);
	SgVarRefExp* arg3 = SageBuilder::buildVarRefExp(nthreads_name);
	SgIntVal* arg4 = SageBuilder::buildIntVal(num_partitionable_groups);
	SgIntVal* arg5 = SageBuilder::buildIntVal(num_data_arrays);
	SgIntVal* arg6 = SageBuilder::buildIntVal(num_partitionable_loops);
	SgIntVal* arg7 = SageBuilder::buildIntVal(indirection_arrays.size());
	SgVarRefExp* arg8 = SageBuilder::buildVarRefExp(iter_num_count);
	SgVarRefExp* arg9 = SageBuilder::buildVarRefExp(data_num_count);
	SgVarRefExp* arg10 = SageBuilder::buildVarRefExp(ro_mask);
	SgExprListExp* args = SageBuilder::buildExprListExp(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10);
	SgExprStatement* create_inspector_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(create_inspector_fn,args));
	SageInterface::insertStatementBefore(graph_gen,create_inspector_stmt);
  
	// //Set indirection array details
	SgName set_access_param_fn_name("set_access_array_param");
	indarray_details::set_access_param_fn = global_scope->lookup_function_symbol(set_access_param_fn_name);
	assert(indarray_details::set_access_param_fn);
	for( deque< indarray_details* >::iterator it = indirection_arrays.begin() ; it != indirection_arrays.end() ; it++)
		(*it)->GenerateSetAccessParam(graph_gen);

	// Generate data array stride info
	SgName set_array_stride_fn_name("set_array_stride");
	array_details::set_array_stride_fn = global_scope->lookup_function_symbol(set_array_stride_fn_name);
	assert(array_details::set_array_stride_fn);
	for( deque< array_details* >::const_iterator it = data_arrays.begin() ; it != data_arrays.end() ; it++ )
		(*it)->GenerateStrideInfo(graph_gen);

	SgName malloc_fn_name("malloc");
	array_details::malloc_fn = global_scope->lookup_function_symbol(malloc_fn_name);
	assert(array_details::malloc_fn);

	SgName memset_fn_name("memset");
	partitionable_loop::memset_fn = global_scope->lookup_function_symbol(memset_fn_name);
	assert(partitionable_loop::memset_fn);

	SgName free_fn_name("free");
	array_details::free_fn = global_scope->lookup_function_symbol(free_fn_name);
	assert(array_details::free_fn);

	SgName get_proc_iter_size_fn_name("get_proc_iter_size");
	partitionable_loop::get_proc_iter_size_fn = global_scope->lookup_function_symbol(get_proc_iter_size_fn_name);
	assert(partitionable_loop::get_proc_iter_size_fn);

	//INitialize the iteration status arrays
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		curr_part_loop->GenerateLocalBlock(graph_gen,array_details::malloc_fn);
		curr_part_loop->ResetStatusArray(init_counter,array_details::malloc_fn,array_details::free_fn);
	}

	//Initialize the shadow scalars
	for( deque<safe_scalar_details*>::iterator it = safe_scalars.begin() ; it != safe_scalars.end() ; it++ )
		(*it)->InitShadowScalars(graph_gen);

	//Initialize Loop and branch Counters
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		curr_part_loop->InitCounters(init_thread_counter);
	}

	//Enclose graph generator within outer loop
	SgName done_graph_gen_name("done_graph_gen");
	SgFunctionSymbol* done_graph_gen_fn = global_scope->lookup_function_symbol(done_graph_gen_name);
	SgNotOp* done_graph_gen_condn = SageBuilder::buildUnaryExpression<SgNotOp>(SageBuilder::buildFunctionCallExp(done_graph_gen_fn));
	SgNullStatement* done_graph_gen_body = SageBuilder::buildNullStatement();
	SgDoWhileStmt* do_while_graph_gen = SageBuilder::buildDoWhileStmt(done_graph_gen_body,SageBuilder::buildExprStatement(done_graph_gen_condn));
	SageInterface::insertStatementBefore(graph_gen,do_while_graph_gen);

	//Enclose counters within outer loop
	SgNotOp* done_init_counter_condn = SageBuilder::buildUnaryExpression<SgNotOp>(SageBuilder::buildFunctionCallExp(done_graph_gen_fn));
	SgNullStatement* done_init_counter_body = SageBuilder::buildNullStatement();
	SgDoWhileStmt* do_while_init_counter = SageBuilder::buildDoWhileStmt(done_init_counter_body,SageBuilder::buildExprStatement(done_init_counter_condn));
	SageInterface::insertStatementBefore(init_counter,do_while_init_counter);

	//Setup Access arrays;
	for( deque<array_details*>::iterator it = data_arrays.begin() ; it != data_arrays.end() ; it++ )
		(*it)->InitAccessArrays(init_access_array);
  
	//Setup Bounds/Conditional arrays
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		curr_part_loop->InitBoundsArrays(init_access_array,array_details::malloc_fn);
	}

	//Reset counters at appropriate locations
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		curr_part_loop->ResetAllCounters(init_access_array);
	}
 

	SgName partition_fn_name("partition_hypergraph");
	SgFunctionSymbol* partition_fn = global_scope->lookup_function_symbol(partition_fn_name);
	SgExprStatement* partition_fn_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(partition_fn,SageBuilder::buildExprListExp(SageBuilder::buildIntVal(CompilerOpts::PartitionerType()))));
	SageInterface::insertStatementBefore(graph_gen,partition_fn_stmt);


	SgName is_known_fn_name("is_known");
	partitionable_loop::is_known_fn = global_scope->lookup_function_symbol(is_known_fn_name);
	assert(partitionable_loop::is_known_fn);
	SgName get_elem_fn_name("get_elem");
	partitionable_loop::get_elem_fn = global_scope->lookup_function_symbol(get_elem_fn_name);
	assert(partitionable_loop::get_elem_fn);
	SgName add_vertex_fn_name("add_vertex");
	partitionable_loop::add_vertex_fn = global_scope->lookup_function_symbol(add_vertex_fn_name);
	assert(partitionable_loop::add_vertex_fn);
	SgName get_vertex_home_fn_name("get_vertex_home");
	partitionable_loop::get_vertex_home_fn = global_scope->lookup_function_symbol(get_vertex_home_fn_name);
	assert(partitionable_loop::get_vertex_home_fn);
	SgName add_pin_fn_name("add_pin_to_net");
	access_details::add_pin_fn = global_scope->lookup_function_symbol(add_pin_fn_name);
	assert(access_details::add_pin_fn);

	//Generate inspector parts
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		//Initialize temporary counters for the prefetch part of inspector
		curr_part_loop->GenerateCInspector(done_graph_gen_body,done_init_counter_body,init_thread_counter,init_access_array);
		curr_part_loop->FreeStatusArray(init_counter,array_details::free_fn);
	}
    
	SgName renumber_fn_name("renumber_access_array");
	access_details::renumber_fn = global_scope->lookup_function_symbol(renumber_fn_name);
	assert(access_details::renumber_fn);
	SgName renumber_offset_fn_name("renumber_offset_array");
	access_details::renumber_offset_fn = global_scope->lookup_function_symbol(renumber_offset_fn_name);
	assert(access_details::renumber_offset_fn);
	SgName renumber_const_offset_fn_name("renumber_const_offset_array");
	access_details::renumber_const_offset_fn = global_scope->lookup_function_symbol(renumber_const_offset_fn_name);
	assert(access_details::renumber_const_offset_fn);
	SgName malloc_2d_double_fn_name("malloc_2d_double");
	array_details::malloc_2d_double_fn = global_scope->lookup_function_symbol(malloc_2d_double_fn_name);
	assert(array_details::malloc_2d_double_fn);
	SgName malloc_2d_float_fn_name("malloc_2d_float");
	array_details::malloc_2d_float_fn = global_scope->lookup_function_symbol(malloc_2d_float_fn_name);
	assert(array_details::malloc_2d_float_fn);
	SgName get_size_fn_name("get_local_size");
	array_details::get_size_fn = global_scope->lookup_function_symbol(get_size_fn_name);
	assert(array_details::get_size_fn);
	SgName populate_fn_name("populate_local_array");
	array_details::populate_fn = global_scope->lookup_function_symbol(populate_fn_name);
	assert(array_details::populate_fn);

	SgName communicate_reads_for_name("communicate_reads_for");
	array_details::communicate_reads_for_fn = global_scope->lookup_function_symbol(communicate_reads_for_name);
	assert(array_details::communicate_reads_for_fn);
	SgName communicate_writes_for_name("communicate_writes_for");
	array_details::communicate_writes_for_fn = global_scope->lookup_function_symbol(communicate_writes_for_name);
	assert(array_details::communicate_writes_for_fn);
  
	//Allocate local arrays and renumber accesses. Also GetComunication Info
	for( deque< array_details*>::iterator it = data_arrays.begin() ; it != data_arrays.end() ; it++ ){
		(*it)->InitLocalArrays(init_access_array);
		(*it)->GetCommunicationInfo(init_access_array);
	}

	//Setup Executor
	SgName setup_executor_name("setup_executor");
	SgFunctionSymbol* setup_executor_fn = global_scope->lookup_function_symbol(setup_executor_name);
	SgExprListExp* setup_executor_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name));
	SgExprStatement* setup_executor_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(setup_executor_fn,setup_executor_args));
	SageInterface::insertStatementBefore(init_access_array,setup_executor_stmt);
  
	SgName communicate_reads_name("communicate_reads");
	partitionable_loop::communicate_reads_fn = global_scope->lookup_function_symbol(communicate_reads_name);
	assert(partitionable_loop::communicate_reads_fn);
	SgName communicate_writes_name("communicate_writes");
	partitionable_loop::communicate_writes_fn = global_scope->lookup_function_symbol(communicate_writes_name);
	assert(partitionable_loop::communicate_writes_fn);
	SgName init_write_ghosts_name("init_write_ghosts");
	partitionable_loop::init_write_ghosts_fn = global_scope->lookup_function_symbol(init_write_ghosts_name);
	assert(partitionable_loop::init_write_ghosts_fn);
	SgName reduce_scalar_name("reduce_scalar");
	partitionable_loop::reduce_scalar_fn = global_scope->lookup_function_symbol(reduce_scalar_name);
	assert(partitionable_loop::reduce_scalar_fn);
  

	//Generate Executor
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		(*it)->GenerateCExecutor((*it)->orig_for);
	}
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		SageInterface::removeStatement((*it)->orig_for);
	}
  
	SgNullStatement* inspector_end = SageBuilder::buildNullStatement();
	if( CompilerOpts::IsOpenMPMode() ){
		SageInterface::removeStatement(root_node);
		omp_block->append_statement(root_node);
		SageInterface::insertStatementAfter(root_node,inspector_end);
	}
	else{
		SageInterface::insertStatementAfter(root_node,inspector_end);
	}
  
	SgName free_2d_double_fn_name("free_2d_double");
	array_details::free_2d_double_fn = global_scope->lookup_function_symbol(free_2d_double_fn_name);
	assert(array_details::free_2d_double_fn);
	SgName free_2d_float_fn_name("free_2d_float");
	array_details::free_2d_float_fn = global_scope->lookup_function_symbol(free_2d_float_fn_name);
	assert(array_details::free_2d_float_fn);

	//CleanUp Code
  
	SgName populate_global_fn_name("populate_global_arrays");
	SgFunctionSymbol* populate_global_fn = global_scope->lookup_function_symbol(populate_global_fn_name);
	assert(populate_global_fn);
	SgFunctionCallExp* populate_global_exp = SageBuilder::buildFunctionCallExp(populate_global_fn,SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name)));
	SgExprStatement* populate_global_stmt = SageBuilder::buildExprStatement(populate_global_exp);
	SageInterface::insertStatementBefore(inspector_end,populate_global_stmt);

	//Deallocate access arrays
	for( deque<array_details*>::iterator it = data_arrays.begin() ; it != data_arrays.end() ; it++ )
		(*it)->DeleteArrays(inspector_end);
  
	for( deque<partitionable_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
		partitionable_loop* curr_part_loop = *it;
		curr_part_loop->DeleteArrays(inspector_end,array_details::free_fn);
	}
      
	SgName delete_insepector_fn_name("delete_inspector");
	SgFunctionSymbol* delete_inspector_fn = global_scope->lookup_function_symbol(delete_insepector_fn_name);
	SgExprListExp* delete_inspector_args = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(thread_id_name));
	SgExprStatement* delete_inspector_stmt = SageBuilder::buildExprStatement(SageBuilder::buildFunctionCallExp(delete_inspector_fn,delete_inspector_args));
	SageInterface::insertStatementBefore(inspector_end,delete_inspector_stmt);

	if( CompilerOpts::IsOpenMPMode() )
		InsertOMPParallel(omp_block);
}


void codegen::GenerateDataSizeArray(SgStatement* insert_before)
{
	assert(data_arrays.size() >0 );
	assert(data_num_count==NULL && ro_mask == NULL && num_data_arrays == 0);
	array_details** ordered_data_arrays = new array_details*[data_arrays.size()];
	for( deque<array_details*>::iterator it = data_arrays.begin(); it != data_arrays.end() ; it++ ){
		ordered_data_arrays[(*it)->GetArrayNum()] = (*it);
		num_data_arrays++;
	}

	SgExprListExp* arg1 = SageBuilder::buildExprListExp(SageBuilder::buildVarRefExp(ordered_data_arrays[0]->GetArraySize()));
	SgExprListExp* ro_arg1 = SageBuilder::buildExprListExp(SageBuilder::buildIntVal(ordered_data_arrays[0]->GetROMask()));
  
	for( int i = 1 ; i < num_data_arrays ; i++ ){
		arg1 = SageBuilder::buildExprListExp(arg1,SageBuilder::buildVarRefExp(ordered_data_arrays[i]->GetArraySize()));
		ro_arg1 = SageBuilder::buildExprListExp(ro_arg1,SageBuilder::buildIntVal(ordered_data_arrays[i]->GetROMask()));
	}
	string var_string("data_num_count");
	SgArrayType* array_type = SageBuilder::buildArrayType(SageBuilder::buildIntType(),SageBuilder::buildIntVal(num_data_arrays));
	data_num_count = SageBuilder::buildVariableDeclaration(var_string,array_type,SageBuilder::buildAggregateInitializer(arg1));
	SageInterface::insertStatementBefore(insert_before,data_num_count);
  
	string ro_var_string("ro_mask");
	SgArrayType* ro_array_type = SageBuilder::buildArrayType(SageBuilder::buildIntType(),SageBuilder::buildIntVal(num_data_arrays));
	ro_mask = SageBuilder::buildVariableDeclaration(ro_var_string,ro_array_type,SageBuilder::buildAggregateInitializer(ro_arg1));
	SageInterface::insertStatementBefore(insert_before,ro_mask);
	delete[] ordered_data_arrays;
}



void codegen::InsertOMPParallel(SgBasicBlock* private_bb)
{
	set<SgVariableSymbol*> omp_private_vars;
	omp_private_vars.insert(outer_vars.begin(),outer_vars.end());
	omp_private_vars.insert(iterators.begin(),iterators.end());
	omp_private_vars.insert(privatizable.begin(),privatizable.end());
	for( deque<data_scalar*>::iterator it = data_scalars.begin() ; it != data_scalars.end() ; it++ )
		omp_private_vars.insert((*it)->getOriginalVar());
	for( deque<safe_scalar_details*>::iterator it = safe_scalars.begin() ; it != safe_scalars.end() ; it++ )
		omp_private_vars.insert((*it)->getOriginalVar());
	OmpSupport::OmpAttribute* bb_omp_attribute = OmpSupport::buildOmpAttribute(OmpSupport::e_parallel,NULL,false);
	for( set<SgVariableSymbol*>::iterator it = omp_private_vars.begin() ; it != omp_private_vars.end() ; it++ ){
		bb_omp_attribute->addVariable(OmpSupport::e_firstprivate,(*it)->get_name().getString());
	}
	bb_omp_attribute->addVariable(OmpSupport::e_firstprivate,thread_id_name.getString());
	bb_omp_attribute->addVariable(OmpSupport::e_firstprivate,myid_name.getString());
	bb_omp_attribute->addExpression(OmpSupport::e_default,"shared");
	OmpSupport::addOmpAttribute(bb_omp_attribute,private_bb);
	OmpSupport::generatePragmaFromOmpAttribute(private_bb);
}
