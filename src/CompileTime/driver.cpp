/*
 * driver.cpp: This file is part of the IEC project.
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
 * @file: driver.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/driver.hpp"
#include <boost/regex.hpp>
#include "ScopExtractor/ScopExtractor.hpp"
#include "ScopExtractor/SageTools.hpp"

using namespace std;

#define MAX_DEPTH 100
#define MIN(a,b) ( a < b ? a : b )

driver* driver::singleton_driver = NULL;
CompilerOpts* CompilerOpts::singleton_instance = NULL;

ie_init::ie_init(SgGlobal* gs):
  global_scope(gs)
{}


void ie_init::visit(SgNode* curr_node)
{
  if( isSgFunctionDeclaration(curr_node) ){
    if( isSgFunctionDeclaration(curr_node)->get_name().getString().compare("main") == 0 ){
      SgFunctionDefinition* main_defn = isSgFunctionDeclaration(curr_node)->get_definition();
            
      SgBasicBlock* main_body = main_defn->get_body();
      SgStatementPtrList& main_stmt_list = main_body->get_statements();
      SgStatementPtrList::const_iterator main_stmt_iter = main_stmt_list.begin();
      SgStatement* first_stmt = *main_stmt_iter;

      //MPI_Init statement
      SgName mpi_init_fn_name("MPI_Init");
      SgFunctionSymbol* mpi_init_symbol= global_scope->lookup_function_symbol(mpi_init_fn_name);
      assert(mpi_init_symbol);
      SgName mpi_init_arg1_name("argc");
      SgVarRefExp* arg1_var = SageBuilder::buildVarRefExp(mpi_init_arg1_name);
      SgAddressOfOp* arg1 = SageBuilder::buildAddressOfOp(isSgExpression(arg1_var));
      SgName mpi_init_arg2_name("argv");
      SgVarRefExp* arg2_var = SageBuilder::buildVarRefExp(mpi_init_arg2_name);
      SgAddressOfOp* arg2 = SageBuilder::buildAddressOfOp(isSgExpression(arg2_var));
      SgExprListExp* init_args = SageBuilder::buildExprListExp(isSgExpression(arg1),isSgExpression(arg2));
      SgFunctionCallExp* mpi_init_fn = SageBuilder::buildFunctionCallExp(mpi_init_symbol,init_args);
      SgExprStatement* mpi_init_stmt = SageBuilder::buildExprStatement(isSgExpression(mpi_init_fn));
      SageInterface::insertStatementBefore(first_stmt,mpi_init_stmt);

      //Inserting ARMCI init statement
      SgName armci_init_fn_name("ARMCI_Init");
      SgFunctionSymbol* armci_init_fn = global_scope->lookup_function_symbol(armci_init_fn_name);
      assert(armci_init_fn);
      SgFunctionCallExp* armci_init_call  = SageBuilder::buildFunctionCallExp(armci_init_fn);
      SgExprStatement* armci_init_stmt = SageBuilder::buildExprStatement(armci_init_call);
      SageInterface::insertStatementBefore(first_stmt,armci_init_stmt);

      //Inserting call to populate rank 
      string mpi_comm_world = "MPI_COMM_WORLD";
      SgVarRefExp* rank_arg1 = SageBuilder::buildVarRefExp(mpi_comm_world);
      SgName rank_decl("__proc_id__");
      // SgVariableSymbol* rank_decl = global_scope->lookup_variable_symbol("__proc_id__");
      // assert(rank_decl);
      SgAddressOfOp* rank_arg2 = SageBuilder::buildAddressOfOp(isSgExpression(SageBuilder::buildVarRefExp(rank_decl)));
      SgName rank_fn_name("MPI_Comm_rank");
      SgFunctionSymbol* rank_fn = global_scope->lookup_function_symbol(rank_fn_name);
      assert(rank_fn);
      SgExprListExp* rank_fn_args = SageBuilder::buildExprListExp(isSgExpression(rank_arg1),isSgExpression(rank_arg2));
      SgFunctionCallExp* rank_fn_call = SageBuilder::buildFunctionCallExp(rank_fn,rank_fn_args);
      SgExprStatement* rank_stmt = SageBuilder::buildExprStatement(isSgExpression(rank_fn_call));
      SageInterface::insertStatementBefore(first_stmt,rank_stmt);

      //Inserting call to populate nprocs
      SgName nprocs_decl("__nprocs__");
      // SgVariableSymbol* nprocs_decl = global_scope->lookup_variable_symbol("__nprocs__");
      // assert(nprocs_decl);
      SgAddressOfOp* nprocs_arg2 = SageBuilder::buildAddressOfOp(isSgExpression(SageBuilder::buildVarRefExp(nprocs_decl)));
      SgName nprocs_fn_name("MPI_Comm_size");
      SgFunctionSymbol* nprocs_fn = global_scope->lookup_function_symbol(nprocs_fn_name);
      assert(nprocs_fn);
      SgExprListExp* nprocs_fn_args = SageBuilder::buildExprListExp(isSgExpression(rank_arg1),isSgExpression(nprocs_arg2));
      SgFunctionCallExp* nprocs_fn_call = SageBuilder::buildFunctionCallExp(nprocs_fn,nprocs_fn_args);
      SgExprStatement*  nprocs_stmt = SageBuilder::buildExprStatement(isSgExpression(nprocs_fn_call));
      SageInterface::insertStatementBefore(first_stmt,nprocs_stmt);

      SgStatementPtrList& main_new_stmt_list = main_body->get_statements();
      SgStatementPtrList::const_iterator main_new_stmt_iter = main_new_stmt_list.end();
      
      main_new_stmt_iter--;
      SgReturnStmt* last_statement = isSgReturnStmt(*main_new_stmt_iter);
      assert(last_statement);

      //Insert MPI_Barrier;
//       SgName mpi_barrier_name("MPI_Barrier");
//       SgFunctionSymbol* mpi_barrier_fn  = global_scope->lookup_function_symbol(mpi_barrier_name);
//       assert(mpi_barrier_fn);
//       SgExpression* barrier_arg1 = SageBuilder::buildVarRefExp(mpi_comm_world);
//       SgExprListExp* barrier_args = SageBuilder::buildExprListExp(barrier_arg1);
//       SgFunctionCallExp* mpi_barrier_call = SageBuilder::buildFunctionCallExp(mpi_barrier_fn,barrier_args);
//       SgExprStatement* mpi_barrier_stmt = SageBuilder::buildExprStatement(mpi_barrier_call);
//       SageInterface::insertStatementBefore(last_statement,mpi_barrier_stmt);
  
      //Insert ARMCI Finalize
      SgName armci_finalize_name("ARMCI_Finalize");
      SgFunctionSymbol* armci_finalize_fn = global_scope->lookup_function_symbol(armci_finalize_name);
      assert(armci_finalize_fn);
      SgFunctionCallExp* armci_finalize_call = SageBuilder::buildFunctionCallExp(armci_finalize_fn);
      SgExprStatement* armci_finalize_stmt = SageBuilder::buildExprStatement(armci_finalize_call);
      SageInterface::insertStatementBefore(last_statement,armci_finalize_stmt);

      //Insert MPI_Finalize
      SgName mpi_finalize_name("MPI_Finalize");
      SgFunctionSymbol* mpi_finalize_fn = global_scope->lookup_function_symbol(mpi_finalize_name);
      assert(mpi_finalize_fn);
      SgFunctionCallExp* mpi_finalize_call = SageBuilder::buildFunctionCallExp(mpi_finalize_fn);
      SgExprStatement* mpi_finalize_stmt = SageBuilder::buildExprStatement(mpi_finalize_call);
      SageInterface::insertStatementBefore(last_statement,mpi_finalize_stmt);
    }
  }
}



void driver::visit(SgNode* curr_node)
{
  SgPragmaDeclaration* curr_pragma_decl = isSgPragmaDeclaration(curr_node);
  if( curr_pragma_decl ){
    string pragma_string = curr_pragma_decl->get_pragma()->get_pragma();
    if( pragma_string.compare(0,6,"arrays") == 0 ){
      string array_names = pragma_string.substr(7,pragma_string.size()-7);

      boost::regex array_name_regexp("\\s*([^,]*)(.*)");
      boost::cmatch all_array_matches;
  
      if( boost::regex_match(array_names.c_str(),all_array_matches,array_name_regexp) ){
	string curr_array_match(all_array_matches[1].first,all_array_matches[1].second);
	// printf("Array subexp : %s\n",curr_array_match.c_str());
	//array_details* new_array = new array_details(curr_array_match);
	array_string* new_string = new array_string(curr_array_match);
	//data_arrays.push_back(pair<array_string,array_details*>(new_string,NULL));
	array_info.push_back(new_string);
	curr_array_match.assign(all_array_matches[2].first,all_array_matches[2].second);
	string full(all_array_matches[0].first,all_array_matches[0].second);
	// printf("Full : %s ; Remaining : %s; size = %d\n",full.c_str(),curr_array_match.c_str(),all_array_matches.size()); 
	bool flag = true;
	while( flag ){
	  boost::regex rest_array_regexp("\\s*,\\s*([^,]*)(.*)");
	  boost::cmatch rest_array_matches;
	  flag = boost::regex_match(curr_array_match.c_str(),rest_array_matches,rest_array_regexp);
	  if( flag ){
	    string full(rest_array_matches[0].first,rest_array_matches[0].second);
	    string array_name(rest_array_matches[1].first,rest_array_matches[1].second);
	    // printf("Array subexp : %s\n",array_name.c_str());
	      //array_details* new_array = new array_details(array_name);
	    array_string* new_string = new array_string(array_name);
	    array_info.push_back(new_string);
	    curr_array_match.assign(rest_array_matches[2].first,rest_array_matches[2].second);
	    // printf("Full : %s; Remaining : %s; size = %d\n",full.c_str(),curr_array_match.c_str(),rest_array_matches.size());
	  }
	}
      }

    }
    else if( !CompilerOpts::IsDetectOnly() &&  pragma_string.compare(0,9,"orig_loop") == 0 ){
      if( CompilerOpts::IsPragmaDriven() ){
	printf("[IE] : Cant handle two orig_loop pragma regions. Aborting\n");
	exit(1);
      }
      CompilerOpts::SetPragmaDriven();
      SgStatement* root_node = SageInterface::getNextStatement(isSgStatement(curr_node));
      vector<SgPragmaDeclaration*> part_loop_pragmas = SageInterface::querySubTree<SgPragmaDeclaration>(root_node);
      deque<SgForStatement*> target_loops;
      set<SgVariableSymbol*> data_arrays,indirection_arrays,indirection_scalars,parameters,iterators,data_scalars,outer_vars,privatizable;
      for( vector<SgPragmaDeclaration*>::iterator it = part_loop_pragmas.begin() ; it != part_loop_pragmas.end() ; it++ ){
	string part_loop_string = (*it)->get_pragma()->get_pragma();
	if( part_loop_string.compare(0,9,"part_loop") == 0 ){
	  SgForStatement* curr_loop = isSgForStatement(SageInterface::getNextStatement(*it));
	  if(curr_loop){
	    printf("Statement %s is not a loop and follows the '#pragma part_loop' - ignored\n",curr_loop->unparseToString().c_str());
	  }
	  else{
	    partitionable_loop* curr_part_loop = new partitionable_loop(curr_loop);
	    if( curr_part_loop->CheckProperties()){
	      target_loops.push_back(curr_loop);
	      curr_loop->setAttribute("partitionable_loop",curr_part_loop);
	    }
	    else{
	      printf("Loop %s wrongly marked as partitionable loop - ignored",curr_loop->unparseToString().c_str());
	      delete curr_part_loop;
	    }
	  }
	}
	else if( part_loop_string.compare(0,6,"orig_loop")){
	  printf("Target computations cant be nested within each other\n");
	  exit(1);
	}
      }
      if( AnalysePartitionableLoops(root_node,target_loops,data_arrays,data_scalars,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars,privatizable ) ){
	codegen* new_codegen = new codegen(isSgStatement(root_node),target_loops,data_arrays,data_scalars,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars,privatizable,array_info);
	target_sequence.push_back(new_codegen);
      }
      else{
	printf("Orig_loop %s doesnt satisfy all requirements -- ignored\n",root_node->unparseToString().c_str());
      }
    }
  }
}


// driver::driver(int argc, char** argv)
// {
//   openmp_mode = true;
//   CompilerOpts::IsDetectOnly() = false;
//   CompilerOpts::IsPragmaDriven() = false;
//   for( int i = 1 ; i < argc ; i++ ){
//     if( !strcmp(argv[i],"-mode-mpi") ){
//       openmp_mode = false;
//       continue;
//     }
//     if( !strcmp(argv[i],"-detect-only") ){
//       CompilerOpts::IsDetectOnly() = true;
//       continue;
//     }
//   }
// }



driver::~driver()
{
  for( list<codegen*>::iterator it = target_sequence.begin() ; it != target_sequence.end() ; it++ )
    delete (*it);
  target_sequence.clear();
}


bool driver::CheckArray(SgVarRefExp* array_exp)
{
  //THis is to nullify the privatizable arrays trick. All arrays are treated as indirection or data arrays.
  if( singleton_driver->array_info.size() == 0 ){
    return true;
  }
  else{

    assert(singleton_driver);
    bool found_array = false;
    string curr_array_name = array_exp->get_symbol()->get_name().getString();
    for( deque<array_string*>::iterator it = singleton_driver->array_info.begin() ; it != singleton_driver->array_info.end() ; it++ )
      if( (*it)->array_name.compare(curr_array_name) == 0 ){
	found_array = true;
	break;
      }
    return found_array;  
  }
}


void driver::GenerateCCode(SgProject* project)
{
  SgFilePtrList& file_list = project->get_fileList();
  SgFilePtrList::const_iterator file_iter = file_list.begin();
  SgSourceFile* file = isSgSourceFile(*file_iter);
  SgGlobal* global_scope = file->get_globalScope();
  assert(global_scope);

  //Setup MPI and openmp
  ie_init init_ie(global_scope);
  init_ie.traverseInputFiles(project,preorder);

  for( list<codegen*>::iterator it = target_sequence.begin(); it != target_sequence.end() ; it++ ){
    (*it)->GenerateCCode(global_scope);
  }
}


bool driver::CheckTargets(SgProject* project)
{
  if( CompilerOpts::IsPragmaDriven() )
    return true;
  
  ScopExtractor extractor(project);

  for( vector<SgNode*>::iterator it = extractor.getScopRoots().begin() ; it != extractor.getScopRoots().end();  it++ ){
    if( CompilerOpts::IsDetectOnly()){
      printf("[IE] : Scope : %s\n",(*it)->unparseToString().c_str());
    }
    deque<SgForStatement*> target_loops;
    int highest_depth = FindPartitionableLoops(isSgStatement(*it),target_loops);
    set<SgVariableSymbol*> data_arrays,indirection_arrays,indirection_scalars,parameters,iterators,data_scalars,outer_vars,privatizable;
    if( highest_depth != -1 &&  AnalysePartitionableLoops(*it,target_loops,data_arrays,data_scalars,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars,privatizable ) ){
      SgNode* codegen_initial = *it;
      SgNode* codegen_root = *it;
      //Add any surrounding loops to the codegen body
      //Outer loops might be for or while loops but must not contain references to data arrays, indirection arrays or indirection scalars
      while( isSgForStatement(codegen_root->get_parent()) || isSgWhileStmt(codegen_root->get_parent()) ){
	SgNode* parent = codegen_root->get_parent();
	if( isSgForStatement(parent) ){
	  SgForStatement* surrounding_for = isSgForStatement(parent);
	  if( CheckOuterVars(surrounding_for->get_for_init_stmt(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) 
	      && CheckOuterVars(surrounding_for->get_test(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) 
	      && CheckOuterVars(surrounding_for->get_increment(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) )
	    codegen_root = parent;
	  else
	    break;
	}
	else if( isSgWhileStmt(parent) ){
	  SgWhileStmt* surrounding_while = isSgWhileStmt(parent);
	  if( CheckOuterVars(surrounding_while->get_condition(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars))
	    codegen_root = parent;
	  else
	    break;
	}
      }
      //Ignore scope if there is no surrounding loop
      if( highest_depth == 0 && codegen_initial == codegen_root ){
	if( CompilerOpts::IsDetectOnly() ){
	  printf("[IE] : Ignoring Scope \n") ;
	}
      }
      else{
	if( CompilerOpts::IsDetectOnly() ){
	  printf("[IE] : Codegen : %s\n",codegen_root->unparseToString().c_str());
	}
	else{
	  codegen* new_codegen = new codegen(isSgStatement(codegen_root),target_loops,data_arrays,data_scalars,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars,privatizable,array_info);
	  target_sequence.push_back(new_codegen);
	}
      }
    }
  }
  if(  target_sequence.size() > 0 )
    return true;
  else 
    return false;
}



int driver::FindPartitionableLoops(SgStatement* curr_stmt, deque<SgForStatement*>& target_loops)
{
  int highest_depth = -1;
  if( isSgBasicBlock(curr_stmt)) {
    SgBasicBlock* curr_bb = isSgBasicBlock(curr_stmt);
    SgStatementPtrList& curr_list = curr_bb->get_statements();
    SgStatementPtrList::iterator curr_iter;
    for( curr_iter = curr_list.begin() ; curr_iter != curr_list.end() ; curr_iter++) {
      int new_depth = FindPartitionableLoops(*curr_iter,target_loops);
      if( new_depth != -1 )
	highest_depth = ( highest_depth == -1 ? new_depth : MIN(highest_depth,new_depth) );
    }
  }
  else if( isSgForStatement(curr_stmt) ){
    SgForStatement* curr_for = isSgForStatement(curr_stmt);
    partitionable_loop* curr_part_loop = new partitionable_loop(curr_for);
    if(!curr_part_loop->CheckProperties()){
      delete curr_part_loop;
      int new_depth = FindPartitionableLoops(curr_for->get_loop_body(),target_loops);
      highest_depth = ( new_depth == -1 ? -1 : new_depth + 1);
    }
    else{
      //printf("[IE] : Partitionable Loop: %s\n",curr_for->unparseToString().c_str());
      target_loops.push_back(curr_for);
      curr_for->setAttribute("partitionable_loop",curr_part_loop);
      highest_depth = 0;
    }
  }
  else if( isSgWhileStmt(curr_stmt) ){
    SgWhileStmt* curr_while = isSgWhileStmt(curr_stmt);
    int new_depth = FindPartitionableLoops(curr_while->get_body(),target_loops) ;
    highest_depth = ( new_depth == -1 ? -1 : new_depth + 1);
  }
  return highest_depth;
}


bool driver::AnalysePartitionableLoops(SgNode* root_node, deque<SgForStatement*>& target_loops, set<SgVariableSymbol*>& data_arrays,  set<SgVariableSymbol*>& data_scalars,  set<SgVariableSymbol*>& indirection_arrays, set<SgVariableSymbol*>& indirection_scalars, set<SgVariableSymbol*>& iterators, set<SgVariableSymbol*>& parameters, set<SgVariableSymbol*>& outer_vars, set<SgVariableSymbol*>& privatizable)
{
  bool check_property = true;
  for( deque<SgForStatement*>::iterator it = target_loops.begin() ; it != target_loops.end() ; it++ ){
    SgForStatement* curr_for = *it;
    assert(curr_for->attributeExists("partitionable_loop"));
    partitionable_loop* curr_part_loop = static_cast<partitionable_loop*>(curr_for->getAttribute("partitionable_loop"));

    //Check that no indirection array is used as a data array
    const set<SgVariableSymbol*>& curr_ind_arrays = curr_part_loop->getIndirectionArrays();
    for( set<SgVariableSymbol*>::const_iterator jt = curr_ind_arrays.begin() ; jt != curr_ind_arrays.end() ; jt++ )
      if( data_arrays.find(*jt) == data_arrays.end() ){
	indirection_arrays.insert(*jt);
      } 
      else{
	check_property = false;
	break;
      }
    if( !check_property )
      break;

    //Check that no data scalar is used as an indirection scalar or loop iterator or parameter
    const set<SgVariableSymbol*>& curr_read_data_scalars = curr_part_loop->getReadDataScalars();
    for( set<SgVariableSymbol*>::const_iterator jt = curr_read_data_scalars.begin() ; jt != curr_read_data_scalars.end() ; jt++ )
      if( indirection_scalars.find(*jt) == indirection_scalars.end() && iterators.find(*jt) == iterators.end() && parameters.find(*jt) == parameters.end() ){
	data_scalars.insert(*jt);
      } 
      else{
	check_property = false;
	break;
      }    
    if( !check_property )
      break;
    const set<SgVariableSymbol*>& curr_write_data_scalars = curr_part_loop->getWriteDataScalars();
    for( set<SgVariableSymbol*>::const_iterator jt = curr_write_data_scalars.begin() ; jt != curr_write_data_scalars.end() ; jt++ )
      if( indirection_scalars.find(*jt) == indirection_scalars.end() && iterators.find(*jt) == iterators.end() && parameters.find(*jt) == parameters.end() ){
	data_scalars.insert(*jt);
      } 
      else{
	check_property = false;
	break;
      }
    if( !check_property )
      break;

    // Check that no indirection scalar or loop iterator is used as a parameter.
    const set<SgVariableSymbol*>& curr_params = curr_part_loop->getParameters();
    for( set<SgVariableSymbol*>::const_iterator jt = curr_params.begin() ; jt != curr_params.end() ; jt++ )
      if( indirection_scalars.find(*jt) == indirection_scalars.end() && iterators.find(*jt) == iterators.end() ){
	parameters.insert(*jt);
      } 
      else{
	check_property = false;
	break;
      }    
    if( !check_property )
      break;
    else{
      data_arrays.insert(curr_part_loop->getReadDataArrays().begin(),curr_part_loop->getReadDataArrays().end());
      data_arrays.insert(curr_part_loop->getWriteDataArrays().begin(),curr_part_loop->getWriteDataArrays().end());
      data_arrays.insert(curr_part_loop->getLocalReadDataArrays().begin(),curr_part_loop->getLocalReadDataArrays().end());
      data_arrays.insert(curr_part_loop->getLocalWriteDataArrays().begin(),curr_part_loop->getLocalWriteDataArrays().end());
      indirection_scalars.insert(curr_part_loop->getIndirectionScalars().begin(),curr_part_loop->getIndirectionScalars().end());
      iterators.insert(curr_part_loop->getIterators().begin(),curr_part_loop->getIterators().end());
      privatizable.insert(curr_part_loop->getPrivatizable().begin(),curr_part_loop->getPrivatizable().end());
    }
  }

  if( !check_property)
    return false;

  if( !CompilerOpts::IsPragmaDriven() ){
    //Check that no array or indirection scalar is used in any statement outside of partitionable loops
    vector<SgStatement*> all_stmts_vect = SageInterface::querySubTree<SgStatement>(root_node);
    list<SgStatement*> all_stmts;
    all_stmts.insert(all_stmts.begin(),all_stmts_vect.begin(),all_stmts_vect.end());
    all_stmts_vect.clear();
    while( all_stmts.size() != 0 ){
      SgStatement* curr_stmt = all_stmts.front();
      //printf("Nstmts: %d, Stmt: %s\n",all_stmts.size(),curr_stmt->unparseToString().c_str());
      if( isSgExprStatement(curr_stmt) ){
	if( !(curr_stmt)->attributeExists("assignment_type") ){
	  if( !CheckOuterVars((curr_stmt),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) ){
	    check_property = false;
	    break;
	  }
	}
	vector<SgStatement*> inner_stmts = SageInterface::querySubTree<SgStatement>(curr_stmt);
	for( vector<SgStatement*>::iterator it = inner_stmts.begin() ; it != inner_stmts.end() ; it++ )
	  all_stmts.remove(*it);
      }
      else if( isSgForStatement(curr_stmt) ){
	SgForStatement* curr_for = isSgForStatement(curr_stmt);
	if( !(curr_stmt)->attributeExists("partitionable_loop") && !(curr_stmt)->attributeExists("inner_loop") ){
	  if( !CheckOuterVars(curr_for->get_for_init_stmt(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) ){
	    check_property = false;
	    break;
	  }
	  else if( !CheckOuterVars(curr_for->get_test(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) ){
	    check_property = false;
	    break;
	  }
	  else if( !CheckOuterVars(curr_for->get_increment(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) ){
	    check_property = false;
	    break;
	  }
	}
	vector<SgStatement*> init_stmts = SageInterface::querySubTree<SgStatement>(curr_for->get_for_init_stmt());
	for( vector<SgStatement*>::iterator it = init_stmts.begin() ; it != init_stmts.end() ; it++ )
	  all_stmts.remove(*it);
	vector<SgStatement*> test_stmts = SageInterface::querySubTree<SgStatement>(curr_for->get_test());
	for( vector<SgStatement*>::iterator it = test_stmts.begin() ; it != test_stmts.end() ; it++ )
	  all_stmts.remove(*it);
	vector<SgStatement*> inc_stmts = SageInterface::querySubTree<SgStatement>(curr_for->get_increment());
	for( vector<SgStatement*>::iterator it = inc_stmts.begin() ; it != inc_stmts.end() ; it++ )
	  all_stmts.remove(*it);
      }
      else if( isSgWhileStmt(curr_stmt) ){
	SgWhileStmt* curr_while = isSgWhileStmt(curr_stmt);
	if( !CheckOuterVars(curr_while->get_condition(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) ){
	  check_property = false;
	  break;
	}
	vector<SgStatement*> condn_stmts = SageInterface::querySubTree<SgStatement>(curr_while->get_condition());
	for( vector<SgStatement*>::iterator it = condn_stmts.begin() ; it != condn_stmts.end() ; it++ )
	  all_stmts.remove(*it);
      }
      else if( isSgIfStmt(curr_stmt) ){
	if( !(curr_stmt)->attributeExists("then_branch") ){
	  if( !CheckOuterVars(isSgIfStmt(curr_stmt)->get_conditional(),data_arrays,indirection_arrays,indirection_scalars,iterators,parameters,outer_vars) ){
	    check_property = false;
	    break;
	  }
	}
	vector<SgStatement*> condn_stmts = SageInterface::querySubTree<SgStatement>(isSgIfStmt(curr_stmt)->get_conditional());
	for( vector<SgStatement*>::iterator it = condn_stmts.begin() ; it != condn_stmts.end() ; it++ )
	  all_stmts.remove(*it);
      }
      all_stmts.remove(curr_stmt);
    }
  }
  return check_property;
}


bool driver::CheckOuterVars(SgStatement* root_node, set<SgVariableSymbol*>& data_arrays, set<SgVariableSymbol*>& indirection_arrays, set<SgVariableSymbol*>& indirection_scalars, set<SgVariableSymbol*>& iterators, set<SgVariableSymbol*>& parameters, set<SgVariableSymbol*>& outer_vars)
{
  bool check_property = true;
  vector<SgNode*> read_refs;
  vector<SgNode*> write_refs;
  if( !SageTools::collectReadWriteRefs(root_node,read_refs,write_refs,true) )
    return false;
  
  //No read ref must be an indirection array/scalar or data array or loop iterator
  for( vector<SgNode*>::iterator it = read_refs.begin() ; it != read_refs.end() ; it++ ){
    if( isSgVarRefExp(*it) ){
      SgVariableSymbol* curr_symbol = isSgVarRefExp(*it)->get_symbol();
      if( indirection_scalars.find(curr_symbol) != indirection_scalars.end() || iterators.find(curr_symbol) != iterators.end() ){
	check_property = false;
	break;
      }
      else
	outer_vars.insert(curr_symbol);
    }
    else if( isSgPntrArrRefExp(*it) ){
      //In OpenMP mode, no array references here.
      if( CompilerOpts::IsOpenMPMode()){
	check_property = false;
	break;
      }
      SgPntrArrRefExp* arrref = isSgPntrArrRefExp(*it);
      SgExpression* lhs_exp = arrref->get_lhs_operand();
      while( isSgPntrArrRefExp(lhs_exp) ){
	lhs_exp = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
      }
      if( isSgVarRefExp(lhs_exp) ){
	SgVariableSymbol* curr_symbol = isSgVarRefExp(lhs_exp)->get_symbol();
	if( indirection_arrays.find(curr_symbol) != indirection_arrays.end() || data_arrays.find(curr_symbol) != data_arrays.end() ){
	  check_property = false;
	  break;
	}
	else
	  outer_vars.insert(curr_symbol);
      }
      else{
	check_property = false;
	break;
      }
    }
  }
  if( !check_property ) 
    return false;

  //No write ref must be an indirection array/scalar or data array or parameter or loop iterator
  for( vector<SgNode*>::iterator it = write_refs.begin() ; it != write_refs.end() ; it++ ){
    if( isSgVarRefExp(*it) ){
      SgVariableSymbol* curr_symbol = isSgVarRefExp(*it)->get_symbol();
      if( indirection_scalars.find(curr_symbol) != indirection_scalars.end() || iterators.find(curr_symbol) != iterators.end() || parameters.find(curr_symbol) != parameters.end() ){
	check_property = false;
	break;
      }
      else
	outer_vars.insert(curr_symbol);
    }
    else if( isSgPntrArrRefExp(*it) ){
      //In OpenMP mode no array references
      if( CompilerOpts::IsOpenMPMode() ){
	check_property = false;
	break;
      }
      SgPntrArrRefExp* arrref = isSgPntrArrRefExp(*it);
      SgExpression* lhs_exp = arrref->get_lhs_operand();
      while( isSgPntrArrRefExp(lhs_exp) ){
	lhs_exp = isSgPntrArrRefExp(lhs_exp)->get_lhs_operand();
      }
      if( isSgVarRefExp(lhs_exp) ){
	SgVariableSymbol* curr_symbol = isSgVarRefExp(lhs_exp)->get_symbol();
	if( indirection_arrays.find(curr_symbol) != indirection_arrays.end() || data_arrays.find(curr_symbol) != data_arrays.end() ){
	  check_property = false;
	  break;
	}
	else
	  outer_vars.insert(curr_symbol);
      }
      else{
	check_property = false;
	break;
      }
    }
  }
  // if( !check_property ) 
  //   return false;
  return check_property;
}

bool driver::CheckOuterVars(SgExpression* curr_exp, set<SgVariableSymbol*>& data_arrays, set<SgVariableSymbol*>& indirection_arrays, set<SgVariableSymbol*>& indirection_scalars, set<SgVariableSymbol*>& iterators, set<SgVariableSymbol*>& parameters, set<SgVariableSymbol*>& outer_vars)  {
  bool check_property = true;
  if( CompilerOpts::IsOpenMPMode() ){
    //No array references allowed in OpenMP mode.
    vector<SgPntrArrRefExp*> array_refs = SageInterface::querySubTree<SgPntrArrRefExp>(curr_exp);
    if( array_refs.size() != 0 ){
      return false;
    }
  }
  vector<SgVarRefExp*> ref_vars = SageInterface::querySubTree<SgVarRefExp>(curr_exp);
  for( vector<SgVarRefExp*>::iterator jt = ref_vars.begin() ; jt != ref_vars.end() ; jt++ ){
    SgVariableSymbol* curr_symbol = (*jt)->get_symbol();
    //THis is too restrictive. Parametes in expressions should be allowed, but currently this doesnt check for side-effect free expressions.
    //But works out since all expressions apart from increment of for are treated as statements within ROSE.
    if( data_arrays.find(curr_symbol) != data_arrays.end() || indirection_arrays.find(curr_symbol) != indirection_arrays.end() || indirection_scalars.find(curr_symbol) != indirection_scalars.end() || iterators.find(curr_symbol) != iterators.end() || parameters.find(curr_symbol) != parameters.end() ){
      check_property = false;
      break;
    }
    else
      outer_vars.insert(curr_symbol);
  }
  return check_property;
}
