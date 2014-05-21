/*
 * local_inspector.hpp: This file is part of the IEC project.
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
 * @file: local_inspector.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __LOCAL_INSPECTOR_HPP__
#define __LOCAL_INSPECTOR_HPP__

#include <deque>
#include <cassert>
#include <cstdio>

#include "RunTime/local_data.hpp"
#include "RunTime/local_comm.hpp"
class inspector;

class local_inspector{

 private:
  
  const int nprocs;
  
  const int nthreads;

  const int proc_id;
  
  const int thread_id;

  const int myid;

#ifndef NDEBUG
  FILE* data_file;

  FILE* comm_file;
#endif
  
  std::deque<local_data*> all_data;

  std::deque<local_comm*> all_comm;

  local_inspector(int,int,int,int,int);

  static local_inspector** all_local_inspectors;

 public:
  
  ~local_inspector();
  
  static local_inspector* instance(int tid){
    return all_local_inspectors[tid];
  }

  inline void SetupLocalArray(int dn) {
    all_data[dn]->SetupLocalArray();
  }
  
  inline int GetLocalDataSize(int dn) const{
    return all_data[dn]->GetLocalSize();
  }
  
  inline void AddLocalData(int mn, int stride_size, const net** ddni, int oas, bool iro, bool ic){
    local_data* new_data = new local_data_double(mn, nprocs*nthreads,myid,proc_id,nthreads,stride_size,ddni,oas,iro,ic); 
    all_data.push_back(new_data);
  }

  inline void AddReadArray(int in, int an){
    all_comm[in]->read_arrays.push_back(all_data[an]);
  }

  inline void AddWriteArray(int in, int an){
    all_comm[in]->write_arrays.push_back(all_data[an]);
  }
  
  inline void AddIndexAccessed(int dn, int ind, int at){
    all_data[dn]->AddIndexAccessed(ind,at);
  }

  inline void PopulateLocalArray(int an, double* lb, double* oa, int st){
    all_data[an]->PopulateLocalArray(lb,oa,st);
  }

  inline void InsertDirectAccess(int an, int* v, int n){
    all_data[an]->InsertDirectAccess(v,n);
  }
  
  inline void InsertIndirectAccess(int an, int* v, int n){
    all_data[an]->InsertIndirectAccess(v,n);
  }

  inline void RenumberAccessArray(int an, int as, int* aa){
    all_data[an]->RenumberAccessArray(as,aa);
  }

  inline void RenumberOffsetArray(int an, int as, int* aa, int* la){
    all_data[an]->RenumberOffsetArray(as,aa,la);
  }
    
  inline void GenerateGhosts(){
    int counter =0;
    for( std::deque<local_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ ,counter++){
      (*it)->GenerateGhosts();
    }
  }

  inline void GenerateOwned(){
    int counter = 0;
    for( std::deque<local_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ ,counter++){
	(*it)->GenerateOwned();
      }
  }

  inline void InitWriteGhosts(int cn){
    all_comm[cn]->InitWriteGhosts();
  }
  
  void PopulateGlobalArrays();

  inline void SetLocalArray(int an, void* la) { all_data[an]->SetLocalArray(la); }

  inline int GetLocalBlockOffset(int an) { return all_data[an]->GetLocalBlockOffset() ; }

  void print_data(){
#ifndef NDEBUG
    for( std::deque<local_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
      (*it)->print_data(data_file);
    for( std::deque<local_comm*>::iterator it = all_comm.begin() ; it != all_comm.end() ; it++ )
      (*it)->print_comm(comm_file);
#endif
  }

  friend class inspector;
};



#endif
