/*
 * loops.hpp: This file is part of the IEC project.
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
 * @file: loops.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __LOOPS_HPP__
#define __LOOPS_HPP__

#include <rose.h>
#include <cassert>
#include <deque>
#include <set>

class array_details;
class safe_scalar_details;
class codegen;

namespace OmpSupport{
  class OmpAttribute;
};


enum descriptor_type{
  IA_DESCRIPTOR,
  OA_DESCRIPTOR,
  CONSTANT
};


class IAssignmentAttr : public AstAttribute{
private:
  const bool is_i_assignment;

public:
  IAssignmentAttr(bool val): is_i_assignment(val) {}
 
  inline bool IsIAssignment() const  { return is_i_assignment; }
};


class branch : public AstAttribute{

protected:

  bool is_used;
  
  std::set<SgVariableSymbol*> written_iscalars;

public:

  virtual bool IsLoop() const=0;

  inline void AddWriteIScalar(SgVariableSymbol* a){ written_iscalars.insert(a); }

  void GenerateSetGaurds(SgStatement*,SgExpression*) const;

  inline bool IsUsed() const { return is_used; }

  inline void SetUsed() { is_used = true; }

  branch(bool);

  virtual ~branch() { }

  friend class codegen;

  friend class partitionable_loop;
};




class scope : public AstAttribute{

public:

  scope(std::deque<branch*>& curr_scope){ scope_stack.insert(scope_stack.begin(),curr_scope.begin(),curr_scope.end()); }

  std::deque<branch*> const& getScope() const { return scope_stack; }

private:
 
  std::deque<branch*> scope_stack;

};



class loop : public branch{

 protected:
  
  SgVariableSymbol* orig_iterator;

  SgForStatement* orig_for;

  SgExpression* lb_exp;

  SgExpression* ub_exp;

  SgForStatement* executor_for;

  const int depth;

 public:

  virtual ~loop() { }
  
  loop(SgForStatement* , const int);

  inline int GetDepth() const { return depth; }

  inline SgVariableSymbol* GetIterator() const { return orig_iterator; }

  inline SgExpression* GetLBExp() const { return lb_exp; }

  inline SgExpression* GetUBExp() const { return ub_exp; }

  inline void SetExecutorFor(SgForStatement* new_for) { assert(executor_for == NULL) ; executor_for = new_for; }

  inline void InsertBeforeExecutor(SgStatement* new_stmt) const { SageInterface::insertStatementBefore(executor_for,new_stmt); }

  inline bool IsLoop() const { return true; } 

  virtual bool IsPartitionableLoop() const=0;

  virtual int GetLoopNum() const=0;
  
  virtual int GetGroupNum() const=0;

  virtual void InitCounters(SgStatement*)=0;

  virtual void InitBoundsArrays(SgStatement*,SgFunctionSymbol*)=0;

  virtual SgVariableDeclaration* GetBodyCounter() const=0;

  virtual SgVariableDeclaration* GetLoopCounter() const=0;

  virtual std::pair<SgStatement*,SgStatement*> GenerateCArrayStore() const=0;

  virtual SgIfStmt* GenerateIfFirstIterStmt(SgStatement*) const=0;

  virtual void GenerateLoopCounterIncrement(SgStatement*) const=0;

  virtual void GenerateBodyCounterIncrement(SgStatement*) const=0;

  virtual void GenerateLoopCounterIncrementAfter(SgStatement*) const=0;

  virtual void GenerateBodyCounterIncrementAfter(SgStatement*) const=0;
  
  virtual void GenerateResetCounters(SgStatement*) const=0;

  virtual void GenerateResetCountersAfter(SgStatement*&) const=0;

  virtual SgVariableDeclaration* GetLBArray() const=0;

  void AddPrivateVariables(OmpSupport::OmpAttribute*) const;

  inline std::string unparseToString() {
    return orig_for->unparseToString();
  }
		       
  friend class codegen;

  friend class partitionable_loop;
};



class partitionable_loop: public loop{

 private:
  
  int my_num;

  int my_group;

  bool write_communication;

  bool read_communication;

  std::set<SgVariableSymbol*> read_data_arrays;

  std::set<SgVariableSymbol*> write_data_arrays;

  std::set<SgVariableSymbol*> local_read_arrays;

  std::set<SgVariableSymbol*> local_write_arrays;

  std::set<SgVariableSymbol*> read_data_scalars;

  std::set<SgVariableSymbol*> write_data_scalars;

  std::set<SgVariableSymbol*> indirection_scalars;

  std::set<SgVariableSymbol*> indirection_arrays;

  std::set<SgVariableSymbol*> iterators;

  std::set<SgVariableSymbol*> parameters;

  std::set<SgVariableSymbol*> privatizable;

  std::deque<SgExprStatement*> assignment_stmts;

  std::deque<SgExprStatement*> i_assignment_stmts;

  SgNodePtrList inner_loop_nodes;

  SgNodePtrList inner_condn_nodes;

  SgVariableDeclaration* local_block_start;

  SgVariableDeclaration* local_block_end;

  SgVariableDeclaration* body_counter;

  SgVariableDeclaration* nlocal_iters;

  SgVariableDeclaration* status_array;

  static SgVariableDeclaration* iter_num_count;

  bool GetIScalarArrays(SgStatement*,std::deque<branch*>&);

  bool GetIScalarArraysBasicBlock(SgBasicBlock*,std::deque<branch*>&);
  
  bool GetIScalarArraysFor(SgForStatement*,std::deque<branch*>&);

  bool GetIScalarArraysIf(SgIfStmt*,std::deque<branch*>&);

  bool GetIScalarArraysExprStatement(SgExprStatement*,std::deque<branch*>&);

  bool GetIScalarArraysExpression(SgExpression*,std::deque<branch*>&);

  bool CheckPrivatizable(SgStatement*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&);

  bool CheckPrivatizableBasicBlock(SgBasicBlock*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&);
  
  bool CheckPrivatizableFor(SgForStatement*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&);

  bool CheckPrivatizableIf(SgIfStmt*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&);

  bool CheckPrivatizableExprStatement(SgExprStatement*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&);

  bool CheckPrivatizableExpression(SgExpression*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,scope*);

  bool AddIArraysScalars(SgExpression*,const std::deque<branch*>&);

  bool CheckVariables(SgExpression*,std::set<SgVariableSymbol*>&, scope*);

  bool CheckInnerDimensions() const;

  bool CheckInnerExpression(SgExpression*,scope*) const;

  void AnalyseTargetLoop();

  void AnalyseTargetStatement(SgStatement*);

  void AnalyseTargetBasicBlock(SgBasicBlock*);

  void AnalyseTargetFor(SgForStatement*);

  void AnalyseTargetIfStatement(SgIfStmt*);

  void AnalyseTargetExprStatement(SgExprStatement*);

  void AnalyseTargetExpression(SgExpression*,scope*, bool, VariantT);

  std::pair<branch*,descriptor_type> AnalyseSafeExpression(SgExpression*,scope*);


  void GenerateCInspector(SgStatement*,SgStatement*,SgStatement*,SgStatement*);

  void GenerateCStatement(SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgExpression*,SgExpression*);

  void GenerateCBasicBlock(SgBasicBlock*,SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgExpression*,SgExpression*);

  void GenerateCIfStatement(SgIfStmt*,SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgExpression*,SgExpression*);

  void GenerateCForStatement(SgForStatement*,SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgExpression*,SgExpression*);

  void GenerateCExprStatement(SgExprStatement*,SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgExpression*,SgExpression*);

  void GenerateCExpression(SgExpression*,SgStatement*,SgStatement*,SgStatement*,SgStatement*,SgExpression*,SgExpression*);
  
  void ReplaceIndirectionArrays(SgExpression*);

  SgExpression* GenerateCSafeStatement(SgExpression*);

  void GenerateCExecutor(SgStatement*);

  void GenerateCExecutorStatement(SgStatement*,SgStatement*);

  void GenerateCExecutorBasicBlock(SgBasicBlock*,SgStatement*);

  void GenerateCExecutorForStatement(SgForStatement*, SgStatement*);

  void GenerateCExecutorIfStatement(SgIfStmt*,SgStatement*);
  
  void GenerateCExecutorExprStatement(SgExprStatement*,SgStatement*);

  void GenerateCExecutorExpression(SgExpression*,SgExpression*);

  static std::deque<partitionable_loop*> all_loops;

  void DeleteAttributes();

 public:

  partitionable_loop(SgForStatement*);

  ~partitionable_loop();

  inline bool IsPartitionableLoop() const{ return true; }

  bool CheckProperties(); 

  inline const std::set<SgVariableSymbol*>& getReadDataArrays() { return read_data_arrays; }
  inline const std::set<SgVariableSymbol*>& getWriteDataArrays() { return write_data_arrays; }
  inline const std::set<SgVariableSymbol*>& getReadDataScalars() { return read_data_scalars; }
  inline const std::set<SgVariableSymbol*>& getWriteDataScalars() { return write_data_scalars; }
  inline const std::set<SgVariableSymbol*>& getLocalReadDataArrays() { return local_read_arrays; }
  inline const std::set<SgVariableSymbol*>& getLocalWriteDataArrays() { return local_write_arrays; }
  inline const std::set<SgVariableSymbol*>& getIndirectionArrays() { return indirection_arrays; }
  inline const std::set<SgVariableSymbol*>& getIndirectionScalars() { return indirection_scalars; }
  inline const std::set<SgVariableSymbol*>& getIterators() { return iterators; }
  inline const std::set<SgVariableSymbol*>& getParameters() { return parameters; }
  inline const std::set<SgVariableSymbol*>& getPrivatizable() { return privatizable; }

  inline void CommunicateReads() { read_communication = true; }

  inline void CommunicateWrites() { write_communication = true; }

  inline void SetGroupLoopNum( int ln, int gn) { my_num = ln ; my_group = gn; }
  
  inline int GetLoopNum() const { return my_num; }
  
  inline int GetGroupNum() const { return my_group; }

  void InitStatusArray(SgStatement*);

  inline SgVariableDeclaration* GetStatusArray() const { return status_array; }

  void GenerateLocalBlock(SgStatement*,SgFunctionSymbol*);

  void ResetStatusArray(SgStatement*,SgFunctionSymbol*,SgFunctionSymbol*) const;

  void FreeStatusArray(SgStatement*,SgFunctionSymbol*) const;

  void InitCounters(SgStatement*);

  void InitBoundsArrays(SgStatement* , SgFunctionSymbol* );

  SgVariableDeclaration* GetBodyCounter() const { return body_counter; }

  SgVariableDeclaration* GetLoopCounter() const { assert(0); return NULL;  }

  SgVariableDeclaration* GetLocalIters() const { assert(nlocal_iters); return nlocal_iters; }

  inline SgVariableDeclaration* GetLocalBlockStart() const { return local_block_start; }

  inline SgVariableDeclaration* GetLocalBlockEnd() const { return local_block_end; }

  inline std::pair<SgStatement*,SgStatement*> GenerateCArrayStore() const{ assert(0); return std::pair<SgStatement*,SgStatement*>(NULL,NULL);}

  inline SgIfStmt* GenerateIfFirstIterStmt(SgStatement* a) const { assert(0); return NULL;} 

  void GenerateLoopCounterIncrement(SgStatement*a) const { assert(0); }

  void GenerateBodyCounterIncrement(SgStatement*) const;

  void GenerateLoopCounterIncrementAfter(SgStatement*a) const { assert(0); }

  void GenerateBodyCounterIncrementAfter(SgStatement*a) const { assert(0); }

  void GenerateResetCounters(SgStatement*) const;

  inline void GenerateResetCountersAfter(SgStatement*& a) const { }

  void DeleteArrays(SgStatement* , SgFunctionSymbol*) const;

  void ResetShadowScalars(SgStatement*) const;
  
  void ResetAllCounters(SgStatement*) const ;

  inline SgVariableDeclaration* GetLBArray() const { assert(0); return NULL; }

  static SgFunctionSymbol* init_write_ghosts_fn;

  static SgFunctionSymbol* communicate_reads_fn;
  
  static SgFunctionSymbol* communicate_writes_fn;

  static SgFunctionSymbol* reduce_scalar_fn;

  static SgFunctionSymbol* memset_fn;

  static SgFunctionSymbol* get_proc_iter_size_fn;

  static SgFunctionSymbol* add_vertex_fn;

  static SgFunctionSymbol* get_vertex_home_fn;

  static SgFunctionSymbol* is_known_fn;

  static SgFunctionSymbol* get_elem_fn;

  friend class codegen;
};


class inner_loop: public loop{
  
  SgVariableDeclaration* body_counter;

  SgVariableDeclaration* loop_counter;

  SgVariableDeclaration* lb_array;

  SgVariableDeclaration* ub_array;

  const int my_num;

  static int total_inner;

 public:
  
  inner_loop(SgForStatement*, const int);

  ~inner_loop() { }

  inline bool IsPartitionableLoop() const { return false; }

  inline int GetLoopNum() const { assert(0); }

  inline int GetGroupNum() const {assert(0); }

  void InitCounters(SgStatement*);

  void InitBoundsArrays(SgStatement*,SgFunctionSymbol*);

  inline SgVariableDeclaration* GetBodyCounter() const { return body_counter; }

  inline SgVariableDeclaration* GetLoopCounter() const { return loop_counter; }

  std::pair<SgStatement*,SgStatement*> GenerateCArrayStore() const;

  SgIfStmt* GenerateIfFirstIterStmt(SgStatement*) const;

  void GenerateLoopCounterIncrement(SgStatement*) const;

  void GenerateBodyCounterIncrement(SgStatement*) const;

  void GenerateLoopCounterIncrementAfter(SgStatement*) const;

  void GenerateBodyCounterIncrementAfter(SgStatement*) const;

  void GenerateResetCounters(SgStatement*) const;

  void GenerateResetCountersAfter(SgStatement*&) const;

  void DeleteBoundsArrays(SgStatement*,SgFunctionSymbol*) const;

  inline SgVariableDeclaration* GetLBArray() const { return lb_array; }

  inline SgVariableDeclaration* GetUBArray() const { return ub_array; }

  //void GenerateSetGaurds(SgStatement*,SgExpression*) const;

  static SgFunctionSymbol* is_not_known_scalar_fn;
};


#endif
