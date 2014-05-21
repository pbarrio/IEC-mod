/*
 * driver.hpp: This file is part of the IEC project.
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
 * @file: driver.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __DRIVER_HPP__
#define __DRIVER_HPP__

#include <rose.h>
#include "CompileTime/codegen.hpp"
#include <set>
#include <deque>
#include <string.h>
/* #include <CXX_Grammar.h> */


class ie_init : public AstSimpleProcessing{
 private:
  SgGlobal* global_scope;

 protected:
  void visit(SgNode*);

 public:

  ie_init(SgGlobal*);
};


class driver: public AstSimpleProcessing {
  
private:

  std::list<codegen*> target_sequence;

  std::deque<array_string*> array_info;

  static driver* singleton_driver;

  // driver(int,char**) ;

  ~driver();

  int FindPartitionableLoops(SgStatement*,std::deque<SgForStatement*>&);

  bool AnalysePartitionableLoops(SgNode*,std::deque<SgForStatement*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&); 
  bool CheckOuterVars(SgStatement*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&); 
  bool CheckOuterVars(SgExpression*,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&,std::set<SgVariableSymbol*>&); 

  // bool openmp_mode;

  // bool detect_only;

  // bool pragma_driven;

protected:

  void visit(SgNode*);

public:

  static driver* create_driver(){
    if (singleton_driver == NULL )
      singleton_driver = new driver;
    return singleton_driver;
  }

  static driver* get_driver(){
    return singleton_driver;
  }

  static void delete_driver(){
    if( singleton_driver )
      delete singleton_driver;
  }
  
  static bool CheckArray(SgVarRefExp*);

  bool CheckTargets(SgProject*);

  void GenerateCCode(SgProject*);
  
};

#endif
