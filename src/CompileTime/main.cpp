/*
 * main.cpp: This file is part of the IEC project.
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
 * @file: main.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "CompileTime/opts.hpp"
#include "CompileTime/driver.hpp"

int main(int argc, char** argv){

  SgProject* project = frontend(argc,argv);
  //generateDOT(*project);
  
  CompilerOpts::Init(argc,argv);

  driver* main_driver = driver::create_driver();

  main_driver->traverseInputFiles(project,preorder);

  if( main_driver->CheckTargets(project) ){
    main_driver->GenerateCCode(project);
    project->unparse();
  }
  
  driver::delete_driver();

  // ScopExtractor extractor(project);

  return 0;
}
