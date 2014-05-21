/*
 * opts.hpp: This file is part of the IEC project.
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
 * @file: opts.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __OPTS_HPP__
#define __OPTS_HPP__

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>

class CompilerOpts{
  
private:
  bool openmp_mode;

  bool detect_only;

  bool pragma_driven;

  int partitioner_type;
  
  static CompilerOpts* singleton_instance;

  CompilerOpts(int argc, char** argv):
    openmp_mode(true),
    detect_only(false),
    pragma_driven(false),
    partitioner_type(0)
  {
    bool solver_specified = false;
    for( int i = 0 ; i < argc ; i++ ){
      if( strcmp(argv[i],"--mode-mpi") == 0 ){
	openmp_mode = false;
	continue;
      }
      if( strcmp(argv[i],"--detect-only") == 0 ){
	detect_only = true;
	continue;
      }
      if( strcmp(argv[i],"--use-metis") == 0 ){
	if( solver_specified ){
	  fprintf(stderr,"Command Line Error! Multiple partitioner specified\n");
	  exit(1);
	}
	partitioner_type = 1;
	solver_specified = true;
	continue;
      }
      if( strcmp(argv[i],"--use-block") == 0 ){
	if( solver_specified ){
	  fprintf(stderr,"Command Line Error! Multiple partitioner specified\n");
	  exit(1);
	}
	partitioner_type = 2;
	solver_specified = true;
	continue;
      }
      if( strcmp(argv[i],"--use-patoh") == 0 ){
	solver_specified = true;
	continue;
      }
    }
  }

public:
  
  static bool IsOpenMPMode() {
    return singleton_instance->openmp_mode;
  }

  static bool IsDetectOnly() {
    return singleton_instance->detect_only;
  }

  static bool IsPragmaDriven() {
    return singleton_instance->pragma_driven;
  }

  static void SetPragmaDriven() {
    singleton_instance->pragma_driven = true;
  }

  static int PartitionerType() {
    return singleton_instance->partitioner_type;
  }
  
  static void Init(int argc,char** argv){
    assert(singleton_instance == NULL );
    singleton_instance = new CompilerOpts(argc,argv);
  }
};



#endif
