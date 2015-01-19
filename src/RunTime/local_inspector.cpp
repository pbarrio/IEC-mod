/*
 * local_inspector.cpp: This file is part of the IEC project.
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
 * @file: local_inspector.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/local_inspector.hpp"

using namespace std;

/**
 * \brief Inspector bits specific to each thread
 *
 * I should change this to one local inspector per pipeline workgroup
 */
local_inspector** local_inspector::all_local_inspectors = NULL;

/**
 * /brief Constructor
 *
 * /param np Number of processes
 * /param pid This process' id
 * /param nc Number of communicators (one per loop, get rid of it).
 */
local_inspector::local_inspector(int np, int pid, int nc):
	nprocs(np),
	proc_id(pid),
	myid(pid)
{ 
#ifndef NDEBUG
	printf("GID:%d,Local_Inspector:%p\n",myid,this);
	fflush(stdout);
#endif
	for( int i = 0 ; i < nc ; i++ ){
		local_comm* new_comm = new local_comm(i, np, pid);
		all_comm.push_back(new_comm);
	}
#ifndef NDEBUG
	char df_name[18];
	sprintf(df_name,"local_data_%d.dat",myid);
	data_file = fopen(df_name,"w");
	char cf_name[18];
	sprintf(cf_name,"local_comm_%d.dat",myid);
	comm_file = fopen(cf_name,"w");
#endif  
}

local_inspector::~local_inspector()
{
	for( deque<local_data*>::iterator it = all_data.begin() ; it!= all_data.end() ; it++ )
		delete (*it);
	for( deque<local_comm*>::iterator it = all_comm.begin() ; it!=all_comm.end() ; it++ )
		delete(*it);
#ifndef NDEBUG
	fclose(data_file);
	fclose(comm_file);
#endif
}

void local_inspector::PopulateGlobalArrays()
{
	for( deque<local_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
		(*it)->PopulateGlobalArray();
}
