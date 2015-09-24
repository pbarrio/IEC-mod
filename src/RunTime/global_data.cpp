/*
 * global_data.cpp: This file is part of the IEC project.
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
 * @file: global_data.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include <cassert>
#include <cstdio>

#include "RunTime/global_data.hpp"

using namespace std;

/**
 * \brief Constructor
 *
 * \param mn Identifier for this data
 * \param oas Size of the original array
 * \param of Offset: address of the first elem if we put all global data in a
 *           single buffer in ID order.
 * \param iro True if read-only array
 */
global_data::global_data(int mn, int oas, int of, bool iro):
	id(mn), orig_array_size(oas), offset(of), is_read_only(iro){

	is_constrained = false;
}


global_data::~global_data(){

	for (LoopNets::iterator netIt = data_net_info.begin(),
		     netEnd = data_net_info.begin();
	     netIt != netEnd; ++netIt){

		net** nets = netIt->second;
		for (int i = 0; i < orig_array_size; i++)
			delete nets[i];
		delete[] nets;
	}
}


/**
 * \brief Mark the array as used in a loop and initialize some required info
 */
void global_data::use_in_loop(int loopID){

	// Create a net for each element of the array to be used in partitioning
	data_net_info[loopID] = new net*[orig_array_size];
	for (int i = 0; i < orig_array_size; i++)
		data_net_info[loopID][i] = new net(id, i, stride_size, offset + i);
}

/**
 * \brief Constructor
 *
 * \param mn Identifier for this data
 * \param oas Size of the original array
 * \param of Offset: address of the first elem if we put all global data in a
 *           single buffer in ID order.
 * \param iro True if read-only array
 */
global_data_double::global_data_double(int mn, int oas, int of, bool iro):
	global_data(mn, oas, of, iro){

	// Set default stride between elements
	stride_size = sizeof(double);
}

void global_data_double::SetStride(int st){

	stride_size = st * sizeof(double);

	for (LoopNets::iterator netIt = data_net_info.begin(),
		     netEnd = data_net_info.end();
	     netIt != netEnd; ++netIt){

		net** nets = netIt->second;
		for (int i = 0; i < orig_array_size; i++)
			nets[i]->weight = stride_size;
	}
}

global_data_double::~global_data_double(){}
