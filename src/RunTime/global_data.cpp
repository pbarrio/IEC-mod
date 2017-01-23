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
#include <iostream>

#include "RunTime/global_data.hpp"

using namespace std;

/**
 * \brief Constructor
 *
 * <param proc Our process ID
 * \param mn Identifier for this data
 * \param oas Size of the original array
 * \param of Offset: address of the first elem if we put all global data in a
 *           single buffer in ID order
 * \param iro True if read-only array
 */
global_data::global_data(int proc, int mn, int oas, int of, bool iro):
	procId(proc), id(mn), orig_array_size(oas), offset(of), is_read_only(iro){

	is_constrained = false;

	// By default, assume that this array will be last written in the current
	// loop. When we initialize the loops, we'll find out if it's not.
	lastWriteInPipeline = true;
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
 *
 * Note that there is another version of this function below with a slightly
 * different functionality.
 *
 * \param loopID ID of the loop where this array is used
 * \param isRead True if this array is read in the loop
 * \param isWrite True if this array is written in the loop
 */
void global_data::use_in_loop(int loopID, bool isRead, bool isWrite){

	isReadInLoop[loopID] = isRead;
	isWriteInLoop[loopID] = isWrite;

	// Create a net for each element of the array to be used in partitioning
	data_net_info[loopID] = new net*[orig_array_size];
	for (int i = 0; i < orig_array_size; i++)
		data_net_info[loopID][i] = new net(id, i, stride_size, offset + i);
}


/**
 * \brief Mark a position in the array as used in a loop iteration
 */
void global_data::use_in_loop(int loopID, int iter, int index){

	finalUse[loopID][index] = iter;

	// If this is the first time that we visit this index, mark it as first use
	if (!initialUse[loopID].count(index) || (initialUse[loopID][index] > index))
		initialUse[loopID][index] = iter;
}


/**
 * \brief Populate the pipeline send information for this array
 *
 * This is later used to know which parts of the array must be communicated in
 * which iteration.
 */
void global_data::pipe_calc_sends(int myLoop){

	// For all loops that use this array
	for (global_data::LoopNets::iterator netIt = data_net_info.begin(),
		     netEnd = data_net_info.end();
	     netIt != netEnd;
	     ++netIt){

		int consumerId = netIt->first;
		if (consumerId <= procId)
			continue;

		if (isReadInLoop[consumerId])
			for (int index = 0; index < orig_array_size; ++index){
				pipeSendIndexes[finalUse[myLoop][index]][consumerId]
					.push_back(index);
			}

		// If the consumer writes into the array, stop adding more consumers.
		// The next consumers should get the values from him (not from us),
		// because it is the most up-to-date value.
		if (isWriteInLoop[consumerId])
			break;
	}
}


/**
 * \brief Populate the pipeline receive information
 */
void global_data::pipe_calc_recvs(){

	producer = -1;

	// Find the latest loop that writes to this array before our loop
	for (global_data::LoopNets::iterator netIt = data_net_info.begin(),
		     netEnd = data_net_info.end();
	     netIt != netEnd;
	     ++netIt){

		int pId = netIt->first;
		if (pId >= procId)
			continue;

		if (isWriteInLoop[pId] && (pId > producer))
			producer = pId;
	}

	// For each index and iteration, find when the producer will send it to us
	for (int index = 0; index < orig_array_size; ++index){

		int lastUse = finalUse[producer][index];
		int firstUse = initialUse[procId][index];

		pipeRecvIndexes[lastUse][producer].push_back(index);

		// Mark the (producer) iteration where we will receive the data
		if (producerIter.size() <= firstUse){
			producerIter.resize(firstUse + 1, -1);
		}
		if (producerIter[firstUse] < lastUse)
			producerIter[firstUse] = lastUse;
	}
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
global_data_double::global_data_double
(int proc, int mn, int oas, int of, bool iro):
	global_data(proc, mn, oas, of, iro){

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
