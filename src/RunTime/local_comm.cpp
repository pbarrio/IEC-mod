/*
 * local_comm.cpp: This file is part of the IEC project.
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
 * @file: local_comm.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/local_comm.hpp"
#include "RunTime/local_data.hpp"
#include <cassert>
#include <string.h>

using namespace std;

local_comm::local_comm(int np, int pid):
	nprocs(np),
	proc_id(pid),
	myid(pid),
	read_send_offset(new int[np]),
	read_recv_offset(new int[np]),
	write_send_offset(new int[np]),
	write_recv_offset(new int[np]),
	read_send_count(new int[np]),
	read_recv_count(new int[np]),
	write_send_count(new int[np]),
	write_recv_count(new int[np]){

	memset(read_send_count, 0, sizeof(int) * nprocs);
	memset(read_send_offset, 0, sizeof(int) * nprocs);
	memset(read_recv_count, 0, sizeof(int) * nprocs);
	memset(read_recv_offset, 0, sizeof(int) * nprocs);
	memset(write_send_count, 0, sizeof(int) * nprocs);
	memset(write_send_offset, 0, sizeof(int) * nprocs);
	memset(write_recv_count, 0, sizeof(int) * nprocs);
	memset(write_recv_offset, 0, sizeof(int) * nprocs);
#ifdef COMM_TIME
	read_comm_time = 0.0;
	write_comm_time = 0.0;
	read_time1 = 0.0;
	write_time1 = 0.0;
	read_time2 = 0.0;
	write_time2 = 0.0;
	read_time3 = 0.0;
	write_time3 = 0.0;
#endif

}

local_comm::~local_comm(){

#ifdef COMM_TIME
	if (thread_id == 0 && read_comm_time > 0.0)
		printf("MXC:GID:%d,ReadCommTime:%5.4le\n", myid, read_comm_time);
	if (thread_id == 0 && write_comm_time > 0.0)
		printf("MXC:GID:%d,WriteCommTime:%5.4le\n", myid, read_comm_time);
#endif
	delete[] read_send_offset;
	delete[] read_recv_offset;
	delete[] write_send_offset;
	delete[] write_recv_offset;
	delete[] read_send_count;
	delete[] read_recv_count;
	delete[] write_send_count;
	delete[] write_recv_count;
}

/**
 * \brief Unused in quake
 */
int local_comm::GetReadSendCount(const int dest, const int curr_offset){

	assert(read_send_count[dest] == 0);
	read_send_offset[dest] = curr_offset;
	int send_count = 0;
	for (vector<local_data*>::iterator it = read_arrays.begin();
	     it != read_arrays.end(); it++)
		send_count += (*it)->GetOwnedCount(dest);
	read_send_count[dest] = send_count;

	if (dest == proc_id)
		return 0;
	else
		return send_count;
}

int local_comm::GetReadRecvCount(const int source, const int curr_offset){

	assert(read_recv_count[source] == 0);
	read_recv_offset[source] = curr_offset;
	int recv_count = 0;
	for (vector<local_data*>::iterator it = read_arrays.begin();
	     it != read_arrays.end(); it++)
		recv_count += (*it)->GetGhostsCount(source);
	read_recv_count[source] = recv_count;

	if (source == proc_id)
		return 0;
	else
		return recv_count;
}

/**
 * \brief Calculate size of the ghosts to be sent to one process
 *
 * The modified ghosts (those that belong to a "write" array)
 * must be sent back to their owners for aggregation of results.
 * This function also updates "write_send_offset", which tracks the
 * starting position of the ghosts for each destination process
 * inside the array (bad design... bad luck!).
 *
 * \param dest Destination process (owner of the ghosts)
 * \param curr_offset Size of the array so far -> start of this batch
 */
int local_comm::GetWriteSendCount(const int dest, const int curr_offset){

	assert(write_send_count[dest] == 0);
	write_send_offset[dest] = curr_offset;
	int send_count = 0;
	for (vector<local_data*>::iterator it = write_arrays.begin();
	     it != write_arrays.end(); it++){
		send_count += (*it)->GetGhostsCount(dest);
	}
	write_send_count[dest] = send_count;

	if (dest == proc_id)
		return 0;
	else
		return send_count;
}


int local_comm::GetWriteRecvCount(const int source, const int curr_offset){

	assert(write_recv_count[source] == 0);
	write_recv_offset[source] = curr_offset;
	int recv_count = 0;
	for (vector<local_data*>::iterator it = write_arrays.begin();
	     it != write_arrays.end(); it++)
		recv_count += (*it)->GetOwnedCount(source);
	write_recv_count[source] = recv_count;

	if (source == proc_id)
		return 0;
	else
		return recv_count;
}


/**
 * \brief Populate send buffer with data owned by this processor
 *
 * This data will be communicated to other processes
 * that use this data.
 *
 * \param send_buffer The buffer to be populated.
 * \param buffer_size
 */
void local_comm::PopulateReadSendBuffer(char* send_buffer){

	for( vector<local_data*>::iterator it = read_arrays.begin();
	     it != read_arrays.end() ; it++)

		(*it)->SendOwnedData(send_buffer, read_send_offset);
}


/**
 * \brief Move ghosts received in a receive buffer to the local copy
 *
 * This is data that we need to perform our computations, but that
 * we don't own. It has been previously communicated by their owner
 * processes to us.
 *
 * \param recv_buffer The receive buffer to be copied to our local data
 */
void local_comm::ExtractReadRecvBuffer(char* recv_buffer){

	for( vector<local_data*>::iterator it = read_arrays.begin();
	     it  != read_arrays.end() ; it++)

		(*it)->RecvGhostData(recv_buffer, read_recv_offset);
}


/**
 * \brief Write ghosts from all local data to a send buffer.
 *
 * This data must be communicated to their owners so that they
 * can update the data and broadcast it to other user processors.
 *
 * \param send_buffer The send buffer to be populated
 */
void local_comm::PopulateWriteSendBuffer(char* send_buffer){

	for (vector<local_data*>::iterator it = write_arrays.begin();
	     it != write_arrays.end(); it++)

		(*it)->SendGhostData(send_buffer, write_send_offset);
}


/**
 * \brief Get received data owned by us in order to update local copy
 *
 * Since we are the owners of all that data, our copy is the important
 * one!
 *
 * \param send_buffer The receive buffer where preliminary results
 *        are located
 */
void local_comm::ExtractWriteRecvBuffer(char* recv_buffer){

	for( vector<local_data*>::iterator it = write_arrays.begin();
	     it  != write_arrays.end() ; it++)

		(*it)->RecvOwnedData(recv_buffer, write_recv_offset);
}


void local_comm::InitWriteGhosts(){

	for (vector<local_data*>::iterator it = write_arrays.begin();
	     it != write_arrays.end(); it++)

		(*it)->InitWriteGhosts();
}
