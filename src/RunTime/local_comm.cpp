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

local_comm::local_comm(int mn, int np/*, int nt*/, int pid/*, int tid*/):
	my_num(mn),
	nprocs(np),
	//   nthreads(nt),
	proc_id(pid),
	//   thread_id(tid),
	myid(pid/**nt+tid*/),
	read_send_offset(new int[np/**nt*/]),
	read_recv_offset(new int[np/**nt*/]),
	write_send_offset(new int[np/**nt*/]),
	write_recv_offset(new int[np/**nt*/]),
	read_send_count(new int[np/**nt*/]),
	read_recv_count(new int[np/**nt*/]),
	write_send_count(new int[np/**nt*/]),
	write_recv_count(new int[np/**nt*/]),
	offset_array(new int[np/**nt*/])
{ 
	memset(read_send_count,0,sizeof(int)*nprocs/**nthreads*/);
	memset(read_send_offset,0,sizeof(int)*nprocs/**nthreads*/);
	memset(read_recv_count,0,sizeof(int)*nprocs/**nthreads*/);
	memset(read_recv_offset,0,sizeof(int)*nprocs/**nthreads*/);
	memset(write_send_count,0,sizeof(int)*nprocs/**nthreads*/);
	memset(write_send_offset,0,sizeof(int)*nprocs/**nthreads*/);
	memset(write_recv_count,0,sizeof(int)*nprocs/**nthreads*/);
	memset(write_recv_offset,0,sizeof(int)*nprocs/**nthreads*/);
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

local_comm::~local_comm()
{
#ifdef COMM_TIME
	if( thread_id == 0 && read_comm_time > 0.0 )
		printf("MXC:GID:%d,CommNum:%d,ReadCommTime:%5.4le\n",myid,my_num,read_comm_time);
	if( thread_id == 0 && write_comm_time > 0.0 )
		printf("MXC:GID:%d,CommNum:%d,WriteCommTime:%5.4le\n",myid,my_num,read_comm_time);
#endif
	delete[] read_send_offset;
	delete[] read_recv_offset;
	delete[] write_send_offset;
	delete[] write_recv_offset;
	delete[] read_send_count;
	delete[] read_recv_count;
	delete[] write_send_count;
	delete[] write_recv_count;
	delete[] offset_array;
}


int local_comm::GetReadSendCount(const int dest, const int curr_offset)
{
	assert(read_send_count[dest] == 0);
	read_send_offset[dest] = curr_offset;
	int send_count = 0;
	for( vector<local_data*>::iterator it = read_arrays.begin() ; it != read_arrays.end() ; it++ )
		send_count += (*it)->GetOwnedCount(dest);
	read_send_count[dest] = send_count;
  
	if( dest/* / nthreads*/ == proc_id )
		return 0;
	else
		return send_count;
}

int local_comm::GetReadRecvCount(const int source, const int curr_offset )
{
	assert(read_recv_count[source] == 0);
	read_recv_offset[source] = curr_offset;
	int recv_count = 0;
	for( vector<local_data*>::iterator it = read_arrays.begin() ; it != read_arrays.end() ; it++ )
		recv_count += (*it)->GetGhostsCount(source);
	read_recv_count[source] = recv_count;

	if( source /*/ nthreads*/ == proc_id )
		return 0;
	else
		return recv_count;
}

int local_comm::GetWriteSendCount(const int dest, const int curr_offset)
{
	assert(write_send_count[dest] == 0);
	write_send_offset[dest] = curr_offset;
	int send_count = 0;
	for( vector<local_data*>::iterator it = write_arrays.begin() ; it != write_arrays.end() ; it++){
		send_count += (*it)->GetGhostsCount(dest);
	}
	write_send_count[dest] = send_count;

	if( dest /*/ nthreads*/ == proc_id )
		return 0;
	else
		return send_count;
}


int local_comm::GetWriteRecvCount(const int source, const int curr_offset)
{
	assert(write_recv_count[source] == 0);
	write_recv_offset[source] = curr_offset;
	int recv_count = 0;
	for( vector<local_data*>::iterator it = write_arrays.begin(); it != write_arrays.end() ; it++ )
		recv_count += (*it)->GetOwnedCount(source);
	write_recv_count[source] = recv_count;

	if( source /*/ nthreads*/ == proc_id )
		return 0;
	else
		return recv_count;
}


void local_comm::PopulateReadSendBuffer(char* send_buffer, int buffer_size)
{
	for( int i = 0 ; i < nprocs/**nthreads*/ ; i++ )
		offset_array[i]  = read_send_offset[i];

	for( vector<local_data*>::iterator it = read_arrays.begin() ; it  != read_arrays.end() ; it++){
		(*it)->SendOwnedData(send_buffer,offset_array);
	}
#ifndef NDEBUG
	for( int i = 0 ; i < nprocs/* * nthreads*/ ; i++ )
		assert(read_send_offset[i] <= buffer_size );
#endif
}


void local_comm::ExtractReadRecvBuffer(char* recv_buffer, int buffer_size)
{
	for( int i = 0 ; i < nprocs/**nthreads*/ ; i++ )
		offset_array[i]  = read_recv_offset[i];

	for( vector<local_data*>::iterator it = read_arrays.begin() ; it  != read_arrays.end() ; it++){
		(*it)->RecvGhostData(recv_buffer,offset_array);
	}
#ifndef NDEBUG
	for( int i = 0 ; i < nprocs/* * nthreads*/ ; i++ ){
		assert(read_recv_offset[i] <= buffer_size );
	}
#endif

}


/**
 * \brief Populates a send buffer with all the data that must be sent to other processes
 *
 * \param send_buffer
 */
void local_comm::PopulateWriteSendBuffer(char* send_buffer)
{
	// offset_array must be filled in with the positions of the send buffer where each
	// ghost must reside
	for (int i = 0 ; i < nprocs/**nthreads*/ ; i++)
		offset_array[i]  = write_send_offset[i];

	for (vector<local_data*>::iterator it = write_arrays.begin() ; it  != write_arrays.end() ; it++){
		(*it)->SendGhostData(send_buffer,offset_array);
	}
}


void local_comm::ExtractWriteRecvBuffer(char* recv_buffer)
{
	for( int i = 0 ; i < nprocs/**nthreads*/ ; i++ )
		offset_array[i]  = write_recv_offset[i];

	for( vector<local_data*>::iterator it = write_arrays.begin() ; it  != write_arrays.end() ; it++){
		(*it)->RecvOwnedData(recv_buffer,offset_array);
	}
}


void local_comm::InitWriteGhosts()
{
	for( vector<local_data*>::iterator it = write_arrays.begin() ; it  != write_arrays.end() ; it++){
		(*it)->InitWriteGhosts();
	}
}

void local_comm::print_comm(FILE* comm_file)
{
#ifndef NDEBUG
	assert(comm_file);
	fprintf(comm_file,"Local Communicator %d\n",my_num);
	if( read_arrays.size() != 0 ) {
		fprintf(comm_file,"Read Arrays:");
		for( vector<local_data*>::iterator it = read_arrays.begin() ; it != read_arrays.end() ; it++ )
			fprintf(comm_file," %d",(*it)->GetMyNum());
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tReadSendCount :");
		for(int i = 0 ; i < nprocs/* * nthreads*/ ; i++ )
			fprintf(comm_file," %d",read_send_count[i]);
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tReadSendOffset :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",read_send_offset[i]);
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tReadRecvCount :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",read_recv_count[i]);
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tReadRecvOffset :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",read_recv_offset[i]);
		fprintf(comm_file,"\n");
	}
	if( write_arrays.size() != 0 ){
		fprintf(comm_file,"Write Arrays:");
		for( vector<local_data*>::iterator it = write_arrays.begin() ; it != write_arrays.end() ; it++ )
			fprintf(comm_file," %d",(*it)->GetMyNum());
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tWriteSendCount :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",write_send_count[i]);
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tWriteSendOffset :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",write_send_offset[i]);
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tWriteRecvCount :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",write_recv_count[i]);
		fprintf(comm_file,"\n");
		fprintf(comm_file,"\tWriteRecvOffset :");
		for(int i = 0 ; i < nprocs/* * nthreads */; i++ )
			fprintf(comm_file," %d",write_recv_offset[i]);
		fprintf(comm_file,"\n");
	}
	fflush(comm_file);
#endif
}
