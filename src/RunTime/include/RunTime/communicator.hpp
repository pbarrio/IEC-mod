/*
 * communicator.hpp: This file is part of the IEC project.
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
 * @file: communicator.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __COMMUNICATOR_HPP__
#define __COMMUNICATOR_HPP__

#include <cstdio>
#include "mpi.h"

/* The pattern for the buffers. Case 16 process, 4 nodes, example on node 1
   SEnd :
   | 4->0 | 5->0 | 6->0 | 7->0 || 4->1 | 5->1 | 6->1 | 7->1 || 4->2 | 5->2
   | 6->2 | 7->2 || 4->3 | 5->3 | 6->3 | 7->3 ||| 4->4 | 5->4 | 6->4 | 7->4
   || 4->5 | 5->5 | 6->5 | 7->5 || 4->6 | 5->6 | 6->6 | 7->6 || 4->7 | 5->7
   | 6->7 | 7->7 ||| 4->8 | 5->8 | 6->8 | 7->8 || 4->9 | 5->9 | 6->9 | 7->9
   || 4->10 | 5->10 | 6->10 | 7->10 || 4->11 | 5->11 | 6->11 | 7->11 ||| 4->12
   | 5->12 | 6->12 | 7->12 || 4->13 | 5->13 | 6->13 | 7->13 || 4->14 | 5->14
   | 6->14 | 7->14 || 4->15 | 5->15 | 6->15 | 7->15 |||
   Recv:
   | 0->4 | 1->4 | 2->4 | 3->4 || 0->5 | 1->5 | 2->5 | 3->5 || 0->6 | 1->6
   | 2->6 | 3->6 || 0->7 | 1->7 | 2->7 | 3->7 ||| 4->4 | 5->4 | 6->4 | 7->4
   || 4->5 | 5->5 | 6->5 | 7->5 || 4->6 | 5->6 | 6->6 | 7->6 || 4->7 | 5->7
   | 6->7 | 7->7 ||| 8->4 | 9->4 | 10->4 | 11->4 || 8->5 | 9->5 | 10->5 | 11->5
   || 8->6 | 9->6 | 10->6 | 11->6 || 8->7 | 9->7 | 10->7 | 11->7 ||| 12->4
   | 13->4 | 14->4 | 15->4 || 12->5 | 13->5 | 14->5 | 15->5 || 12->6 | 13->6
   | 14->6 | 15->6 || 12->7 | 13->7 | 14->7 | 15->7 |||
*/


class inspector;

/**
 * \brief Info to communicate data to other processes
 *
 * There could be multiple global communicators, one per loop. In the past,
 * different loops required to use different comm details. When we move to
 * one local comm per team, remove this (we can use only one global comm because
 * each process will only work in one loop).
 */
class global_comm{
private:

	int* const read_send_offset;

	int* const read_recv_offset;

	int* const write_send_offset;

	int* const write_recv_offset;

	int* const read_send_count;

	int* const read_recv_count;

	int* const write_send_count;

	int* const write_recv_count;

	int* const read_put_displ;

	int* const write_put_displ;

	static char* read_send_signal;

	static char* read_recv_signal;

	static char* write_send_signal;

	static char* write_recv_signal;

	int nprocs_read_send;

	int nprocs_read_recv;

	int nprocs_write_send;

	int nprocs_write_recv;

	int* proc_id_read_send;

	int* proc_id_read_recv;

	int* proc_id_write_send;

	int* proc_id_write_recv;

	static int max_nprocs_read_send;

	static int max_nprocs_read_recv;

	static int max_nprocs_write_send;

	static int max_nprocs_write_recv;

	static MPI_Request* read_send_end_request;

	static MPI_Status* read_recv_end_status;

	static MPI_Status* read_send_end_status;

	static MPI_Request* read_recv_start_request;

	static MPI_Status* read_send_start_status;

	static MPI_Status* read_recv_start_status;

	static MPI_Request* write_send_end_request;

	static MPI_Status* write_recv_end_status;

	static MPI_Status* write_send_end_status;

	static MPI_Request* write_recv_start_request;

	static MPI_Status* write_send_start_status;

	static MPI_Status* write_recv_start_status;

	const int my_num;

	const int nprocs;

	const int proc_id;

	static char** put_buffer;

	static char* send_buffer;

	static char* recv_buffer;

	static int max_send_size;

	static int max_recv_size;

#ifndef NDEBUG
	static FILE* comm_file;
#endif

	global_comm(int,int,int);

	~global_comm();

	void CommunicateReads();

	void CommunicateWrites();

	void InitWriteGhosts();

public:

	static MPI_Comm global_iec_communicator;
	static MPI_Comm team_communicator;

	void print_comm();

	friend class Inspector;
};


#endif
