/*
 * local_comm.hpp: This file is part of the IEC project.
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
 * @file: local_comm.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __LOCAL_COMM_HPP__
#define __LOCAL_COMM_HPP__

#include <cstdio>
#include <vector>

class inspector;
class local_inspector;
class local_data;

class local_comm{
  
private:

	/// Read-only arrays in the target code
	std::vector<local_data*> read_arrays;
  
	/// Writable arrays in the target code
	std::vector<local_data*> write_arrays;

	const int proc_id;

	const int nprocs;

	const int myid;

	const int my_num;
  
	/*
	 * read =  data that we use.
	 * write = data that we modify after computations.
	 * ghost = data that we do not own
	 * owned = data that we own
	 */

	/// Start index in the send buffer, for each client process,
	/// of the OWNED data sent FROM THIS process.
	int* const read_send_offset;
  
	/// Start index in the recv buffer, for each owner process,
	/// of the GHOST data sent TO THIS process.
	int* const read_recv_offset;
  
	/// Start index in the send buffer, for each owner process,
	/// of the GHOST to be sent FROM THIS process.
	int* const write_send_offset;
  
	/// Start index in the recv buffer, for each client process,
	/// of the updated OWNED data sent TO THIS process.
	int* const write_recv_offset;

	int* const read_send_count;
  
	int* const read_recv_count;
  
	int* const write_send_count;
  
	int* const write_recv_count;
  
public:
  
	local_comm(int,int,int);
  
	~local_comm();

	///Unused in quake
	int GetReadSendCount(const int, const int);

	int GetReadRecvCount(const int, const int);

	int GetWriteSendCount(const int, const int);

	int GetWriteRecvCount(const int, const int);
  
	void PopulateReadSendBuffer(char*,int);
  
	void ExtractReadRecvBuffer(char*,int);

	void PopulateWriteSendBuffer(char*);
  
	void ExtractWriteRecvBuffer(char*);

	void InitWriteGhosts();

	void print_comm(FILE*);

#ifdef COMM_TIME
	double read_comm_time;
	double write_comm_time;
	double read_time1,write_time1,read_time2,write_time2,read_time3,write_time3;
#endif
  
	friend class local_inspector;

	friend class inspector;
};



#endif
