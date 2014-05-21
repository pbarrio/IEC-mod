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
  std::vector<local_data*> read_arrays;
  
  std::vector<local_data*> write_arrays;
  
  const int thread_id;
  
  const int proc_id;

  const int nthreads;
  
  const int nprocs;

  const int myid;

  const int my_num;
  
  int* const read_send_offset;
  
  int* const read_recv_offset;
  
  int* const write_send_offset;
  
  int* const write_recv_offset;

  int* const read_send_count;
  
  int* const read_recv_count;
  
  int* const write_send_count;
  
  int* const write_recv_count;

  int* const offset_array;
  
 public:
  
  local_comm(int,int,int,int,int);
  
  ~local_comm();

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
