/*
 * communicator.cpp: This file is part of the IEC project.
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
 * @file: communicator.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/communicator.hpp"
#include <cassert>

#ifndef NDEBUG
FILE* global_comm::comm_file = NULL;
#endif

global_comm::global_comm(int mn,int np, int md):
  nprocs(np),
  proc_id(md),
  my_num(mn),
  read_send_count(new int[np]),
  read_send_offset(new int[np+1]),
  read_recv_count(new int[np]),
  read_recv_offset(new int[np+1]),
  read_put_displ(new int[np]),
  nprocs_read_send(0),
  nprocs_read_recv(0),
  proc_id_read_send(NULL),
  proc_id_read_recv(NULL),
  write_send_count(new int[np]),
  write_send_offset(new int[np+1]),
  write_recv_count(new int[np]),
  write_recv_offset(new int[np+1]),
  write_put_displ(new int[np]),
  nprocs_write_send(0),
  nprocs_write_recv(0),
  proc_id_write_send(NULL),
  proc_id_write_recv(NULL)
{
  for( int i = 0; i < np ; i++ ){
    read_send_count[i] = 0;
    read_send_offset[i] = 0;
    read_recv_count[i] = 0;
    read_recv_offset[i] = 0;
    read_put_displ[i] = 0;
    write_send_count[i] = 0;
    write_send_offset[i] = 0;
    write_recv_count[i] = 0;
    write_recv_offset[i] = 0;
    write_put_displ[i] = 0;
  }
  read_send_offset[np] = 0;
  read_recv_offset[np] = 0;
  write_send_offset[np] = 0;
  write_recv_offset[np] = 0;

#ifndef NDEBUG
  if( comm_file == NULL ){
    char cf_name[20];
    sprintf(cf_name,"global_comm_%d.dat",proc_id);
    comm_file = fopen(cf_name,"w");
  }
#endif
}

global_comm::~global_comm()
{
  delete[] read_send_count;
  delete[] read_send_offset;
  delete[] read_recv_count;
  delete[] read_recv_offset;
  delete[] read_put_displ;
  delete[] write_send_count;
  delete[] write_send_offset;
  delete[] write_recv_count;
  delete[] write_recv_offset;
  delete[] write_put_displ;
#ifndef NDEBUG
  if( comm_file){
    fclose(comm_file);
    comm_file = NULL;
  }
#endif
}






void global_comm::print_comm()
{
#ifndef NDEBUG
  assert(comm_file);
  fprintf(comm_file,"\tRead Send:");
  for( int i = 0 ; i < nprocs ; i++ )
    fprintf(comm_file," %d(%d)",read_send_count[i],read_send_offset[i]);
  fprintf(comm_file," (%d)\n",read_send_offset[nprocs]);
  fprintf(comm_file,"\tRead Recv:");
  for( int i = 0 ; i < nprocs ; i++ )
    fprintf(comm_file," %d(%d)",read_recv_count[i],read_recv_offset[i]);
  fprintf(comm_file," (%d)\n",read_recv_offset[nprocs]);
  fprintf(comm_file,"Read Put Displ:");
  for( int i = 0 ; i < nprocs ; i++ )
    fprintf(comm_file," %d",read_put_displ[i]);
  fprintf(comm_file,"\n");
  fprintf(comm_file,"\tWrite Send:");
  for( int i = 0 ; i < nprocs ; i++ )
    fprintf(comm_file," %d(%d)",write_send_count[i],write_send_offset[i]);
  fprintf(comm_file," (%d)\n",write_send_offset[nprocs]);
  fprintf(comm_file,"\tWrite Recv:");
  for( int i = 0 ; i < nprocs ; i++ )
    fprintf(comm_file," %d(%d)",write_recv_count[i],write_recv_offset[i]);
  fprintf(comm_file," (%d)\n",write_recv_offset[nprocs]);
  fprintf(comm_file,"Write Put Displ:");
  for( int i = 0 ; i < nprocs ; i++ )
    fprintf(comm_file," %d",write_put_displ[i]);
  fprintf(comm_file,"\n");
#endif

}
