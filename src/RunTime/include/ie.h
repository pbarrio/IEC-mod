/*
 * ie.h: This file is part of the IEC project.
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
 * @file: ie.h
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __IE_H__
#define __IE_H__

#include "mpi.h"
#include "armci.h"
#include "omp.h"
#include "string.h"
#include "stdlib.h"

static int __nprocs__ = 1;
/* static int __nthreads__ = 1; */
static int __nprocs_y__ = 1;
static int __myid__ = 0;
static int __proc_id__ = 0;
/* static int __thread_id__ = 0; */
static int __proc_coordy__ = 1;

enum partition_type{
  PARTITION_PATOH,
  PARTITION_METIS,
  PARTITION_BLOCK
};

extern void create_inspector(int md, int np, int team, int pid_team,
                             int teamsize, int nloops, int ndata, int nc,
                             int nac, int *nic, int* ndc, int* ro);
extern void set_access_array_param(int, int, int, int*);
extern int done_graph_gen();
extern int is_known(int, int);
extern int get_elem(int, int);
extern void set_array_stride(int, int);
extern void partition_hypergraph(enum partition_type);
extern int add_vertex(int, int);
extern void add_pin_to_net(int, int, int, int, int);
extern void add_index_from_proc(int, int, int, int);
extern int get_vertex_home(int, int);
extern int get_local_size(int);
extern int populate_local_array(int, double*, double*, int);
extern int renumber_access_array(int, int, int*);
extern int renumber_offset_array(int, int, int* ,int*);
extern int renumber_const_offset_array(int, int, int, int* ,int);
extern int get_proc_iter_size(int);
extern void setup_executor();
extern void communicate_reads_for(int, int);
extern void communicate_writes_for(int, int);
extern void communicate_reads(int);
extern void communicate_reads_start(int, int);
extern void communicate_reads_end(int, int);
extern void communicate_writes(int);
extern void init_write_ghosts(int);
extern void reduce_scalar(double*);
extern void print_hypergraph();
extern void print_data(int);
extern void print_access();
extern void populate_global_arrays();
extern double rtclock();
extern void delete_inspector();
extern double** malloc_2d_double(int, int);
extern void free_2d_double(double**);
extern float** malloc_2d_float(int, int);
extern void free_2d_float(float**);


/*
 * NEW FUNCTIONS FOR PIPELINING
 */
extern void pipe_endExternalIter();


#endif
