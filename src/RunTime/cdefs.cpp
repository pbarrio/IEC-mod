/*
 * cdefs.cpp: This file is part of the IEC project.
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
 * @file: cdefs.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/inspector.hpp"
#include "sys/time.h"
#include "mpi.h"
#include "armci.h"

using namespace std;

static double* scalar_holder = NULL;
static double final_val = 0.0;
static ret_data_access data_access_info;
static double inspector_start = 0.0;
static double inspector_stop = 0.0;
static double executor_start = 0.0;
static double executor_stop = 0.0;


extern "C" {

	void armci_init_(){
		ARMCI_Init();
	}


	void armci_finalize_(){

		ARMCI_Finalize();
	}


	double rtclock() { 
		struct timezone Tzp; 
		struct timeval Tp; 
		int stat; 
		stat = gettimeofday (&Tp, &Tzp); 
		if (stat != 0) printf("Error return from gettimeofday: %d",stat); 
		return(Tp.tv_sec + Tp.tv_usec*1.0e-6); 
	}


	void create_inspector(int md, int np, int team, int pid_team, int teamsize,
	                      int nloops, int ndata, int nc, int nac, int* nic,
	                      int* ndc, int* ro){

		inspector_start = rtclock();
		inspector* new_inspector =
			inspector::instance(md, np, team, pid_team, teamsize, nloops, ndata,
			                    nc, nac, nic, ndc, ro);
		scalar_holder = new double[1];
#ifndef NDEBUG
		fprintf(stdout, "PID:%d,Creating inspector:%p\n", md, new_inspector);
		fflush(stdout);
#endif
	}


	void create_inspector_(int* md, int* np, int team, int pid_team,
	                       int teamsize, int* nloops, int* ndata, int* nc,
	                       int *nac, int *nic, int* ndc, int* ro){

		create_inspector(*md, *np, team, pid_team, teamsize, *nloops, *ndata,
		                 *nc, *nac, nic, ndc, ro);
	}


	void set_access_array_param(int an,int as,int st,int* oa){
		inspector::instance()->SetAccessArrayParam(an,as,st,oa);
	}

	void set_access_array_param_(int* an, int *as, int* st, int *oa){
		set_access_array_param(*an,*as,*st,oa);
	}
  

	/**
	 * \brief Ask if we have the value of a position in an indirection array.
	 *
	 * \param an ID of the indirection array
	 * \param indx Position of the array that we are looking for.
	 *
	 * \return 0 if value not found, 1 if found
	 */
	int is_known(int an, int indx){
		if( inspector::instance()->HaveIndex(an,indx) )
			return 1;
		else
			return 0;
	}

	void is_known_(int* an, int* indx, int* val){
		*val = is_known(*an,*indx);
	}


	/**
	 * \brief Get value of a position in an indirection array.
	 *
	 * \param an ID of the indirection array
	 * \param index Position in the array
	 *
	 * \return Value
	 */
	int get_elem(int an, int index){
		return inspector::instance()->GetIndex(an,index);
	}

	void get_elem_(int* an, int* index, int* val){
		*val = get_elem(*an,*index);
	}

	int done_graph_gen(){
		if( inspector::instance()->DoneGraphGeneration() )
			return 1;
		else
			return 0;
	}

	void done_graph_gen_(int* val){
		*val = done_graph_gen();
	}

	/**
	 * \brief Add a vertex corresponding to a loop iteration.
	 *
	 * \param iter_num Loop identifier
	 * \param iter_value Iteration number
	 */
	void add_vertex(int iter_num, int iter_value) {
		inspector* my_inspect = inspector::instance();
		my_inspect->AddVertex(iter_num,iter_value);
	}

  
	void add_vertex_(int* iter_num, int* iter_value){
		inspector* my_inspect = inspector::instance();
		my_inspect->AddVertex(*iter_num,*iter_value);
	}

	/**
	 * \brief Connect a vertex to a net.
	 *
	 * Vertices are loop iterations, nets represent data.
	 * When we detect that some data is used in an iteration,
	 * we must connect the iteration vertex to the data net.
	 * This allows us to know which processes will use which data,
	 * since each vertex will be mapped to a single process.
	 *
	 * \param data_num Identifier of the data array
	 * \param index Position in the array
	 * \param isdirect !=0 if addressing is affine; =0 if it depends on indirection array.
	 * \param isploop !=0 if the access comes from a partitionable loop. What is a p. loop?
	 */
	void add_pin_to_net(int data_num, int index, int isdirect, int isploop) {
		inspector* my_inspect = inspector::instance();
		my_inspect->AddNet(data_num,index,isdirect,isploop);
	}


	void add_pin_to_net_(int* data_num, int* index, int* id, int* isploop){
		inspector* my_inspect = inspector::instance();
		my_inspect->AddNet(*data_num,*index,*id,*isploop);
	}

	void partition_hypergraph(int use_type)  {
		switch(use_type){
		case 0:
			inspector::instance()->PatohPartition();
			return;
		case 1:
			inspector::instance()->MetisPartition();
			return;
		case 2:
			inspector::instance()->BlockPartition();
			return;
		default:
			assert(0);
		}
	}

	void partition_hypergraph_(int *type){
		partition_hypergraph(*type);
	}

  
	void print_hypergraph() {  
#ifndef NDEBUG
		printf("ID:%d,Here:\n",inspector::instance()->get_proc_id());
		fflush(stdout);
		char graph_file[18];
		sprintf(graph_file,"hypergraph_%d.dat",inspector::instance()->get_proc_id());
		FILE* outfile = fopen(graph_file,"w");
		inspector::instance()->print_hypergraph(outfile);
		fclose(outfile);
#endif
	}

	void print_hypergraph_(){
		print_hypergraph();
	}
  
	int get_proc_iter_size(int in){
		return inspector::instance()->GetProcLocal(in);
	}

	void get_proc_iter_size_(int *in, int *out){
		*out = inspector::instance()->GetProcLocal(*in);
	}

	void set_array_stride(int an, int st) {
		inspector::instance()->SetStride(an,st);
	}


	void set_array_stride_(int* an, int* st) {
		inspector::instance()->SetStride(*an,*st);
	}

	int get_vertex_home(int in ,int iv){
		return inspector::instance()->GetVertexHome(in,iv);
	}


	void get_vertex_home_(int* in ,int* iv, int* out){
		*out = inspector::instance()->GetVertexHome(*in,*iv);
	}

	int get_local_size(int an){
			data_access_info = inspector::instance()->GetLocalAccesses(an);

		int * proc_recv_buffer = data_access_info.recvbuffer;
		int * proc_recv_displ = data_access_info.recvdispl;
		int * proc_recv_count = data_access_info.recvcount;
		local_inspector* curr_local_inspector = local_inspector::instance(/*thread_id*/);
		int nprocs = inspector::instance()->GetNProcs();
		for( int j = 0 ; j < nprocs ; j++ ){
			curr_local_inspector->InsertDirectAccess(an,proc_recv_buffer+proc_recv_displ[j*2],proc_recv_count[j*2]);
			curr_local_inspector->InsertIndirectAccess(an,proc_recv_buffer+proc_recv_displ[j*2+1],proc_recv_count[j*2+1]);
		}
			delete[] data_access_info.recvbuffer;
			delete[] data_access_info.recvcount;
			delete[] data_access_info.recvdispl;
		curr_local_inspector->SetupLocalArray(an);
		return curr_local_inspector->GetLocalDataSize(an);
	}


	void get_local_size_(int* an, int* out){
		*out = get_local_size(*an);
	}

	/**
	 * \brief Unused in quake
	 */
	void communicate_reads_for(int in, int an){
		local_inspector::instance()->AddReadArray(in,an);
	}

	/**
	 * \brief Unused in quake
	 */
	void communicate_reads_for_(int *in, int *an){
		local_inspector::instance()->AddReadArray(*in,*an);
	}

	void communicate_writes_for(int in, int an){
		local_inspector::instance()->AddWriteArray(in,an);
	}

	void communicate_writes_for_(int *in, int *an){
		local_inspector::instance()->AddWriteArray(*in,*an);
	}

	void add_index_from_proc(int data_num, int index, int access_type){
		local_inspector::instance()->AddIndexAccessed(data_num,index,access_type);
	}

	void add_index_from_proc_(int *data_num, int *index, int *access_type){
		local_inspector::instance()->AddIndexAccessed(*data_num,*index,*access_type);
	}


	void populate_local_array(int an, double*lb, double*oa, int st){
		local_inspector::instance()->PopulateLocalArray(an,lb,oa,st);
	}

	void populate_local_array_(int *an, double*lb, double*oa, int *st){
		local_inspector::instance()->PopulateLocalArray(*an,lb,oa,*st);
	}


	void renumber_access_array(int an, int as, int* aa){
		local_inspector::instance()->RenumberAccessArray(an,as,aa);
	}

	void renumber_access_array_(int *an, int *as, int* aa){
		local_inspector::instance()->RenumberAccessArray(*an,*as,aa);
	}

	void renumber_offset_array(int an, int as, int* aa, int* la){
		local_inspector::instance()->RenumberOffsetArray(an,as,aa,la);
	}

	void renumber_offset_array_(int *an, int *as, int* aa, int* la){
		local_inspector::instance()->RenumberOffsetArray(*an,*as,aa,la);
	}

	/**
	 * \brief Start executor
	 *
	 * This function must run immediately before the start of the computations.
	 */
	void setup_executor(){
		local_inspector::instance()->GenerateGhosts();
    
#ifndef NDEBUG
		printf("LocalGhostsDone\n");
		fflush(stdout);
#endif

		inspector::instance()->CommunicateGhosts();
#ifndef NDEBUG
		printf("Global ghosts done\n");
		fflush(stdout);
#endif

		local_inspector::instance()->GenerateOwned();

		inspector::instance()->GetBufferSize();
		inspector_stop = rtclock();
		if( inspector::instance()->get_proc_id() == 0 )
			fprintf(stderr,"[IEC]:InspectorTime:%lf\n",inspector_stop-inspector_start);

		executor_start = rtclock();
	}

  
	void setup_executor_(){
		setup_executor();
	}

	void communicate_reads(int cn){
		inspector::instance()->CommunicateReads(cn);
	}

	void communicate_reads_(int *cn){
		inspector::instance()->CommunicateReads(*cn);
	}

	void communicate_to_next(){
		inspector::instance()->CommunicateToNext();
	}

	// void communicate_writes(int cn){
	// 	inspector::instance()->CommunicateWrites(cn);
	// }

	// void communicate_writes_(int *cn){
	// 	inspector::instance()->CommunicateWrites(*cn);
	// }

	void init_write_ghosts(int cn){
		local_inspector::instance()->InitWriteGhosts(cn);
	}

	void init_write_ghosts_(int *cn){
		local_inspector::instance()->InitWriteGhosts(*cn);
	}
	void reduce_scalar(double *val){
		scalar_holder[0] = *val;
		final_val = scalar_holder[0];      
		MPI_Allreduce(MPI_IN_PLACE,&final_val,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		*val = final_val;
	}

	void reduce_scalar_(double* val){
		reduce_scalar(val);
	}

	void populate_global_arrays(){
			executor_stop = rtclock();
			if( inspector::instance()->get_proc_id() == 0 )
				fprintf(stderr,"[IEC]:ExecutorTime:%lf\n",executor_stop-executor_start);
 
			local_inspector::instance()->PopulateGlobalArrays();
	}

	void populate_global_arrays_(){
		local_inspector::instance()->PopulateGlobalArrays();
	}


	void delete_inspector(){

		local_inspector* li = local_inspector::instance();
		if (li)
			delete li;

		{
			delete inspector::instance();
			delete[] scalar_holder;
		}

	}

	void delete_inspector_(){
		delete_inspector();
	}

	void print_data(){
#ifndef NDEBUG
		local_inspector::instance()->print_data();
			inspector::instance()->print_comm();
#endif
	}

	void print_data_(){
		print_data();
	}


	void print_access(){
#ifndef NDEBUG
		printf("PID:%d,PrinitingAccess\n",inspector::instance()->get_proc_id());
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		inspector::instance()->print_access();
#endif
	}

	void print_access_(){
		print_access();
	}

  
	//dim1 = number of rows, dim2 number of columns. Row major (colums is fastest varying)
	double** malloc_2d_double(int dim1, int dim2){
		double* base_array = (double*)malloc(dim1*dim2*sizeof(double));
		double** actual_array = (double**)malloc(dim1*sizeof(double*));
		for( int i = 0 ; i < dim1 ; i++ ){
			actual_array[i] = base_array + i*dim2;
		}
		return actual_array;
	}

	float** malloc_2d_float(int dim1, int dim2){
		float* base_array = (float*)malloc(dim1*dim2*sizeof(float));
		float** actual_array = (float**)malloc(dim1*sizeof(float*));
		for( int i = 0 ; i < dim1 ; i++ ){
			actual_array[i] = base_array + i*dim2;
		}
		return actual_array;
	}

	void free_2d_double(double** actual_array){
		free(actual_array[0]);
		free(actual_array);
	}

	void free_2d_float(float** actual_array){
		free(actual_array[0]);
		free(actual_array);
	}


	int  init_solver(int s){
		return inspector::instance()->InitSolver(s);
	}

	void init_solver_(int *s, int *sn){
		*sn = inspector::instance()->InitSolver(*s);
	}
  
	void add_unknown(int sn, int an, int idx, int rn){
		assert(sn == 0);
		inspector::instance()->AddUnknown(sn,an,idx,rn);
	}

	void add_unknown_(int *sn, int *an, int *idx, int *rn){
		assert(*sn == 0);
		inspector::instance()->AddUnknown(*sn,*an,*idx,*rn);
	}

	void renumber_solver_rows(int sn, int* oa, int as){
		assert(sn == 0 );
		inspector::instance()->RenumberGlobalRows(sn,oa,as);
	}

	void renumber_solver_rows_(int *sn, int* oa, int *as){
		assert(*sn == 0 );
		inspector::instance()->RenumberGlobalRows(*sn,oa,*as);
	}

	void print_solver(){
		inspector::instance()->print_solver();
	}

	void print_solver_(){
		inspector::instance()->print_solver();
	}

	int get_local_rows(int sn/*, int tid*/){
		return inspector::instance()->GetLocalRows(sn/*,tid*/);
	}

	void get_local_rows_(int *sn, /*int *tid,*/ int *lr){
		*lr =  inspector::instance()->GetLocalRows(*sn/*,*tid*/);
	}


	/**
	 * \brief Unused in quake
	 */
	void set_constraint(int an){
		inspector::instance()->SetConstraint(an);
	}


	/**
	 * \brief Unused in quake
	 */
	void set_constraint_(int *an){
		inspector::instance()->SetConstraint(*an);
	}

	void set_local_array(/*int tid,*/ int an, void* la){
		local_inspector::instance(/*tid*/)->SetLocalArray(an,la);
	}

	void set_local_array_(/*int* tid,*/ int* an, void* la){
		local_inspector::instance(/**tid*/)->SetLocalArray(*an,la);
	}

	int get_local_block_offset(/*int tid,*/ int an){
		return local_inspector::instance(/*tid*/)->GetLocalBlockOffset(an);
	}

	void get_local_block_offset_(/*int* tid,*/ int* an, int *bo){
		*bo = local_inspector::instance(/**tid*/)->GetLocalBlockOffset(*an);
	}

	void fflush_(){
		fflush(stdout);
		fflush(stderr);
	}


	/*
	 * NEW FUNCTIONS FOR PIPELINING
	 */

	void pipe_endExternalIter(){

		inspector::instance()->pipe_endExternalIter();
	}

}
