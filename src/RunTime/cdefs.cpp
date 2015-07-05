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
#include "RunTime/Inspector.hpp"
#include "sys/time.h"
#include "mpi.h"
#include "armci.h"

using namespace std;

static double* scalar_holder = NULL;
static double final_val = 0.0;
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


	/**
	 * \brief For more info, see the Inspector constructor
	 */
	void create_inspector(int md, int np, int team, int pid_team, int teamsize,
	                      int nloops, int ndata, int nc, int nac, int* nic,
	                      int* ndc, int* ro){

		inspector_start = rtclock();
		Inspector* new_inspector =
			Inspector::instance(md, np, team, pid_team, teamsize, nloops, ndata,
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
		Inspector::instance()->SetAccessArrayParam(an,as,st,oa);
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
		if( Inspector::instance()->HaveIndex(an,indx) )
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
		return Inspector::instance()->GetIndex(an,index);
	}

	void get_elem_(int* an, int* index, int* val){
		*val = get_elem(*an,*index);
	}

	int done_graph_gen(){
		if( Inspector::instance()->DoneGraphGeneration() )
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
	 * \param loop Loop identifier
	 * \param iter_value Iteration number
	 */
	void add_vertex(int loop, int iter_value) {
		Inspector* my_inspect = Inspector::instance();
		my_inspect->AddVertex(loop, iter_value);
	}

	void add_vertex_(int* iter_num, int* iter_value){
		Inspector* my_inspect = Inspector::instance();
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
	 * \param loop The pin is meaningful in the context of this loop.
	 * \param isdirect !=0 if the addressing is affine; =0 if depends on
	 *        indirection array.
	 * \param isploop !=0 if the access comes from a partitionable loop.
	 */
	void add_pin_to_net
	(int data_num, int index, int loop, int isdirect, int isploop){

		Inspector* my_inspect = Inspector::instance();
		my_inspect->AddPinToNet(data_num, index, loop, isdirect, isploop);
	}

	void add_pin_to_net_
	(int* data_num, int* index, int* loop, int* id, int* isploop){

		add_pin_to_net(*data_num, *index, *loop, *id, *isploop);
	}

	/**
	 * \brief Partition the hypergraphs for all loops
	 *
	 * \param useType Partition algorithm: 0 = Patoh, 1 = Metis, 2 = block
	 */
	void partition_hypergraph(int useType){
		switch(useType){
		case 0:
			Inspector::instance()->PatohPartitionAll();
			return;
		case 1:
			Inspector::instance()->MetisPartitionAll();
			return;
		case 2:
			Inspector::instance()->BlockPartition();
			return;
		default:
			assert(0);
		}
	}

	void partition_hypergraph_(int *type){
		partition_hypergraph(*type);
	}

	int get_proc_iter_size(int in){
		return Inspector::instance()->GetProcLocal(in);
	}

	void get_proc_iter_size_(int *in, int *out){
		*out = Inspector::instance()->GetProcLocal(*in);
	}

	void set_array_stride(int an, int st) {
		Inspector::instance()->SetStride(an, st);
	}


	void set_array_stride_(int* an, int* st) {
		Inspector::instance()->SetStride(*an, *st);
	}

	int get_vertex_home(int in ,int iv){
		return Inspector::instance()->GetVertexHome(in, iv);
	}

	void get_vertex_home_(int* in ,int* iv, int* out){
		*out = Inspector::instance()->GetVertexHome(*in, *iv);
	}


	/**
	 * \brief Get the local size of an array in this process
	 *
	 * This is the size of the array portion that is locally used
	 * by this process.
	 *
	 * \param an Array ID
	 */
	int get_local_size(int an){

		int *buf, *displ, *count;
		Inspector* inspector = Inspector::instance();
		inspector->GetLocalAccesses(an, &buf, &displ, &count);

		int nprocs = Inspector::instance()->GetNProcs();
		for (int j = 0; j < nprocs; j++){
			inspector->InsertDirectAccess(an, buf + displ[j * 2], count[j * 2]);
			inspector->InsertIndirectAccess(an, buf + displ[j * 2 + 1],
			                                count[j * 2 + 1]);
		}

		delete[] buf;
		delete[] count;
		delete[] displ;

		inspector->SetupLocalArray(an);
		return inspector->GetLocalDataSize(an);
	}


	void get_local_size_(int* an, int* out){
		*out = get_local_size(*an);
	}

	/**
	 * \brief Unused in quake
	 */
	void communicate_reads_for(int in, int an){
		Inspector::instance()->AddReadArray(in,an);
	}
	void communicate_reads_for_(int *in, int *an){
		Inspector::instance()->AddReadArray(*in,*an);
	}

	void communicate_writes_for(int in, int an){
		Inspector::instance()->AddWriteArray(in,an);
	}
	void communicate_writes_for_(int *in, int *an){
		Inspector::instance()->AddWriteArray(*in,*an);
	}

	void add_index_from_proc(int data_num, int index, int access_type){
		Inspector::instance()->
			AddIndexAccessed(data_num, index, access_type);
	}
	void add_index_from_proc_(int *data_num, int *index, int *access_type){
		Inspector::instance()->
			AddIndexAccessed(*data_num, *index, *access_type);
	}

	/**
	 * \brief Populates a local array from its corresponding global array
	 *
	 * \param an ID of the local array to be populated
	 * \param lb Allocated clean array to be populated
	 * \param oa Original array
	 * \param st Stride. e.g. A bidimensional array containing 100 coordinates
	 *           (x, y, z) would be a[100][3] and the stride would be 3. Knowing
	 *           the stride allows the program to use contiguous memory for
	 *           N-dimensional arrays. The stride comes from the IEC pragmas.
	 */
	void populate_local_array(int an, double* lb, double* oa, int st){
		Inspector::instance()->PopulateLocalArray(an, lb, oa, st);
	}
	void populate_local_array_(int* an, double* lb, double* oa, int* st){
		Inspector::instance()->PopulateLocalArray(*an, lb, oa, *st);
	}

	void renumber_access_array(int an, int as, int* aa){
		Inspector::instance()->RenumberAccessArray(an, as, aa);
	}
	void renumber_access_array_(int *an, int *as, int* aa){
		Inspector::instance()->RenumberAccessArray(*an, *as, aa);
	}

	void renumber_offset_array(int an, int as, int* aa, int* la){
		Inspector::instance()->RenumberOffsetArray(an, as, aa, la);
	}
	void renumber_offset_array_(int *an, int *as, int* aa, int* la){
		Inspector::instance()->RenumberOffsetArray(*an, *as, aa, la);
	}

	/**
	 * \brief Start executor
	 *
	 * This function must run immediately before the start of the computations.
	 */
	void setup_executor(){
		Inspector::instance()->GenerateGhosts();

#ifndef NDEBUG
		printf("LocalGhostsDone\n");
		fflush(stdout);
#endif

		Inspector::instance()->CommunicateGhosts();
#ifndef NDEBUG
		printf("Global ghosts done\n");
		fflush(stdout);
#endif

		Inspector::instance()->GenerateOwned();

		Inspector::instance()->GetBufferSize();
		inspector_stop = rtclock();
		if (Inspector::instance()->get_proc_id() == 0)
			fprintf(stderr,"[IEC]:InspectorTime:%lf\n",
			        inspector_stop-inspector_start);

		executor_start = rtclock();
	}

	void setup_executor_(){
		setup_executor();
	}

	void communicate_reads(int cn){
		Inspector::instance()->CommunicateReads(cn);
	}
	void communicate_reads_(int *cn){
		Inspector::instance()->CommunicateReads(*cn);
	}

	void communicate_to_next(){
		Inspector::instance()->CommunicateToNext();
	}

	void get_from_previous(){
		Inspector::instance()->GetFromPrevious();
	}

	void init_write_ghosts(int cn){
		Inspector::instance()->InitWriteGhosts(cn);
	}

	void init_write_ghosts_(int *cn){
		Inspector::instance()->InitWriteGhosts(*cn);
	}

	void reduce_scalar(double *val){
		scalar_holder[0] = *val;
		final_val = scalar_holder[0];
		MPI_Allreduce(MPI_IN_PLACE, &final_val, 1, MPI_DOUBLE, MPI_SUM,
		              MPI_COMM_WORLD);

		*val = final_val;
	}

	void reduce_scalar_(double* val){
		reduce_scalar(val);
	}

	void populate_global_arrays(){
		executor_stop = rtclock();
		if (Inspector::instance()->get_proc_id() == 0)
			fprintf(stderr, "[IEC]:ExecutorTime:%lf\n",
			        executor_stop-executor_start);

		Inspector::instance()->PopulateGlobalArrays();
	}

	void populate_global_arrays_(){
		Inspector::instance()->PopulateGlobalArrays();
	}


	void delete_inspector(){

		Inspector* li = Inspector::instance();
		if (li)
			delete li;

		delete[] scalar_holder;
	}

	void delete_inspector_(){
		delete_inspector();
	}

	// dim1 = number of rows, dim2 number of columns.
	// Row major (colums is fastest varying)
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
		return Inspector::instance()->InitSolver(s);
	}

	void init_solver_(int *s, int *sn){
		*sn = Inspector::instance()->InitSolver(*s);
	}
  
	void add_unknown(int sn, int an, int idx, int rn, int l){
		assert(sn == 0);
		Inspector::instance()->AddUnknown(sn, an, idx, rn, l);
	}

	void add_unknown_(int *sn, int *an, int *idx, int *rn, int *l){
		add_unknown(*sn, *an, *idx, *rn, *l);
	}

	void renumber_solver_rows(int sn, int* oa, int as){
		assert(sn == 0 );
		Inspector::instance()->RenumberGlobalRows(sn, oa, as);
	}

	void renumber_solver_rows_(int *sn, int* oa, int *as){
		assert(*sn == 0 );
		Inspector::instance()->RenumberGlobalRows(*sn, oa, *as);
	}

	void print_solver(){
		Inspector::instance()->print_solver();
	}

	void print_solver_(){
		Inspector::instance()->print_solver();
	}

	int get_local_rows(int sn){
		return Inspector::instance()->GetLocalRows(sn);
	}

	void get_local_rows_(int *sn, int *lr){
		*lr =  Inspector::instance()->GetLocalRows(*sn);
	}


	/**
	 * \brief Unused in quake
	 */
	void set_constraint(int an){
		Inspector::instance()->SetConstraint(an);
	}


	/**
	 * \brief Unused in quake
	 */
	void set_constraint_(int *an){
		Inspector::instance()->SetConstraint(*an);
	}

	void set_local_array(int an, void* la){
		Inspector::instance()->SetLocalArray(an,la);
	}

	void set_local_array_(int* an, void* la){
		Inspector::instance()->SetLocalArray(*an,la);
	}

	int get_local_block_offset(int an){
		return Inspector::instance()->GetLocalBlockOffset(an);
	}

	void get_local_block_offset_(int* an, int *bo){
		*bo = Inspector::instance()->GetLocalBlockOffset(*an);
	}

	void fflush_(){
		fflush(stdout);
		fflush(stderr);
	}


	/*
	 * NEW FUNCTIONS FOR PIPELINING
	 */

	void pipe_initLoop(int loop, int usedArrays[], int usedArraysSize){
		vector<int> usedArraysV(usedArrays, usedArrays + usedArraysSize);
		Inspector::instance()->init_loop(loop, usedArraysV);
	}

	void pipe_endExternalIter(){

		Inspector::instance()->pipe_endExternalIter();
	}

}
