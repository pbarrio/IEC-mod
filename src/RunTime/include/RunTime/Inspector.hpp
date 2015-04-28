/*
 * Inspector.hpp: This file is part of the IEC project.
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
 * @file: Inspector.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __INSPECTOR_HPP__
#define __INSPECTOR_HPP__

#include "RunTime/hypergraph.hpp"
#include "RunTime/global_data.hpp"
#include "RunTime/global_loop.hpp"
#include "RunTime/communicator.hpp"
#include "RunTime/access_data.hpp"
#include "RunTime/external.hpp"
#include "RunTime/local_comm.hpp"
#include "RunTime/local_data.hpp"
#include <deque>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>

#define PIPE_TAG 9813 // For instance :D

//Class for the inspector to generate and partition hypergraph

extern "C" double rtclock();

#define MAX(A,B) ( A > B ? A : B)

struct ret_data_access{
	int* recvbuffer;
	int* recvcount;
	int* recvdispl;
	ret_data_access() : recvbuffer(NULL),recvdispl(NULL),recvcount(NULL) { }
	ret_data_access(int* a, int* b , int * c) :
		recvbuffer(a), recvcount(b), recvdispl(c) {}
};

class Inspector{
  
private:
  
	const int proc_id, nprocs, team_num, id_in_team, team_size;

	///Used by the partitioner
	int pins_size;

	///All arrays in pragma, except for indirection arrays.
	std::deque<global_data*> all_data;

	/// Same as above but for local data. I wonder if I could get rid of these!
	std::deque<local_data*> all_local_data;

	///All loops specified as parallel
	std::deque<global_loop*> all_loops;

	/// Loop that we process
	global_loop* my_loop;

	/// Loops that give us data (we are consumers of them)
	std::map<int, global_loop*> producer_loops;

	/// Loops that wait for our data (we are producers for them)
	std::map<int, global_loop*> consumer_loops;

	///All solvers
	std::deque<petsc_solve*> all_solvers;

	int* const iter_num_offset;

	int* const data_num_offset;

	/// Node level communicator
	std::deque<global_comm*> all_comm; 

	std::deque<local_comm*> all_local_comm;

	/// Indirection arrays are tracked here
	std::deque<access_data*> all_access_data;

	Inspector(int pid, int np, int team, int pid_team, int teamsize,
	          int nl, int nd, int nc, int nad, int* iter_num_count,
	          int* data_num_count, int* ro);
  
	///Singleton inspector object. There is one inspector per process.
	static Inspector* singleton_inspector;

	///For a given iteration, the vertex value is set first, all subsequent
	///calls to add nets uses the curr_vertex value
	vertex* curr_vertex;

#ifndef NDEBUG
	FILE* access_file;
#endif

	std::set<int>** send_info;


	void AfterPartition();


public:
	~Inspector();

	static Inspector* instance(int pid, int np, int team, int pid_team,
	                           int teamsize, int nl, int nd, int nc, int nad,
	                           int* iter_num_count, int* data_num_count,
	                           int* ro){

		if( singleton_inspector == NULL )
			singleton_inspector =
				new Inspector(pid, np, team, pid_team, teamsize, nl, nd, nc, nad,
				              iter_num_count, data_num_count, ro);
		return singleton_inspector;
	}
    
	static Inspector* instance() {
		assert(singleton_inspector);
		return singleton_inspector;
	}

	inline int get_proc_id() const{return proc_id;}

	inline int GetVertexHome(int loop_id, int iter) const{
		return all_loops[loop_id]->GetVertexHome(iter);
	}

	inline int GetNProcs() const { return nprocs; }

	inline void SetStride(int an, int st){
		all_data[an]->SetStride(st);
	}

	inline int GetProcLocal(int in) const { return all_loops[in]->nproc_local; }

	///Add a new vertex to the hypergraph
	///Arguments :
	///  1) the loop_number
	///  2) the iterator value
	void AddVertex(int,int);
  
	///Add a new pin to a net to the hypergraph
	///Arguments :
	///  1) the array_number
	///  2) the index of the element
	///  3) if the access is through a direct access ( 1=true, 0=false)
	///  4) if the access is direct from the partitionable loops ( 1 = true, 0 = false)
	void AddNet(int,int,int,int);

	void print_hypergraph(FILE*);

	void print_comm(){
#ifndef NDEBUG
		fprintf(global_comm::comm_file,
		        "Max Send size = %d, Max Recv Size = %d\n",
		        global_comm::max_send_size, global_comm::max_recv_size);
		int counter = 0;
		for (std::deque<global_comm*>::iterator it = all_comm.begin();
		     it != all_comm.end(); it++, counter++ ){

			fprintf(global_comm::comm_file,"Communicator %d\n",counter);
			(*it)->print_comm();
		}
#endif
	}

	///Combine the pieces of the hypergraph from different processes
	void PatohPrePartition();
	void PatohPartition();
	void PatohAfterPartition();

	void MetisReplicateHypergraph();
	void BlockPartition();
	void MetisPartition();
	void MetisPrePartition();
	void MetisAfterPartition();

	void CommunicateGhosts();

	void GetBufferSize();

	void CommunicateReads(int);
  
	void CommunicateToNext();


	/**
	 * \param a Identifier of the indirection array
	 * \param b Size of the array
	 * \param c Stride
	 * \param d Pointer to the start of the array piece assigned to this process
	 */
	inline void SetAccessArrayParam(int a,int b,int c ,int* d){
		all_access_data[a]->SetParams(b,c,d);
	}

	/**
	 * \brief Ask an indirection array if it knows one of its values.
	 *
	 * The value can be unknown if it's owned by another process. Remember that
	 * indirection arrays are partitioned among all the processes in a team.
	 *
	 * \param an ID of the indirection array.
	 * \param indx Position of the array to look for.
	 */
	inline bool HaveIndex(int an, int indx){
		return all_access_data[an]->HaveIndex(indx);
	}


	/**
	 * \brief Get value of a position in an indirection array.
	 *
	 * \param an ID of the indirection array
	 * \param indx Position in the array
	 *
	 * \return int Value
	 */
	inline int GetIndex(int an, int indx) const{
		return all_access_data[an]->GetIndex(indx);
	}

	bool DoneGraphGeneration();

	///Retrieve values of indirection array elements from other processes
	void GetDontHave();

	ret_data_access GetLocalAccesses(int);

	inline int InitSolver(int s){
		assert(all_solvers.size() == 0);
		all_solvers.push_back(new petsc_solve(proc_id,nprocs/*,nthreads*/,s));
		return all_solvers.size() - 1;
	}
  
	void AddUnknown(int,int,int,int);

	void RenumberGlobalRows(int sn, int* oa, int as) const{
		assert(sn == 0);
		all_solvers[sn]->RenumberGlobalRows(oa,as);
	}

	inline int GetLocalRows(int sn) const{
		assert(sn==0);
		return all_solvers[sn]->GetLocalRows();
	}

	/**
	 * \brief Unimportant for the quake benchmark
	 */
	void SetConstraint(int an) {
		all_data[an]->SetConstraint();
	}

#ifndef NDEBUG
	inline void print_access(){
		for (std::deque<access_data*>::iterator it = all_access_data.begin();
		     it != all_access_data.end(); it++)

			(*it)->print_access(access_file);
	}
#endif

	inline void print_solver(){
		char file_name[15];
		sprintf(file_name,"solver_%d.dat", proc_id);
		FILE* outfile = fopen(file_name, "w");
		all_solvers[0]->print(outfile);
		fclose(outfile);
	}


	/*
	 * STOLEN FROM LOCAL INSPECTOR
	 */

	inline void SetupLocalArray(int dn) {
		all_local_data[dn]->SetupLocalArray();
	}

	inline int GetLocalDataSize(int dn) const{
		return all_local_data[dn]->GetLocalSize();
	}

	/**
	 * \brief Add local array corresponding to a global_data
	 *
	 * \param mn ID of the global_data
	 * \param stride_size Size of a position (e.g. int = 4)
	 * \param ddni Nets for all the positions of the global data
	 * \param oas Size of the original array
	 * \param iro True if the array is read-only
	 * \param ic Unused in the quake benchmark
	 */
	inline void add_local_data(int mn, int stride_size, const net** ddni,
	                           int oas, bool iro, bool ic){

		local_data* new_data =
			new local_data_double(mn, team_size, team_num, proc_id, stride_size, ddni,
			                      oas, iro, ic);
		all_local_data.push_back(new_data);
	}

	inline void AddReadArray(int in, int an){
		all_local_comm[in]->read_arrays.push_back(all_local_data[an]);
	}

	inline void AddWriteArray(int in, int an){
		all_local_comm[in]->write_arrays.push_back(all_local_data[an]);
	}

	inline void AddIndexAccessed(int dn, int ind, int at){
		all_local_data[dn]->AddIndexAccessed(ind,at);
	}

	inline void PopulateLocalArray(int an, double* lb, double* oa, int st){
		all_local_data[an]->PopulateLocalArray(lb,oa,st);
	}

	inline void InsertDirectAccess(int an, int* v, int n){
		all_local_data[an]->InsertDirectAccess(v,n);
	}

	inline void InsertIndirectAccess(int an, int* v, int n){
		all_local_data[an]->InsertIndirectAccess(v,n);
	}

	inline void RenumberAccessArray(int an, int as, int* aa){
		all_local_data[an]->RenumberAccessArray(as,aa);
	}

	inline void RenumberOffsetArray(int an, int as, int* aa, int* la){
		all_local_data[an]->RenumberOffsetArray(as,aa,la);
	}

	inline void GenerateGhosts(){
		int counter =0;
		for (std::deque<local_data*>::iterator it = all_local_data.begin();
		     it != all_local_data.end(); it++, counter++)

			(*it)->GenerateGhosts();
	}

	inline void GenerateOwned(){
		int counter = 0;
		for (std::deque<local_data*>::iterator it = all_local_data.begin();
		     it != all_local_data.end(); it++, counter++)

			(*it)->GenerateOwned();
	}

	inline void InitWriteGhosts(int cn){
		all_local_comm[cn]->InitWriteGhosts();
	}

	inline void SetLocalArray(int an, void* la) {
		all_local_data[an]->SetLocalArray(la);
	}

	inline int GetLocalBlockOffset(int an) {
		return all_local_data[an]->GetLocalBlockOffset();
	}

	void print_local_data();
	void PopulateGlobalArrays();
	void print_local_inspector_data();


	/*
	 * NEW FUNCTIONS FOR PIPELINING
	 */

	inline void pipe_endExternalIter(){

		MPI_Barrier(global_comm::global_iec_communicator);
	}

};


#endif
