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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <vector>

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

	/// Used by the partitioner
	int pins_size;

	/// All arrays in pragma, except for indirection arrays.
	std::deque<global_data*> allData;

	/// Same as above but for local data. I wonder if I could get rid of these!
	std::deque<local_data*> allLocalData;


	/*
	 * LOOP INFO
	 */

	/// All loops specified as parallel
	std::deque<global_loop*> allLoops;

	/// Loop that we process
	global_loop* myLoop;

	/// Loops that give us data (we are consumers of them)
	std::map<int, global_loop*> producerLoops;

	/// Loops that wait for our data (we are producers for them)
	std::map<int, global_loop*> consumerLoops;

	/*-----------------*/


	/// All solvers
	std::deque<petsc_solve*> all_solvers;

	int* const iter_num_offset;

	int* const data_num_offset;

	/// Node level communicator
	std::deque<global_comm*> all_comm;

	/// Communication info useful in the team. There used to be multiple of them
	/// (one per loop), but now that each process only processes one loop, there
	/// is only one.
	local_comm* team_comm;

	/// Indirection arrays are tracked here
	std::deque<access_data*> all_access_data;

	Inspector(int pid, int np, int team, int pid_team, int teamsize,
	          int nl, int nd, int nc, int nad, int* iter_num_count,
	          int* data_num_count, int* ro);

	/// Singleton inspector object. There is one inspector per process.
	static Inspector* singleton_inspector;

	/// For a given iteration, the vertex value is set first, all subsequent
	/// calls to add nets uses the curr_vertex value
	vertex* curr_vertex;

	std::set<int>** send_info;

	void AfterPartition(int loop);


public:
	~Inspector();

	static Inspector* instance(int pid, int np, int team, int pid_team,
	                           int teamsize, int nl, int nd, int nc, int nad,
	                           int* iter_num_count, int* data_num_count,
	                           int* ro){

		if( singleton_inspector == NULL )
			singleton_inspector =
				new Inspector(pid, np, team, pid_team, teamsize, nl, nd, nc,
				              nad, iter_num_count, data_num_count, ro);
		return singleton_inspector;
	}

	static Inspector* instance() {
		assert(singleton_inspector);
		return singleton_inspector;
	}

	inline int get_proc_id() const{return proc_id;}

	inline int GetVertexHome(int loop_id, int iter) const{
		return allLoops[loop_id]->GetVertexHome(iter);
	}

	inline int GetNProcs() const {return nprocs;}
	inline int GetTeamSize() const {return team_size;}

	inline void SetStride(int an, int st){
		allData[an]->SetStride(st);
	}

	inline int GetProcLocal(int in) const {return allLoops[in]->nproc_local;}

	void init_loop(int, std::vector<int>);
	void AddVertex(int, int);
	void AddPinToNet(int, int, int, int, int);

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

	// Combine the pieces of the hypergraph from different processes
	void PatohPartitionAll();
	void PatohPrePartition(int loop);
	void PatohPartition(int loop);
	void PatohAfterPartition(int loop);

	void MetisReplicateHypergraph(int loop);
	void BlockPartition();
	void MetisPartitionAll();
	void MetisPartition(int loop);
	void MetisPrePartition(int loop);
	void MetisAfterPartition(int loop);

	void GetBufferSize();

	/*
	 * These are the communication functions. The first two communicate
	 * data among processes in the same team. The last two communicate
	 * data between producers and consumers.
	 */
	void CommunicateReads(int);
	void CommunicateGhosts();
	void CommunicateToNext();
	void GetFromPrevious();


	/**
	 * \param a Identifier of the indirection array
	 * \param b Size of the array
	 * \param c Stride
	 * \param d Pointer to the start of the array piece assigned to this process
	 */
	inline void SetAccessArrayParam(int a, int b, int c, int* d){
		all_access_data[a]->SetParams(b, c, d);
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

	void GetLocalAccesses(int, int**, int**, int**);

	inline int InitSolver(int s){
		assert(all_solvers.size() == 0);
		all_solvers.push_back(new petsc_solve(proc_id, nprocs, s));
		return all_solvers.size() - 1;
	}

	void AddUnknown(int, int, int, int, int);

	void RenumberGlobalRows(int sn, int* oa, int as) const{
		assert(sn == 0);
		all_solvers[sn]->RenumberGlobalRows(oa, as);
	}

	inline int GetLocalRows(int sn) const{
		assert(sn==0);
		return all_solvers[sn]->GetLocalRows();
	}

	/**
	 * \brief Unimportant for the quake benchmark
	 */
	void SetConstraint(int an){
		allData[an]->SetConstraint();
	}

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

	inline void SetupLocalArray(int dn){
		allLocalData[dn]->SetupLocalArray();
	}

	inline int GetLocalDataSize(int dn) const{
		return allLocalData[dn]->GetLocalSize();
	}

	/**
	 * \brief Add local array corresponding to a global_data
	 *
	 * \param mn ID of the global_data
	 * \param stride_size Number of array positions in the "pragma" stride
	 * \param ddni Nets for all the positions of the global data
	 * \param oas Size of the original array
	 * \param iro True if the array is read-only
	 * \param ic Unused in the quake benchmark
	 */
	inline void add_local_data(int mn, int stride_size, const net** ddni,
	                           int oas, bool iro, bool ic){

		local_data* new_data =
			new local_data_double(mn, team_size, proc_id, stride_size, ddni,
			                      oas, iro, ic);
		allLocalData.push_back(new_data);
	}

	/*
	 * \brief Add array to the set of data that is read in this process
	 *
	 * \param an Id of the array
	 */
	inline void AddReadArray(int an){
		team_comm->read_arrays.push_back(allLocalData[an]);
	}

	/*
	 * \brief Add array to the set of data that is written in this process
	 *
	 * \param an Id of the array
	 */
	inline void AddWriteArray(int an){
		team_comm->write_arrays.push_back(allLocalData[an]);
	}


	inline void AddIndexAccessed(int dn, int ind, int at){
		allLocalData[dn]->AddIndexAccessed(ind, at);
	}

	/**
	 * \brief Populates a local array from its corresponding global array
	 *
	 * \param an ID of the local array to be populated
	 * \param lb Allocated clean array to be populated
	 * \param oa Original array
	 * \param st Stride. Not sure what this is.
	 */
	inline void PopulateLocalArray(int an, double* lb, double* oa, int st){
		all_local_data[an]->PopulateLocalArray(lb, oa, st);
	}

	inline void InsertDirectAccess(int an, int* v, int n){
		allLocalData[an]->InsertDirectAccess(v, n);
	}


	inline void InsertIndirectAccess(int an, int* v, int n){
		allLocalData[an]->InsertIndirectAccess(v, n);
	}


	/**
	 * \param an ID of the local array to be populated
	 * \param as Size of this local array
	 * \param aa Array that will be used to index this local array
	 */
	inline void RenumberAccessArray(int an, int as, int* aa){
		allLocalData[an]->RenumberAccessArray(as, aa);
	}


	inline void RenumberOffsetArray(int an, int as, int* aa, int* la){
		allLocalData[an]->RenumberOffsetArray(as, aa, la);
	}


	inline void GenerateGhosts(){
		int counter = 0;
		for (std::deque<local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++, counter++)

			(*it)->GenerateGhosts();
	}


	inline void GenerateOwned(){
		int counter = 0;
		for (std::deque<local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++, counter++)

			(*it)->GenerateOwned();
	}

	inline void InitWriteGhosts(int cn){
		team_comm->InitWriteGhosts();
	}

	inline void SetLocalArray(int an, void* la) {
		allLocalData[an]->SetLocalArray(la);
	}

	inline int GetLocalBlockOffset(int an) {
		return allLocalData[an]->GetLocalBlockOffset();
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
