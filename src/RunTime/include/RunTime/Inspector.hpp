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

	const int procId, nProcs, teamNum, idInTeam, teamSize;

	/// Used by the partitioner
	int pins_size;

	/// All arrays in pragma, except for indirection arrays
	std::deque<global_data*> allData;

	/// Same as global arrays but for local data. Maps array ids to arrays.
	std::map<int, local_data*> allLocalData;


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

	int* const iterNumOffset;

	int* const dataNumOffset;

	/// Node level communicator
	std::deque<global_comm*> all_comm;

	/// Communication info useful in the team. There used to be multiple of them
	/// (one per loop), but now that each process only processes one loop, there
	/// is only one.
	local_comm* team_comm;

	/// Indirection arrays are tracked here
	std::deque<access_data*> all_access_data;

	/// Singleton inspector object. There is one inspector per process.
	static Inspector* singleton_inspector;

	std::set<int>** send_info;

	/// Temporary buffers for pipeline comms. They belong to the Inspector
	/// so that they are only allocated once.
	int *pipeSendCounts, *pipeSendDispls, *pipeRecvCounts, *pipeRecvDispls;
	char *pipeSendBuf, *pipeRecvBuf;

	/// Contains, for each iteration here, the iteration in the producers after
	/// which it is safe to start computations (because we have all required
	/// data).
	std::vector<int> safeIter;

	/**
	 * \brief Returns the team ID of a process
	 *
	 * \param globalId Global ID of the process
	 */
	int get_team_id(int globalId){return globalId % teamSize;}

	/**
	 * \brief Get first global process ID in the team
	 */
	int get_first_id_in_team(){return teamNum * teamSize;}

	/**
	 * \brief Get last global process ID in the team
	 */
	int get_last_id_in_team(){return (teamNum + 1) * teamSize - 1;}

	void AfterPartition(int loop);

	void pipe_reset_counts_and_displs();

	Inspector(int pid, int np, int team, int pid_team, int teamsize,
	          int nl, int nd, int nad,
	          int* iter_num_count, int* data_num_count, int* ro);

public:
	~Inspector();

	static Inspector* instance(int pid, int np, int team, int pid_team,
	                           int teamsize, int nl, int nd, int nad,
	                           int* iter_num_count, int* data_num_count,
	                           int* ro){

		if( singleton_inspector == NULL )
			singleton_inspector =
				new Inspector(pid, np, team, pid_team, teamsize, nl, nd,
				              nad, iter_num_count, data_num_count, ro);
		return singleton_inspector;
	}

	static Inspector* instance() {
		assert(singleton_inspector);
		return singleton_inspector;
	}

	inline int get_proc_id() const{return procId;}

	inline int GetVertexHome(int loop_id, int iter) const{
		return allLoops[loop_id]->GetVertexHome(iter);
	}

	inline int GetNProcs() const {return nProcs;}
	inline int GetTeamSize() const {return teamSize;}

	inline void SetStride(int an, int st){
		allData[an]->SetStride(st);
	}

	inline int GetProcLocal(int in) const {return allLoops[in]->nproc_local;}

	void AddVertex(int, int);
	void AddPinToNet(int, int, int, int, int, int);

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

	// These are the communication functions among processes in the same team
	void CommunicateReads(int);
	void CommunicateGhosts();


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
		all_solvers.push_back(new petsc_solve(procId, nProcs, s));
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


	/*
	 * STOLEN FROM LOCAL INSPECTOR
	 */

	inline void SetupLocalArray(int dn){
		allLocalData[dn]->SetupLocalArray();
	}

	inline int GetLocalDataSize(int dn) const{
		return allLocalData.at(dn)->GetLocalSize();
	}

	/**
	 * \brief Add local array corresponding to a global_data
	 *
	 * \param mn ID of the global_data
	 * \param strideSize Number of array positions in the "pragma" stride
	 * \param ddni Nets for all the positions of the global data
	 * \param oas Size of the original array
	 * \param iro True if the array is read-only
	 * \param ic Unused in the quake benchmark
	 */
	inline void add_local_data(int mn, int strideSize, const net** ddni,
	                           int oas, bool iro, bool ic){

		local_data* new_data =
			new local_data_double(mn, teamSize, procId, strideSize, ddni,
			                      oas, iro, ic);
		allLocalData[mn] = new_data;
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
		for (std::map<int, local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++, counter++)

			it->second->GenerateGhosts();
	}


	inline void GenerateOwned(){
		int counter = 0;
		for (std::map<int, local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++, counter++)

			it->second->GenerateOwned();
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
	void PopulateLocalArray(int, double*, double*, int);


	/*
	 * NEW FUNCTIONS FOR PIPELINING
	 */

	void pipe_init_loop(const int, const int[], const bool[], const bool[],
	                    const bool[], const int);

	void pipe_calculate_comm_info();

	void pipe_init_comm_structs();

	inline void pipe_endExternalIter(){

		MPI_Barrier(global_comm::global_iec_communicator);
	}

	// Communication functions between producers and consumers
	void pipe_receive(int iter);
	void pipe_send(int iter);
};


#endif
