/*
 * local_data.hpp: This file is part of the IEC project.
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
 * @file: local_data.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __LOCAL_DATA_HPP__
#define __LOCAL_DATA_HPP__

#include <cstdio>
#include <cassert>
#include <vector>

#include "RunTime/hypergraph.hpp"
#include "RunTime/global_data.hpp"

int binary_search(int* const, const int, const int);


/**
 * \brief Local handler of a global array (global_data)
 *
 * I *think* that all "local_*" classes were introduced to allow
 * OpenMP threads. Now that we don't have threads, it might be
 * possible to get rid of them.
 */
class local_data{

protected:

	/// Number of parts in which this array will be divided
	const int nparts;

	/// The id of this process
	const int proc_id;

	const int stride;

	/// Id of this local array
	const int my_num;

	const net** const data_net_info;

	const int orig_array_size;

	const bool is_read_only;

	const bool is_constrained;

	int local_array_size;

	/// Temp structure used to calculate the direct accesses
	std::set<int> direct_access;

	/// Number of direct accesses
	long direct_access_size;

	/// Mapping from local to global indexes for direct accesses. If
	/// direct_access_array[a] = b, then index 'b' in the global array
	/// corresponds to a local index 'a'.
	int* direct_access_array;

	/// Temp structure used to calculate the indirect accesses
	std::set<int> indirect_access;

	/// Number of indirect accesses
	long indirect_access_size;

	/// Same as direct_access_array but for indirect accesses
	int* indirect_access_array;

	int* const ghosts_offset;

	/// Indices to the positions of this local data owned by other processes.
	int* ghosts;

	int* const owned_offset;

	/// Indices to the positions of this local data owned by this process.
	int* owned;

	int* l_to_g;

	int block_owned_offset;

	/// Maximum number of sent and received items (in bytes)
	int maxSendCounts, maxRecvCounts;

	/// For each iteration and producer, the receive size in bytes.
	global_data::CountsPerProcPerIter pipeRecvCounts;

	/// For each iteration and consumer, the send size in bytes.
	global_data::CountsPerProcPerIter pipeSendCounts;

	/// For each iteration and producer, a list of positions in the local array
	/// that we need to receive from the the producer.
	global_data::IdxsPerProcPerIter pipeRecvIndexes;

	/// For each iteration and consumer, a list of positions in the local array
	/// that we need to send to the consumer.
	global_data::IdxsPerProcPerIter pipeSendIndexes;

	int GetLocalIndex(int global_index) const;

public:

	std::set<int> *global_ghosts;

	std::set<int> *global_owned;

	local_data(int, int, int, int, const net**, int, bool, bool);

	virtual ~local_data();

	virtual int GetGhostsCount(int) const = 0;

	virtual int GetOwnedCount(int) const = 0;

	virtual void PopulateLocalArray
	(const global_data*, double*, double*, int) = 0;

	virtual void PopulateGlobalArray() = 0;

	inline int GetLocalSize() const {
		return direct_access_size + indirect_access_size;
	}

	void AddIndexAccessed(int, int);

	void InsertDirectAccess(const int*, const int);

	void InsertIndirectAccess(const int*, const int);

	void RenumberAccessArray(int, int*);

	void RenumberOffsetArray(int, int*, int*);

	void GenerateGhosts();

	void GenerateOwned();

	void SetupLocalArray();

	virtual void SendOwnedData(char*, const int*) = 0;

	virtual void RecvGhostData(char*, const int*) = 0;

	virtual void SendGhostData(char*, const int*) = 0;

	virtual void RecvOwnedData(char*, const int*) = 0;

	virtual void PopulateLocalGhosts(local_data*, int) = 0;

	virtual void UpdateLocalOwned(local_data*, int) = 0;

	virtual void InitWriteGhosts() = 0;

	virtual int get_size() const = 0;

	virtual int get_stride_size() const = 0;

	inline int GetMyNum() const{return my_num;}

	inline int GetLocalBlockOffset() const {
		assert(is_constrained);
		return block_owned_offset;
	}

	virtual void SetLocalArray(void*) = 0;

	/**
	 * \brief Get the send size in bytes for an iteration and process
	 */
	virtual int pipe_get_sendcounts(int iter, int proc) = 0;

	/**
	 * \brief Get the receive size in bytes for an iteration and process
	 */
	virtual int pipe_get_recvcounts(int iter, int proc) = 0;

	virtual int pipe_get_max_sendcounts() = 0;

	virtual int pipe_get_max_recvcounts() = 0;

	virtual int pipe_populate_send_buf(int, int, char*) = 0;
	virtual void pipe_update(int, int, char*) = 0;

	friend class Inspector;
};


class local_data_double: public local_data{

private:

	double* orig_array;

	/// Contains the local copy of the array
	double* local_array;

public:

	local_data_double(int, int, int, int, const net** const, int, bool, bool);

	~local_data_double();

	void PopulateLocalArray(const global_data*, double*, double*, int);

	inline int get_size() const{return sizeof(double);}

	inline int get_stride_size() const{return sizeof(double) * stride;}

	void print_data();

	int GetGhostsCount(int source) const{
		return (ghosts_offset[source + 1] - ghosts_offset[source]) * stride *
			sizeof(double);
	}

	int GetOwnedCount (int dest) const{
		return (owned_offset[dest + 1] - owned_offset[dest]) * stride *
			sizeof(double);
	}

	void SendOwnedData(char*, const int*);

	void RecvGhostData(char*, const int*);

	void SendGhostData(char*, const int*);

	void RecvOwnedData(char*, const int*);

	void InitWriteGhosts();

	void PopulateLocalGhosts(local_data*,int);

	void UpdateLocalOwned(local_data*,int);

	void PopulateGlobalArray();

	inline void SetLocalArray(void* la){
		assert( is_constrained && local_array == NULL );
		local_array = static_cast<double*>(la);
	}

	int pipe_get_sendcounts(int iter, int proc){
		return pipeSendIndexes[iter][proc].size() * sizeof(double) * stride;
	}

	int pipe_get_recvcounts(int iter, int proc){
		return pipeRecvIndexes[iter][proc].size() * sizeof(double) * stride;
	}

	int pipe_get_max_sendcounts(){
		return maxSendCounts * sizeof(double) * stride;
	}

	int pipe_get_max_recvcounts(){
		return maxRecvCounts * sizeof(double) * stride;
	}

	int pipe_populate_send_buf(int, int, char*);

	void pipe_update(int, int, char*);

	friend class Inspector;
};


#endif
