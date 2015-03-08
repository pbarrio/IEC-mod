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
#include "RunTime/hypergraph.hpp"

int binary_search(int * const, const int , const int );


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

	/// Identifier of this local data
	const int myid;

	/// The id of this process
	const int proc_id;

	const int stride;

	const int my_num;

	const net** const data_net_info;

	const int orig_array_size;

	const bool is_read_only;

	const bool is_constrained;
  
	int local_array_size;

	std::set<int> direct_access;

	std::set<int> indirect_access;

	long direct_access_size;
 
	int* direct_access_array;

	int indirect_access_size; 
 
	int* indirect_access_array;

	int* const ghosts_offset;

	/// Indices to the positions of this local data owned by other processes.
	int* ghosts;

	int* const owned_offset;
  
	/// Indices to the positions of this local data owned by this process.
	int* owned;

	int* l_to_g;
  
	int GetLocalIndex(int global_index) const;

	int block_owned_offset;

public:

	std::set<int> *global_ghosts;

	std::set<int> *global_owned;

	local_data(int,int,int,int,int,const net**,int,bool,bool);

	virtual ~local_data();

	virtual int GetGhostsCount(int) const=0;

	virtual int GetOwnedCount(int) const=0;

	virtual void PopulateLocalArray(double*,double*,int)=0;

	virtual void PopulateGlobalArray()=0; 

	inline int GetLocalSize() const{ return direct_access_size + indirect_access_size; }

	void AddIndexAccessed(int, int);

	void InsertDirectAccess(const int*,const int);

	void InsertIndirectAccess(const int*,const int);
  
	void RenumberAccessArray(int,int*);

	void RenumberOffsetArray(int,int*,int*);

	void GenerateGhosts();

	void GenerateOwned();

	void SetupLocalArray();

	virtual void SendOwnedData(char*, const int*) = 0;

	virtual void RecvGhostData(char*, const int*) = 0;

	virtual void SendGhostData(char*, const int*) = 0;

	virtual void RecvOwnedData(char*, const int*) = 0;

	virtual void PopulateLocalGhosts(local_data*,int)=0;

	virtual void UpdateLocalOwned(local_data*,int)=0;

	virtual void InitWriteGhosts()=0;

	virtual void print_data(FILE*);

	virtual int get_size() const=0;

	virtual int get_stride_size() const=0;

	inline int GetMyNum() const{return my_num;}

	inline int GetLocalBlockOffset() const { assert(is_constrained) ; return block_owned_offset; }

	virtual void SetLocalArray(void*)=0;

	friend class inspector;

};


class local_data_double: public local_data{

private:

	double* orig_array;

	/// Contains the local copy of the array
	double* local_array;

public:

	local_data_double(int, int, int, int, int, const net** const, int, bool, bool);
  
	~local_data_double();
  
	void PopulateLocalArray(double*,double*,int);

	inline int get_size() const{return sizeof(double);}

	inline int get_stride_size() const{return sizeof(double)*stride;}

	void print_data();

	int GetGhostsCount(int source) const{
		return (ghosts_offset[source+1]-ghosts_offset[source])*stride*sizeof(double);
	}

	int GetOwnedCount (int dest) const{
		return (owned_offset[dest+1]-owned_offset[dest])*stride*sizeof(double);
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

	void print_data(FILE*);
  
	friend class inspector;

};


#endif
