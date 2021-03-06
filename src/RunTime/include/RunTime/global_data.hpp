/*
 * global_data.hpp: This file is part of the IEC project.
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
 * @file: global_data.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __GLOBAL_DATA_HPP__
#define __GLOBAL_DATA_HPP__

#include <map>
#include <vector>
#include "RunTime/hypergraph.hpp"

class inspector;

/**
 * \brief An array of data used in the code's computations.
 *
 * Indirection arrays are explicitly excluded.
 */
class global_data{

public:
	/// For each position of an array, the iteration of the:
	///     - last write, if this is a producer loop. We need to know when the
	///       loop can safely communicate the final value.
	///     - first use, if this is my loop. We need to know when we can start
	///       an iteration because we already received all data from producers.
	typedef std::map<int, int> ArrayUseMapInLoop;

	/// A map to write/use information for all indices, for all loops
	typedef std::map<int,  ArrayUseMapInLoop> ArrayUseMap;

	/// Pointers to nets for each array value, for each team
	typedef std::map<int, net**> LoopNets;

	typedef std::map<int, std::map<int, int> > CountsPerProcPerIter;

	// map < iterations, map < processes, vector<indexes> > >
	typedef std::map<int, std::map<int, std::vector<int> > > IdxsPerProcPerIter;

	// map < processes, vector<indexes> >
	typedef std::map<int, std::vector<int> > IdxsPerProc;

private:
	/// For each loop, true if the array is read in the loop.
	std::map<int, bool> isReadInLoop;

	/// For each loop, true if the array is updated (written).
	/// This is only meaningful for pipelining. The original IEC only required
	/// to know if the array was RO/RW throughout the entire code, not in a
	/// per-loop basis.
	std::map<int, bool> isWriteInLoop;

	/// Id of the loop that is the last in the pipeline to write to the array.
	int lastWriteInPipeline;

	/// List of loop IDs that require data to be transferred from the previous
	/// outer iteration.
	std::set<int> needsInterIterReceive;

	/// Map to the write/use information for this array for all loops.
	/// Consumers don't need to keep track of these values. Note that these
	/// specification is not enforced in this class (TODO) but the class gives
	/// methods use_in_loop() and is_already_used_in_loop() to enforce it.
	ArrayUseMap finalUse;

	/// Track the first time that a position is used. This allows consumers to
	/// contrast their first uses with the data received from producers and
	/// decide whether they can start the current iteration or they need to wait
	/// for more data.
	ArrayUseMap initialUse;

	/// The producer that gives us the data (in case we use this array)
	int producer;

protected:

	/// Our process Id
	int procId;

	/// The identifier assigned to this global_data.
	const int id;

	/// All the nets, one for each position of this global array, for each team.
	/// e.g. data_net_info[2][6] is a net* with the net for position 6, loop 2.
	LoopNets data_net_info;

	/// Size of the original array
	const int orig_array_size;

	/// Size of each group of elements of the array in bytes. e.g. for an array
	/// of 3D integer coordinates where each element = (x, y, z), we will have
	/// 4 bytes/component * 3 components = 12. The stride is 12 bytes.
	int stride_size;

	/// Offset: address of the first elem if we put all global data in a single
	/// buffer in ID order.
	int offset;

	/// True if the array is read-only in the original code
	const bool is_read_only;

	/// Unimportant for the quake benchmark
	bool is_constrained;

	/// For each iteration and producer, the receive size in bytes.
	CountsPerProcPerIter pipeRecvCounts;

	/// For each iteration and consumer, the send size in bytes.
	CountsPerProcPerIter pipeSendCounts;

	/// For each iteration and producer, a list of positions in the global array
	/// that we need to receive from the producer.
	IdxsPerProcPerIter pipeRecvIndexes;

	/// For each iteration and consumer, a list of positions in the global array
	/// that we need to send to the consumer.
	IdxsPerProcPerIter pipeSendIndexes;

	/// For each process that needs data from us in the next iteration, a list
	/// of positions that we need to send.
	IdxsPerProc pipeInterIterSendIndexes;

	/// For each process from the previous iteration that needs to send us data,
	/// a list of positions that we need to receive.
	IdxsPerProc pipeInterIterRecvIndexes;

	/// For each iteration, the corresponding iteration in the producer that
	/// ensures that we have all the required data (in this array) to start
	/// doing the computations. THIS IS PROBABLY WRONG, because we may have more
	/// than one producer per array (e.g. different processes in the same team).
	std::vector<int> producerIter;

public:
	global_data(int, int, int, int, bool);

	virtual ~global_data();

	virtual void SetStride(int) = 0;

	virtual int get_stride_size() const = 0;

	/**
	 * \brief Get this array ID
	 */
	int getId(){return id;}

	/**
	 * \brief Unimportant for the quake benchmark
	 */
	inline void SetConstraint(){is_constrained = true;}

	void use_in_loop(int, bool, bool);

	void use_in_loop(int, int, int);

	/**
	 * \brief Check if this array is used in the current loop
	 *
	 * \param loop Loop ID
	 */
	bool is_used_in_loop(int loop){
		if (data_net_info.find(loop) != data_net_info.end())
			return true;
		return false;
	}

	/**
	 * \brief See if a specific index has been already marked used in the loop
	 *
	 * \param index Index in this array to be checked
	 * \param loopID Loop identifier
	 */
	bool is_already_used_in_loop(int index, int loopID){
		return (finalUse[loopID].find(index) != finalUse[loopID].end());
	}

	int get_final_use(int loop, int index){
		return finalUse[loop][index];
	}

	bool is_read(int loop){return isReadInLoop[loop];}
	bool is_write(int loop){return isWriteInLoop[loop];}
	bool is_last_write_in_pipeline(int p){return lastWriteInPipeline == p;}
	unsigned get_last_write_in_pipeline(){return lastWriteInPipeline;}
	void set_last_write_in_pipeline(int p){lastWriteInPipeline = p;}

	bool needs_interiter_receive(int p){
		return needsInterIterReceive.count(p) != 0;
	}

	void set_needs_interiter_receive(int p){
		needsInterIterReceive.insert(p);
	}

	void pipe_calc_sends(int myLoop);
	void pipe_calc_recvs();
	void pipe_calc_interiter_sends();
	void pipe_calc_interiter_recvs();
	int pipe_safe_iteration(int index) {
		return producerIter.size() > index ? producerIter[index] : -1;
	}

	friend class Inspector;
	friend class local_data;
	friend class local_data_double;
};

class global_data_double: public global_data{
public:

	global_data_double(int, int, int, int, bool);

	~global_data_double();

	void SetStride(int);

	inline int get_stride_size() const{ return stride_size; }
};


#endif
