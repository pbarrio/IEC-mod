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

private:
	/// For each loop, true if the array is read in the loop.
	std::map<int, bool> isReadInLoop;

	/// For each loop, true if the array is updated (written).
	/// This is only meaningful for pipelining. The original IEC only required
	/// to know if the array was RO/RW throughout the entire code, not in a
	/// per-loop basis.
	std::map<int, bool> isWriteInLoop;

	/// Map to the write/use information for this array for all loops.
	/// Consumers don't need to keep track of these values. Note that these
	/// specification is not enforced in this class (TODO) but the class gives
	/// methods use_in_loop() and is_already_used_in_loop() to enforce it.
	ArrayUseMap finalUse;

protected:

	/// Pointers to nets for each array value, for each team
	typedef std::map<int, net**> LoopNets;

	/// The identifier assigned to this global_data.
	const int id;

	/// All the nets, one for each position of this global array, for each team.
	/// e.g. data_net_info[2][6] is a net* with the net for position 6, loop 2.
	LoopNets data_net_info;

	/// Size of the original array
	const int orig_array_size;

	/// Size of each element of the array in bytes (e.g. for int. stride_size=4)
	int stride_size;

	/// Offset: address of the first elem if we put all global data in a single
	/// buffer in ID order.
	int offset;

	/// True if the array is read-only in the original code
	const bool is_read_only;

	/// Unimportant for the quake benchmark
	bool is_constrained;

public:
	global_data(int, int, int, bool);

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

	bool is_read(int loop){return isReadInLoop[loop];}
	bool is_write(int loop){return isWriteInLoop[loop];}

	friend class Inspector;
};

class global_data_double: public global_data{
public:

	global_data_double(int, int, int, bool);

	~global_data_double();

	void SetStride(int);

	inline int get_stride_size() const{ return stride_size; }
};


#endif
