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

protected:

	/// The identifier assigned to this global_data.
	const int id;

	/// All the nets, one for each position of this global array.
	net**  data_net_info;

	/// Size of the original array
	const int orig_array_size;

	/// Size of each element of the array in bytes (e.g. for int. stride_size=4)
	int stride_size;

	/// True if the array is read-only in the original code
	const bool is_read_only;

	/// Unimportant for the quake benchmark
	bool is_constrained;

public:
	global_data(int,int,bool);

	virtual ~global_data();

	virtual void SetStride(int) = 0;

	virtual int get_stride_size() const = 0;

	/**
	 * \brief Unimportant for the quake benchmark
	 */
	inline void SetConstraint(){is_constrained = true;}

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
