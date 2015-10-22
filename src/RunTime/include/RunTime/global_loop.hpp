/*
 * global_loop.hpp: This file is part of the IEC project.
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
 * @file: global_loop.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __GLOBAL_LOOP_HPP__
#define __GLOBAL_LOOP_HPP__

#include <set>
#include "RunTime/hypergraph.hpp"

class inspector;

class global_loop{

public:
	typedef std::set<int> ArrayIDList;
	typedef std::map<int, ArrayIDList> ArraysPerProcess;

private:

	typedef enum {MY_LOOP, PRODUCER, CONSUMER} LoopType;
	LoopType type;

	vertex** iter_vertex;

	const int my_num;

	const int num_iters;

	int nproc_local;

	/// Arrays needed in this loop (per producer)
	ArraysPerProcess usedArrays;

	/// Arrays calculated in this loop (per consumer)
	ArraysPerProcess computedArrays;

public:

	global_loop(int, int, int);

	~global_loop();

	/**
	 * \param iter_value Iteration of the loop
	 */
	inline int GetVertexHome(int iter_value) const{
		return iter_vertex[iter_value]->home;
	}

	void set_as_producer(){type = PRODUCER;}
	void set_as_consumer(){type = CONSUMER;}
	void set_as_my_loop(){type = MY_LOOP;}

	bool is_producer(){return type == PRODUCER;}
	bool is_consumer(){return type == CONSUMER;}
	bool is_my_loop(){return type == MY_LOOP;}

	void add_used_array(int proc, int arrayID){
		usedArrays[proc].insert(arrayID);
	}

	void add_computed_array(int proc, int arrayID){
		computedArrays[proc].insert(arrayID);
	}

	ArrayIDList::iterator computed_arrays_begin(int proc){
		return computedArrays[proc].begin();
	}

	ArrayIDList::iterator computed_arrays_end(int proc){
		return computedArrays[proc].end();
	}

	ArrayIDList::iterator used_arrays_begin(int proc){
		return usedArrays[proc].begin();
	}

	ArrayIDList::iterator used_arrays_end(int proc){
		return usedArrays[proc].end();
	}

	friend class Inspector;
};

#endif
