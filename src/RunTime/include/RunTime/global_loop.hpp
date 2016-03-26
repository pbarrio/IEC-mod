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
	typedef std::map<int, ArrayIDList> ArraysPerIter;

private:

	typedef enum {MY_LOOP, PRODUCER, CONSUMER} LoopType;
	LoopType type;

	vertex** iter_vertex;

	const int my_num;

	const int num_iters;

	int nproc_local;

	/// Arrays needed in each iteration of this loop
	ArraysPerIter usedArrays;

	/// Arrays calculated in this loop (per consumer)
	ArraysPerIter computedArrays;

public:

	global_loop(int, int, int);

	~global_loop();

	int get_loop_id(){return my_num;}

	/**
	 * \brief Get the owner process of this loop iteration
	 *
	 * Each iteration of the loop is run in only one process. This method
	 * returns the process that will execute the iteration (the owner).
	 *
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

	void add_used_array(int iter, int arrayID){
		usedArrays[iter].insert(arrayID);
	}

	void add_computed_array(int iter, int arrayID){
		computedArrays[iter].insert(arrayID);
	}

	ArrayIDList::iterator computed_arrays_begin(int iter){
		return computedArrays[iter].begin();
	}

	ArrayIDList::iterator computed_arrays_end(int iter){
		return computedArrays[iter].end();
	}

	ArrayIDList::iterator used_arrays_begin(int iter){
		return usedArrays[iter].begin();
	}

	ArrayIDList::iterator used_arrays_end(int iter){
		return usedArrays[iter].end();
	}

	friend class Inspector;
};

#endif
