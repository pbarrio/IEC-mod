/*
 * hypergraph.hpp: This file is part of the IEC project.
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
 * @file: hypergraph.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#ifndef __HYPERGRAPH_HPP__
#define __HYPERGRAPH_HPP__

#include <set>
#include <map>

#ifdef USE_QPHASH
#include "DataStructures/qphash.hpp"
#elif USE_KSHASH
#include "DataStructures/kshash.hpp"
#endif

/** 
 * \brief A vertex in the hypergraph
 *
 * Each vertex represents one iteration of a loop
 */
struct vertex{
  const short iter_num;
  const int iter_value;
  const int my_num;
  short home;
  vertex(short,int,int);
  ~vertex();
  double clear_adjvertex();
#ifdef USE_QPHASH
  qphash* adjvertex;
#elif USE_KSHASH
  kshash* adjvertex;
#elif defined USE_SORT
  int* adjvertex;
  int nadjvertex;
#else
  std::map<int,int> adjvertex;
#endif
};


struct pin_info{
  vertex* pin;
  bool is_direct;
pin_info(vertex* a, bool b): pin(a),is_direct(b) {}
};


struct pin_comparator{
  bool operator()( const pin_info lhs , const pin_info rhs ) const{
    return lhs.pin->my_num < rhs.pin->my_num;
  }
};

/**
 * \brief A net in the hypergraph
 *
 * Each net represents one position of an array.
 */
struct net{

  /// Identifier of the array
  const short data_num;

  /// Index of the data (relative to the array) that this net represents
  const int data_index;

  /// Unique net identifier
  const int my_num;

  short home;
  std::set<pin_info,pin_comparator> pins;
  short weight;
  vertex* direct_vertex;
  net(short,int,short,int);
  ~net();
};

#endif
