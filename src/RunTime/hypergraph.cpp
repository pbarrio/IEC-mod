/*
 * hypergraph.cpp: This file is part of the IEC project.
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
 * @file: hypergraph.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include <iostream>
#include "RunTime/hypergraph.hpp"

using namespace std;

vertex::vertex(short in, int iv, int mn):
  iter_num(in),
  iter_value(iv),
  my_num(mn)
{
  home = -1;
#if defined USE_QPHASH 
  adjvertex = new qphash(197);
#elif defined USE_KSHASH
  adjvertex = new kshash(197);
#endif
}

vertex::~vertex()
{
#if defined USE_QPHASH || defined USE_KSHASH
  if( adjvertex )
    delete adjvertex;
#else
  adjvertex.clear();
#endif
}


double vertex::clear_adjvertex()
{
#if defined USE_QPHASH || defined USE_KSHASH
  double ret_val = 0.0 ;//adjvertex->metric();
  delete adjvertex;
  adjvertex = NULL;
#else
  double ret_val = 0.0;
  adjvertex.clear();
#endif  

  return ret_val;
}


/**
 * \brief Constructor
 *
 * Every position of an array has a corresponding net.
 * Hence, a net is defined by the array ID + index.
 *
 * \param dn ID of the array
 * \param di Index of the array position
 * \param wt ??? Looks like it's initialized to the size of the data (e.g. int = 4)
 * \param mn ID of the net
 */
net::net(short dn, int di, short wt, int mn):
  data_num(dn),
  data_index(di),
  weight(wt),
  my_num(mn),
  direct_vertex(NULL)
{
  home = -1;
}

net::~net()
{
  pins.clear();
}

