/*
 * access_data.hpp: This file is part of the IEC project.
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
 * @file: access_data.hpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/access_data.hpp"
#include <cassert>

using namespace std;

/**
 * \param mn ID for this indirection array
 * \param md ID of this process
 * \param np Number of processors
 * \param teamID ID of this processor in the team
 * \param teamSize Size of the team
 */
access_data::access_data(int mn, int md, int np, int teamID, int teamSize):
  my_num(mn),
  myid(md),
  nprocs(np),
  myTeamID(teamID),
  myTeamSize(teamSize)
{}


/**
 * \param as Size of the entire original array
 * \param st Stride
 * \param oa Pointer to the start of the array piece assigned to this process.
 */
void access_data::SetParams(int as, int st, int* oa)
{
  array_size = as;
  stride = st;
  orig_array = oa;
  assert(stride == 1);
  int split = array_size / myTeamSize;
  local_start_index = split * myTeamID;
  local_end_index =  ( myTeamID == myTeamSize - 1 ? array_size : split * (myTeamID + 1));
#ifdef USE_HSQPHASH
  have_set = new qphash(401);
#endif
}


access_data::~access_data() 
{
#ifdef USE_HSQPHASH
  delete have_set;
#else
  have_set.clear();
#endif
  dont_have_set.clear();
}

/**
 * \brief Checks if this process knows the value of the array position.
 *
 * This process has two possibilities to know the value of this array position.
 * The first is that the position falls within the array chunk initially assigned
 * to this process (between local_start_index and local_end_index). The second is
 * that another process has shared that information with the current process.
 *
 * \param index Position in the indirection array
 * \return true if the process knows the value, false otherwise
 */
bool access_data::HaveIndex(int index)
{
  bool does_have = false;

  // Do I have the information from the beginning?
  if( index >= local_start_index && index < local_end_index )
    does_have = true; //value = orig_array[index-local_start_index];

  // Did I receive the information from another process?
  else{
#ifdef USE_HSQPHASH
    int value = have_set->return_key_val(index);
    if( value == -1 )
      dont_have_set.insert(index);
    else
      does_have = true;
#else
    map<int,int>::iterator ret = have_set.find(index);
    if( ret != have_set.end() ){
      does_have = true; //value = (*ret).second;
    }
    else
      dont_have_set.insert(index);
#endif
  }
  //printf("MD:ID:%d,IndArray:%d,Index:%d,DoesHave:%d\n",myid,my_num,index,does_have);
  return does_have;
}


/**
 * \brief Get the value of a position in an indirection array
 *
 * If the position is owned by this process from the very beginning, return it
 * immediately. If not, it should have been communicated to us by another process,
 * and we can find it in the "have_set" map.
 *
 * \param index Position of the array
 *
 * \return int Value of the position
 */
int access_data::GetIndex(int index) const
{
  int value = -1;
  if( index >= local_start_index && index < local_end_index)
    value = orig_array[index-local_start_index];
  else{
#ifdef USE_HSQPHASH
    value = have_set->return_key_val(index);
    //assert(value != -1);
#else
    map<int,int>::const_iterator ret = have_set.find(index);
    if( ret == have_set.end() )
      printf("ID=%d,Didnt Find %d of array %d\n",myid,index,my_num);
    assert(ret != have_set.end());
    value = (*ret).second;
#endif
  }
  return value;
}


/**
 * \brief Calculates how much data to receive from each other process
 *
 * \param n_elems After the function, the i-th position of this array must contain the #elems to receive from process i
 * \param buffer_stride Equal to the number of indirection arrays 
 */
void access_data::GetSendCounts(int *n_elems, int buffer_stride) const
{
  int split = array_size / myTeamSize;
  for( set<int>::iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ ){
    int owner_process = (*it) / split;
    if( owner_process > myTeamSize - 1 ) owner_process = myTeamSize - 1;
    assert(owner_process != myid );
    n_elems[owner_process*buffer_stride+my_num]++;
  }
}  


void access_data::PopulateBuffer(int * buffer, int buffer_size, int* curr_offset) const
{
  int split = array_size / myTeamSize ; 
  //if( buffer )
    for( set<int>::const_iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ ){
      int owner_process = (*it) / split;
      if( owner_process > myTeamSize - 1 ) owner_process = myTeamSize - 1;
      buffer[curr_offset[owner_process]++] = *it;
      assert(curr_offset[owner_process] <= buffer_size );
    }
}


void access_data::GetRequestedValue(int* buffer, int buffer_size, int * curr_offset, int * send_n_elems, int n_elems_stride) const
{
  int split = array_size / myTeamSize;
  for( int j = 0 ; j < myTeamSize ; j++ )
    for( int i = 0  ; i < send_n_elems[j*n_elems_stride + my_num] ; i++ ){
      assert(buffer[curr_offset[j]] >= myid*split);
#ifndef NDEBUG
      int max_index = (myid == myTeamSize -1 ? array_size : (myid+1)*split);
      assert(buffer[curr_offset[j]] < max_index );
#endif
      int access_index = buffer[curr_offset[j]]- myid*split;
      buffer[curr_offset[j]] =  orig_array[access_index*stride];
      curr_offset[j]++;
    }
}


void access_data::AddToHaveSet(int* buffer, int buffer_size, int* curr_offset)
{
  int split = array_size / myTeamSize;
  for( set<int>::iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ ){
	  int owner_process = (*it) / split;
	  if( owner_process > myTeamSize - 1 ) owner_process = myTeamSize - 1;
	  int value = buffer[curr_offset[owner_process]++];
#ifdef USE_HSQPHASH
	  have_set->hash_insert_only((*it),value);
#else
      have_set.insert(pair<int,int>((*it),value));
#endif
      assert(curr_offset[owner_process] <= buffer_size );
    }  
  dont_have_set.clear();
}


#ifndef NDEBUG
void access_data::print_access(FILE* outfile)
{
  fprintf(outfile,"Access Array = %d, Array_size = %d, stride = %d, orig_array = %p\n",my_num,array_size,stride,orig_array);
#ifndef USE_HSQPHASH
  fprintf(outfile,"\tHave Set (%d):",(int)have_set.size());
#ifdef HIGH_DETAILS
  for( map<int,int>::iterator it = have_set.begin() ; it != have_set.end() ; it++ )
    fprintf(outfile," (%d->%d)",(*it).first,(*it).second);
#endif
  fprintf(outfile,"\n");
  fprintf(outfile,"\tDont Have Set (%d):",(int)dont_have_set.size());
#ifdef HIGH_DETAILS
  for( set<int>::iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ )
    fprintf(outfile," %d",*it);
#endif
  fprintf(outfile,"\n\n");
  fflush(outfile);
#endif
}
#endif
