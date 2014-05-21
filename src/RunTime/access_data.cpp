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

access_data::access_data(int mn, int md, int np):
  my_num(mn),
  myid(md),
  nprocs(np)
{}


void access_data::SetParams(int as, int st, int* oa)
{
  array_size = as;
  stride = st;
  orig_array = oa;
  assert(stride == 1);
  int split = array_size / nprocs;
  local_start_index = split * myid;
  local_end_index =  ( myid == nprocs - 1 ? array_size : split*(myid+1));
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


bool access_data::HaveIndex(int index)
{
  bool does_have = false;
  if( index >= local_start_index && index < local_end_index )
    does_have = true; //value = orig_array[index-local_start_index];
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


void access_data::GetSendCounts(int *n_elems, int buffer_stride) const
{
  int split = array_size / nprocs;
  for( set<int>::iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ ){
    int owner_process = (*it) / split;
    if( owner_process > nprocs - 1 ) owner_process = nprocs - 1;
    assert(owner_process != myid );
    n_elems[owner_process*buffer_stride+my_num]++;
  }
}  


void access_data::PopulateBuffer(int * buffer, int buffer_size, int* curr_offset) const
{
  int split = array_size / nprocs ; 
  //if( buffer )
    for( set<int>::const_iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ ){
      int owner_process = (*it) / split;
      if( owner_process > nprocs - 1 ) owner_process = nprocs - 1;
      buffer[curr_offset[owner_process]++] = *it;
      assert(curr_offset[owner_process] <= buffer_size );
    }
}


void access_data::GetRequestedValue(int* buffer, int buffer_size, int * curr_offset, int * send_n_elems, int n_elems_stride) const
{
  int split = array_size / nprocs;
  for( int j = 0 ; j < nprocs ; j++ )
    for( int i = 0  ; i < send_n_elems[j*n_elems_stride + my_num] ; i++ ){
      //printf("ID = %d, Buffer[i] = %d, split = %d, my_num = %d\n",myid,buffer[curr_offset[j]],split,my_num);
      assert(buffer[curr_offset[j]] >= myid*split);
#ifndef NDEBUG
      int max_index = (myid == nprocs -1 ? array_size : (myid+1)*split);
      assert(buffer[curr_offset[j]] < max_index );
#endif
      int access_index = buffer[curr_offset[j]]- myid*split;
      buffer[curr_offset[j]] =  orig_array[access_index*stride];
      //printf("ID = %d, Buffer[i] = %d, access = %d , my_num = %d\n",myid,buffer[curr_offset[j]],access_index,my_num);
      curr_offset[j]++;
    }
}


void access_data::AddToHaveSet(int* buffer, int buffer_size, int* curr_offset)
{
  int split = array_size/nprocs;
  //if( buffer)
    for( set<int>::iterator it = dont_have_set.begin() ; it != dont_have_set.end() ; it++ ){
      int owner_process = (*it) / split;
      if( owner_process > nprocs - 1 ) owner_process = nprocs - 1;
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
