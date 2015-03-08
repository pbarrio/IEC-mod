/*
 * inspector_metis.cpp: This file is part of the IEC project.
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
 * @file: inspector_metis.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "mpi.h"
#include "RunTime/inspector.hpp"
#include "parmetis.h"
#include <cassert>
//#include "iecxx.h"

using namespace std;

void insertion_sort(int* const array, const int left, const int right, int* const mirror_array)
{
  if( left >= right )
    return;
  for( int i = left + 1 ; i <= right ; i++ ){
    int val = array[i];
    int mirror_val = mirror_array[i];
    int final_pos = i;
    while( final_pos > left && val < array[final_pos-1] ){
      array[final_pos] = array[final_pos-1];
      mirror_array[final_pos] = mirror_array[final_pos-1];
      final_pos--;
    }
    array[final_pos] = val;
    mirror_array[final_pos] = mirror_val;
  }
}


static int top[32];
static int bottom[32];

void quick_sort(int* const array, const int left, const int right, int* const mirror_array)
{
  if( right - left <= 4 ){
    insertion_sort(array,left,right,mirror_array);
    return;
  }
  
  int sp = 0;
  top[sp] = right;
  bottom[sp] = left;
  int iter = 0;
  do{
    iter++;
    const int curr_left = bottom[sp];
    const int curr_right = top[sp];
    const int pivot = (curr_left+curr_right)/2;
    sp--;
    int temp;
    temp = array[curr_right];
    array[curr_right] = array[pivot];
    array[pivot] = temp;
    temp = mirror_array[curr_right];
    mirror_array[curr_right] = mirror_array[pivot];
    mirror_array[pivot] = temp;

    int final_pos = curr_left;
    for( int k = curr_left ; k < curr_right ; k++ )
      if( array[k] < array[curr_right] ){
        temp = array[final_pos];
        array[final_pos] = array[k];
        array[k] = temp;
        temp = mirror_array[final_pos];
        mirror_array[final_pos] = mirror_array[k];
        mirror_array[k] = temp;
        final_pos++;
      }
    temp = array[curr_right];
    array[curr_right] = array[final_pos];
    array[final_pos] = temp;
    temp = mirror_array[curr_right];
    mirror_array[curr_right] = mirror_array[final_pos];
    mirror_array[final_pos] = temp;
    
    if( final_pos -1 > curr_left ){
      sp++;
      top[sp] = final_pos - 1;
      bottom[sp] = curr_left;
    }
    if( curr_right > final_pos+1 ){
      sp++;
      top[sp] = curr_right;
      bottom[sp] = final_pos+1;
    }
  }while( sp >= 0 );

}



int remove_duplicates(int* const array, const int size, int* const mirror_array)
{
  if( size == 0 )
    return 0;

  int actual_posn = 0;
  for( int i = 1 ; i < size ; i++ )
    if( array[actual_posn] != array[i] ){
      array[++actual_posn] = array[i];
      mirror_array[actual_posn] = mirror_array[i];
    }
    else
      mirror_array[actual_posn] += mirror_array[i];

  return actual_posn+1;
}


void inspector::MetisReplicateHypergraph()
{
	if (team_size < 2)
		return;

	const int num_nets = data_num_offset[all_data.size()];

	//First Step , broadcast the number of pins a process has for each net
  
	int sendcount[team_size],senddispl[team_size + 1],recvcount[team_size],recvdispl[team_size + 1];

	int num_local_nets = 0;
	for( int i = 0 ; i < team_size ; i++) {
		sendcount[i] = 0;
		recvcount[i] = 0;
	}
    
	for( int i = 0 ; i < all_data.size() ; i++ ){
		int curr_size = data_num_offset[i+1] - data_num_offset[i];
		int split = curr_size / team_size;
		num_local_nets += ( proc_id == team_size - 1 ? curr_size - split * proc_id : split );
		for( int j = 0 ; j < team_size ; j++ ){
			sendcount[j] += (j == proc_id ? 0 : ( j == team_size - 1 ? curr_size - split * j : split ) );
		}
	}
    
	senddispl[0] = 0; recvdispl[0] = 0;
	for( int i =0 ; i < team_size ; i++ ){
		senddispl[i+1] = senddispl[i] + sendcount[i];
		recvcount[i] = ( i == proc_id ? 0 : num_local_nets ) ;
		recvdispl[i+1] = recvdispl[i] + recvcount[i];
	}

	int* const net_send_npins = new int[senddispl[team_size]];
	int* const net_recv_npins = new int[recvdispl[team_size]];
	net** const local_nets = new net*[num_local_nets];
  
	int curr_displ[team_size];
	for( int i= 0; i < team_size ; i++ )
		curr_displ[i] = senddispl[i];

	int local_count = 0;
	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ ){
		// if( !(*it)->is_read_only ){
		int curr_array_size = (*it)->orig_array_size ;
		int curr_split = curr_array_size / team_size;
		int my_start = curr_split * proc_id;
		int my_end = ( proc_id == team_size - 1 ? curr_array_size : my_start + curr_split );
		for( int j = 0,curr_proc = 0 ; j < curr_array_size ; j++){
			if( j < my_start || j >= my_end )
				net_send_npins[curr_displ[curr_proc]++] = (*it)->data_net_info[j]->pins.size()*2+1;
			else
				local_nets[local_count++] = (*it)->data_net_info[j];
			if( (j+1)%curr_split == 0 )
				curr_proc = ( curr_proc >= team_size - 1 ? team_size - 1 : curr_proc+1 );
		}
	}
#ifndef NDEBUG
	for( int i = 0; i < team_size ; i++ )
		assert(curr_displ[i] == senddispl[i+1]);
	assert(local_count == num_local_nets);
#endif
  
	MPI_Alltoallv(net_send_npins,sendcount,senddispl,MPI_INT,net_recv_npins,recvcount,recvdispl,MPI_INT,global_comm::global_iec_communicator);
	int send_pins_count[team_size],send_pins_displ[team_size + 1],recv_pins_count[team_size],recv_pins_displ[team_size + 1];
  
	recv_pins_displ[0] = 0; send_pins_displ[0] = 0;
	for( int i =0 ; i < team_size ; i++ ){
		send_pins_count[i] = 0;
		for( int j = senddispl[i] ; j < senddispl[i+1] ; j++ )
			send_pins_count[i] += net_send_npins[j];
		send_pins_displ[i+1] = send_pins_displ[i] + send_pins_count[i] ;
		recv_pins_count[i] = 0;
		for( int j = recvdispl[i] ; j < recvdispl[i+1] ; j++ )
			recv_pins_count[i] += net_recv_npins[j];
		recv_pins_displ[i+1] = recv_pins_displ[i] + recv_pins_count[i] ;
	}

	//Second Step : Send the pins of net to the respective process.

	int* const send_pins = new int[send_pins_displ[team_size]];
	int* const recv_pins = new int[recv_pins_displ[team_size]];

	for( int i = 0 ; i < team_size ; i++ )
		curr_displ[i] = send_pins_displ[i];

	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ ){
		int curr_array_size = (*it)->orig_array_size ;
		int curr_split = curr_array_size / team_size;
		int my_start = curr_split * proc_id;
		int my_end = ( proc_id == team_size - 1 ? curr_array_size : my_start + curr_split );
		for( int j = 0 , curr_proc = 0 ; j < curr_array_size ; j++ ){
			if( j < my_start || j >= my_end ){
				net* curr_net = (*it)->data_net_info[j];
				for( set<pin_info,pin_comparator>::iterator jt = curr_net->pins.begin() ; jt != curr_net->pins.end() ; jt++ ){
					send_pins[curr_displ[curr_proc]++] = (*jt).pin->my_num;
					send_pins[curr_displ[curr_proc]++] = ( (*jt).is_direct ? 1 : 0 );
				}
				if( curr_net->direct_vertex ){
					send_pins[curr_displ[curr_proc]++] = curr_net->direct_vertex->my_num;
				}
				else{
					send_pins[curr_displ[curr_proc]++]  = -1;
				}
			}
			if( (j+1)%curr_split == 0 )
				curr_proc = ( curr_proc >= team_size - 1 ? team_size - 1 : curr_proc + 1);
		}
	}
#ifndef NDEBUG
	for( int i = 0; i < team_size ; i++ )
		assert(curr_displ[i] == send_pins_displ[i+1]);
#endif

	MPI_Alltoallv(send_pins,send_pins_count,send_pins_displ,MPI_INT,recv_pins,recv_pins_count,recv_pins_displ,MPI_INT,global_comm::global_iec_communicator);
    
	int vcount[team_size];
	for( int i =0 ; i < team_size ; i++ )
		vcount[i] = recv_pins_displ[i];
	for( int i = 0 ; i < num_local_nets ; i++ ){
		net* curr_net = local_nets[i];
		for( int j = 0 ; j < team_size ; j++ )
			if( j != proc_id ){
				int num_pins = net_recv_npins[recvdispl[j]+i]-1;
				for( int l = 0 ; l < num_pins ; l+=2 ){
					int vertex_num = recv_pins[vcount[j]++];
					int k = 0 ;
					for( k = 0 ; k < all_loops.size() ; k++ )
						if( iter_num_offset[k+1] > vertex_num )
							break;
					assert(k != all_loops.size());
					assert( vertex_num - iter_num_offset[k] >= 0 && vertex_num - iter_num_offset[k] < all_loops[k]->num_iters);
					pin_info new_pin(all_loops[k]->iter_vertex[vertex_num - iter_num_offset[k]],(recv_pins[vcount[j]++] != 0 ? true : false));
					curr_net->pins.insert(new_pin);
				}
				int direct_vertex_num = recv_pins[vcount[j]++];
				if(direct_vertex_num != -1 ){
					int k = 0 ;
					for( k = 0 ; k < all_loops.size() ; k++ )
						if( iter_num_offset[k+1] > direct_vertex_num )
							break;
					assert(k != all_loops.size());
					assert( direct_vertex_num - iter_num_offset[k] >= 0 && direct_vertex_num - iter_num_offset[k] < all_loops[k]->num_iters);
					if( curr_net->direct_vertex == NULL )
						curr_net->direct_vertex = all_loops[k]->iter_vertex[direct_vertex_num-iter_num_offset[k]];
					else
						assert(curr_net->direct_vertex == all_loops[k]->iter_vertex[direct_vertex_num-iter_num_offset[k]]);
				}
			}
	}

	delete[] send_pins;
	delete[] recv_pins;
	delete[] local_nets;
	delete[] net_recv_npins;
	delete[] net_send_npins;
}


void inspector::MetisPrePartition()
{
  MetisReplicateHypergraph();

  //Go through the local nets to generate information about the adjacency list of the graph

  for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
    if( !(*it)->is_read_only ){
      int curr_array_size = (*it)->orig_array_size;
      int curr_split = curr_array_size / nprocs ;
      int curr_start = curr_split * proc_id;
      int curr_end = ( proc_id == nprocs - 1 ? curr_array_size : curr_split * ( proc_id + 1 ) );
      for( int j = curr_start ; j < curr_end ; j++ ){
	net* curr_net = (*it)->data_net_info[j];
	int npins = curr_net->pins.size();
	vertex* curr_pins[npins];
	int counter = 0;
	for( set<pin_info,pin_comparator>::iterator jt = curr_net->pins.begin() ; jt != curr_net->pins.end() ; jt++ )
	  curr_pins[counter++] = (*jt).pin;
	if( npins > 1 ){
	  for( int k = 0 ; k < npins ; k++ )
	    for( int l = k+1 ; l < npins ; l++ ){
	      vertex* vertex1 = curr_pins[k];
	      vertex* vertex2 = curr_pins[l];
#if defined USE_QPHASH || defined USE_KSHASH
	      vertex1->adjvertex->hash_insert(vertex2->my_num,curr_net->weight);
	      vertex2->adjvertex->hash_insert(vertex1->my_num,curr_net->weight);
#else
	      map<int,int>::iterator kt = vertex1->adjvertex.find(vertex2->my_num);
	      if( kt == vertex1->adjvertex.end() ){
		assert( vertex2->adjvertex.find(vertex1->my_num) == vertex2->adjvertex.end() );
		vertex1->adjvertex.insert(pair<int,int>(vertex2->my_num,curr_net->weight));
		vertex2->adjvertex.insert(pair<int,int>(vertex1->my_num,curr_net->weight));
	      }
	      else{
		(*kt).second += curr_net->weight;
		map<int,int>::iterator lt = vertex2->adjvertex.find(vertex1->my_num);
		assert( lt != vertex2->adjvertex.end());
		(*lt).second += curr_net->weight;
	      }
#endif
	    }
	}
      }
    }
}


void inspector::BlockPartition()
{
  MetisReplicateHypergraph();

  int nparts = team_size /** nthreads*/;
  
  for( deque<global_loop*>::iterator it = all_loops.begin() ; it!= all_loops.end() ; it++ ){
    int numiters = (*it)->num_iters;
    int split = numiters / nparts;
    int home = 0;
    for( int j = 0 ; j < numiters ; j++ ){
      vertex* curr_vertex = (*it)->iter_vertex[j];
      if( home != (nparts -1) && j == (home+1)*split )
	home++;
      curr_vertex->home = home;
      if( home/* / nthreads*/ == proc_id )
	(*it)->nproc_local++;
    }
  }
  MetisAfterPartition();
}


void inspector::MetisPartition()
{
  MetisPrePartition();

  const int num_vertex = iter_num_offset[all_loops.size()];
  const int vertex_split = num_vertex / nprocs;
  const int vertex_start = proc_id * vertex_split;
  const int vertex_end = ( proc_id == nprocs - 1 ? num_vertex - 1 : vertex_split*(proc_id+1) - 1 );
  const int num_local_vertex = vertex_end - vertex_start +1;
  
  vertex** metis_vertex = new vertex*[num_local_vertex];
  int counter = 0;
  if ( nprocs > 1 ){  
    int recvcount[nprocs],recvdispl[nprocs+1];
    int sendcount[nprocs],senddispl[nprocs+1];
    recvdispl[0] = 0 ; senddispl[0] = 0;
    for( int i =0 ; i < nprocs ; i++ ){
      sendcount[i] = ( i == proc_id ? 0 : ( i != nprocs - 1 ? vertex_split : num_vertex - i*vertex_split ) );
      recvcount[i] = ( i == proc_id ? 0 : num_local_vertex );
      senddispl[i+1] = senddispl[i] + sendcount[i];
      recvdispl[i+1] = recvdispl[i] + recvcount[i];
    }
    assert(senddispl[nprocs] == num_vertex - num_local_vertex );
    assert(recvdispl[nprocs] == num_local_vertex * (nprocs-1) );
    
    int * send_adj_count = new int[senddispl[nprocs]];
    int * recv_adj_count = new int[recvdispl[nprocs]];
  
    counter = 0;int countv = 0;
    for( deque<global_loop*>::iterator it = all_loops.begin() ; it!= all_loops.end() ; it++ ){
      for( int j = 0 ; j < (*it)->num_iters ; j++ ){
	vertex* curr_vertex = (*it)->iter_vertex[j];
	if( counter >= vertex_start && counter <= vertex_end )
	  metis_vertex[counter - vertex_start] = curr_vertex;
	else{
#if defined USE_QPHASH || defined USE_KSHASH
	  send_adj_count[countv++] = curr_vertex->adjvertex->num_elems();
#else
	  send_adj_count[countv++] = curr_vertex->adjvertex.size();
#endif
	}
	counter++;
      }
    }
  
    MPI_Alltoallv(send_adj_count,sendcount,senddispl,MPI_INT,recv_adj_count,recvcount,recvdispl,MPI_INT,MPI_COMM_WORLD);


    int sendcount_ja[nprocs],recvcount_ja[nprocs],senddispl_ja[nprocs+1],recvdispl_ja[nprocs+1];

    senddispl_ja[0] = 0 ; recvdispl_ja[0] =0;
    for( int i = 0 ; i < nprocs ; i++ ){
      sendcount_ja[i] = 0;
      for( int j = senddispl[i] ; j < senddispl[i+1] ; j++ )
	sendcount_ja[i] += send_adj_count[j] * 2;
      senddispl_ja[i+1] = senddispl_ja[i] + sendcount_ja[i];
      recvcount_ja[i] = 0;
      for( int j = recvdispl[i] ; j < recvdispl[i+1] ; j++ )
	recvcount_ja[i] += recv_adj_count[j] * 2;
      recvdispl_ja[i+1] = recvdispl_ja[i] + recvcount_ja[i];
    }

    int* send_adj_list = new int[senddispl_ja[nprocs]];
    int* recv_adj_list = new int[recvdispl_ja[nprocs]];

    counter =0;countv=0;
    for( deque<global_loop*>::iterator it = all_loops.begin() ; it!= all_loops.end() ; it++ )
      for( int j = 0 ; j < (*it)->num_iters ; j++ ){
	if( countv < vertex_start || countv > vertex_end ){
	  vertex* curr_vertex = (*it)->iter_vertex[j];
#if defined USE_QPHASH || defined USE_KSHASH
	  curr_vertex->adjvertex->return_keys_vals(send_adj_list+counter);
	  counter += curr_vertex->adjvertex->num_elems()*2;
	  curr_vertex->clear_adjvertex();
#else
	  for( map<int,int>::iterator jt = curr_vertex->adjvertex.begin() ; jt != curr_vertex->adjvertex.end() ; jt++ ){
	    send_adj_list[counter++] = (*jt).first;
	    send_adj_list[counter++] = (*jt).second;
	  }
#endif
	}
	countv++;
      }
    assert(counter == senddispl_ja[nprocs]);
  
    MPI_Alltoallv(send_adj_list,sendcount_ja,senddispl_ja,MPI_INT,recv_adj_list,recvcount_ja,recvdispl_ja,MPI_INT,MPI_COMM_WORLD);
  
    for( int i = 0 ; i < nprocs ; i++ )
      if( i != proc_id ){
	counter = recvdispl_ja[i];
	for( int j = 0 ; j < num_local_vertex ; j++ ){
	  int nadj = recv_adj_count[j+recvdispl[i]];
	  vertex* curr_vertex = metis_vertex[j];
	  for( int k = 0 ; k < nadj ; k++ ){
	    int vertex_num = recv_adj_list[counter++];
	    int edge_weight = recv_adj_list[counter++];
#if defined USE_QPHASH || defined USE_KSHASH
	    curr_vertex->adjvertex->hash_insert(vertex_num,edge_weight);
#else
	    map<int,int>::iterator kt = curr_vertex->adjvertex.find(vertex_num);
	    if( kt != curr_vertex->adjvertex.end() )
	      (*kt).second += edge_weight;
	    else
	      curr_vertex->adjvertex.insert(pair<int,int>(vertex_num,edge_weight));
#endif
	  }
	}
	assert(counter == recvdispl_ja[i+1]);
      }

    delete[] send_adj_count;
    delete[] recv_adj_count;
    delete[] send_adj_list;
    delete[] recv_adj_list;
  }
  else{
    counter = 0;
    for( deque<global_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ )
      for( int  j =0 ; j < (*it)->num_iters ; j++ )
	metis_vertex[counter++] = (*it)->iter_vertex[j];
  }


  int *xadj = new int[num_local_vertex+1];
  int *vwgt = new int[num_local_vertex*all_loops.size()];
  int *vtxdist = new int[nprocs+1];
  int ncon = all_loops.size();
  int nparts = nprocs/**nthreads*/;
  int wgtflag = 3;
  int numflag = 0;

  vtxdist[0] = 0;
  for( int i = 0 ; i < nprocs - 1 ; i++ )
    vtxdist[i+1] = vtxdist[i] + vertex_split;
  vtxdist[nprocs] = num_vertex;

  xadj[0] = 0;
  for( int i =  0 ; i < num_local_vertex ; i++ ){
    vertex* curr_vertex = metis_vertex[i];
#if defined USE_QPHASH || defined USE_KSHASH
    xadj[i+1] = xadj[i] + curr_vertex->adjvertex->num_elems();
#else
    xadj[i+1] = xadj[i] + curr_vertex->adjvertex.size();
#endif
    for( int k = 0 ; k < ncon ; k++ )
      vwgt[i*ncon+k] = 0;
    vwgt[i*ncon+curr_vertex->iter_num] = 1;
  }


  int* adjncy = new int[xadj[num_local_vertex]];
  int* adjwgt = new int[xadj[num_local_vertex]];

  counter = 0;
  for(int i = 0 ; i < num_local_vertex ; i++ ){
#if defined USE_QPHASH || defined USE_KSHASH
    metis_vertex[i]->adjvertex->return_keys_vals(adjncy+counter,adjwgt+counter);
    counter += metis_vertex[i]->adjvertex->num_elems();
    metis_vertex[i]->clear_adjvertex();
#else
    for( map<int,int>::iterator kt = metis_vertex[i]->adjvertex.begin() ; kt != metis_vertex[i]->adjvertex.end() ; kt++ ){
      adjncy[counter] = (*kt).first;
      adjwgt[counter] = (*kt).second;
      counter++;
    }
#endif
  }


  float*  tpwgts = new float[ncon*nparts];
  
  for( int i = 0 ; i < nparts ; i++ )
    for( int j = 0 ; j < ncon ; j++ )
      tpwgts[i*ncon+j] = 1.0 / nparts;

  float* ubvec = new float[ncon];
  for( int i =0 ; i < ncon ; i++ )
    ubvec[i] = 1.05;

  int options[3];
  options[0] = 1;
  options[1] = 0;
  options[2] = 42;


  int nedgecuts;
  int *vertex_home = new int[num_local_vertex];
  MPI_Comm world_comm = MPI_COMM_WORLD;
  
  ParMETIS_V3_PartKway(vtxdist,xadj,adjncy,vwgt,adjwgt,&wgtflag,&numflag,&ncon,&nparts,tpwgts,ubvec,options,&nedgecuts,vertex_home,&world_comm);
  
#ifndef NDEBUG
  printf("ID=%d, Done graphpartition\n",proc_id);
#endif

  if( nprocs > 1 ) {
    int recvcount[nprocs],recvdispl[nprocs+1];
    
    recvdispl[0] = 0;
    for( int i = 0; i < nprocs -1 ; i++ ){
      recvcount[i] = vertex_split;
      recvdispl[i+1] = recvdispl[i] + recvcount[i];
    }
    recvdispl[nprocs] = num_vertex;
    recvcount[nprocs-1] = recvdispl[nprocs] - recvdispl[nprocs-1];
    assert( recvcount[proc_id] == num_local_vertex);
    
    int* recvbuffer = new int[recvdispl[nprocs]];

    MPI_Allgatherv(vertex_home,num_local_vertex,MPI_INT,recvbuffer,recvcount,recvdispl,MPI_INT,MPI_COMM_WORLD);
    
    counter = 0;
    for( deque<global_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ )
      for( int j = 0; j < (*it)->num_iters ; j++ ){
	      if( recvbuffer[counter] /*/ nthreads*/ == proc_id )
	  (*it)->nproc_local++;
	(*it)->iter_vertex[j]->home = recvbuffer[counter++];
      }
    delete[] recvbuffer;
  }
  else{
    counter = 0;
    for( deque<global_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ ){
      (*it)->nproc_local = (*it)->num_iters;
      for( int j = 0; j < (*it)->num_iters ; j++ )
	(*it)->iter_vertex[j]->home = vertex_home[counter++];
    }
  }
    

#ifndef NDEBUG  
  printf("ID=%d,PrintingGraphFile\n",proc_id);
  char graph_file_name[14];
  sprintf(graph_file_name,"graph_%d.dat",proc_id);
  FILE* graph_file = fopen(graph_file_name,"w");
  fprintf(graph_file,"vwgt(%d): ",num_local_vertex);
  for( int i = 0; i < num_local_vertex ; i++ )
    for( int j = 0 ; j < ncon ; j++ )
      fprintf(graph_file," %d",vwgt[i*ncon+j]);
  fprintf(graph_file,"\n");
  fprintf(graph_file,"xadj(%d): ",num_local_vertex+1);
  for( int i = 0; i < num_local_vertex + 1 ; i++ )
    fprintf(graph_file," %d",xadj[i]);
  fprintf(graph_file,"\n");
  fprintf(graph_file,"adjncy(%d): ",xadj[num_local_vertex]);
  for( int i = 0; i < xadj[num_local_vertex] ; i++ )
    fprintf(graph_file," %d",adjncy[i]);
  fprintf(graph_file,"\n");
  fprintf(graph_file,"adjwgt(%d): ",xadj[num_local_vertex]);
  for( int i = 0; i < xadj[num_local_vertex] ; i++ )
    fprintf(graph_file," %d",adjwgt[i]);
  fprintf(graph_file,"\n");
  fprintf(graph_file,"vtxdist(%d): ",nprocs+1);
  for( int i = 0; i < nprocs +1; i++ )
    fprintf(graph_file," %d",vtxdist[i]);
  fprintf(graph_file,"\n");
  fprintf(graph_file,"home(%d): ",num_local_vertex);
  for( int i = 0; i < num_local_vertex; i++ )
    fprintf(graph_file," %d",vertex_home[i]);
  fprintf(graph_file,"\n");
  fclose(graph_file);
#endif
  
  delete[] xadj;
  delete[] vwgt;
  delete[] vtxdist;
  delete[] adjncy;
  delete[] adjwgt;
  delete[] tpwgts;
  delete[] ubvec;
  delete[] vertex_home;
  delete[] metis_vertex;

  MetisAfterPartition();
}

void inspector::MetisAfterPartition()
{
	int num_local_nets = 0;
	int* const recvcount= new int[nprocs];
	for( int i = 0 ; i < nprocs ; i++ ){
		recvcount[i] = 0;
	}
  
	for( int i = 0 ; i < all_data.size() ; i++ )
		if( !all_data[i]->is_read_only ){
			int curr_array_size = all_data[i]->orig_array_size;
			int curr_split = curr_array_size / nprocs;
			num_local_nets += ( proc_id == nprocs - 1 ? curr_array_size - curr_split * proc_id : curr_split );
			for( int j = 0 ; j < nprocs -1 ; j++ )
				recvcount[j] += curr_split;
			recvcount[nprocs-1] += curr_array_size - curr_split * (nprocs - 1) ;
		}
	int* const send_home = new int[num_local_nets];
  
	int* const recvdispl = new int[nprocs+1];
	recvdispl[0]=0;
	for( int i = 0 ; i < nprocs ; i++ )
		recvdispl[i+1] = recvdispl[i] + recvcount[i];
  
	int* const recv_home = new int[recvdispl[nprocs]];
	for( int i =0 ; i < recvdispl[nprocs] ; i++ )
		recv_home[i] = -1;
  
	int* const possible = new int[nprocs/**nthreads*/];
	int countv = 0;
	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )

		if( !(*it)->is_read_only ){

			int curr_array_size = (*it)->orig_array_size;
			int curr_split = curr_array_size / nprocs ; 
			int curr_start = curr_split * proc_id ; 
			int curr_end = ( proc_id == nprocs - 1 ? curr_array_size : curr_start + curr_split );
			if( (*it)->is_constrained ){
				int thread_split = curr_array_size / (nprocs/* * nthreads*/);
				for( int j = curr_start ; j < curr_end ; j++ ){
					net* curr_net = (*it)->data_net_info[j];
					int home = (j / thread_split > nprocs/**nthreads*/ - 1 ? nprocs/**nthreads*/ - 1 : j / thread_split ) ;
					curr_net->home = home;
					send_home[countv++] = curr_net->home;
				}	
			}

			else{

				for( int j = curr_start ; j < curr_end ; j++ ){

					net* curr_net = (*it)->data_net_info[j];
					int home = -1;
					if( curr_net->direct_vertex ){
						home = curr_net->direct_vertex->home;
					}

					else{

						for( int i = 0 ; i < nprocs/**nthreads*/ ; i++ )
							possible[i] = 0;
						for( set<pin_info,pin_comparator>::iterator jt = curr_net->pins.begin() ; jt != curr_net->pins.end() ; jt++ ){
							assert((*jt).pin->home >= 0 && (*jt).pin->home < (nprocs/**nthreads*/));
							possible[(*jt).pin->home]++;
						}

						int maxval = -1;
						int counter = 0 , i = 0 ;
						while( counter < nprocs/* * nthreads */){
							if( possible[i] > maxval ) {
								maxval = possible[i];
								home = i;
							}
							counter++;
							i = (i+1)%(nprocs/**nthreads*/);
						}
					}
					assert( curr_net->home == -1 && home >= 0 && home <= (nprocs/* * nthreads*/));
					curr_net->home = home; 
	
					assert(curr_net->home >= 0 && curr_net->home < (nprocs/* * nthreads*/) );
					assert(countv < num_local_nets);
					send_home[countv++] = curr_net->home;
				}
			}
		}
	assert(countv == num_local_nets);

	MPI_Allgatherv(send_home,num_local_nets,MPI_INT,recv_home,recvcount,recvdispl,MPI_INT,global_comm::global_iec_communicator);

	int curr_displ[nprocs];
	for( int i = 0; i < nprocs ; i++ )
		curr_displ[i] = recvdispl[i];

	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
		if( !(*it)->is_read_only ){
			int curr_array_size = (*it)->orig_array_size;
			int curr_split = curr_array_size / nprocs;
			int my_start = curr_split*proc_id;
			int my_end = ( proc_id == nprocs - 1 ? curr_array_size : my_start + curr_split );
			for( int j = 0, curr_proc = 0 ; j < curr_array_size ; j++ ){
				if( j < my_start || j >= my_end ){
					net* curr_net = (*it)->data_net_info[j];
					assert(curr_net->home == -1 );
					assert(recv_home[curr_displ[curr_proc]] >= 0 && recv_home[curr_displ[curr_proc]] < nprocs/**nthreads*/ );
					curr_net->home = recv_home[curr_displ[curr_proc]++];
				}
				if( (j+1)%curr_split == 0 )
					curr_proc = ( curr_proc >= nprocs - 1 ? nprocs - 1 : curr_proc + 1 ) ;
			}
		}

	delete[] possible;
	delete[] recvcount;
	delete[] recvdispl;
	delete[] send_home;
	delete[] recv_home;

	AfterPartition();
}

