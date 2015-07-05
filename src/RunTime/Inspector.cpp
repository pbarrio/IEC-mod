/*
 * inspector.cpp: This file is part of the IEC project.
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
 * @file: Inspector.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/Inspector.hpp"
#include "armci.h"
#include "mpi.h"
#include "omp.h"
#include "patoh.h"

using namespace std;

Inspector* Inspector::singleton_inspector = NULL;
char** global_comm::put_buffer = NULL;
char* global_comm::send_buffer = NULL;
char* global_comm::recv_buffer = NULL;
int global_comm::max_send_size = 0;
int global_comm::max_recv_size = 0;
int global_comm::max_nprocs_read_send = 0;
int global_comm::max_nprocs_read_recv = 0;
int global_comm::max_nprocs_write_send = 0;
int global_comm::max_nprocs_write_recv = 0;
MPI_Comm global_comm::global_iec_communicator = NULL;
MPI_Comm global_comm::team_communicator = NULL;
MPI_Status* global_comm::read_send_start_status = NULL;
MPI_Status* global_comm::read_recv_start_status = NULL;
MPI_Request* global_comm::read_recv_start_request = NULL;
MPI_Status* global_comm::read_send_end_status = NULL;
MPI_Status* global_comm::read_recv_end_status = NULL;
MPI_Request* global_comm::read_send_end_request = NULL;
MPI_Status* global_comm::write_send_start_status = NULL;
MPI_Status* global_comm::write_recv_start_status = NULL;
MPI_Request* global_comm::write_recv_start_request = NULL;
MPI_Status* global_comm::write_send_end_status = NULL;
MPI_Status* global_comm::write_recv_end_status = NULL;
MPI_Request* global_comm::write_send_end_request = NULL;
char* global_comm::read_send_signal = NULL;
char* global_comm::read_recv_signal = NULL;
char* global_comm::write_send_signal = NULL;
char* global_comm::write_recv_signal = NULL;


/**
 * \brief Sets up singleton inspector
 *
 * \param pid Process identifier in the distributed program
 * \param np Number of processes in the distributed program
 * \param team Identifier of this process' team
 * \param pid_team Identifier of this processor private to the team
 * \param teamsize Size (in number of processors) of each team
 * \param nl Number of loops
 * \param nd Number of direct access arrays
 * \param nc Number of communicators. There is one per loop.
 * \param nad Number of indirection arrays.
 * \param iter_num_count Array of iteration limits (must be of size nl)
 * \param data_num_count Array of direct array sizes (must be of size nd)
 * \param ro Array of read-only flags for the direct arrays (must be of size nd)
 */
Inspector::Inspector(int pid, int np, int team, int pid_team, int teamsize,
                     int nl, int nd, int nc, int nad, int* iter_num_count,
                     int* data_num_count, int* ro):
	proc_id(pid),
	team_num(team),
	id_in_team(pid_team),
	team_size(teamsize),
	nprocs(np),
	pins_size(-1),
	my_loop(NULL),
	iter_num_offset(new int[nl + 1]), // The first '1' is for my loop
	data_num_offset(new int[nd + 1]){

	// Create new communicator associated to the inspector. This avoids having
	// to use MPI_COMM_WORLD.
	// Also create a communicator for the team.
	MPI_Comm_dup(MPI_COMM_WORLD, &global_comm::global_iec_communicator);
	MPI_Comm_split(MPI_COMM_WORLD, team_num, 0,
	               &global_comm::team_communicator);
	MPI_Barrier(MPI_COMM_WORLD);

	if (pid >= np)
		return;

	/// Create all loops
	iter_num_offset[0] = 0;
	for (int iLoop = 0 ; iLoop < nl ; ++iLoop){

		global_loop* new_loop =
			new global_loop(iLoop,
			                iter_num_count[iLoop],
			                iter_num_offset[iLoop]);
		all_loops.push_back(new_loop);
		iter_num_offset[iLoop + 1] =
			iter_num_offset[iLoop] + iter_num_count[iLoop];
	}

	/// Create data arrays (excludes indirection arrays)
	data_num_offset[0] = 0;
	for (int i = 0; i < nd; i++){
		global_data* new_data;
		new_data = new global_data_double(i, data_num_count[i],
		                                  data_num_offset[i],
		                                  ro[i] == 1? true : false);
		all_data.push_back(new_data);

		data_num_offset[i+1] = data_num_offset[i] + data_num_count[i];
	}

	// Set up the global communicator for each loop
	// TODO: we should get rid of multiple global communicators (use only one).
	for( int i = 0 ; i < nc ; i++ ){
		global_comm* new_comm = new global_comm(i, nprocs, proc_id);
		all_comm.push_back(new_comm);
	}

	/// Indirection arrays (called access data).
	for( int i = 0 ; i < nad ; i++ ){
		access_data* new_access_data =
			new access_data(i, proc_id, nprocs, pid_team, teamsize);
		all_access_data.push_back(new_access_data);
	}

	send_info = new set<int>*[nprocs*2];
	for( int i = 0; i < nprocs*2 ; i++ )
		send_info[i] = new set<int>;

#ifndef NDEBUG
	char access_file_name[20];
	sprintf(access_file_name,"access_data_%d.dat",proc_id);
	access_file = fopen(access_file_name,"w");
#endif

	MPI_Barrier(global_comm::global_iec_communicator);
}


Inspector::~Inspector(){

	for (deque<global_loop*>::iterator it = all_loops.begin();
	     it != all_loops.end() ; it++)

		delete (*it);

	delete[] iter_num_offset;
	all_loops.clear();
	for (deque<global_data*>::iterator it = all_data.begin();
	     it != all_data.end(); it++ )

		delete (*it);

	delete[] data_num_offset;
	all_data.clear();
	int counter = 0;
	for (deque<global_comm*>::iterator it = all_comm.begin();
	     it != all_comm.end(); it++ , counter++){

		delete (*it);

	}
	all_comm.clear();
	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		delete (*it);

	all_access_data.clear();

	if( global_comm::max_send_size > 0 )
		ARMCI_Free_local(global_comm::send_buffer);
	if (global_comm::put_buffer){
		ARMCI_Free(global_comm::put_buffer[proc_id]);
		delete[] global_comm::put_buffer;
	}
	if( global_comm::max_nprocs_read_send > 0 ){
		delete[] global_comm::read_send_start_status;
		delete[] global_comm::read_send_end_status;
		delete[] global_comm::read_send_end_request;
		delete[] global_comm::read_send_signal;
	}
	if( global_comm::max_nprocs_read_recv > 0 ){
		delete[] global_comm::read_recv_start_status;
		delete[] global_comm::read_recv_end_status;
		delete[] global_comm::read_recv_start_request;
		delete[] global_comm::read_recv_signal;
	}
	if( global_comm::max_nprocs_write_send > 0 ){
		delete[] global_comm::write_send_start_status;
		delete[] global_comm::write_send_end_status;
		delete[] global_comm::write_send_end_request;
		delete[] global_comm::write_send_signal;
	}
	if( global_comm::max_nprocs_write_recv > 0 ){
		delete[] global_comm::write_recv_start_status;
		delete[] global_comm::write_recv_end_status;
		delete[] global_comm::write_recv_start_request;
		delete[] global_comm::write_recv_signal;
	}

	if (send_info){
		for( int i = 0 ; i < nprocs * 2 ; i++ ){
			send_info[i]->clear();
					delete send_info[i];
		}
		delete[] send_info;
	}

	// Deallocate MPI communicator common to all participating processes
	// in the inspector/executor.
	if (global_comm::team_communicator != MPI_COMM_NULL)
		MPI_Comm_free(&global_comm::team_communicator);
	if (global_comm::global_iec_communicator != MPI_COMM_NULL)
		MPI_Comm_free(&global_comm::global_iec_communicator);
}


void Inspector::GetDontHave(){

	int stride = all_access_data.size(); // Number of indirection arrays
	int* get_n_elems = new int[team_size * stride];
	int* send_n_elems = new int[team_size * stride];

	/// Get the number of elements to get from other processes
	for( int i = 0 ; i < team_size * stride ; i++ )
		get_n_elems[i] = 0;
	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->GetSendCounts(get_n_elems,stride);

	///Initial communication to setup the buffer
	MPI_Alltoall(get_n_elems, stride, MPI_INT, send_n_elems, stride, MPI_INT,
	             global_comm::global_iec_communicator);

	int send_count[team_size];
	int recv_count[team_size];
	int send_offset[team_size + 1];
	int recv_offset[team_size + 1];

	// Compute the size of message sent and received by each process
	send_offset[0] = 0; recv_offset[0] = 0;
	for (int i = 0; i < team_size; i++){

		send_count[i] = 0; recv_count[i] = 0;
		for (int j = 0; j < stride; j++){
			send_count[i] += get_n_elems[i * stride + j];
			recv_count[i] += send_n_elems[i * stride + j];
		}
		send_offset[i + 1] = send_offset[i] + send_count[i];
		recv_offset[i + 1] = recv_offset[i] + recv_count[i];
	}

	int *sendbuf, *recvbuf;

	// Allocate send buffer
	if (send_offset[team_size] > 0)
		sendbuf = new int[send_offset[team_size]];
	else
		sendbuf = NULL;

	// Allocate recv buffer
	if (recv_offset[team_size] > 0)
		recvbuf = new int[recv_offset[team_size]];
	else
		recvbuf = NULL;

	// Populate the send buffer with the indices requested by other processes
	int* curr_offset = new int[team_size];
	for (int i = 0; i < team_size; i++)
		curr_offset[i] = send_offset[i];
	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->PopulateBuffer(sendbuf, send_offset[team_size], curr_offset);

	MPI_Alltoallv(sendbuf, send_count, send_offset, MPI_INT, recvbuf,
	              recv_count, recv_offset, MPI_INT,
	              global_comm::global_iec_communicator);
  
	// Populate the receive buffer with the values of all indices requested
	// by this process.
	for (int i = 0; i < team_size; i++)
		curr_offset[i] = recv_offset[i];

	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->GetRequestedValue(recvbuf, recv_offset[team_size], curr_offset,
		                         send_n_elems, stride);
  
	MPI_Alltoallv(recvbuf, recv_count, recv_offset, MPI_INT, sendbuf,
	              send_count, send_offset, MPI_INT,
	              global_comm::global_iec_communicator);
  
	///Add to set of elements that are known on each process
	for (int i = 0; i < team_size; i++)
		curr_offset[i] = send_offset[i];

	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->AddToHaveSet(sendbuf,send_offset[team_size],curr_offset);

	delete[] get_n_elems;
	delete[] send_n_elems;
	if( sendbuf )
		delete[] sendbuf;
	if( recvbuf )
		delete[] recvbuf;
	delete[] curr_offset;
}


bool Inspector::DoneGraphGeneration(){
	int n_dont = 0;
	for( int i =0 ; i < all_access_data.size() ; i++ )
		n_dont += all_access_data[i]->dont_have_set.size();

	// Compute the number of indirect array positions missing from all processes
	MPI_Allreduce(MPI_IN_PLACE, &n_dont, 1, MPI_INT, MPI_SUM,
	              global_comm::global_iec_communicator);
#ifndef NDEBUG
	printf("MXD:ID:%d,n_dont=%d\n",proc_id,n_dont);
	fflush(stdout);
#endif

	// If there is still one array position missing,
	// the graph generation is not done.
	if (n_dont != 0){
		GetDontHave();
		return false;
	}
	else
		return true;
}


/**
 * \brief Initializes structures required by the Inspector for this loop
 *
 * Mark arrays as used
 */
void Inspector::init_loop(int loopID, vector<int> dataIDs){

	for (vector<int>::iterator dataID = dataIDs.begin(),
		     dataIDEnd = dataIDs.end();
	     dataID != dataIDEnd; dataID++){

		global_data* data = all_data[*dataID];
		data->use_in_loop(loopID);
	}
}


/**
 * \brief Adds a vertex to the iteration/data hypergraph
 *
 * This function also classifies the loop as a "producer", "consumer",
 * or executed loop. The executed loop is the one that will be executed by this
 * process, while other loops are executed by other processes and will either
 * supply data to this process (producer) or wait for data to be communicated by
 * this process (consumers).
 *
 * \param loop_id The vertex will be added to this loop.
 * \param iter_value The vertex represents this iteration number.
 */
void Inspector::AddVertex(int loop_id, int iter_value){

	assert(loop_id < all_loops.size());

	// Loop classification
	global_loop* curr_loop = all_loops[loop_id];
	if (loop_id == team_num)
		my_loop = curr_loop;
	else if (loop_id < team_num)
		producer_loops[loop_id] = curr_loop;
	else if (loop_id > team_num)
		consumer_loops[loop_id] = curr_loop;

	curr_vertex = curr_loop->iter_vertex[iter_value];
}


/**
 * \brief Adds a pin from the current vertex to a network representing some data
 *
 * \param data_num Identifier of the array to be referenced.
 * \param index Index in the array to be referenced.
 * \param loop The pin is meaningful in the context of this loop.
 * \param is_direct !=0 if addressing is affine; =0 if it depends on indirection
 *                  array.
 * \param is_ploop !=0 if the access is direct from the partitionable loops.
 */
void Inspector::AddPinToNet
(int data_num, int index, int loop, int is_direct, int is_ploop){

	assert(data_num < all_data.size());

	///Get the current net
	net* net_num = all_data[data_num]->data_net_info[loop][index];

	pin_info new_pin(curr_vertex, is_direct != 0 ? true : false);

	/// If access is from a partitionable loop and is direct access
	if (is_ploop && is_direct){
		// The compile time analysis should ensure that there are not other
		// "direct accesses" to this array. The direct_vertex field should be
		// NULL or the same vertex as the curr_vertex
		assert(net_num->direct_vertex == NULL ||
		       net_num->direct_vertex == curr_vertex);
		net_num->direct_vertex = curr_vertex;
	}

	///Check for the pin already existing in the set
	set<pin_info,pin_comparator>::iterator it = net_num->pins.find(new_pin);
	if( it != net_num->pins.end() ){
		// If it is not a direct access, and current access is direct,
		// promote it to direct access.
		if( is_direct && !(*it).is_direct ){
			net_num->pins.erase(*it);
			net_num->pins.insert(new_pin);
		}
	}
	else{
		net_num->pins.insert(new_pin);
	}
}

void Inspector::AddUnknown
(int solver_num, int array_num, int index, int row_num, int loop){

	assert(array_num < all_data.size() && solver_num < all_solvers.size());
	assert(index < all_data[array_num]->orig_array_size);
	net* curr_net = all_data[array_num]->data_net_info[index][loop];

	assert(all_solvers[solver_num]->orig_net[row_num] == NULL ||
	       all_solvers[solver_num]->orig_net[row_num] == curr_net);

	all_solvers[solver_num]->orig_net[row_num] = curr_net;
}

void Inspector::PatohPrePartition(int loop){

	if( nprocs > 1 ) {
		const int ia_size = data_num_offset[all_data.size()];

		int** armci_net_ia = new int*[nprocs];
		ARMCI_Malloc((void**)armci_net_ia, (ia_size + 1) * sizeof(int));
		int* const net_ia = (int*)armci_net_ia[proc_id];

		net_ia[0] = 0;
		int counter = 0;

		// Store the pin information for nets in CSR format, compute the ia
		for (int i = 0; i < all_data.size(); i++)
			for (int j = 0; j < all_data[i]->orig_array_size; j++){
				net_ia[counter + 1] = net_ia[counter] +
					all_data[i]->data_net_info[loop][j]->pins.size() * 2 + 1;
				counter++;
			}
		assert(counter == ia_size);

		int net_ja_size[nprocs];

		///Total number of pins on each process
		MPI_Allgather(&(net_ia[ia_size]), 1, MPI_INT, net_ja_size, 1, MPI_INT,
		              global_comm::global_iec_communicator);

		int max_net_ja_size = 0;
		///Maximum size needed for the pin information from each process
		for( int i =0 ; i < nprocs ; i++)
			max_net_ja_size = (net_ja_size[i] > max_net_ja_size ?
			                   net_ja_size[i] : max_net_ja_size);

		int** armci_net_ja = new int*[nprocs];
		ARMCI_Malloc((void**)armci_net_ja, max_net_ja_size * sizeof(int));
		int* const net_ja = (int*)armci_net_ja[proc_id];

		///Populate the pin information on the local process
		counter = 0;
		for (int i = 0; i < all_data.size(); i++)
			for (int j = 0; j < all_data[i]->orig_array_size; j++){
				net* curr_net = all_data[i]->data_net_info[loop][j];

				for (set<pin_info, pin_comparator>::iterator it =
					     curr_net->pins.begin();
				     it != curr_net->pins.end();
				     it++){

					net_ja[counter++] = (*it).pin->my_num;
					net_ja[counter++] = ((*it).is_direct ? 1 : 0);
				}

				if (curr_net->direct_vertex)
					net_ja[counter++] = curr_net->direct_vertex->my_num;
				else
					net_ja[counter++] = -1;
			}
		assert(counter == net_ia[ia_size]);

		// Local buffer to hold the pin information from different processes
		int* const recv_ia =
			(int*)ARMCI_Malloc_local((ia_size + 1) * sizeof(int));
		int* const recv_ja =
			(int*)ARMCI_Malloc_local(max_net_ja_size * sizeof(int));

		MPI_Barrier(global_comm::global_iec_communicator);

		for (int i = 1; i < nprocs; i++){

			// Get the pin information from each process
			int dest_id = (proc_id + i) % nprocs;
			int source_id = (proc_id + nprocs - i) % nprocs;

			ARMCI_Get((void*)armci_net_ia[source_id], (void*)recv_ia,
			          (ia_size + 1) * sizeof(int), source_id);
			ARMCI_Get((void*)armci_net_ja[source_id], (void*)recv_ja,
			          recv_ia[ia_size] * sizeof(int), source_id);

			counter = 0;
			for (deque<global_data*>::iterator of_data = all_data.begin();
			     of_data != all_data.end(); of_data++){

				global_data* curr_data = (*of_data);

				///Update the local information
				for( int i = 0 ; i < curr_data->orig_array_size ; i++ ){

					net* curr_net = curr_data->data_net_info[loop][i];
					for (int j = recv_ia[counter];
					     j < recv_ia[counter + 1] - 1;
					     j+=2){

						int vertex_num = recv_ja[j];
						int k = 0;
						for( k = 0 ; k < all_loops.size() ; k++ )
							if( iter_num_offset[k+1] > vertex_num )
								break;
						assert(k != all_loops.size());
						pin_info new_pin
							(all_loops[k]->iter_vertex[vertex_num -
							                           iter_num_offset[k]],
							 (recv_ja[j + 1] != 0 ? true : false));
						curr_net->pins.insert(new_pin);
					}
					if( recv_ja[recv_ia[counter+1]-1] != -1 ){
						int vertex_num = recv_ja[recv_ia[counter+1]-1];
						int k = 0;
						for( k = 0 ; k < all_loops.size() ; k++ )
							if( iter_num_offset[k+1] > vertex_num )
								break;
						assert(k != all_loops.size());
						curr_net->direct_vertex = all_loops[k]->
							iter_vertex[vertex_num - iter_num_offset[k]];
					}
					counter++;
				}
			}
		}

		ARMCI_Free(armci_net_ja[proc_id]);
		ARMCI_Free(armci_net_ia[proc_id]);
		ARMCI_Free_local(recv_ia);
		ARMCI_Free_local(recv_ja);
		delete[] armci_net_ia;
		delete[] armci_net_ja;

#ifndef NDEBUG
		printf("PID:%d,ReplicationDone\n",proc_id);
		fflush(stdout);
		MPI_Barrier(global_comm::global_iec_communicator);
#endif
	}
}


void Inspector::PatohPartitionAll(){

	for (int loop = 0; loop < all_loops.size(); ++loop)
		PatohPartition(loop);
}


void Inspector::PatohPartition(int loop){

	PatohPrePartition(loop);

	if( nprocs > 1 ){
		int patoh_n,patoh_c,patoh_nc,patoh_cutsize;
		int *patoh_cellwts, *patoh_netwts, *patoh_xpins, *patoh_pins,
			*patoh_partwts, *patoh_partvec;
		float* patoh_targetwts=NULL;
		int i,j,k;

		PaToH_Parameters patoh;

		PaToH_Initialize_Parameters(&patoh, PATOH_CONPART,
		                            PATOH_SUGPARAM_QUALITY);

		//Number of vertices
		patoh_c = iter_num_offset[all_loops.size()];
		//number of nets
		patoh_n = 0;
		for( i = 0 ; i < all_data.size() ; i++ )
			if( !all_data[i]->is_read_only )
				patoh_n += all_data[i]->orig_array_size;
		//Number of constraints
		patoh_nc = all_loops.size();
		//Number of partitions
		patoh._k = nprocs;

		//Weight for cells
		patoh_cellwts = new int[patoh_c*patoh_nc];
		int counter = 0;
		for( i = 0 ; i < all_loops.size() ; i++ )
			for( j = 0 ; j < all_loops[i]->num_iters ; j++ )  {
				for( k = 0 ; k < patoh_nc ; k++ )
					patoh_cellwts[counter*patoh_nc + k] = 0;
				patoh_cellwts[ counter*patoh_nc + i ] = 1;
				counter++;
			}

		//Weight and xpins for nets
		patoh_netwts = new int[patoh_n];
		patoh_xpins = new int[patoh_n+1];
		patoh_xpins[0] = 0;
		counter = 0;
		for( i = 0 ; i < all_data.size() ; i++ )
			if( !all_data[i]->is_read_only)
				for( j = 0 ; j < all_data[i]->orig_array_size ; j++ ) {
					patoh_netwts[counter] =
						all_data[i]->data_net_info[loop][j]->weight;
					patoh_xpins[counter + 1] = patoh_xpins[counter] +
						all_data[i]->data_net_info[loop][j]->pins.size();
					counter++;
				}

		//pins for nets;
		pins_size = patoh_xpins[patoh_n];
		patoh_pins = new int[patoh_xpins[patoh_n]];
		counter = 0;
		for( i = 0 ; i < all_data.size() ; i++ )
			if( !all_data[i]->is_read_only )
				for( j = 0 ; j < all_data[i]->orig_array_size ; j++ ) {
					set<pin_info, pin_comparator>::iterator jt =
						all_data[i]->data_net_info[loop][j]->pins.begin();
					for (;
					     jt != all_data[i]->data_net_info[loop][j]->pins.end();
					     jt++)
						patoh_pins[counter++] = (*jt).pin->my_num;
				}

		PaToH_Alloc(&patoh, patoh_c, patoh_n, patoh_nc, patoh_cellwts,
		            patoh_netwts, patoh_xpins, patoh_pins);

		patoh_partvec = new int[patoh_c];
		patoh_partwts = new int[patoh_nc*patoh._k];
		patoh.seed = 42;

		PaToH_Part(&patoh, patoh_c, patoh_n, patoh_nc, 0, patoh_cellwts,
		           patoh_netwts, patoh_xpins, patoh_pins, patoh_targetwts,
		           patoh_partvec, patoh_partwts, &patoh_cutsize);

		counter = 0;
		for( i = 0 ; i < all_loops.size() ; i++ )
			for( j = 0 ; j < all_loops[i]->num_iters ; j++ ){
				if( patoh_partvec[counter] == proc_id )
					all_loops[i]->nproc_local++;
				all_loops[i]->iter_vertex[j]->home = patoh_partvec[counter++];
			}

		delete[] patoh_cellwts;
		delete[] patoh_netwts;
		delete[] patoh_xpins;
		delete[] patoh_pins;
		delete[] patoh_partvec;
		delete[] patoh_partwts;

		PaToH_Free();
	}
	else{
		for( int i = 0 ; i < all_loops.size() ; i++ )
			for( int j = 0 ; j < all_loops[i]->num_iters ; j++ ){
				all_loops[i]->nproc_local++;
				all_loops[i]->iter_vertex[j]->home = 0;
			}
	}

	PatohAfterPartition(loop);
}

void Inspector::PatohAfterPartition(int loop){

	//Decide homes for the nets
	int possible[nprocs];

	for (deque<global_data*>::iterator it = all_data.begin();
	     it != all_data.end();
	     it++)

		if( !(*it)->is_read_only ){
			if( !(*it)->is_constrained ){
				for( int j = 0 ; j < (*it)->orig_array_size ; j++ ){
					net* curr_net = (*it)->data_net_info[loop][j];
					int home = -1;
					///If there is a direct access use the home of that vertex
					if( curr_net->direct_vertex ){
						home = curr_net->direct_vertex->home;
					}
					else{
						///Add it to the same process as the one that accesses
						///it the most (really doesn't matter)
						for (int i = 0 ; i < nprocs ; i++ )
							possible[i] = 0;
						for (set<pin_info,pin_comparator>::iterator jt =
							     curr_net->pins.begin();
						     jt != curr_net->pins.end();
						     jt++)

							possible[(*jt).pin->home]++;

						int maxval = -1;
						int counter = 0 , i = 0 ;
						while( counter < nprocs ){
							if( possible[i] > maxval ) {
								maxval = possible[i];
								home = i;
							}
							counter++;
							i = (i+1)%(nprocs);
						}
					}
					assert(curr_net->home == -1 && home >= 0 && home <= nprocs);
					curr_net->home = home;
				}
			}
			else{
				///If constrained, then the array is assumed to be
				///block partitioned (ghosts added later)
				int curr_proc = 0;
				const int array_split = (*it)->orig_array_size / (nprocs);
				for( int j = 0 ; j < (*it)->orig_array_size ; j++ ){
					net* curr_net = (*it)->data_net_info[loop][j];
					curr_net->home = curr_proc;
					if( (j+1)%array_split == 0 )
						curr_proc = (curr_proc + 1 > nprocs -1 ?
						             curr_proc : curr_proc + 1);
				}
			}
		}
	AfterPartition(loop);
}

/**
 * \brief Set up solvers and the local inspector.
 *
 * This function is intended to be called after the partitioner has finished.
 * Therefore, at this point we know the mapping of loop iterations to processes.
 * At the moment, I'm not very clear about what "solvers" are.
 * Previously, there was one local inspector per thread in the process.
 * Now, there is only one thread.
 */
void Inspector::AfterPartition(int loop){

	/// What are exactly the solvers?
	for (deque<petsc_solve*>::iterator it = all_solvers.begin();
	    it != all_solvers.end();
	    it++ ){

		//Replicate the net info
		petsc_solve* curr_solver = (*it);
		const int curr_size = (*it)->size;
		int* send_buf = new int[curr_solver->size];
		int* recv_buf = new int[curr_solver->size];

		int source_proc = proc_id -1;
		for (int dest_proc = (proc_id + 1) % nprocs;
		     dest_proc != proc_id;
		     dest_proc = (dest_proc + 1) % nprocs, source_proc--){

			MPI_Request i_send_request;
			MPI_Status i_recv_status,i_send_status;

			if( source_proc < 0 )
				source_proc = nprocs - 1;

			for( int j = 0 ; j < curr_size ; j++ )
				if( curr_solver->orig_net[j] != NULL )
					send_buf[j] = curr_solver->orig_net[j]->my_num;
				else
					send_buf[j] = -1;

			MPI_Isend(send_buf, curr_size, MPI_INT, dest_proc, 0,
			          global_comm::global_iec_communicator, &i_send_request);
			MPI_Recv(recv_buf, curr_size, MPI_INT, source_proc, 0,
			         global_comm::global_iec_communicator,&i_recv_status);

			for (int j = 0; j < curr_size; j++){
				int net_num = recv_buf[j];
				if (net_num != -1){
					int k = 0;
					for (; k < all_data.size(); k++)
						if (data_num_offset[k + 1] > net_num)
							break;
					if (k == all_data.size())
						printf("ID=%d,Source=%d,net_num=%d,k=%d\n", proc_id,
						       source_proc, net_num, k);
					assert(k != all_data.size());
					assert(!all_data[k]->is_read_only);
					assert(curr_solver->orig_net[j] == NULL ||
					       curr_solver->orig_net[j] ==
					       all_data[k]->data_net_info[loop][net_num -
					                                        data_num_offset[k]
					                                        ]);
					curr_solver->orig_net[j] =
						all_data[k]->data_net_info[loop][net_num -
						                                 data_num_offset[k]];
				}
			}
			MPI_Wait(&i_send_request,&i_send_status);
		}
		delete[] send_buf;
		delete[] recv_buf;
		(*it)->FindNewRowNumbers();
	}

	/// Take all global data and intialize the corresponding local data
	// It might be possible to get rid of this
	for (deque<global_data*>::iterator it = all_data.begin();
	     it != all_data.end(); it++){

		if ((*it)->data_net_info[loop] == NULL)
			continue;

		int array_id = (*it)->id;
		int stride_size = (*it)->stride_size;
		int orig_array_size = (*it)->orig_array_size;
		bool read_only_flag =(*it)->is_read_only;
		bool is_constrained = (*it)->is_constrained; // Irrelevant for "quake"
		const net** data_net_info =
			const_cast<const net**>((*it)->data_net_info[loop]);
		add_local_data(array_id, stride_size, data_net_info, orig_array_size,
		               read_only_flag, is_constrained);
	}
}


void Inspector::GetLocalAccesses(int array_num, int** recvbuf, int** displ,
                                 int** count){

	const global_data* curr_array = all_data[array_num];
	net** const curr_nets = curr_array->data_net_info.at(team_num);
	const int curr_array_size = curr_array->orig_array_size;
	const int curr_split = curr_array_size / team_size;
	const int curr_start = curr_split * id_in_team;
	const int curr_end = (id_in_team == team_size - 1 ?
	                      curr_array_size :
	                      curr_split * (id_in_team + 1));

	for (int i = 0; i < team_size * 2; i++)
		send_info[i]->clear();
	int* const sendcount_mpi = new int[team_size];
	bool* is_direct = new bool[team_size];
	bool* flags = new bool[team_size];

	// Find all the processes that access the blocked part of array owned
	// by the process
	for (int i = curr_start; i < curr_end; i++){

		const net* curr_net = curr_nets[i];
		for (int j = 0; j < team_size; j++){
			is_direct[j] = false;
			flags[j] = false;
		}

		for (set<pin_info, pin_comparator>::const_iterator it =
			     curr_net->pins.begin();
		     it != curr_net->pins.end();
		     it++)

			is_direct[(*it).pin->home] =
				is_direct[(*it).pin->home] || (*it).is_direct;

		for (set<pin_info, pin_comparator>::const_iterator it =
			     curr_net->pins.begin();
		     it != curr_net->pins.end();
		     it++ ){

			const int access_proc = (*it).pin->home;
			assert(access_proc != -1 && access_proc < team_size);
			if (!flags[access_proc]){
				if (is_direct[access_proc])
					send_info[access_proc * 2]->insert(curr_net->data_index);
				else
					send_info[access_proc * 2 + 1]->insert(curr_net->data_index);
				flags[access_proc] = true;
			}
		}
	}
	delete[] is_direct;
	delete[] flags;

	int* const sendcount = new int[team_size * 2];
	int* const curr_displ = new int[team_size * 2 + 1];
	int* const senddispl = new int[team_size + 1];
	senddispl[0] = 0;
	curr_displ[0] = 0;
	for (int i = 0; i < team_size; i++){
		sendcount_mpi[i] = 0;
		for (int k = 0; k < 2; k++){
			sendcount[i * 2 + k] = send_info[i * 2 + k]->size();
			curr_displ[i * 2 + k + 1] =
				curr_displ[i * 2 + k] + sendcount[i * 2 + k];
			sendcount_mpi[i] += sendcount[i * 2 + k];
		}
		senddispl[i + 1] = senddispl[i] + sendcount_mpi[i];
	}

	int* const recvcount = new int[team_size * 2];
	int* const recvcount_mpi = new int[team_size];

	// Send the number of elements sent from each process
	MPI_Alltoall(sendcount, 2, MPI_INT, recvcount, 2, MPI_INT,
	             global_comm::team_communicator);

	int* const sendbuffer = new int[senddispl[team_size]];

	// Send the actual elements
	int counter = 0;
	for (int i = 0; i < team_size; i++)
		for (int k = 0; k < 2; k++)
			for (set<int>::iterator it = send_info[i * 2 + k]->begin();
			     it != send_info[i * 2 + k]->end(); it++)

				sendbuffer[counter++] = (*it);

	assert(counter == senddispl[team_size]);

	int* const recvdispl = new int[team_size + 1];
	*displ = new int[team_size * 2];
	*count = new int[team_size * 2];
	recvdispl[0] = 0;
	int curr_recv_displ = 0;
	for (int i = 0; i < team_size; i++){
		recvcount_mpi[i] = 0;
		for (int k = 0; k < 2; k++){
			*count[i * 2 + k] = recvcount[i * 2 + k];
			recvcount_mpi[i] += recvcount[i * 2 + k];
			*displ[i * 2 + k] = curr_recv_displ;
			curr_recv_displ += recvcount[i * 2 + k];
		}
		recvdispl[i + 1] = recvdispl[i] + recvcount_mpi[i];
	}

	*recvbuf = new int[recvdispl[team_size]];

	MPI_Alltoallv(sendbuffer, sendcount_mpi, senddispl, MPI_INT,
	              *recvbuf, recvcount_mpi, recvdispl, MPI_INT,
	              global_comm::team_communicator);

	delete[] sendcount;
	delete[] sendcount_mpi;
	delete[] senddispl;
	delete[] curr_displ;
	delete[] recvcount;
	delete[] recvcount_mpi;
	delete[] sendbuffer;
	delete[] recvdispl;
}


void Inspector::GetBufferSize(){

	global_comm::max_send_size = 0;
	global_comm::max_recv_size = 0;
	global_comm::max_nprocs_read_send = 0;
	global_comm::max_nprocs_read_recv = 0;
	global_comm::max_nprocs_write_send = 0;
	global_comm::max_nprocs_read_send = 0;

	for( int iter_num = 0 ; iter_num < all_comm.size(); iter_num++ ){
		global_comm* curr_global_comm = all_comm[iter_num];
		int send_read_count = 0, recv_read_count = 0;
		int send_write_count = 0, recv_write_count = 0;

		curr_global_comm->read_send_offset[0] = 0;
		curr_global_comm->read_recv_offset[0] = 0;
		curr_global_comm->write_send_offset[0] = 0;
		curr_global_comm->write_recv_offset[0] = 0;

		for( int i = 0 ; i < nprocs ; i++ ){

			local_comm* local_recv_comm = all_local_comm[iter_num];
			local_comm* local_send_comm = all_local_comm[iter_num];
			send_read_count +=
				local_send_comm->GetReadSendCount(i,send_read_count);
			recv_read_count +=
				local_recv_comm->GetReadRecvCount(i,recv_read_count);
			send_write_count +=
				local_send_comm->GetWriteSendCount(i,send_write_count);
			recv_write_count +=
				local_recv_comm->GetWriteRecvCount(i,recv_write_count);
			curr_global_comm->read_send_count[i] =
				send_read_count - curr_global_comm->read_send_offset[i];
			curr_global_comm->read_send_offset[i+1] = send_read_count;
			curr_global_comm->read_recv_count[i] =
				recv_read_count - curr_global_comm->read_recv_offset[i];
			curr_global_comm->read_recv_offset[i+1] = recv_read_count;
			curr_global_comm->write_send_count[i] =
				send_write_count - curr_global_comm->write_send_offset[i];
			curr_global_comm->write_send_offset[i+1] = send_write_count;
			curr_global_comm->write_recv_count[i] = recv_write_count - curr_global_comm->write_recv_offset[i];
			curr_global_comm->write_recv_offset[i+1] = recv_write_count;
		}
    
		global_comm::max_send_size = ( global_comm::max_send_size > curr_global_comm->write_send_offset[nprocs] ? global_comm::max_send_size : curr_global_comm->write_send_offset[nprocs] );
		global_comm::max_send_size = ( global_comm::max_send_size > curr_global_comm->read_send_offset[nprocs] ? global_comm::max_send_size : curr_global_comm->read_send_offset[nprocs] );
		global_comm::max_recv_size = ( global_comm::max_recv_size > curr_global_comm->write_recv_offset[nprocs] ? global_comm::max_recv_size : curr_global_comm->write_recv_offset[nprocs] );
		global_comm::max_recv_size = ( global_comm::max_recv_size > curr_global_comm->read_recv_offset[nprocs] ? global_comm::max_recv_size : curr_global_comm->read_recv_offset[nprocs] );

		MPI_Alltoall(curr_global_comm->read_recv_offset,1,MPI_INT,curr_global_comm->read_put_displ,1,MPI_INT,global_comm::global_iec_communicator);
		MPI_Alltoall(curr_global_comm->write_recv_offset,1,MPI_INT,curr_global_comm->write_put_displ,1,MPI_INT,global_comm::global_iec_communicator);

		for( int i = 0 ; i < nprocs ; i++ ){
			if( curr_global_comm->read_send_count[i] != 0 )
				curr_global_comm->nprocs_read_send++;
			if( curr_global_comm->read_recv_count[i] != 0 )
				curr_global_comm->nprocs_read_recv++;
			if( curr_global_comm->write_send_count[i] != 0 )
				curr_global_comm->nprocs_write_send++;
			if( curr_global_comm->write_recv_count[i] != 0 )
				curr_global_comm->nprocs_write_recv++;
		}

		if( curr_global_comm->nprocs_read_send  != 0 )
			curr_global_comm->proc_id_read_send = new int[curr_global_comm->nprocs_read_send];
		if( curr_global_comm->nprocs_read_recv  != 0 )
			curr_global_comm->proc_id_read_recv = new int[curr_global_comm->nprocs_read_recv];
		if( curr_global_comm->nprocs_write_send  != 0 )
			curr_global_comm->proc_id_write_send = new int[curr_global_comm->nprocs_write_send];
		if( curr_global_comm->nprocs_write_recv  != 0 )
			curr_global_comm->proc_id_write_recv = new int[curr_global_comm->nprocs_write_recv];

		int r_s_count = 0,r_r_count = 0,w_s_count =0,w_r_count = 0;
		for( int i = 0 ; i < nprocs ; i++ ){
			if( curr_global_comm->read_send_count[i] != 0 )
				curr_global_comm->proc_id_read_send[r_s_count++] = i;
			if( curr_global_comm->read_recv_count[i] != 0 )
				curr_global_comm->proc_id_read_recv[r_r_count++] = i;
			if( curr_global_comm->write_send_count[i] != 0 )
				curr_global_comm->proc_id_write_send[w_s_count++] = i;
			if( curr_global_comm->write_recv_count[i] != 0 )
				curr_global_comm->proc_id_write_recv[w_r_count++] = i;
		}
		assert( r_s_count == curr_global_comm->nprocs_read_send );
		assert( r_r_count == curr_global_comm->nprocs_read_recv );
		assert( w_s_count == curr_global_comm->nprocs_write_send );
		assert( w_r_count == curr_global_comm->nprocs_write_recv );

		global_comm::max_nprocs_read_send = MAX( curr_global_comm->nprocs_read_send, global_comm::max_nprocs_read_send);
		global_comm::max_nprocs_read_recv = MAX( curr_global_comm->nprocs_read_recv, global_comm::max_nprocs_read_recv);
		global_comm::max_nprocs_write_send = MAX( curr_global_comm->nprocs_write_send, global_comm::max_nprocs_write_send);
		global_comm::max_nprocs_write_recv = MAX( curr_global_comm->nprocs_write_recv, global_comm::max_nprocs_write_recv);
	}

	if( global_comm::max_send_size > 0 )
		global_comm::send_buffer =
			(char*)ARMCI_Malloc_local(global_comm::max_send_size);
  
	global_comm::put_buffer = new char*[nprocs];
	ARMCI_Malloc((void**)global_comm::put_buffer,global_comm::max_recv_size+1);

	if( global_comm::max_recv_size > 0 )
		global_comm::recv_buffer = (char*)global_comm::put_buffer[proc_id];

	if( global_comm::max_nprocs_read_send > 0 ){
		int array_size = global_comm::max_nprocs_read_send;
		global_comm::read_send_start_status = new MPI_Status[array_size];
		global_comm::read_send_end_status = new MPI_Status[array_size];
		global_comm::read_send_end_request = new MPI_Request[array_size];
		global_comm::read_send_signal = new char[array_size];
	}
	if( global_comm::max_nprocs_read_recv > 0 ){
		int array_size = global_comm::max_nprocs_read_recv;
		global_comm::read_recv_start_status = new MPI_Status[array_size];
		global_comm::read_recv_end_status = new MPI_Status[array_size];
		global_comm::read_recv_start_request = new MPI_Request[array_size];
		global_comm::read_recv_signal = new char[array_size];
	}
	if( global_comm::max_nprocs_write_send > 0 ){
		int array_size = global_comm::max_nprocs_write_send;
		global_comm::write_send_start_status = new MPI_Status[array_size];
		global_comm::write_send_end_status = new MPI_Status[array_size];
		global_comm::write_send_end_request = new MPI_Request[array_size];
		global_comm::write_send_signal = new char[array_size];
	}
	if( global_comm::max_nprocs_write_recv > 0 ){
		int array_size = global_comm::max_nprocs_write_recv;
		global_comm::write_recv_start_status = new MPI_Status[array_size];
		global_comm::write_recv_end_status = new MPI_Status[array_size];
		global_comm::write_recv_start_request = new MPI_Request[array_size];
		global_comm::write_recv_signal = new char[array_size];
	}


#ifndef NDEBUG
	printf("PID:%d,SendBufferSize:%d,", proc_id,global_comm::max_send_size);
	printf("Send_Buffer:%p, ", global_comm::send_buffer);
	printf("RecvBufferSize:%d, ", global_comm::max_recv_size);
	printf("RecvBuffer:%p\n", global_comm::recv_buffer);
	fflush(stdout);
#endif
}


/**
 * \brief Communicate all ghosts from/to members of the team
 */
void Inspector::CommunicateGhosts(){

	// Number of send-ghosts for each process, all arrays
	int sendcount[nprocs];

	// senddispl[i] will contain the index where ghosts start for process i,
	// senddispl[i + 1] the index where ghosts end (not included) for process i
	int senddispl[nprocs + 1];

	// Number of receive-ghosts for each process, all arrays
	int recvcount[nprocs];

	// Same explanation as senddispl but for receive-ghosts
	int recvdispl[nprocs + 1];

	// Number of R/W arrays
	int n_data = 0;

	for (deque<global_data*>::iterator it = all_data.begin();
	     it != all_data.end(); it++ )

		if (!(*it)->is_read_only)
			n_data++;

	// #send-ghosts per target process, per R/W array
	int* sendghosts_count = new int[team_size * n_data];

	for (int i = 0; i < team_size; i++){
		sendcount[i] = 0;
		recvcount[i] = 0;
	}

	// Count #send-ghosts for all receiver processes
	for (int i = 0; i < team_size; i++){

		int dest_proc = i;
		int d = 0;

		// Count #send-ghosts for all arrays
		for (deque<local_data*>::iterator it = all_local_data.begin();
		     it != all_local_data.end(); it++)

			if (!(*it)->is_read_only){
				int nghosts = (*it)->global_ghosts[dest_proc].size();
				sendghosts_count[dest_proc * n_data + d] = nghosts;
				sendcount[i] += nghosts;
				d++;
			}
	}

	// Calculate indices where send-ghosts start & end for all processes
	senddispl[0] = 0;
	for (int i = 0; i < team_size; i++)
		senddispl[i + 1] = senddispl[i] + sendcount[i];

	// #receive-ghosts per target process, per R/W array
	int* recvghosts_count = new int[nprocs*n_data];

	// Send counts of send-ghosts and receive counts of receive-ghosts
	MPI_Alltoall(sendghosts_count, n_data, MPI_INT, recvghosts_count, n_data,
	             MPI_INT, global_comm::global_iec_communicator);

	// Extract #receive-ghosts for all sender processes
	for (int i = 0; i < nprocs; i++){
		int send_proc = i;
		for (int d = 0; d < n_data; d++){
			int nghosts = recvghosts_count[i * n_data + d];
			recvcount[i] += nghosts;
		}
	}

	// Calculate indices where receive-ghosts start & end for all processes
	recvdispl[0] = 0;
	for (int i = 0; i < nprocs; i++)
		recvdispl[i + 1] = recvdispl[i] + recvcount[i];

	int* ghosts_send_val = new int[senddispl[team_size]];
	int* ghosts_recv_val = new int[recvdispl[team_size]];

	// Populate the actual ghosts from data in the local arrays
	int counter = 0;
	for (int i = 0; i < nprocs; i++){

		int d = 0;
		for (deque<local_data*>::iterator it = all_local_data.begin();
		     it != all_local_data.end(); it++)

			if (!(*it)->is_read_only){
				for (set<int>::iterator jt = (*it)->global_ghosts[i].begin();
				     jt != (*it)->global_ghosts[i].end(); jt++){

					ghosts_send_val[counter] = (*jt);
					++counter;
				}
				d++;
			}
	}

	assert(counter == senddispl[nprocs]);

	// Communicate the send- and receive-ghosts
	MPI_Alltoallv(ghosts_send_val, sendcount, senddispl, MPI_INT,
	              ghosts_recv_val, recvcount, recvdispl, MPI_INT,
	              global_comm::global_iec_communicator);

	// Update the owned data with the received ghost values
	counter = 0;
	for (int i = 0; i < nprocs; i++){

		int d = 0;
		for (deque<local_data*>::iterator it = all_local_data.begin();
		     it != all_local_data.end(); it++)

			if (!(*it)->is_read_only){
				int nghosts = recvghosts_count[i * n_data + d];
				for (int g = 0; g < nghosts; g++)
					(*it)->global_owned[i].insert(ghosts_recv_val[counter++]);
				d++;
			}
	}

	delete[] ghosts_send_val;
	delete[] ghosts_recv_val;
	delete[] recvghosts_count;
	delete[] sendghosts_count;
}


void Inspector::CommunicateReads(int comm_num){

	int nparts = nprocs;

#ifdef COMM_TIME
	double start_t,stop_t,start_t1,stop_t1,start_t2,stop_t2,start_t3,stop_t3;
	start_t = rtclock();
#endif

	local_comm* curr_local_comm = all_local_comm[comm_num];

	if( global_comm::max_send_size > 0 )
		curr_local_comm->PopulateReadSendBuffer(global_comm::send_buffer);

	// For all participants in the communication with id "comm_num",
	// send and receive alive signals.
	for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_recv ; i++ )
		MPI_Isend(global_comm::read_recv_signal, 1, MPI_CHAR,
		          all_comm[comm_num]->proc_id_read_recv[i], comm_num,
		          global_comm::global_iec_communicator,
		          global_comm::read_recv_start_request+i);

	for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_send ; i++ )
		MPI_Recv(global_comm::read_send_signal+i, 1, MPI_CHAR,
		         all_comm[comm_num]->proc_id_read_send[i],
		         comm_num,global_comm::global_iec_communicator,
		         global_comm::read_send_start_status+i);

	if( all_comm[comm_num]->nprocs_read_recv > 0 )
		MPI_Waitall(all_comm[comm_num]->nprocs_read_recv,
		            global_comm::read_recv_start_request,
		            global_comm::read_recv_start_status);

	if( global_comm::max_send_size > 0 ){
		global_comm* curr_communicator = all_comm[comm_num];
		int count;
		for( int i =0 ; i < nprocs ; i++ )
			if( ( count = curr_communicator->read_send_count[i]  ) != 0 )
				ARMCI_NbPut(global_comm::send_buffer +
				            curr_communicator->read_send_offset[i],

				            global_comm::put_buffer[i] +
				            curr_communicator->read_put_displ[i],

				            count, i, NULL);
	}

	if( curr_local_comm->read_recv_count[0] > 0 ){

		vector<local_data*>::iterator it =
			curr_local_comm->read_arrays.begin();

		for (; it != curr_local_comm->read_arrays.end(); it++)
			(*it)->PopulateLocalGhosts
				(all_local_data[(*it)->my_num], proc_id);
	}

	ARMCI_WaitAll();
	ARMCI_AllFence();

	for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_send ; i++ )
		MPI_Isend(global_comm::read_send_signal, 1, MPI_CHAR,
		          all_comm[comm_num]->proc_id_read_send[i], comm_num,
		          global_comm::global_iec_communicator,
		          global_comm::read_send_end_request+i);

	for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_recv ; i++ )
		MPI_Recv(global_comm::read_recv_signal+i, 1, MPI_CHAR,
		         all_comm[comm_num]->proc_id_read_recv[i], comm_num,
		         global_comm::global_iec_communicator,
		         global_comm::read_recv_end_status+i);

	if( all_comm[comm_num]->nprocs_read_send > 0 )
		MPI_Waitall(all_comm[comm_num]->nprocs_read_send,
		            global_comm::read_send_end_request,
		            global_comm::read_send_end_status);

	if( global_comm::max_recv_size > 0 )
		curr_local_comm->ExtractReadRecvBuffer(global_comm::recv_buffer);

#ifdef COMM_TIME
	stop_t = rtclock();
	curr_local_comm->read_comm_time += stop_t - start_t;
#endif

}


/**
 * \brief Sends data to all consumer processes.
 */
void Inspector::CommunicateToNext(){

	assert(0 && "CommunicateToNext not implemented yet");
}

/**
 * \brief Receives data from producer processes
 */
void Inspector::GetFromPrevious(){
	assert(0 && "GetFromPrevious not implemented yet");
}


/*
 * FUNCTIONS STOLEN FROM local_inspector
 */

void Inspector::PopulateGlobalArrays(){

	for (deque<local_data*>::iterator it = all_local_data.begin();
	     it != all_local_data.end(); it++)

		(*it)->PopulateGlobalArray();
}
