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

#include "RunTime/global_loop.hpp"
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
 * \param pidTeam Identifier of this processor private to the team
 * \param teamsize Size (in number of processors) of each team
 * \param nl Number of loops
 * \param nd Number of "data" arrays
 * \param nad Number of "indirection" arrays
 * \param iterNumCount Array of iteration limits (must be of size nl)
 * \param dataNumCount Array of direct array sizes (must be of size nd)
 * \param ro Array of read-only flags for the direct arrays (must be of size nd)
 */
Inspector::Inspector(int pid, int np, int team, int pidTeam, int teamsize,
                     int nl, int nd, int nad, int* iterNumCount,
                     int* dataNumCount, int* ro):
	procId(pid),
	teamNum(team),
	idInTeam(pidTeam),
	teamSize(teamsize),
	nProcs(np),
	pins_size(-1),
	myLoop(NULL),
	iterNumOffset(new int[nl + 1]), // The first '1' is for my loop
	dataNumOffset(new int[nd + 1]){

	// Create new communicator associated to the inspector. This avoids having
	// to use MPI_COMM_WORLD.
	// Also create a communicator for the team.
	MPI_Comm_dup(MPI_COMM_WORLD, &global_comm::global_iec_communicator);
	MPI_Comm_split(MPI_COMM_WORLD, teamNum, 0,
	               &global_comm::team_communicator);
	MPI_Barrier(MPI_COMM_WORLD);

	if (pid >= np)
		return;

	/// Create all loops
	iterNumOffset[0] = 0;
	for (int iLoop = 0; iLoop < nl; ++iLoop){

		global_loop* new_loop =
			new global_loop(iLoop,
			                iterNumCount[iLoop],
			                iterNumOffset[iLoop]);
		allLoops.push_back(new_loop);
		iterNumOffset[iLoop + 1] =
			iterNumOffset[iLoop] + iterNumCount[iLoop];
	}

	/// Create data arrays (excludes indirection arrays)
	dataNumOffset[0] = 0;
	for (int i = 0; i < nd; i++){
		global_data* new_data;
		new_data = new global_data_double
			(procId, i, dataNumCount[i], dataNumOffset[i],
			 ro[i] == 1? true : false);
		allData.push_back(new_data);

		dataNumOffset[i + 1] = dataNumOffset[i] + dataNumCount[i];
	}

	// Set up the global IEC communicator for each loop
	// TODO: we should get rid of multiple global communicators (use only one).
	for (int i = 0; i < nl; i++){
		global_comm* new_comm = new global_comm(i, nProcs, procId);
		all_comm.push_back(new_comm);
	}

	// Set up IEC communicator local to the team
	team_comm = new local_comm(teamSize, idInTeam);

	// Indirection arrays (called access data).
	for (int i = 0; i < nad; i++){
		access_data* new_access_data =
			new access_data(i, procId, nProcs, pidTeam, teamsize);
		all_access_data.push_back(new_access_data);
	}

	// Allocate arrays for the pipeline communications
	pipeSendCounts = new int[np];
	pipeSendDispls = new int[np];
	pipeRecvCounts = new int[np];

	send_info = new set<int>*[teamSize * 2];
	for (int i = 0; i < teamSize * 2; i++)
		send_info[i] = new set<int>;

	MPI_Barrier(global_comm::global_iec_communicator);
}


Inspector::~Inspector(){

	for (deque<global_loop*>::iterator it = allLoops.begin();
	     it != allLoops.end(); it++)

		delete (*it);

	delete[] iterNumOffset;
	allLoops.clear();
	for (deque<global_data*>::iterator it = allData.begin();
	     it != allData.end(); it++)

		delete (*it);

	delete[] dataNumOffset;
	allData.clear();
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
		ARMCI_Free(global_comm::put_buffer[procId]);
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
		for (int i = 0; i < teamSize * 2; i++){
			send_info[i]->clear();
			delete send_info[i];
		}
		delete[] send_info;
	}

	// Free pipeline buffers
	delete [] pipeSendCounts;
	delete [] pipeSendDispls;
	delete [] pipeRecvCounts;

	// Deallocate MPI communicator common to all participating processes
	// in the inspector/executor.
	if (global_comm::team_communicator != MPI_COMM_NULL)
		MPI_Comm_free(&global_comm::team_communicator);
	if (global_comm::global_iec_communicator != MPI_COMM_NULL)
		MPI_Comm_free(&global_comm::global_iec_communicator);

	if (pipeSendBuf)
		delete pipeSendBuf;

	for (std::map<int, char*>::iterator
		     buf = pipeRecvBuf.begin(),
		     end = pipeRecvBuf.end();
	     buf != end;
	     ++buf)

		if (buf->second)
			free(pipeRecvBuf[buf->first]);

	if (internalMPIRequest)
		delete [] internalMPIRequest;
}


void Inspector::GetDontHave(){

	int stride = all_access_data.size(); // Number of indirection arrays
	int* get_n_elems = new int[teamSize * stride];
	int* send_n_elems = new int[teamSize * stride];

	/// Get the number of elements to get from other processes
	for( int i = 0 ; i < teamSize * stride ; i++ )
		get_n_elems[i] = 0;
	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->GetSendCounts(get_n_elems, stride);

	///Initial communication to setup the buffer
	MPI_Alltoall(get_n_elems, stride, MPI_INT, send_n_elems, stride, MPI_INT,
	             global_comm::global_iec_communicator);

	int send_count[teamSize];
	int recv_count[teamSize];
	int send_offset[teamSize + 1];
	int recv_offset[teamSize + 1];

	// Compute the size of message sent and received by each process
	send_offset[0] = 0; recv_offset[0] = 0;
	for (int i = 0; i < teamSize; i++){

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
	if (send_offset[teamSize] > 0)
		sendbuf = new int[send_offset[teamSize]];
	else
		sendbuf = NULL;

	// Allocate recv buffer
	if (recv_offset[teamSize] > 0)
		recvbuf = new int[recv_offset[teamSize]];
	else
		recvbuf = NULL;

	// Populate the send buffer with the indices requested by other processes
	int* curr_offset = new int[teamSize];
	for (int i = 0; i < teamSize; i++)
		curr_offset[i] = send_offset[i];
	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->PopulateBuffer(sendbuf, send_offset[teamSize], curr_offset);

	MPI_Alltoallv(sendbuf, send_count, send_offset, MPI_INT, recvbuf,
	              recv_count, recv_offset, MPI_INT,
	              global_comm::global_iec_communicator);

	// Populate the receive buffer with the values of all indices requested
	// by this process.
	for (int i = 0; i < teamSize; i++)
		curr_offset[i] = recv_offset[i];

	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->GetRequestedValue(recvbuf, recv_offset[teamSize], curr_offset,
		                         send_n_elems, stride);

	MPI_Alltoallv(recvbuf, recv_count, recv_offset, MPI_INT, sendbuf,
	              send_count, send_offset, MPI_INT,
	              global_comm::global_iec_communicator);

	///Add to set of elements that are known on each process
	for (int i = 0; i < teamSize; i++)
		curr_offset[i] = send_offset[i];

	for (deque<access_data*>::iterator it = all_access_data.begin();
	     it != all_access_data.end(); it++)

		(*it)->AddToHaveSet(sendbuf, send_offset[teamSize], curr_offset);

	delete[] get_n_elems;
	delete[] send_n_elems;
	if (sendbuf)
		delete[] sendbuf;
	if (recvbuf)
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
 * \brief Adds a vertex to the iteration/data hypergraph
 *
 * This function also classifies the loop as a "producer", "consumer",
 * or executed loop. The executed loop is the one that will be executed by this
 * process, while other loops are executed by other processes and will either
 * supply data to this process (producer) or wait for data to be communicated by
 * this process (consumers).
 *
 * \param loopID The vertex will be added to this loop.
 * \param iterValue The vertex represents this iteration number.
 */
void Inspector::AddVertex(int loopID, int iterValue){

	assert(loopID < allLoops.size());

	// Loop classification
	global_loop* currLoop = allLoops[loopID];
	if (loopID == teamNum)
		myLoop = currLoop;
	else if (loopID < teamNum)
		producerLoops[loopID] = currLoop;
	else if (loopID > teamNum)
		consumerLoops[loopID] = currLoop;
}


/**
 * \brief Adds a pin from the current vertex to a network representing some data
 *
 * \param arrayID Identifier of the array to be referenced.
 * \param index Index in the array to be referenced.
 * \param loop The pin is meaningful in the context of this loop.
 * \param iter Iteration of the loop where this data is used.
 * \param is_direct !=0 if addressing is affine; =0 if it depends on indirection
 *                  array.
 * \param is_ploop !=0 if the access is direct from the partitionable loops.
 */
void Inspector::AddPinToNet
(int arrayID, int index, int loopID, int iter, int is_direct, int is_ploop){

	assert(arrayID < allData.size());

	global_loop* loop = allLoops[loopID];
	global_data* data = allData[arrayID];

	// Only care about read arrays if we are in a loop that we own (because
	// we must wait for a producer to give us the data).
	if (loop->is_my_loop() && data->is_read(loopID))
		loop->add_used_array(iter, arrayID);
	// We care about write arrays in producers and our loops, because that's how
	// we calculate when the producer will communicate data to our loop.
	if ((loop->is_producer() || loop->is_my_loop()) && data->is_write(loopID)){
		loop->add_computed_array(iter, arrayID);
		data->use_in_loop(loopID, iter, index);
	}

	// Get the current net
	net* net_num = data->data_net_info[loopID][index];
	vertex* curr_vertex = loop->iter_vertex[iter];
	pin_info new_pin(curr_vertex, is_direct != 0 ? true : false);

	// If access is from a partitionable loop and is direct access
	if (is_ploop && is_direct){
		// The compile time analysis should ensure that there are not other
		// "direct accesses" to this array. The direct_vertex field should be
		// NULL or the same vertex as the curr_vertex
		assert(net_num->direct_vertex == NULL ||
		       net_num->direct_vertex == curr_vertex);
		net_num->direct_vertex = curr_vertex;
	}

	// Check for the pin already existing in the set
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

	assert(array_num < allData.size() && solver_num < all_solvers.size());
	assert(index < allData[array_num]->orig_array_size);
	net* curr_net = allData[array_num]->data_net_info[index][loop];

	assert(all_solvers[solver_num]->orig_net[row_num] == NULL ||
	       all_solvers[solver_num]->orig_net[row_num] == curr_net);

	all_solvers[solver_num]->orig_net[row_num] = curr_net;
}

void Inspector::PatohPrePartition(int loop){

	if (nProcs > 1){
		const int ia_size = dataNumOffset[allData.size()];

		int** armci_net_ia = new int*[nProcs];
		ARMCI_Malloc((void**)armci_net_ia, (ia_size + 1) * sizeof(int));
		int* const net_ia = (int*)armci_net_ia[procId];

		net_ia[0] = 0;
		int counter = 0;

		// Store the pin information for nets in CSR format, compute the ia
		for (int i = 0; i < allData.size(); i++)
			for (int j = 0; j < allData[i]->orig_array_size; j++){
				net_ia[counter + 1] = net_ia[counter] +
					allData[i]->data_net_info[loop][j]->pins.size() * 2 + 1;
				counter++;
			}
		assert(counter == ia_size);

		int net_ja_size[nProcs];

		///Total number of pins on each process
		MPI_Allgather(&(net_ia[ia_size]), 1, MPI_INT, net_ja_size, 1, MPI_INT,
		              global_comm::global_iec_communicator);

		int max_net_ja_size = 0;
		///Maximum size needed for the pin information from each process
		for (int i = 0; i < nProcs; i++)
			max_net_ja_size = (net_ja_size[i] > max_net_ja_size ?
			                   net_ja_size[i] : max_net_ja_size);

		int** armci_net_ja = new int*[nProcs];
		ARMCI_Malloc((void**)armci_net_ja, max_net_ja_size * sizeof(int));
		int* const net_ja = (int*)armci_net_ja[procId];

		///Populate the pin information on the local process
		counter = 0;
		for (int i = 0; i < allData.size(); i++)
			for (int j = 0; j < allData[i]->orig_array_size; j++){
				net* curr_net = allData[i]->data_net_info[loop][j];

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

		for (int i = 1; i < nProcs; i++){

			// Get the pin information from each process
			int dest_id = (procId + i) % nProcs;
			int source_id = (procId + nProcs - i) % nProcs;

			ARMCI_Get((void*)armci_net_ia[source_id], (void*)recv_ia,
			          (ia_size + 1) * sizeof(int), source_id);
			ARMCI_Get((void*)armci_net_ja[source_id], (void*)recv_ja,
			          recv_ia[ia_size] * sizeof(int), source_id);

			counter = 0;
			for (deque<global_data*>::iterator of_data = allData.begin();
			     of_data != allData.end(); of_data++){

				global_data* curr_data = (*of_data);

				///Update the local information
				for( int i = 0 ; i < curr_data->orig_array_size ; i++ ){

					net* curr_net = curr_data->data_net_info[loop][i];
					for (int j = recv_ia[counter];
					     j < recv_ia[counter + 1] - 1;
					     j+=2){

						int vertex_num = recv_ja[j];
						int k = 0;
						for( k = 0 ; k < allLoops.size() ; k++ )
							if( iterNumOffset[k+1] > vertex_num )
								break;
						assert(k != allLoops.size());
						pin_info new_pin
							(allLoops[k]->iter_vertex[vertex_num -
							                           iterNumOffset[k]],
							 (recv_ja[j + 1] != 0 ? true : false));
						curr_net->pins.insert(new_pin);
					}
					if( recv_ja[recv_ia[counter+1]-1] != -1 ){
						int vertex_num = recv_ja[recv_ia[counter+1]-1];
						int k = 0;
						for (k = 0; k < allLoops.size(); k++)
							if (iterNumOffset[k + 1] > vertex_num)
								break;
						assert(k != allLoops.size());
						curr_net->direct_vertex = allLoops[k]->
							iter_vertex[vertex_num - iterNumOffset[k]];
					}
					counter++;
				}
			}
		}

		ARMCI_Free(armci_net_ja[procId]);
		ARMCI_Free(armci_net_ia[procId]);
		ARMCI_Free_local(recv_ia);
		ARMCI_Free_local(recv_ja);
		delete [] armci_net_ia;
		delete [] armci_net_ja;
	}
}


void Inspector::PatohPartitionAll(){

	for (int loop = 0; loop < allLoops.size(); ++loop)
		PatohPartition(loop);
}


void Inspector::PatohPartition(int loop){

	PatohPrePartition(loop);

	if (nProcs > 1){
		int patoh_n,patoh_c,patoh_nc,patoh_cutsize;
		int *patoh_cellwts, *patoh_netwts, *patoh_xpins, *patoh_pins,
			*patoh_partwts, *patoh_partvec;
		float* patoh_targetwts=NULL;
		int i,j,k;

		PaToH_Parameters patoh;

		PaToH_Initialize_Parameters(&patoh, PATOH_CONPART,
		                            PATOH_SUGPARAM_QUALITY);

		//Number of vertices
		patoh_c = iterNumOffset[allLoops.size()];
		//number of nets
		patoh_n = 0;
		for( i = 0 ; i < allData.size() ; i++ )
			if( !allData[i]->is_read_only )
				patoh_n += allData[i]->orig_array_size;
		//Number of constraints
		patoh_nc = allLoops.size();
		//Number of partitions
		patoh._k = nProcs;

		//Weight for cells
		patoh_cellwts = new int[patoh_c * patoh_nc];
		int counter = 0;
		for (i = 0; i < allLoops.size(); i++)
			for (j = 0; j < allLoops[i]->num_iters; j++){
				for (k = 0; k < patoh_nc; k++)
					patoh_cellwts[counter*patoh_nc + k] = 0;
				patoh_cellwts[ counter*patoh_nc + i ] = 1;
				counter++;
			}

		//Weight and xpins for nets
		patoh_netwts = new int[patoh_n];
		patoh_xpins = new int[patoh_n+1];
		patoh_xpins[0] = 0;
		counter = 0;
		for( i = 0 ; i < allData.size() ; i++ )
			if( !allData[i]->is_read_only)
				for( j = 0 ; j < allData[i]->orig_array_size ; j++ ) {
					patoh_netwts[counter] =
						allData[i]->data_net_info[loop][j]->weight;
					patoh_xpins[counter + 1] = patoh_xpins[counter] +
						allData[i]->data_net_info[loop][j]->pins.size();
					counter++;
				}

		//pins for nets;
		pins_size = patoh_xpins[patoh_n];
		patoh_pins = new int[patoh_xpins[patoh_n]];
		counter = 0;
		for( i = 0 ; i < allData.size() ; i++ )
			if( !allData[i]->is_read_only )
				for( j = 0 ; j < allData[i]->orig_array_size ; j++ ) {
					set<pin_info, pin_comparator>::iterator jt =
						allData[i]->data_net_info[loop][j]->pins.begin();
					for (;
					     jt != allData[i]->data_net_info[loop][j]->pins.end();
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
		for (i = 0; i < allLoops.size(); i++)
			for (j = 0; j < allLoops[i]->num_iters; j++){
				if (patoh_partvec[counter] == procId)
					allLoops[i]->nproc_local++;
				allLoops[i]->iter_vertex[j]->home = patoh_partvec[counter++];
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
		for( int i = 0 ; i < allLoops.size() ; i++ )
			for( int j = 0 ; j < allLoops[i]->num_iters ; j++ ){
				allLoops[i]->nproc_local++;
				allLoops[i]->iter_vertex[j]->home = 0;
			}
	}

	PatohAfterPartition(loop);
}

void Inspector::PatohAfterPartition(int loop){

	//Decide homes for the nets
	int possible[nProcs];

	for (deque<global_data*>::iterator it = allData.begin();
	     it != allData.end();
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
						for (int i = 0; i < nProcs; i++)
							possible[i] = 0;
						for (set<pin_info,pin_comparator>::iterator jt =
							     curr_net->pins.begin();
						     jt != curr_net->pins.end();
						     jt++)

							possible[(*jt).pin->home]++;

						int maxval = -1;
						int counter = 0 , i = 0 ;
						while (counter < nProcs){
							if (possible[i] > maxval){
								maxval = possible[i];
								home = i;
							}
							counter++;
							i = (i + 1) % nProcs;
						}
					}
					assert(curr_net->home == -1 && home >= 0 && home <= nProcs);
					curr_net->home = home;
				}
			}
			else{
				///If constrained, then the array is assumed to be
				///block partitioned (ghosts added later)
				int curr_proc = 0;
				const int array_split = (*it)->orig_array_size / (nProcs);
				for( int j = 0 ; j < (*it)->orig_array_size ; j++ ){
					net* curr_net = (*it)->data_net_info[loop][j];
					curr_net->home = curr_proc;
					if ((j + 1) % array_split == 0)
						curr_proc = (curr_proc + 1 > nProcs - 1 ?
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

		int source_proc = procId - 1;
		for (int dest_proc = (procId + 1) % nProcs;
		     dest_proc != procId;
		     dest_proc = (dest_proc + 1) % nProcs, source_proc--){

			MPI_Request i_send_request;
			MPI_Status i_recv_status,i_send_status;

			if( source_proc < 0 )
				source_proc = nProcs - 1;

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
					for (; k < allData.size(); k++)
						if (dataNumOffset[k + 1] > net_num)
							break;
					if (k == allData.size())
						printf("ID=%d,Source=%d,net_num=%d,k=%d\n", procId,
						       source_proc, net_num, k);
					assert(k != allData.size());
					assert(!allData[k]->is_read_only);
					assert(curr_solver->orig_net[j] == NULL ||
					       curr_solver->orig_net[j] ==
					       allData[k]->data_net_info[loop][net_num -
					                                        dataNumOffset[k]
					                                        ]);
					curr_solver->orig_net[j] =
						allData[k]->data_net_info[loop][net_num -
						                                dataNumOffset[k]];
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
	if (loop == teamNum){
		for (deque<global_data*>::iterator it = allData.begin();
		     it != allData.end(); it++){

			if ((*it)->data_net_info[loop] == NULL)
				continue;

			int array_id = (*it)->id;
			int stride_size = (*it)->stride_size;
			int orig_array_size = (*it)->orig_array_size;
			bool read_only_flag =(*it)->is_read_only;
			bool is_constrained = (*it)->is_constrained; // Irrelevant in quake
			const net** data_net_info =
				const_cast<const net**>((*it)->data_net_info[loop]);
			add_local_data(array_id, stride_size, data_net_info,
			               orig_array_size, read_only_flag, is_constrained);
		}
	}
}


void Inspector::GetLocalAccesses(int array_num, int** recvbuf, int** displ,
                                 int** count){

	const global_data* curr_array = allData[array_num];
	net** const curr_nets = curr_array->data_net_info.at(teamNum);
	const int curr_array_size = curr_array->orig_array_size;
	const int curr_split = curr_array_size / teamSize;
	const int curr_start = curr_split * idInTeam;
	const int curr_end = (idInTeam == teamSize - 1 ?
	                      curr_array_size :
	                      curr_split * (idInTeam + 1));

	for (int i = 0; i < teamSize * 2; i++)
		send_info[i]->clear();
	int* const sendcount_mpi = new int[teamSize];
	bool* is_direct = new bool[teamSize];
	bool* flags = new bool[teamSize];

	// Find all the processes that access the part of array owned by the process
	for (int i = curr_start; i < curr_end; i++){

		const net* curr_net = curr_nets[i];
		for (int j = 0; j < teamSize; j++){
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

			const int access_proc = get_team_id((*it).pin->home);
			assert(access_proc != -1);
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

	// Calculate data counts and displacements
	int* const sendcount = new int[teamSize * 2];
	int* const curr_displ = new int[teamSize * 2 + 1];
	int* const senddispl = new int[teamSize + 1];
	senddispl[0] = 0;
	curr_displ[0] = 0;
	for (int i = 0; i < teamSize; i++){
		sendcount_mpi[i] = 0;
		for (int k = 0; k < 2; k++){
			sendcount[i * 2 + k] = send_info[i * 2 + k]->size();
			curr_displ[i * 2 + k + 1] =
				curr_displ[i * 2 + k] + sendcount[i * 2 + k];
			sendcount_mpi[i] += sendcount[i * 2 + k];
		}
		senddispl[i + 1] = senddispl[i] + sendcount_mpi[i];
	}

	int* const recvcount = new int[teamSize * 2];
	int* const recvcount_mpi = new int[teamSize];

	// Send the number of elements sent from each process
	MPI_Alltoall(sendcount, 2, MPI_INT, recvcount, 2, MPI_INT,
	             global_comm::team_communicator);

	int* const sendbuffer = new int[senddispl[teamSize]];

	// Populate the send buffer with the actual elements
	int counter = 0;
	for (int i = 0; i < teamSize; i++)

		// k takes values 0 (direct accesses) or 1 (indirect accesses)
		for (int k = 0; k < 2; k++){
			set<int>* send_info_item = send_info[i * 2 + k];
			for (set<int>::iterator it = send_info_item->begin();
			     it != send_info_item->end(); it++)
				sendbuffer[counter++] = (*it);
		}

	assert(counter == senddispl[teamSize]);

	int* const recvdispl = new int[teamSize + 1];
	*displ = new int[teamSize * 2];
	*count = new int[teamSize * 2];
	recvdispl[0] = 0;
	int curr_recv_displ = 0;

	for (int i = 0; i < teamSize; i++){
		recvcount_mpi[i] = 0;

		// k=0 for direct accesses, k=1 for indirect accesses
		for (int k = 0; k < 2; k++){
			(*count)[i * 2 + k] = recvcount[i * 2 + k];
			recvcount_mpi[i] += recvcount[i * 2 + k];
			(*displ)[i * 2 + k] = curr_recv_displ;
			curr_recv_displ += recvcount[i * 2 + k];
		}
		recvdispl[i + 1] = recvdispl[i] + recvcount_mpi[i];
	}

	*recvbuf = new int[recvdispl[teamSize]];

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

	for (int iter_num = 0 ; iter_num < all_comm.size(); iter_num++){
		global_comm* curr_global_comm = all_comm[iter_num];
		int send_read_count = 0, recv_read_count = 0;
		int send_write_count = 0, recv_write_count = 0;

		curr_global_comm->read_send_offset[0] = 0;
		curr_global_comm->read_recv_offset[0] = 0;
		curr_global_comm->write_send_offset[0] = 0;
		curr_global_comm->write_recv_offset[0] = 0;

		for (int i = 0; i < teamSize; i++){

			send_read_count +=
				team_comm->GetReadSendCount(i, send_read_count);
			recv_read_count +=
				team_comm->GetReadRecvCount(i, recv_read_count);
			send_write_count +=
				team_comm->GetWriteSendCount(i, send_write_count);
			recv_write_count +=
				team_comm->GetWriteRecvCount(i, recv_write_count);

			curr_global_comm->read_send_count[i] =
				send_read_count - curr_global_comm->read_send_offset[i];
			curr_global_comm->read_send_offset[i + 1] = send_read_count;
			curr_global_comm->read_recv_count[i] =
				recv_read_count - curr_global_comm->read_recv_offset[i];
			curr_global_comm->read_recv_offset[i + 1] = recv_read_count;
			curr_global_comm->write_send_count[i] =
				send_write_count - curr_global_comm->write_send_offset[i];
			curr_global_comm->write_send_offset[i + 1] = send_write_count;
			curr_global_comm->write_recv_count[i] =
				recv_write_count - curr_global_comm->write_recv_offset[i];
			curr_global_comm->write_recv_offset[i + 1] = recv_write_count;
		}

		global_comm::max_send_size =
			MAX(global_comm::max_send_size,
			    MAX(curr_global_comm->write_send_offset[nProcs],
			        curr_global_comm->read_send_offset[nProcs]));
		global_comm::max_recv_size =
			MAX(global_comm::max_recv_size,
			    MAX(curr_global_comm->write_recv_offset[nProcs],
			        curr_global_comm->read_recv_offset[nProcs]));

		MPI_Alltoall(curr_global_comm->read_recv_offset, 1, MPI_INT,
		             curr_global_comm->read_put_displ, 1, MPI_INT,
		             global_comm::global_iec_communicator);
		MPI_Alltoall(curr_global_comm->write_recv_offset, 1, MPI_INT,
		             curr_global_comm->write_put_displ, 1, MPI_INT,
		             global_comm::global_iec_communicator);

		for (int i = 0; i < nProcs; i++){
			if (curr_global_comm->read_send_count[i] != 0)
				curr_global_comm->nprocs_read_send++;
			if (curr_global_comm->read_recv_count[i] != 0)
				curr_global_comm->nprocs_read_recv++;
			if (curr_global_comm->write_send_count[i] != 0)
				curr_global_comm->nprocs_write_send++;
			if (curr_global_comm->write_recv_count[i] != 0)
				curr_global_comm->nprocs_write_recv++;
		}

		if (curr_global_comm->nprocs_read_send != 0)
			curr_global_comm->proc_id_read_send =
				new int[curr_global_comm->nprocs_read_send];
		if (curr_global_comm->nprocs_read_recv != 0)
			curr_global_comm->proc_id_read_recv =
				new int[curr_global_comm->nprocs_read_recv];
		if (curr_global_comm->nprocs_write_send != 0)
			curr_global_comm->proc_id_write_send =
				new int[curr_global_comm->nprocs_write_send];
		if (curr_global_comm->nprocs_write_recv != 0)
			curr_global_comm->proc_id_write_recv =
				new int[curr_global_comm->nprocs_write_recv];

		int r_s_count = 0, r_r_count = 0, w_s_count = 0, w_r_count = 0;
		for (int i = 0; i < nProcs; i++){
			if (curr_global_comm->read_send_count[i] != 0)
				curr_global_comm->proc_id_read_send[r_s_count++] = i;
			if (curr_global_comm->read_recv_count[i] != 0)
				curr_global_comm->proc_id_read_recv[r_r_count++] = i;
			if (curr_global_comm->write_send_count[i] != 0)
				curr_global_comm->proc_id_write_send[w_s_count++] = i;
			if (curr_global_comm->write_recv_count[i] != 0)
				curr_global_comm->proc_id_write_recv[w_r_count++] = i;
		}
		assert(r_s_count == curr_global_comm->nprocs_read_send);
		assert(r_r_count == curr_global_comm->nprocs_read_recv);
		assert(w_s_count == curr_global_comm->nprocs_write_send);
		assert(w_r_count == curr_global_comm->nprocs_write_recv);

		global_comm::max_nprocs_read_send =
			MAX(curr_global_comm->nprocs_read_send,
			    global_comm::max_nprocs_read_send);
		global_comm::max_nprocs_read_recv =
			MAX(curr_global_comm->nprocs_read_recv,
			    global_comm::max_nprocs_read_recv);
		global_comm::max_nprocs_write_send =
			MAX(curr_global_comm->nprocs_write_send,
			    global_comm::max_nprocs_write_send);
		global_comm::max_nprocs_write_recv =
			MAX(curr_global_comm->nprocs_write_recv,
			    global_comm::max_nprocs_write_recv);
	}

	if( global_comm::max_send_size > 0 )
		global_comm::send_buffer =
			(char*)ARMCI_Malloc_local(global_comm::max_send_size);

	global_comm::put_buffer = new char * [nProcs];
	ARMCI_Malloc((void**)global_comm::put_buffer,
	             global_comm::max_recv_size + 1);

	if( global_comm::max_recv_size > 0 )
		global_comm::recv_buffer = (char*)global_comm::put_buffer[procId];

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
}


/**
 * \brief Communicate all ghosts from/to members of the team
 */
void Inspector::CommunicateGhosts(){

	// Number of send-ghosts for each process, all arrays
	int sendcount[teamSize];

	// senddispl[i] will contain the index where ghosts start for process i,
	// senddispl[i + 1] the index where ghosts end (not included) for process i
	int senddispl[teamSize + 1];

	// Number of receive-ghosts for each process, all arrays
	int recvcount[teamSize];

	// Same explanation as senddispl but for receive-ghosts
	int recvdispl[teamSize + 1];

	// Number of R/W arrays used in this team
	int nArrays = 0;

	for (deque<global_data*>::iterator it = allData.begin();
	     it != allData.end(); it++ )

		if ((*it)->is_write(myLoop->get_loop_id()))
			nArrays++;

	// #send-ghosts per target process, per R/W array
	int* sendghosts_count = new int[teamSize * nArrays];

	for (int i = 0; i < teamSize; i++){
		sendcount[i] = 0;
		recvcount[i] = 0;
	}

	// Count #send-ghosts for all receiver processes
	for (int i = 0; i < teamSize; i++){

		int dest_proc = i;
		int d = 0;

		// Count #send-ghosts for all arrays
		for (map<int, local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++){

			local_data* array = it->second;

			if (!array->is_read_only){
				int nghosts = array->global_ghosts[dest_proc].size();
				sendghosts_count[dest_proc * nArrays + d] = nghosts;
				sendcount[i] += nghosts;
				d++;

			}
		}
	}

	// Calculate indexes where send-ghosts start & end for all processes
	senddispl[0] = 0;
	for (int i = 0; i < teamSize; i++)
		senddispl[i + 1] = senddispl[i] + sendcount[i];

	// #receive-ghosts per target process, per R/W array
	int* recvghosts_count = new int[teamSize * nArrays];

	// Send counts of send-ghosts and receive counts of receive-ghosts
	MPI_Alltoall(sendghosts_count, nArrays, MPI_INT, recvghosts_count, nArrays,
	             MPI_INT, global_comm::team_communicator);

	// Extract #receive-ghosts for all sender processes
	for (int i = 0; i < teamSize; i++)
		// For all arrays
		for (int d = 0; d < nArrays; d++)
			recvcount[i] += recvghosts_count[i * nArrays + d];

	// Calculate indexes where receive-ghosts start & end for all processes
	recvdispl[0] = 0;
	for (int i = 0; i < teamSize; i++)
		recvdispl[i + 1] = recvdispl[i] + recvcount[i];

	int* ghosts_send_val = new int[senddispl[teamSize]];
	int* ghosts_recv_val = new int[recvdispl[teamSize]];

	// Populate the actual ghosts from data in the local arrays
	int counter = 0;
	for (int i = 0; i < teamSize; i++){

		int d = 0;
		for (map<int, local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++){

			local_data* array = it->second;

			if (!array->is_read_only){
				for (set<int>::iterator jt = array->global_ghosts[i].begin();
				     jt != array->global_ghosts[i].end(); jt++){

					ghosts_send_val[counter] = *jt;
					++counter;
				}
				d++;
			}
		}
	}

	assert(counter == senddispl[teamSize]);

	// Communicate the send- and receive-ghosts
	MPI_Alltoallv(ghosts_send_val, sendcount, senddispl, MPI_INT,
	              ghosts_recv_val, recvcount, recvdispl, MPI_INT,
	              global_comm::team_communicator);

	// Update the owned data with the received ghost values
	counter = 0;
	for (int i = 0; i < teamSize; i++){

		int d = 0;
		for (map<int, local_data*>::iterator it = allLocalData.begin();
		     it != allLocalData.end(); it++){

			local_data* array = it->second;
			global_data* globalArray = allData[array->GetMyNum()];

			if (globalArray->is_write(myLoop->get_loop_id())){

				int nghosts = recvghosts_count[i * nArrays + d];
				for (int g = 0; g < nghosts; g++)
					array->global_owned[i].insert(ghosts_recv_val[counter++]);
				d++;
			}
		}
	}

	delete[] ghosts_send_val;
	delete[] ghosts_recv_val;
	delete[] recvghosts_count;
	delete[] sendghosts_count;
}


void Inspector::CommunicateReads(int comm_num){

	int nparts = nProcs;

#ifdef COMM_TIME
	double start_t,stop_t,start_t1,stop_t1,start_t2,stop_t2,start_t3,stop_t3;
	start_t = rtclock();
#endif

	if( global_comm::max_send_size > 0 )
		team_comm->PopulateReadSendBuffer(global_comm::send_buffer);

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
		for (int i = 0; i < nProcs; i++)
			if ((count = curr_communicator->read_send_count[i]) != 0)
				ARMCI_NbPut(global_comm::send_buffer +
				            curr_communicator->read_send_offset[i],

				            global_comm::put_buffer[i] +
				            curr_communicator->read_put_displ[i],

				            count, i, NULL);
	}

	if (team_comm->read_recv_count[0] > 0){

		vector<local_data*>::iterator it =
			team_comm->read_arrays.begin();

		for (; it != team_comm->read_arrays.end(); it++)
			(*it)->PopulateLocalGhosts(allLocalData[(*it)->my_num], procId);
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
		team_comm->ExtractReadRecvBuffer(global_comm::recv_buffer);

#ifdef COMM_TIME
	stop_t = rtclock();
	team_comm->read_comm_time += stop_t - start_t;
#endif

}


/*
 * FUNCTIONS STOLEN FROM local_inspector
 */

void Inspector::PopulateGlobalArrays(){

	for (map<int, local_data*>::iterator it = allLocalData.begin();
	     it != allLocalData.end(); it++){

		global_data* global_array = allData[it->first];
		if (global_array->is_last_write_in_pipeline())
			it->second->PopulateGlobalArray();
	}
}


/**
 * \brief Populates a local array from its corresponding global array
 *
 * \param an ID of the array to be populated
 * \param lb Allocated clean array to be populated
 * \param oa Original array
 * \param st Stride. Not sure what this is.
 */
void Inspector::PopulateLocalArray(int an, double* lb, double* oa, int st){

	local_data* local = allLocalData[an];
	local->PopulateLocalArray(allData[an], lb, oa, st);
}


/*
 * NEW PIPELINING FUNCTIONS
 */

/**
 * \brief Initializes structures required by the Inspector for this loop
 *
 * Mark arrays as used
 *
 * \param loopID Loop identifier
 * \param usedArrays List of array identifiers used in this loop
 * \param readInfo Flags for each array in usedArrays; =1 if the array is read
 *                 in the loop
 * \param writeInfo Flags for each array in usedArrays; =1 if the array is
 *                  written in the loop
 * \param lastWrite Flags for each array in usedArrays; =1 if the array is
 *                  written in this loop for the last time in the pipeline
 * \param nArrays Number of arrays used in the loop
 */
void Inspector::pipe_init_loop(const int loopID,
                               const int usedArrays[],
                               const bool readInfo[],
                               const bool writeInfo[],
                               const bool lastWrite[],
                               const int nArrays){

	// Find loop type
	global_loop* loop = allLoops[loopID];
	if (loopID < teamNum){
		loop->set_as_producer();
		producerLoops[loopID] = loop;
	}
	else if (loopID > teamNum){
		loop->set_as_consumer();
		consumerLoops[loopID] = loop;
	}
	else{
		loop->set_as_my_loop();
		myLoop = loop;
	}

	for (int i = 0; i < nArrays; ++i){

		global_data* array = allData[usedArrays[i]];

		if (!lastWrite[i])
			array->set_not_last_write_in_pipeline();

		array->use_in_loop(loopID, readInfo[i], writeInfo[i]);
	}
}


/**
 * \brief Calculate required communications between us and producer/consumers.
 */
void Inspector::pipe_calculate_comm_info(){

	/*
	 * SEND INFO
	 */

	// For all arrays in the code
	for (std::deque<global_data*>::iterator dataIt = allData.begin(),
		     dataEnd = allData.end();
	     dataIt != dataEnd;
	     ++dataIt){

		global_data* array = *dataIt;

		// Only if the array is writable in our loop, calculate send info
		int myLoopId = myLoop->get_loop_id();
		if (array->is_write(myLoopId))
			array->pipe_calc_sends(myLoopId);

		// Only if the array is readable in our loop, calculate receive info
		if (array->is_read(myLoopId))
			array->pipe_calc_recvs();
	}

	/*
	 * RECEIVE INFO
	 */

	// Make a list with all the producers for this process
	std::set<int> allProducers;
	safeIter.resize(myLoop->num_iters);
	for (int iter = 0; iter < myLoop->num_iters; ++iter){

		for (global_loop::ArrayIDList::iterator
			     aIt = myLoop->used_arrays_begin(iter),
			     aEnd = myLoop->used_arrays_end(iter);
		     aIt != aEnd;
		     ++aIt){

			local_data* localArray = allLocalData[*aIt];
			global_data* globalArray = allData[*aIt];
			int producer = globalArray->producer;
			if (producer >= 0)
				allProducers.insert(producer);
			int safe = globalArray->pipe_safe_iteration(iter);
			if (safeIter[iter][producer] < safe)
				safeIter[iter][producer] = safe;
		}
	}

	// Init required MPI structures for unblocking receives
	internalMPIRequest = new MPI_Request[allProducers.size()];
	int iPos = 0;
	for (std::set<int>::iterator iProc = allProducers.begin();
	     iProc != allProducers.end(); ++iProc, iPos++){

		recvWaitStruct[*iProc] = internalMPIRequest + iPos;
		recvWaitStructInv[iPos] = *iProc;
		receivedSoFar[*iProc] = -1;
	}
}


/**
 * \brief Initialize common structures required for communications
 */
void Inspector::pipe_init_comm_structs(){

	int maxItemsSent = 0;
	int maxItemsRecvd = 0;

	for (std::map<int, local_data*>::iterator
		     localIt = allLocalData.begin(),
		     localEnd = allLocalData.end();
	     localIt != localEnd;
	     ++localIt){

		local_data* array = localIt->second;
		maxItemsSent += array->pipe_get_max_sendcounts();

		int arrayMaxRecvd = array->pipe_get_max_recvcounts();
		if (arrayMaxRecvd > maxItemsRecvd)
			maxItemsRecvd = arrayMaxRecvd;
	}

	pipeSendBuf = (char*)malloc(maxItemsSent * sizeof(char));

	for (std::map<int, global_loop*>::iterator
		     pIt = producerLoops.begin(),
		     pEnd = producerLoops.end();
	     pIt != pEnd;
	     ++pIt)

		pipeRecvBuf[pIt->first] = (char*)malloc(maxItemsRecvd * sizeof(char));
}


/**
 * \brief Zeroes all send pipeline communication control arrays (counts, displs)
 */
void Inspector::pipe_reset_counts_and_displs(){

	for (int i = 0; i < nProcs; ++i){
		pipeSendCounts[i] = pipeSendDispls[i] = 0;
	}
}


/**
 * \brief Communicates to consumers
 *
 * This function must be called at the end of the iteration. It will
 * take care of sending data calculated in the previous iteration to the
 * consumers. This function assumes for now that each process only computes
 * one loop. All other loops are consumers, producers or unrelated loops.
 *
 * \param iter The iteration of the loop that just finished.
 */
void Inspector::pipe_send(int iter){

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Sending" << endl;
#endif

	pipe_reset_counts_and_displs();

	int nReqs = 0;
	MPI_Request reqs[nProcs];

	// Prepare structures for sending
	for (int iProc = 0; iProc < nProcs; ++iProc){

		// This is the "send" part
		for (global_loop::ArrayIDList::iterator
			     arrayIt = myLoop->computed_arrays_begin(iter),
			     arrayEnd = myLoop->computed_arrays_end(iter);
		     arrayIt != arrayEnd;
		     ++arrayIt){

			local_data* localArray = allLocalData[*arrayIt];
			pipeSendCounts[iProc] +=
				localArray->pipe_get_sendcounts(iter, iProc);

			localArray->pipe_populate_send_buf
				(iter,
				 iProc,
				 pipeSendBuf + pipeSendDispls[iProc]);
		}

		// Update the send displacements for the next process
		if (iProc + 1 < nProcs)
			pipeSendDispls[iProc + 1] =
				pipeSendDispls[iProc] + pipeSendCounts[iProc];

		if (pipeSendCounts[iProc] != 0)
			MPI_Issend(pipeSendBuf + pipeSendDispls[iProc],
			           pipeSendCounts[iProc],
			           MPI_BYTE,
			           iProc,
			           PIPE_TAG,
			           global_comm::global_iec_communicator,
			           &reqs[nReqs++]);
	}

	if (nReqs != 0)
		MPI_Waitall(nReqs, reqs, MPI_STATUS_IGNORE);

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Sending DONE" << endl;
#endif
}


void Inspector::pipe_initial_receive(){

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Issuing initial receives" << endl;
#endif

	for (std::map<int, MPI_Request*>::iterator
		     prodIt = recvWaitStruct.begin(), prodEnd = recvWaitStruct.end();
	     prodIt != prodEnd;
	     ++prodIt){

		// Try to receive from the producer in that iteration until we find the
		// first iteration where we really have to receive something from it.
		int iter;
		int prod = prodIt->first;
		for (iter = 0; !internal_issue_recv(iter, prodIt->first); ++iter);
	}

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Issuing initial receives DONE" << endl;
#endif
}


/**
 * \brief Gets data from producers
 *
 * This function must be called before the start of the iteration. It will
 * take care of receiving data needed for the next iteration from the
 * producers. This function assumes that each process only computes
 * one loop. All other loops are consumers, producers or unrelated loops.
 */
void Inspector::pipe_receive(int iter){

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Receiving for local iter " << iter << endl;
#endif

	while (!pipe_ready(iter)){

		int receivedIdx;
		MPI_Waitany(recvWaitStruct.size(), internalMPIRequest, &receivedIdx,
		            MPI_STATUS_IGNORE);
		assert(receivedIdx != MPI_UNDEFINED);
		int receivedProd = recvWaitStructInv[receivedIdx];
		int& iterReceivedFromProd = receivedSoFar[receivedProd];
		iterReceivedFromProd++;

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Received iter " << iterReceivedFromProd
	     << " for producer " << receivedProd << " DONE" << endl;
#endif

		// Move received data to the local array
		char* receivedData = pipeRecvBuf[receivedProd];
		for (global_loop::ArrayIDList::iterator
			     arrayIt = myLoop->used_arrays_begin(iterReceivedFromProd),
			     arrayEnd = myLoop->used_arrays_end(iterReceivedFromProd);
		     arrayIt != arrayEnd;
		     ++arrayIt){

			local_data* localArray = allLocalData[*arrayIt];
			localArray->pipe_update(iterReceivedFromProd,
			                        receivedProd,
			                        receivedData);
			receivedData += localArray->pipe_get_recvcounts
				(iterReceivedFromProd, receivedProd);
		}

		internal_issue_recv(iterReceivedFromProd + 1, receivedProd);
	}

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Receiving DONE" << endl;
#endif
}


bool Inspector::internal_issue_recv(int iter, int iProc){

#ifdef INSPECTOR_DBG
	cout << "[Proc " << procId << "] Issuing receive " << (iter)
	     << " for producer " << iProc << endl;
#endif

	pipeRecvCounts[iProc] = 0;

	// Prepare structures for receiving
	for (global_loop::ArrayIDList::iterator
		     arrayIt = myLoop->used_arrays_begin(iter),
		     arrayEnd = myLoop->used_arrays_end(iter);
	     arrayIt != arrayEnd;
	     ++arrayIt){

		local_data* localArray = allLocalData[*arrayIt];
		pipeRecvCounts[iProc] +=
			localArray->pipe_get_recvcounts(iter, iProc);
	}

	if (pipeRecvCounts[iProc] != 0){
		MPI_Irecv(pipeRecvBuf[iProc],
		          pipeRecvCounts[iProc],
		          MPI_BYTE,
		          iProc,
		          PIPE_TAG,
		          global_comm::global_iec_communicator,
		          recvWaitStruct[iProc]);
		return true;
	}

	return false;
}


/*
 * \brief Check that this process received everything needed for this iteration.
 *
 * \param iter The local iteration that we are receiving for.
 */
bool Inspector::pipe_ready(int iter){

	for (std::map<int, MPI_Request*>::iterator
		     procIt = recvWaitStruct.begin(),
		     procEnd = recvWaitStruct.end();
	     procIt != procEnd;
	     ++procIt){

		int producer = procIt->first;
		if (receivedSoFar[producer] < safeIter[iter][producer]){

#ifdef INSPECTOR_DBG
			cout << "[Proc " << procId << "] Have iter "
			     << receivedSoFar[producer] << " from producer " << producer
			     << ", need at least iter " << safeIter[iter][producer]
			     << ". NOT READY." << endl;
#endif
			return false;
		}

#ifdef INSPECTOR_DBG
		cout << "[Proc " << procId << "] Have iter "
		     << receivedSoFar[producer] << " from producer " << producer
		     << ", need at least iter " << safeIter[iter][producer]
		     << ". CONTINUE." << endl;
#endif
	}
	return true;
}
