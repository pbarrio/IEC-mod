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
 * @file: inspector.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */
#include "RunTime/inspector.hpp"
#include "armci.h"
#include "mpi.h"
#include "omp.h"
#include "patoh.h"

using namespace std;

inspector* inspector::singleton_inspector = NULL;
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
 * \param team Identifier of this processor's team
 * \param pid_team Identifier of this processor private to the team
 * \param teamsize Size (in number of processors) of each team
 * \param nl Number of loops
 * \param nd Number of direct access arrays
 * \param nc Number of... communicators?
 * \param nad Number of indirection arrays.
 * \param iter_num_count Array of iteration limits (must be of size nl)
 * \param data_num_count Array of direct array sizes (must be of size nd)
 * \param ro Array of read-only flags for the direct arrays (must be of size nd)
 */
inspector::inspector(int pid, int np, int team, int pid_team, int teamsize/*, int nt*/, int nl, int nd, int nc, int nad, int* iter_num_count, int* data_num_count, int* ro):
	proc_id(pid),
	team_num(team),
	id_in_team(pid_team),
	team_size(teamsize),
	nprocs(np),
	// 	nthreads(nt),
	pins_size(-1),
	my_loop(NULL),
	iter_num_offset(new int[nl + 1]), // The first '1' is for my loop
	data_num_offset(new int[nd+1])
{

	// Create new communicator associated to the inspector. This avoids having to use MPI_COMM_WORLD.
	//--- Pablo{
	MPI_Group worldGroup;
	MPI_Group iecGroup;
	int ranges[1][3] = {{0, np-1, 1}};
	MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
	MPI_Group_range_incl(worldGroup, 1, ranges, &iecGroup);
	MPI_Comm_create(MPI_COMM_WORLD, iecGroup, &global_comm::global_iec_communicator);
	MPI_Barrier(MPI_COMM_WORLD);

	if (pid >= np)
		return;
	//--- }

	// Create all loops
	iter_num_offset[0] = 0;
	for (int iLoop = 0 ; iLoop < nl ; ++iLoop){

		global_loop* new_loop =
			new global_loop(iLoop, iter_num_count[iLoop], iter_num_offset[iLoop]);
		all_loops.push_back(new_loop);
		iter_num_offset[iLoop + 1] = iter_num_offset[iLoop] + iter_num_count[iLoop];
	}

	data_num_offset[0] = 0;
	for( int i = 0 ; i < nd ; i++ ){
		global_data* new_data;
		if( ro[i] == 1 )
			new_data = new global_data_double(i,data_num_count[i],data_num_offset[i],true);
		else
			new_data = new global_data_double(i,data_num_count[i],data_num_offset[i],false);
		all_data.push_back(new_data);
    
		data_num_offset[i+1] = data_num_offset[i] + data_num_count[i];
	}

	for( int i = 0 ; i < nc ; i++ ){
		global_comm* new_comm = new global_comm(i,nprocs,proc_id);
		all_comm.push_back(new_comm);
	}

	for( int i = 0 ; i < nad ; i++ ){
		access_data* new_access_data = new access_data(i,proc_id,nprocs);
		all_access_data.push_back(new_access_data);
	}

	// 	send_info = new set<int>*[nprocs*nthreads*2];
	// 	for( int i = 0; i < nprocs*nthreads*2 ; i++ )
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


inspector::~inspector()
{
	for( deque<global_loop*>::iterator it = all_loops.begin() ; it != all_loops.end() ; it++ )
		delete (*it);
	delete[] iter_num_offset;
	all_loops.clear();
	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
		delete (*it);
	delete[] data_num_offset;
	all_data.clear();
	int counter = 0;
	for( deque<global_comm*>::iterator it = all_comm.begin() ; it != all_comm.end() ; it++ , counter++){
		delete (*it);
	}
	all_comm.clear();
	for( deque<access_data*>::iterator it = all_access_data.begin() ; it != all_access_data.end() ; it++ )
		delete (*it);
	all_access_data.clear();
  
	delete[] local_inspector::all_local_inspectors;
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

	// 	for( int i = 0 ; i < nprocs * nthreads * 2 ; i++ ){
	if (send_info){
		for( int i = 0 ; i < nprocs * 2 ; i++ ){
			send_info[i]->clear();
					delete send_info[i];
		}
		delete[] send_info;
	}

	// Deallocate MPI communicator common to all participating processes in the inspector/executor.
	if (global_comm::global_iec_communicator != MPI_COMM_NULL)
		MPI_Comm_free(&global_comm::global_iec_communicator);
}


void inspector::GetDontHave()
{
	int stride = all_access_data.size();
	int* get_n_elems = new int[nprocs*stride];
	int* send_n_elems = new int[nprocs*stride];

	for( int i = 0 ; i < nprocs * stride ; i++ )
		get_n_elems[i] = 0;
	for( deque<access_data*>::iterator it = all_access_data.begin() ; it != all_access_data.end() ; it++ )
		(*it)->GetSendCounts(get_n_elems,stride);

	///Initial communication to setup the buffer
	MPI_Alltoall(get_n_elems,stride,MPI_INT,send_n_elems,stride,MPI_INT,global_comm::global_iec_communicator);

	int send_count[nprocs];
	int recv_count[nprocs];
	int send_offset[nprocs+1];
	int recv_offset[nprocs+1];
  
	///Compute the size of message sent and recieved by each process
	send_offset[0] = 0; recv_offset[0] = 0;
	for( int i =0 ; i < nprocs ; i++ ){
		send_count[i] = 0; recv_count[i] = 0;
		for( int j = 0 ; j < stride ; j++ ){
			send_count[i] += get_n_elems[i*stride+j];
			recv_count[i] += send_n_elems[i*stride+j];
		}
		send_offset[i+1] = send_offset[i] + send_count[i];
		recv_offset[i+1] = recv_offset[i] + recv_count[i];
	}

	int *sendbuf,*recvbuf;
#ifndef NDEBUG
	printf("ID = %d, SendBufferSize = %d, recvbuffersize = %d\n",proc_id,send_offset[nprocs],recv_offset[nprocs]);
#endif
	///Allocate send buffer
	if( send_offset[nprocs] > 0 )
		sendbuf = new int[send_offset[nprocs]];
	else
		sendbuf = NULL;

	///Allocate recv buffer
	if( recv_offset[nprocs] > 0 )
		recvbuf = new int[recv_offset[nprocs]];
	else
		recvbuf = NULL;

	///Populate the buffer with the indices requested by each process
	int* curr_offset = new int[nprocs];
	for( int i =0 ; i < nprocs ; i++ )
		curr_offset[i] = send_offset[i];
	for( deque<access_data*>::iterator it = all_access_data.begin() ; it != all_access_data.end() ; it++ )
		(*it)->PopulateBuffer(sendbuf,send_offset[nprocs],curr_offset);

	MPI_Alltoallv(sendbuf,send_count,send_offset,MPI_INT,recvbuf,recv_count,recv_offset,MPI_INT,global_comm::global_iec_communicator);
  
	///Populate the buffer with the values of all requested indices
	for( int i = 0; i < nprocs ; i++ )
		curr_offset[i] = recv_offset[i];
	for( deque<access_data*>::iterator it = all_access_data.begin() ; it != all_access_data.end() ; it++ )
		(*it)->GetRequestedValue(recvbuf,recv_offset[nprocs],curr_offset,send_n_elems,stride);
  
	MPI_Alltoallv(recvbuf,recv_count,recv_offset,MPI_INT,sendbuf,send_count,send_offset,MPI_INT,global_comm::global_iec_communicator);
  
	///Add to set of elements that are known on each process
	for( int i =0 ; i < nprocs ; i++ )
		curr_offset[i] = send_offset[i];
	for( deque<access_data*>::iterator it = all_access_data.begin() ; it != all_access_data.end() ; it++ )
		(*it)->AddToHaveSet(sendbuf,send_offset[nprocs],curr_offset);

	delete[] get_n_elems;
	delete[] send_n_elems;
	if( sendbuf )
		delete[] sendbuf;
	if( recvbuf )
		delete[] recvbuf;
	delete[] curr_offset;
}

  
bool inspector::DoneGraphGeneration(){
	int n_dont = 0;
	for( int i =0 ; i < all_access_data.size() ; i++ )
		n_dont += all_access_data[i]->dont_have_set.size();
	MPI_Allreduce(MPI_IN_PLACE,&n_dont,1,MPI_INT,MPI_SUM,global_comm::global_iec_communicator);
#ifndef NDEBUG
	printf("MXD:ID:%d,n_dont=%d\n",proc_id,n_dont);
	fflush(stdout);
#endif
	if( n_dont != 0 ){
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
 * \param iter_num The vertex will be added to this loop.
 * \param iter_value The vertex represents this iteration number.
 */
void inspector::AddVertex(int iter_num, int iter_value)
{
	assert(iter_num < all_loops.size());

	global_loop* curr_loop = all_loops[iter_num];
	if (iter_num == team_num)
		my_loop = curr_loop;
	else if (iter_num < team_num)
		producer_loops[iter_num] = curr_loop;
	else if (iter_num > team_num)
		consumer_loops[iter_num] = curr_loop;

	curr_vertex = curr_loop->iter_vertex[iter_value];
}


/**
 * \brief Adds a pin from the current vertex to a network representing some data.
 *
 * \param data_num Identifier of the array to be referenced
 * \param index Identifier of the array position to be referenced
 * \param is_direct =0 if addressing is affine; !=0 if it depends on indirection array.
 * \param is_ploop !=0 if the access comes from a partitionable loop
 */
void inspector::AddNet(int data_num, int index, int is_direct, int is_ploop)
{
	assert(data_num < all_data.size());
	global_data* curr_data = all_data[data_num];

	///Get the current net
	net* net_num = all_data[data_num]->data_net_info[index];

	pin_info new_pin(curr_vertex,is_direct != 0 ? true : false);

	///If access is from a partitionable loops and is direct access
	if( is_ploop && is_direct ){
		///The compile time analysis should ensure that there are not other "direct accesses" to this array. The direct_vertex field should be NULL or the same vertex as the curr_vertex
		assert(net_num->direct_vertex == NULL || net_num->direct_vertex == curr_vertex);
		net_num->direct_vertex = curr_vertex;
	}

	///Check for the pin already existing in the set
	set<pin_info,pin_comparator>::iterator it = net_num->pins.find(new_pin);
	if( it != net_num->pins.end() ){
		///If it is not a direct access, and current access is direct, promote it to directo access.
		if( is_direct && !(*it).is_direct ){
			net_num->pins.erase(*it);
			net_num->pins.insert(new_pin);
			//printf("ID=%d,ReplacingIndirectWithDirect,array=%d,index=%d\n",proc_id,data_num,index);
		}
	}
	else{
		net_num->pins.insert(new_pin);
		//printf("ID=%d,AddingAccess,array=%d,index=%d,direct=%d\n",proc_id,data_num,index,is_direct);
	}
}

void inspector::AddUnknown(int solver_num, int array_num, int index, int row_num)
{
	assert(array_num < all_data.size() && solver_num < all_solvers.size() );
	assert(index < all_data[array_num]->orig_array_size);
	net* curr_net = all_data[array_num]->data_net_info[index];
	assert(all_solvers[solver_num]->orig_net[row_num] == NULL || all_solvers[solver_num]->orig_net[row_num] == curr_net );
	all_solvers[solver_num]->orig_net[row_num] = curr_net;
}

void inspector::PatohPrePartition()
{
	if( nprocs > 1 ) {
		const int ia_size = data_num_offset[all_data.size()];

		int** armci_net_ia = new int*[nprocs];
		ARMCI_Malloc((void**)armci_net_ia,(ia_size+1)*sizeof(int));
		int* const net_ia = (int*)armci_net_ia[proc_id];

		net_ia[0] = 0;
		int counter = 0;
		///Store the pin information for nets in CSR format, compute the ia
		for( int i = 0 ; i < all_data.size() ; i++ )
			for( int j = 0 ; j < all_data[i]->orig_array_size ; j++ ){
				net_ia[counter+1] = net_ia[counter] + all_data[i]->data_net_info[j]->pins.size() * 2 + 1;
				counter++;
			}
		assert( counter == ia_size );

		int net_ja_size[nprocs];
    
		///Total number of pins on each process
		MPI_Allgather(&(net_ia[ia_size]),1,MPI_INT,net_ja_size,1,MPI_INT,global_comm::global_iec_communicator);
    
		int max_net_ja_size = 0;
		///Maximum size needed for the pin information from each process
		for( int i =0 ; i < nprocs ; i++)
			max_net_ja_size = ( net_ja_size[i] > max_net_ja_size ? net_ja_size[i] : max_net_ja_size );
    
		int** armci_net_ja = new int*[nprocs];
		ARMCI_Malloc((void**)armci_net_ja,max_net_ja_size*sizeof(int));
		int* const net_ja = (int*)armci_net_ja[proc_id];

		///Populate the pin information on the local process
		counter = 0;
		for( int i = 0 ; i < all_data.size() ; i++ )
			for( int j = 0 ; j < all_data[i]->orig_array_size ; j++ ){
				net* curr_net = all_data[i]->data_net_info[j];
				for( set<pin_info,pin_comparator>::iterator it = curr_net->pins.begin() ; it != curr_net->pins.end() ; it++ ){
					net_ja[counter++] = (*it).pin->my_num;
					net_ja[counter++] = ( (*it).is_direct ? 1 : 0 );
				}
				if( curr_net->direct_vertex )
					net_ja[counter++] = curr_net->direct_vertex->my_num;
				else
					net_ja[counter++] = -1;
			}
		assert(counter == net_ia[ia_size]);

		///Local buffer to hold the pin information from different processes
		int* const recv_ia = (int*)ARMCI_Malloc_local((ia_size+1)*sizeof(int));
		int* const recv_ja = (int*)ARMCI_Malloc_local(max_net_ja_size*sizeof(int));

		MPI_Barrier(global_comm::global_iec_communicator);
    
		for( int i = 1 ; i < nprocs ; i++ ){
			///Get the pin information from each process
			int dest_id = (proc_id+i)%nprocs;
			int source_id = (proc_id+nprocs-i)%nprocs;
			ARMCI_Get((void*)armci_net_ia[source_id],(void*)recv_ia,(ia_size+1)*sizeof(int),source_id);
			ARMCI_Get((void*)armci_net_ja[source_id],(void*)recv_ja,recv_ia[ia_size]*sizeof(int),source_id);

			counter = 0;
			for( deque<global_data*>::iterator of_data = all_data.begin() ; of_data != all_data.end() ; of_data++ ){
				global_data* curr_data = (*of_data);

				///Update the local information
				for( int i = 0 ; i < curr_data->orig_array_size ; i++ ){
					net* curr_net = curr_data->data_net_info[i];
					for( int j = recv_ia[counter] ; j < recv_ia[counter+1] - 1 ; j+=2 ){
						int vertex_num = recv_ja[j];
						int k = 0;
						for( k = 0 ; k < all_loops.size() ; k++ )
							if( iter_num_offset[k+1] > vertex_num )
								break;
						assert(k != all_loops.size());
						pin_info new_pin(all_loops[k]->iter_vertex[vertex_num-iter_num_offset[k]],(recv_ja[j+1] != 0 ? true : false ));
						//curr_net->pins.insert(all_loops[k]->iter_vertex[vertex_num-iter_num_offset[k]]);
						curr_net->pins.insert(new_pin);
					}
					if( recv_ja[recv_ia[counter+1]-1] != -1 ){
						int vertex_num = recv_ja[recv_ia[counter+1]-1];
						int k = 0;
						for( k = 0 ; k < all_loops.size() ; k++ )
							if( iter_num_offset[k+1] > vertex_num )
								break;
						assert(k != all_loops.size());
						curr_net->direct_vertex = all_loops[k]->iter_vertex[vertex_num-iter_num_offset[k]];
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

void inspector::PatohPartition()
{
	PatohPrePartition();

	if( nprocs > 1 ){
		int patoh_n,patoh_c,patoh_nc,patoh_cutsize;
		int *patoh_cellwts,*patoh_netwts,*patoh_xpins,*patoh_pins,*patoh_partwts,*patoh_partvec;
		float* patoh_targetwts=NULL;  
		int i,j,k;

		PaToH_Parameters patoh;

		PaToH_Initialize_Parameters(&patoh,PATOH_CONPART,PATOH_SUGPARAM_QUALITY);

		//Number of vertices
		patoh_c = iter_num_offset[all_loops.size()]; //all_vertex.size();
		//number of nets
		patoh_n = 0;
		for( i = 0 ; i < all_data.size() ; i++ )
			if( !all_data[i]->is_read_only )
				patoh_n += all_data[i]->orig_array_size;
		//Number of constraints 
		patoh_nc = all_loops.size();
		//Number of partitions
		// 		patoh._k = nprocs*nthreads;
		patoh._k = nprocs;
#ifndef NDEBUG
		printf("PID:%d,Nvertices:%d,NNets:%d,Nconstraints:%d,Npartitions:%d\n",proc_id,patoh_c,patoh_n,patoh_nc,patoh._k);
		fflush(stdout);
		MPI_Barrier(global_comm::global_iec_communicator);
#endif

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
					patoh_netwts[counter] = all_data[i]->data_net_info[j]->weight;
					patoh_xpins[counter+1] = patoh_xpins[counter] + all_data[i]->data_net_info[j]->pins.size();
					counter++;
				}
  

		//pins for nets;
		pins_size = patoh_xpins[patoh_n];
		patoh_pins = new int[patoh_xpins[patoh_n]];
		counter = 0;
		for( i = 0 ; i < all_data.size() ; i++ )
			if( !all_data[i]->is_read_only )
				for( j = 0 ; j < all_data[i]->orig_array_size ; j++ ) {
					set<pin_info,pin_comparator>::iterator jt = all_data[i]->data_net_info[j]->pins.begin();
					for( ; jt != all_data[i]->data_net_info[j]->pins.end() ; jt++ )
						patoh_pins[counter++] = (*jt).pin->my_num;
				}

		PaToH_Alloc(&patoh,patoh_c,patoh_n,patoh_nc,patoh_cellwts,patoh_netwts,patoh_xpins,patoh_pins);

		patoh_partvec = new int[patoh_c];
		patoh_partwts = new int[patoh_nc*patoh._k];
		patoh.seed = 42;

		PaToH_Part(&patoh,patoh_c,patoh_n,patoh_nc,0,patoh_cellwts,patoh_netwts,patoh_xpins,patoh_pins,patoh_targetwts,patoh_partvec,patoh_partwts,&patoh_cutsize);

		counter = 0;
		for( i = 0 ; i < all_loops.size() ; i++ )
			for( j = 0 ; j < all_loops[i]->num_iters ; j++ ){
				// 				if( patoh_partvec[counter]/nthreads == proc_id )
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
	PatohAfterPartition();

}

void inspector::PatohAfterPartition()
{
	//Decide homes for the nets
	// 	int possible[nprocs*nthreads];
	int possible[nprocs];

	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
		if( !(*it)->is_read_only ){
			if( !(*it)->is_constrained ){
				for( int j = 0 ; j < (*it)->orig_array_size ; j++ ){
					net* curr_net = (*it)->data_net_info[j];
					int home = -1;
					///If there is a direct access use the home of that vertex
					if( curr_net->direct_vertex ){
						home = curr_net->direct_vertex->home;
					}
					else{
						///Add it to the same process as the one that accesses it the most (really doesnt matter)
						// 						for( int i = 0 ; i < nprocs*nthreads ; i++ )
						for( int i = 0 ; i < nprocs ; i++ )
							possible[i] = 0;
						for( set<pin_info,pin_comparator>::iterator jt = curr_net->pins.begin() ; jt != curr_net->pins.end() ; jt++ )
							possible[(*jt).pin->home]++;
						int maxval = -1;
						int counter = 0 , i = 0 ;
						// 						while( counter < nprocs * nthreads ){
						while( counter < nprocs ){
							if( possible[i] > maxval ) {
								maxval = possible[i];
								home = i;
							}
							counter++;
							// 							i = (i+1)%(nprocs*nthreads);
							i = (i+1)%(nprocs);
						}
					}
					// 					assert( curr_net->home == -1 && home >= 0 && home <= (nprocs * nthreads));
					assert( curr_net->home == -1 && home >= 0 && home <= (nprocs));
					curr_net->home = home;      
				}
			}
			else{
				///If constrained, then the array is assumed to be block partitioned (ghosts added later)
				int curr_proc = 0;
				// 				const int array_split = (*it)->orig_array_size / (nprocs*nthreads);
				const int array_split = (*it)->orig_array_size / (nprocs);
				for( int j = 0 ; j < (*it)->orig_array_size ; j++ ){
					net* curr_net = (*it)->data_net_info[j];
					curr_net->home = curr_proc;
					if( (j+1)%array_split == 0 )
						// 						curr_proc = ( curr_proc + 1 > nprocs*nthreads -1 ? curr_proc : curr_proc + 1);
						curr_proc = ( curr_proc + 1 > nprocs -1 ? curr_proc : curr_proc + 1);
				}
			}
		}
	AfterPartition();
}


void inspector::AfterPartition()
{
	for(deque<petsc_solve*>::iterator it = all_solvers.begin() ; it != all_solvers.end() ; it++ ){
		//Replicate the net info
		petsc_solve* curr_solver = (*it);
		const int curr_size = (*it)->size;
		int* send_buf = new int[curr_solver->size];
		int* recv_buf = new int[curr_solver->size];
    
		int source_proc = proc_id -1;
		for( int dest_proc = (proc_id + 1)%nprocs ; dest_proc != proc_id ; dest_proc = (dest_proc+1) % nprocs , source_proc--){
			MPI_Request i_send_request;
			MPI_Status i_recv_status,i_send_status;
      
			if( source_proc < 0 )
				source_proc = nprocs - 1;

			for( int j = 0 ; j < curr_size ; j++ )
				if( curr_solver->orig_net[j] != NULL )
					send_buf[j] = curr_solver->orig_net[j]->my_num;
				else
					send_buf[j] = -1;

			MPI_Isend(send_buf,curr_size,MPI_INT,dest_proc,0,global_comm::global_iec_communicator,&i_send_request);
			MPI_Recv(recv_buf,curr_size,MPI_INT,source_proc,0,global_comm::global_iec_communicator,&i_recv_status);
    
			for( int j = 0 ; j < curr_size ; j++ ){
				int net_num = recv_buf[j];
				if( net_num != -1 ){
					int k = 0;
					for( ; k < all_data.size() ; k++ )
						if( data_num_offset[k+1] > net_num)
							break;
					if( k == all_data.size() )
						printf("ID=%d,Source=%d,net_num=%d,k=%d\n",proc_id,source_proc,net_num,k);
					assert( k != all_data.size());
					assert(!all_data[k]->is_read_only );
					assert( curr_solver->orig_net[j] == NULL || curr_solver->orig_net[j] == all_data[k]->data_net_info[net_num - data_num_offset[k]] );
					curr_solver->orig_net[j] = all_data[k]->data_net_info[net_num - data_num_offset[k]];
				}
			}
			MPI_Wait(&i_send_request,&i_send_status);
		}
		delete[] send_buf;
		delete[] recv_buf;
		(*it)->FindNewRowNumbers();
	}


	///Setup the local inspectors (as many as number of threads)
	// TODO: Change into one local inspector per "loop" process in the chain of pipelined loops
	// 	local_inspector::all_local_inspectors = new local_inspector*[nthreads];
	// 	for( int i = 0 ; i < nthreads ; i++ )
	// 		local_inspector::all_local_inspectors[i] = new local_inspector(nprocs,nthreads,proc_id,i,all_comm.size());
	local_inspector::all_local_inspectors = new local_inspector*[1];
	local_inspector::all_local_inspectors[0] = new local_inspector(nprocs,proc_id,all_comm.size());
  
	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ ){
		int data_num = (*it)->my_num;
		int stride_size = (*it)->stride_size;
		int orig_array_size = (*it)->orig_array_size;
		bool read_only_flag =(*it)->is_read_only;
		bool is_constrained = (*it)->is_constrained;
		assert( (*it)->data_net_info != NULL );
		const net** data_net_info = const_cast<const net**>((*it)->data_net_info);
		// 		for( int i = 0 ; i < nthreads ; i++ ){
		// 			local_inspector::all_local_inspectors[i]->AddLocalData(data_num,stride_size,data_net_info,orig_array_size,read_only_flag,is_constrained);
		// 		}
		local_inspector::all_local_inspectors[0]->AddLocalData(data_num,stride_size,data_net_info,orig_array_size,read_only_flag,is_constrained);
	}
}


ret_data_access inspector::GetLocalAccesses(int array_num)
{
	const global_data* curr_array = all_data[array_num];
	net** curr_nets = curr_array->data_net_info;
	const int curr_array_size = curr_array->orig_array_size;
	const int curr_split = curr_array_size / nprocs;
	const int curr_start = curr_split * proc_id;
	const int curr_end = ( proc_id == nprocs - 1 ? curr_array_size : curr_split * ( proc_id + 1 ) );

	// 	for( int i =0 ; i < nprocs * nthreads * 2 ;i++  )
	for( int i = 0 ; i < nprocs * 2 ;i++  )
		send_info[i]->clear();
	int* const sendcount_mpi = new int[nprocs];
	// 	bool* is_direct = new bool[nprocs*nthreads];
	// 	bool* flags = new bool[nprocs*nthreads];
	bool* is_direct = new bool[nprocs];
	bool* flags = new bool[nprocs];

	///Find all the processes that access the blocked part of array owned by the process
	for( int i = curr_start ; i < curr_end ; i++ ){
		const net* curr_net = curr_nets[i];
		// 		for( int j = 0 ; j < nprocs*nthreads ; j++ ){
		for( int j = 0 ; j < nprocs ; j++ ){
			is_direct[j] = false;
			flags[j] = false;
		}

		for( set<pin_info,pin_comparator>::const_iterator it = curr_net->pins.begin() ;
		     it != curr_net->pins.end();
		     it++ )

			is_direct[(*it).pin->home] = is_direct[(*it).pin->home] || (*it).is_direct;

		for( set<pin_info,pin_comparator>::const_iterator it = curr_net->pins.begin();
		     it != curr_net->pins.end();
		     it++ ){

			const int access_proc = (*it).pin->home;
			assert(access_proc != -1 && access_proc < nprocs );
			if( !flags[access_proc] ){
				if( is_direct[access_proc] )
					send_info[access_proc*2]->insert(curr_net->data_index);
				else
					send_info[access_proc*2+1]->insert(curr_net->data_index);
				flags[access_proc] = true;
			}
		}
	}
	delete[] is_direct;
	delete[] flags;

	// 	int* const sendcount = new int[nprocs*nthreads*2];
	// 	int* const curr_displ = new int[nprocs*nthreads*2+1];
	int* const sendcount = new int[nprocs*2];
	int* const curr_displ = new int[nprocs*2+1];
	int* const senddispl = new int[nprocs+1];
	senddispl[0] = 0;
	curr_displ[0]=0;
	for( int i = 0 ; i < nprocs ; i++ ){
		sendcount_mpi[i] = 0;
		// 		for( int j = 0 ; j < nthreads ; j++ )
		for( int k = 0 ; k < 2 ; k++ ){
			// 				sendcount[i*nthreads*2+j*2+k] = send_info[(i*nthreads+j)*2+k]->size();
			// 				curr_displ[i*nthreads*2+j*2+k+1] = curr_displ[(i*nthreads+j)*2+k] + sendcount[i*nthreads*2+j*2+k];
			// 				sendcount_mpi[i] += sendcount[i*nthreads*2+j*2+k];
			sendcount[i*2+k] = send_info[i*2+k]->size();
			curr_displ[i*2+k+1] = curr_displ[i*2+k] + sendcount[i*2+k];
			sendcount_mpi[i] += sendcount[i*2+k];
		}
		senddispl[i+1] = senddispl[i] + sendcount_mpi[i];
	}

	// 	int* const recvcount = new int[nprocs*nthreads*2];
	int* const recvcount = new int[nprocs*2];
	int* const recvcount_mpi = new int[nprocs];

	///Send the number of elements sent from each process
	MPI_Alltoall(sendcount,/*nthreads**/2,MPI_INT,recvcount,/*nthreads**/2,MPI_INT,global_comm::global_iec_communicator);

	int* const sendbuffer = new int[senddispl[nprocs]];

	///Send the actual elements 
	int counter = 0;
	for( int i = 0 ; i < nprocs/**nthreads*/ ; i++ )
		for( int k = 0 ; k < 2 ; k++ )
			for( set<int>::iterator it = send_info[i*2+k]->begin() ; it != send_info[i*2+k]->end() ; it++ )
				sendbuffer[counter++] = (*it);
	assert(counter == senddispl[nprocs]);
      
	int* const recvdispl = new int[nprocs+1];
	int* const recv_thread_displ = new int[nprocs/**nthreads*/*2];
	int* const recv_thread_count = new int[nprocs/**nthreads*/*2];
	recvdispl[0] = 0; 
	int curr_recv_displ = 0;
	for( int i = 0; i < nprocs; i++ ){
		recvcount_mpi[i] = 0;
// 		for( int j = 0 ; j < nthreads ; j++ )
		for( int k = 0 ; k < 2 ; k++ ){
			recv_thread_count[/*j*nprocs*2+*/i*2+k] = recvcount[i/* * nthreads*/ *2 + /*2*j +*/ k];
			recvcount_mpi[i] += recvcount[i/* * nthreads*/ * 2 + /*2 * j +*/ k];
			recv_thread_displ[/*j*nprocs*2+*/i*2+k] = curr_recv_displ;
			curr_recv_displ += recvcount[i/* * nthreads*/ * 2 + /*2*j +*/ k ];
		}
		recvdispl[i+1] = recvdispl[i] + recvcount_mpi[i];
	}

 	int* const recvbuffer = new int[recvdispl[nprocs]];

	MPI_Alltoallv(sendbuffer,sendcount_mpi,senddispl,MPI_INT,
	              recvbuffer,recvcount_mpi,recvdispl,MPI_INT,
	              global_comm::global_iec_communicator);

	delete[] sendcount;
	delete[] sendcount_mpi;
	delete[] senddispl;
	delete[] curr_displ;
	delete[] recvcount;
	delete[] recvcount_mpi;
	delete[] sendbuffer;
	delete[] recvdispl;

	ret_data_access ret_val(recvbuffer,recv_thread_count,recv_thread_displ);
	return ret_val;
}



void inspector::GetBufferSize()
{
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
    
		curr_global_comm->read_send_offset[0] = 0, curr_global_comm->read_recv_offset[0] = 0;
		curr_global_comm->write_send_offset[0] = 0, curr_global_comm->write_recv_offset[0] = 0;
		for( int i = 0 ; i < nprocs ; i++ ){
      
// 			for( int j = 0 ; j < nthreads ; j++ ){
			local_comm* local_recv_comm = local_inspector::all_local_inspectors[0/*j*/]->all_comm[iter_num];
// 				for( int k = 0 ; k < nthreads ; k++ ){
			local_comm* local_send_comm = local_inspector::all_local_inspectors[0/*k*/]->all_comm[iter_num];
			send_read_count += local_send_comm->GetReadSendCount(i/**nthreads+j*/,send_read_count);
			recv_read_count += local_recv_comm->GetReadRecvCount(i/**nthreads+k*/,recv_read_count);
			send_write_count += local_send_comm->GetWriteSendCount(i/**nthreads+j*/,send_write_count);
			recv_write_count += local_recv_comm->GetWriteRecvCount(i/**nthreads+k*/,recv_write_count);
// 				}
// 			}
			curr_global_comm->read_send_count[i] = send_read_count - curr_global_comm->read_send_offset[i];
			curr_global_comm->read_send_offset[i+1] = send_read_count;
			curr_global_comm->read_recv_count[i] = recv_read_count - curr_global_comm->read_recv_offset[i];
			curr_global_comm->read_recv_offset[i+1] = recv_read_count;
      
			curr_global_comm->write_send_count[i] = send_write_count - curr_global_comm->write_send_offset[i];
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
		global_comm::send_buffer = (char*)ARMCI_Malloc_local(global_comm::max_send_size);
  
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
	printf("PID:%d,SendBufferSize:%d,Send_Buffer:%p,RecvBufferSize:%d,RecvBuffer:%p\n",proc_id,global_comm::max_send_size,global_comm::send_buffer,global_comm::max_recv_size,global_comm::recv_buffer);
	fflush(stdout);
#endif
}


void inspector::CommunicateGhosts()
{
	int sendcount[nprocs],senddispl[nprocs+1];
	int recvcount[nprocs],recvdispl[nprocs+1];

	int n_data = 0;
	for( deque<global_data*>::iterator it = all_data.begin() ; it != all_data.end() ; it++ )
		if( !(*it)->is_read_only )
			n_data++;

	int* sendghosts_count = new int[nprocs/**nthreads*nthreads*/*n_data];

	for( int i = 0 ; i < nprocs ; i++ ){
		sendcount[i] = 0;
		recvcount[i] = 0;
	}
  
	for( int i = 0 ; i < nprocs ; i++ ){
		// 		for( int j = 0 ; j < nthreads ; j++ ){
		int dest_proc = i/**nthreads+j*/;
		// 			for( int k = 0 ; k < nthreads ;k ++ ){
		// 				int source_thread = k;
		local_inspector* source_inspector = local_inspector::instance(/*k*/);
		int d=0;
		for( deque<local_data*>::iterator it = source_inspector->all_data.begin() ; it != source_inspector->all_data.end() ; it++ )
			if( !(*it)->is_read_only ){
				int nghosts = (*it)->global_ghosts[dest_proc].size();

				sendghosts_count[dest_proc/**nthreads*/*n_data+/*source_thread*n_data*/+d] = nghosts;
				sendcount[i] += nghosts;
				d++;
			}
		// 			}
		// 		}
	}
	senddispl[0] = 0;
	for(int i = 0 ; i < nprocs ; i++ )
		senddispl[i+1] = senddispl[i] + sendcount[i];
  
	int* recvghosts_count = new int[nprocs/**nthreads*nthreads*/*n_data];

	MPI_Alltoall(sendghosts_count,/*nthreads*nthreads**/n_data,MPI_INT,recvghosts_count,/*nthreads*nthreads**/n_data,MPI_INT,global_comm::global_iec_communicator);

	for( int i = 0 ; i < nprocs ; i++ ){
		// 		for( int j = 0 ; j < nthreads ; j++ ){
		// 			int recv_thread = j;
		local_inspector* recv_inspector = local_inspector::instance(/*j*/);
		// 			for( int k = 0 ; k < nthreads ; k++ ){
		int send_proc = i/**nthreads+k*/;
		for( int d = 0 ; d < n_data ; d++ ){
			int nghosts = recvghosts_count[i/**nthreads*nthreads*/*n_data+/*recv_thread*nthreads*n_data+k*n_data*/+d];
			recvcount[i] += nghosts;
		}
		// 	        }
		// 		}
	}
	recvdispl[0] = 0;
	for( int i = 0 ; i < nprocs ; i++ )
		recvdispl[i+1] = recvdispl[i] + recvcount[i];

  
	int* ghosts_send_val = new int[senddispl[nprocs]];
	int* ghosts_recv_val = new int[recvdispl[nprocs]];

	int counter = 0;
	for( int i = 0 ; i < nprocs ; i++ ){
// 		for( int j = 0 ; j < nthreads ; j++ ){
// 			int dest_proc = i*nthreads+j;
// 			for( int k = 0 ; k < nthreads ;k ++ ){
// 				int source_thread = k;
// 				local_inspector* source_inspector = local_inspector::instance(source_thread);

				local_inspector* source_inspector = local_inspector::instance();
				int d = 0; 

				for( deque<local_data*>::iterator it = source_inspector->all_data.begin() ; it != source_inspector->all_data.end() ; it++ )
					if( !(*it)->is_read_only ){
// 						for( set<int>::iterator jt = (*it)->global_ghosts[dest_proc].begin() ; jt != (*it)->global_ghosts[dest_proc].end() ; jt++ )

						for( set<int>::iterator jt = (*it)->global_ghosts[i].begin() ; jt != (*it)->global_ghosts[i].end() ; jt++ ){
							ghosts_send_val[counter/*++*/] = (*jt);
							++counter;
						}

						d++;
					}
// 			}
// 		}
	}

	assert(counter==senddispl[nprocs]);

	MPI_Alltoallv(ghosts_send_val,sendcount,senddispl,MPI_INT,ghosts_recv_val,recvcount,recvdispl,MPI_INT, global_comm::global_iec_communicator);

	counter = 0;
	for( int i = 0 ; i < nprocs ; i++ ){
// 		for( int j = 0 ; j < nthreads ; j++ ){
// 			int recv_thread = j;
		local_inspector* recv_inspector = local_inspector::instance(/*j*/);
// 			for( int k = 0 ; k < nthreads ; k++ ){
// 				int send_proc = i*nthreads+k;
				int d =0 ;
				for( deque<local_data*>::iterator it = recv_inspector->all_data.begin() ; it!= recv_inspector->all_data.end() ; it++ )
					if( !(*it)->is_read_only ){	  
						int nghosts = recvghosts_count[i/**nthreads*nthreads*/*n_data+/*recv_thread*nthreads*n_data+k*n_data*/+d];
						for( int g = 0 ; g < nghosts ; g++ )
// 							(*it)->global_owned[send_proc].insert(ghosts_recv_val[counter++]);
							(*it)->global_owned[i].insert(ghosts_recv_val[counter++]);
						d++;
					}
// 			}
// 		}
	}

	delete[] ghosts_send_val;
	delete[] ghosts_recv_val;
	delete[] recvghosts_count;
	delete[] sendghosts_count;

}


void inspector::CommunicateReads(/*int thread_id,*/ int comm_num)
{
	int nparts = nprocs/**nthreads*/;

#ifdef COMM_TIME
	double start_t,stop_t,start_t1,stop_t1,start_t2,stop_t2,start_t3,stop_t3;
	start_t = rtclock();
#endif

// 	local_comm* curr_local_comm = local_inspector::all_local_inspectors[thread_id]->all_comm[comm_num];
	local_comm* curr_local_comm = local_inspector::all_local_inspectors[0]->all_comm[comm_num];
	if( global_comm::max_send_size > 0 )
		curr_local_comm->PopulateReadSendBuffer(global_comm::send_buffer,global_comm::max_send_size);

// #pragma omp barrier

// #pragma omp master
	{
		//MPI_Barrier(MPI_COMM_WORLD);

		// For all participants in the communication with id "comm_num", send and receive alive signals
		for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_recv ; i++ )
			MPI_Isend(global_comm::read_recv_signal,1,MPI_CHAR,all_comm[comm_num]->proc_id_read_recv[i],comm_num,global_comm::global_iec_communicator,global_comm::read_recv_start_request+i);
		for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_send ; i++ )
			MPI_Recv(global_comm::read_send_signal+i,1,MPI_CHAR,all_comm[comm_num]->proc_id_read_send[i],comm_num,global_comm::global_iec_communicator,global_comm::read_send_start_status+i);
		if( all_comm[comm_num]->nprocs_read_recv > 0 )
			MPI_Waitall(all_comm[comm_num]->nprocs_read_recv,global_comm::read_recv_start_request,global_comm::read_recv_start_status);

		if( global_comm::max_send_size > 0 ){
			global_comm* curr_communicator = all_comm[comm_num];
			int count;
			for( int i =0 ; i < nprocs ; i++ )
				if( ( count = curr_communicator->read_send_count[i]  ) != 0 )
					ARMCI_NbPut(global_comm::send_buffer+curr_communicator->read_send_offset[i],global_comm::put_buffer[i]+curr_communicator->read_put_displ[i],count,i,NULL);
		}
	}

// 	for( int i = 1; i < nthreads ; i++ ){
// 		int source_thread = (thread_id+i)%nthreads;
// 		int source_id = proc_id*nthreads+source_thread;
// 		if( curr_local_comm->read_recv_count[source_id] > 0 ){
		if( curr_local_comm->read_recv_count[0] > 0 ){
// 			local_inspector* source_inspector = local_inspector::all_local_inspectors[source_thread];
			local_inspector* source_inspector = local_inspector::all_local_inspectors[0];
			vector<local_data*>::iterator it = curr_local_comm->read_arrays.begin();
			for( ; it != curr_local_comm->read_arrays.end() ; it++ )
// 				(*it)->PopulateLocalGhosts(source_inspector->all_data[(*it)->my_num],source_id);
				(*it)->PopulateLocalGhosts(source_inspector->all_data[(*it)->my_num],proc_id);
		}
// 	}

// #pragma omp master
 	{
		ARMCI_WaitAll();
		//ARMCI_Barrier();
		ARMCI_AllFence();

		for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_send ; i++ )
			MPI_Isend(global_comm::read_send_signal,1,MPI_CHAR,all_comm[comm_num]->proc_id_read_send[i],
			          comm_num,global_comm::global_iec_communicator,global_comm::read_send_end_request+i);

		for( int i = 0 ; i < all_comm[comm_num]->nprocs_read_recv ; i++ )
			MPI_Recv(global_comm::read_recv_signal+i,1,MPI_CHAR,all_comm[comm_num]->proc_id_read_recv[i],
			         comm_num,global_comm::global_iec_communicator,global_comm::read_recv_end_status+i);

		if( all_comm[comm_num]->nprocs_read_send > 0 )
			MPI_Waitall(all_comm[comm_num]->nprocs_read_send,global_comm::read_send_end_request,global_comm::read_send_end_status);
 	}

// #pragma omp barrier

	if( global_comm::max_recv_size > 0 )
		curr_local_comm->ExtractReadRecvBuffer(global_comm::recv_buffer,global_comm::max_recv_size);

#ifdef COMM_TIME
	stop_t = rtclock();
	curr_local_comm->read_comm_time += stop_t - start_t;
#endif

}


void inspector::CommunicateWrites(/*int thread_id,*/ int comm_num)
{

#ifdef COMM_TIME
	double start_t = rtclock();
#endif

// 	local_comm* curr_local_comm = local_inspector::all_local_inspectors[thread_id]->all_comm[comm_num];
	local_comm* curr_local_comm = local_inspector::all_local_inspectors[0]->all_comm[comm_num];

	// Fills in the buffer that will be sent to all neighbour processes. Here's all the meat!
	if( global_comm::max_send_size > 0 )
		curr_local_comm->PopulateWriteSendBuffer(global_comm::send_buffer);

// #pragma omp barrier

 
// #pragma omp master
	{

		// Sync with the receivers and senders
		for( int i = 0 ; i < all_comm[comm_num]->nprocs_write_recv ; i++ )
			MPI_Isend(global_comm::write_recv_signal,1,MPI_CHAR,all_comm[comm_num]->proc_id_write_recv[i],comm_num,global_comm::global_iec_communicator,global_comm::write_recv_start_request+i);
		for( int i = 0 ; i < all_comm[comm_num]->nprocs_write_send ; i++ )
			MPI_Recv(global_comm::write_send_signal+i,1,MPI_CHAR,all_comm[comm_num]->proc_id_write_send[i],comm_num,global_comm::global_iec_communicator,global_comm::write_send_start_status+i);
		if( all_comm[comm_num]->nprocs_write_recv > 0 )
			MPI_Waitall(all_comm[comm_num]->nprocs_write_recv,global_comm::write_recv_start_request,global_comm::write_recv_start_status);

		if( global_comm::max_send_size > 0 ){
			global_comm* curr_communicator = all_comm[comm_num];
			int count;
			for(int i = 0; i < nprocs; i++){

				if( ( count = curr_communicator->write_send_count[i]  ) != 0 )

					ARMCI_NbPut(global_comm::send_buffer + curr_communicator->write_send_offset[i],
					            global_comm::put_buffer[i] + curr_communicator->write_put_displ[i],
					            count,
					            i,
					            NULL);
			}
		}
	}

// 	for( int i = 1; i < nthreads ; i++ ){
// 		int source_thread = (thread_id + i)%nthreads;
// 		int source_id = proc_id*nthreads+source_thread;
// 		if( curr_local_comm->write_recv_count[source_id] > 0 ){
// 			local_inspector* source_inspector = local_inspector::all_local_inspectors[source_thread];
// 			vector<local_data*>::iterator it = curr_local_comm->write_arrays.begin();
      
// 			for(; it != curr_local_comm->write_arrays.end() ; it++ ){
// 				(*it)->UpdateLocalOwned(source_inspector->all_data[(*it)->my_num],source_id);
// 			}
// 		}
// 	}  

// #pragma omp master
	{
		ARMCI_WaitAll();
		//ARMCI_Barrier();
		ARMCI_AllFence();
		for( int i = 0 ; i < all_comm[comm_num]->nprocs_write_send ; i++ )
			MPI_Isend(global_comm::write_send_signal,1,MPI_CHAR,all_comm[comm_num]->proc_id_write_send[i],comm_num,global_comm::global_iec_communicator,global_comm::write_send_end_request+i);
		for( int i = 0 ; i < all_comm[comm_num]->nprocs_write_recv ; i++ )
			MPI_Recv(global_comm::write_recv_signal+i,1,MPI_CHAR,all_comm[comm_num]->proc_id_write_recv[i],comm_num,global_comm::global_iec_communicator,global_comm::write_recv_end_status+i);
		if( all_comm[comm_num]->nprocs_write_send > 0 )
			MPI_Waitall(all_comm[comm_num]->nprocs_write_send,global_comm::write_send_end_request,global_comm::write_send_end_status);
	}
  
// #pragma omp barrier

	if( global_comm::max_recv_size > 0 )
		curr_local_comm->ExtractWriteRecvBuffer(global_comm::recv_buffer);

#ifdef COMM_TIME
	double stop_t = rtclock();
	curr_local_comm->write_comm_time += (stop_t - start_t);
#endif

}


void inspector::print_hypergraph(FILE* outfile)
{
	printf("PID:%d,PrintingHypergraph\n",proc_id);
	fflush(stdout);
	fprintf(outfile,"%% Input file to be used with PaToH\n");
	int pins_size = 0;
	fprintf(outfile,"0 %d %d %d 3 %d\n",iter_num_offset[all_loops.size()],data_num_offset[all_data.size()],pins_size,(int)all_loops.size());
	fflush(outfile);
	for( int i = 0 ; i  < all_data.size() ; i++ )
		for( int j = 0 ; j <  all_data[i]->orig_array_size ; j++ )
			{
				fprintf(outfile,"%d ",all_data[i]->data_net_info[j]->weight);
				set<pin_info>::iterator it;
				for( it = all_data[i]->data_net_info[j]->pins.begin() ; it != all_data[i]->data_net_info[j]->pins.end() ; it++ )
					fprintf(outfile,"%d ",(*it).pin->my_num);
#ifndef NO_COMMENTS
				fprintf(outfile,"%% Array : %d, Index : %d, Home : %d, type = ",all_data[i]->data_net_info[j]->data_num,all_data[i]->data_net_info[j]->data_index,all_data[i]->data_net_info[j]->home);
				for( it = all_data[i]->data_net_info[j]->pins.begin() ; it != all_data[i]->data_net_info[j]->pins.end() ; it++ )
					fprintf(outfile," %d",(*it).is_direct);
#endif
				fprintf(outfile,"\n");
			}
	fprintf(outfile,"%% Printing weights for vertices\n");
	for( int i = 0 ; i < all_loops.size() ; i++ )
		for( int j = 0 ; j < all_loops[i]->num_iters ; j++ )
			{
				for( int k = 0 ; k < all_loops.size() ; k++ )
					if( k != i )
						fprintf(outfile,"0 ");
					else
						fprintf(outfile,"1 ");
#ifndef NO_COMMENTS
				fprintf(outfile,"%% Loop : %d, Value : %d, Home = %d",all_loops[i]->iter_vertex[j]->iter_num,all_loops[i]->iter_vertex[j]->iter_value,all_loops[i]->iter_vertex[j]->home);
#endif
				fprintf(outfile,"\n");
			}
	fflush(outfile);
	MPI_Barrier(global_comm::global_iec_communicator);
}
