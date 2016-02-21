/*
 * local_data.cpp: This file is part of the IEC project.
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
 * @file: local_data.cpp
 * @author: Mahesh Ravishankar <ravishan@cse.ohio-state.edu>
 */

#include "RunTime/communicator.hpp"
#include "RunTime/local_data.hpp"

#include "mpi.h"

#include <string.h>

using namespace std;

int binary_search(int* const array, const int n, const int value){

	int left = 0, right = n - 1;
	int mid;
	while (left <= right){
		mid = (left + right) / 2;
		if (array[mid] == value)
			break;
		else if (array[mid] < value){
			left = mid + 1;
		}
		else
			right = mid - 1;
	}
	if (left > right)
		return -1;
	else
		return mid;
}


/**
 * \brief Constructor
 *
 * \param mn Identifier of the global array
 * \param np The original array is split into these many parts for comms
 * \param pid Identifier of the current process
 * \param st Size of each position in the array (in array units, not bytes)
 * \param dni Nets corresponding to all the positions in the array
 * \param oas Size of the original array
 * \param iro True if the data is read-only
 * \param ic Unused in quake
 */
local_data::local_data(int mn, int np, int pid, int st, const net** dni,
                       int oas, bool iro, bool ic):
	my_num(mn),
	nparts(np),
	stride(st),
	proc_id(pid),
	data_net_info(dni),
	orig_array_size(oas),
	ghosts_offset(new int[np + 1]),
	owned_offset(new int[np + 1]),
	local_array_size(-1),
	is_read_only(iro),
	is_constrained(ic){

	ghosts = NULL;
	owned = NULL;
	l_to_g = NULL;
	direct_access_array = NULL;
	indirect_access_array = NULL;
	direct_access_size = 0;
	indirect_access_size = 0;
	block_owned_offset = -1;
}


local_data::~local_data(){

	delete[] ghosts_offset;
	delete[] owned_offset;
	if(ghosts)
		delete[] ghosts;
	if( owned )
		delete[] owned;
	if( l_to_g )
		delete[] l_to_g;
	delete[] direct_access_array;
	delete[] indirect_access_array;
}


/**
 * \brief Mark an access to a position of the array as "direct" or "indirect"
 *
 * A direct access is that done with a number as the index (e.g. a[1]). An
 * indirect access is that done by means of another array (e.g. a[b[3]], where
 * "a" is the accessed array and "b" is the indirection array).
 *
 * \param index Index accessed in the array
 * \param access_type 1 if direct access, else indirect
 */
void local_data::AddIndexAccessed(int index, int access_type){

	if (access_type == 1){
		set<int>::iterator it = indirect_access.find(index);
		if (it != indirect_access.end())
			indirect_access.erase(it);
		direct_access.insert(index);
	}
	else{
		set<int>::iterator it = direct_access.find(index);
		if (it == direct_access.end())
			indirect_access.insert(index);
	}
}


/**
 * \brief Insert direct access from this process to a region of data
 *
 * \param read_buffer Pointer to the start of the data
 * \param num_vals Size of the data to be directly accessed
 */
void local_data::InsertDirectAccess(const int* read_buffer, const int num_vals){
	direct_access.insert(read_buffer, read_buffer + num_vals);
}


/**
 * \brief Insert indirect access from this process to a region of data
 *
 * \param read_buffer Pointer to the start of the data
 * \param num_vals Size of the data to be indirectly accessed
 */
void local_data::InsertIndirectAccess
(const int* read_buffer, const int num_vals){

	indirect_access.insert(read_buffer, read_buffer + num_vals);
}


/**
 * \brief Initializes important attributes of the local array
 *
 * This function DOES NOT populate the local array from the global array.
 * There is a function PopulateLocalArray for that.
 */
void local_data::SetupLocalArray(){

	if (is_constrained){
		const int block_start = 0;
		const int block_stop = orig_array_size;
		assert(direct_access.size() == 0);
		for (int i = block_start; i < block_stop; i++)
			indirect_access.insert(i);
	}
	direct_access_size = direct_access.size();
	direct_access_array = new int[direct_access_size];
	indirect_access_size = indirect_access.size();
	indirect_access_array = new int[indirect_access_size];

	int i = 0;
	for (set<int>::const_iterator it = direct_access.begin();
	     it != direct_access.end(); it++, i++)
		direct_access_array[i] = *it;
	i = 0;
	for (set<int>::const_iterator it = indirect_access.begin();
	     it != indirect_access.end(); it++, i++)
		indirect_access_array[i] = *it;

	direct_access.clear();
	indirect_access.clear();
}

/**
 * \brief Find the local index corresponding to a global index
 *
 * \param global_index The global index to look for
 */
int local_data::GetLocalIndex(int global_index) const{

	int local_index = -1;

	if (direct_access_size > 0)
		local_index = binary_search(direct_access_array, direct_access_size,
		                            global_index);

	if (local_index == -1){
		assert(indirect_access_size > 0);
		local_index = binary_search(indirect_access_array, indirect_access_size,
		                            global_index);
		assert(local_index != -1);
		local_index += direct_access_size;
	}

	return local_index;
}


/**
 * \brief Renumber an access array from the global indexes to the local ones
 *
 * An access array contains the precomputed access indexes for arrays that are
 * accessed through indirection arrays. These indexes computed in the inspector
 * refer to the global array. This function translates them to the corresponding
 * local indexes.
 *
 * \param array_size Size of this local array
 * \param access_array Array that will be used to index this local array
 */
void local_data::RenumberAccessArray(int array_size, int* access_array){

	for (int i = 0; i < array_size; i++)
		access_array[i] = GetLocalIndex(access_array[i]);
}


void local_data::RenumberOffsetArray(int array_size, int* offset_array,
                                     int* lower_bound){

	for (int i = 0; i < array_size; i++){

		int curr_index = offset_array[i];
		int curr_posn = binary_search(direct_access_array, direct_access_size,
		                              curr_index);
		if (curr_posn == -1){
			offset_array[i] = -1;
		}
		else
			offset_array[i] = curr_posn - lower_bound[i];
	}
}


/**
 * \brief Calculate the ghosts of this local data on this process
 *
 * A ghost is a position in a local data that is owned by another process.
 * Therefore, if we want to use it locally, it must be communicated by the
 * owner first.
 */
void local_data::GenerateGhosts(){

	assert(local_array_size == -1);
	if (!is_read_only){
		assert(data_net_info);

		// Various initializations
		global_ghosts = new set<int>[nparts];
		global_owned = new set<int>[nparts];
		local_array_size = direct_access_size + indirect_access_size;
		assert(l_to_g == NULL);
		l_to_g = new int[local_array_size];
		for (int i = 0; i < local_array_size; i++)
			l_to_g[i] = -1;

		int total_counter = 0;
		for (int i = 0; i < direct_access_size; i++, total_counter++){
			int curr_index = direct_access_array[i];
			const net* curr_net = data_net_info[curr_index];
			l_to_g[total_counter] = curr_index;
		}
		for (int i = 0; i < indirect_access_size; i++, total_counter++){
			int curr_index = indirect_access_array[i];
			const net* curr_net = data_net_info[curr_index];
			l_to_g[total_counter] = curr_index;
		}

		if (is_constrained){
			block_owned_offset = 0;
			for (int i = 0; i < local_array_size; i++)
				if (l_to_g[i] == -1)
					block_owned_offset++;
				else
					break;
		}
	}
}

void local_data::GenerateOwned(){

	if (is_read_only)
		return;

	ghosts_offset[0] = 0;
	owned_offset[0] = 0;
	for (int i = 0; i < nparts; i++){
		ghosts_offset[i + 1] = ghosts_offset[i] + global_ghosts[i].size();
		owned_offset[i + 1] = owned_offset[i] + global_owned[i].size();
	}

	ghosts = new int[ghosts_offset[nparts]];
	owned = new int[owned_offset[nparts]];

	int counter = 0;
	for (int i = 0; i < nparts; i++)
		for (set<int>::iterator it = global_ghosts[i].begin();
		     it != global_ghosts[i].end(); it++)
			ghosts[counter++] = GetLocalIndex(*it);

	assert(counter == ghosts_offset[nparts]);

	counter = 0;
	for (int i = 0; i < nparts; i++)
		for (set<int>::iterator it = global_owned[i].begin();
		     it != global_owned[i].end(); it++)
			owned[counter++] = GetLocalIndex(*it);

	assert(counter == owned_offset[nparts]);

	delete[] global_owned;
	delete[] global_ghosts;
}


/**
 * \brief Constructor
 *
 * \param mn Identifier of the global array
 * \param np The original array is split into these many parts for comms
 * \param pid Identifier of the current process
 * \param st Size of each position in the array (e.g. int = 4 bytes)
 * \param dni Nets corresponding to all the positions in the array
 * \param oas Size of the original array
 * \param iro True if the data is read-only
 * \param ic Unused in quake
 */
local_data_double::local_data_double(int mn, int np, int pid, int st,
                                     const net** dni, int oas, bool iro,
                                     bool ic):

	local_data(mn, np, pid, st / sizeof(double), dni, oas, iro, ic){

	orig_array = NULL;
	local_array = NULL;
}


local_data_double::~local_data_double(){}


/**
 * \brief Populate local array with the contents of the global array
 *
 * \param local_base Allocated clean array to be populated
 * \param orig Original array
 * \param st Stride. See populate_local_array() in cdefs.cpp for more info.
 */
void local_data_double::PopulateLocalArray(const global_data* globalArray,
                                           double* local_base, double* orig,
                                           int st){

	assert(stride == st && orig_array == NULL && local_array == NULL);
	orig_array = orig;
	local_array = local_base;

	std::map<int, int> posMap; // Mapping from global index to local index
	int counter = 0;

	if (stride != 1){

		for (int j = 0; j < direct_access_size; j++, counter++){
			int orig_offset = direct_access_array[j] * stride;
			for (int i = 0; i < stride; i++){
				int localIndex = counter * stride + i;
				int globalIndex = orig_offset + i;
				local_array[localIndex] = orig_array[globalIndex];
				posMap[globalIndex] = localIndex;
			}
		}

		for (int j = 0; j < indirect_access_size; j++, counter++){
			int orig_offset = indirect_access_array[j] * stride;
			for (int i = 0; i < stride; i++){
				int localIndex = counter * stride + i;
				int globalIndex = orig_offset + i;
				local_array[localIndex] = orig_array[globalIndex];
				posMap[globalIndex] = localIndex;
			}
		}
	}
	else{
		for (int i = 0; i < direct_access_size; i++, counter++){
			local_array[counter] = orig_array[direct_access_array[i]];
			posMap[direct_access_array[i]] = counter;
		}
		for (int i = 0; i < indirect_access_size; i++, counter++){
			local_array[counter] = orig_array[indirect_access_array[i]];
			posMap[indirect_access_array[i]] = counter;
		}
	}

	// Convert comm info from global to local indexes; register in local array.

	// Send info
	maxSendCounts = 0;

	// For all iterations in our loop
	for (global_data::IdxsPerProcPerIter::const_iterator
		     iterIt = globalArray->pipeSendIndexes.begin(),
		     iterEnd = globalArray->pipeSendIndexes.end();
	     iterIt != iterEnd;
	     ++iterIt){

		int iter = iterIt->first;
		const std::map<int, std::vector<int> > idxsPerProc = iterIt->second;

		// For all consumer processes
		int count = 0;
		for (std::map<int, std::vector<int> >::const_iterator
			     procIt = idxsPerProc.begin(),
			     procEnd = idxsPerProc.end();
		     procIt != procEnd;
		     ++procIt){

			int proc = procIt->first;
			std::vector<int> idxs = procIt->second;
			for (std::vector<int>::const_iterator idxIt = idxs.begin(),
				     idxEnd = idxs.end();
			     idxIt != idxEnd;
			     ++idxIt)

				if (posMap.find(*idxIt) != posMap.end()){
					++count;
					pipeSendIndexes[iter][proc].push_back(posMap[*idxIt]);
				}
		}

		if (count > maxSendCounts)
			maxSendCounts = count;
	}

	// Receive info
	maxRecvCounts = 0;

	// For all iterations in our loop
	for (global_data::IdxsPerProcPerIter::const_iterator
		     iterIt = globalArray->pipeRecvIndexes.begin(),
		     iterEnd = globalArray->pipeRecvIndexes.end();
	     iterIt != iterEnd;
	     ++iterIt){

		int iter = iterIt->first;
		const std::map<int, std::vector<int> > idxsPerProc = iterIt->second;

		// For all producer processes
		int count = 0;
		for (std::map<int, std::vector<int> >::const_iterator
			     procIt = idxsPerProc.begin(),
			     procEnd = idxsPerProc.end();
		     procIt != procEnd;
		     ++procIt){

			int proc = procIt->first;
			std::vector<int> idxs = procIt->second;
			for (std::vector<int>::const_iterator idxIt = idxs.begin(),
				     idxEnd = idxs.end();
			     idxIt != idxEnd;
			     ++idxIt)

				if (posMap.find(*idxIt) != posMap.end()){
					++count;
					pipeRecvIndexes[iter][proc].push_back(posMap[*idxIt]);
				}
		}

		if (count > maxRecvCounts)
			maxRecvCounts = count;
	}
}


/**
 * \brief Populate send buffer with the data owned by this processor
 *
 * This data will be communicated to other user processes.
 *
 * \param send_buffer The buffer to be populated.
 * \param offset Array containing the start address, in the send
 *        buffer, of the chunk going to each user process.
 */
void local_data_double::SendOwnedData(char* send_buffer,
                                      const int* offset){

	for (int i = 0; i < nparts; i++)
		if (i != proc_id){

			int curr_offset = offset[i];
			double* buffer =
				reinterpret_cast<double*>(send_buffer + curr_offset);
			int counter = 0;
			if (stride == 1)
				for (int j = owned_offset[i]; j < owned_offset[i + 1];
				     j++, counter++)

					buffer[counter] = local_array[owned[j]];
			else
				for (int j = owned_offset[i]; j < owned_offset[i + 1];
				     j++, counter++)

					for (int k = 0; k < stride; k++)
						buffer[counter * stride + k] =
							local_array[owned[j] * stride + k];

			curr_offset += counter * stride * sizeof(double);
		}
}


/**
 * \brief Write ghosts from a receive buffer to the local data.
 *
 * Updated data owned by other processes is received in the receive
 * buffer and copied to the local data.
 *
 * \param recv_buffer The receive buffer to be copied to our local data
 * \param offset Array containing the start address, in the receive
 *        buffer, of the chunk coming from each owner process.
 */
void local_data_double::RecvGhostData(char* recv_buffer,
                                      const int* offset){

	for( int i = 0 ; i < nparts ; i++ )
		if( i != proc_id ){

			int curr_offset = offset[i];
			double* buffer =
				reinterpret_cast<double*>(recv_buffer + curr_offset);
			int counter = 0;
			if( stride == 1 )
				for( int j = ghosts_offset[i] ; j < ghosts_offset[i+1] ; j++ , counter++){
					local_array[ghosts[j]] = buffer[counter];
				}
			else
				for( int j = ghosts_offset[i] ; j < ghosts_offset[i+1] ; j++ , counter++ ){
					for( int k = 0 ; k < stride ; k++ )
						local_array[ghosts[j]*stride+k] = buffer[counter*stride+k];
				}
			curr_offset += counter * stride * sizeof(double);
		}
}


/**
 * \brief Write ghosts to a send buffer.
 *
 * After this function, the send buffer contains the positions of this
 * local data owned by other processess but computed (maybe partially)
 * by this process. This data must be communicated to their owners so
 * that these can update the data and broadcast it to other user
 * processors.
 *
 * \param send_buffer The send buffer to be populated
 * \param offset Contains the start position of the ghosts to be
 *        sent to each owner process
 */
void local_data_double::SendGhostData(char* send_buffer,
                                      const int* offset){

	// For all the parts of this local data, each of which is owned
	// by one process...
	for( int i = 0 ; i < nparts ; i++ )
		if( i != proc_id ){

			int curr_offset = offset[i];
			double* buffer =
				reinterpret_cast<double*>(send_buffer + curr_offset);
			int counter = 0;
			if( stride == 1 )
				for( int j = ghosts_offset[i] ; j < ghosts_offset[i+1] ; j++ , counter++){
					buffer[counter] = local_array[ghosts[j]];
				}
			else
				for( int j = ghosts_offset[i] ; j < ghosts_offset[i+1] ; j++ , counter++ ){
					for( int k = 0 ; k < stride ; k++ )
						buffer[counter*stride+k] = local_array[ghosts[j]*stride+k];
				}
			curr_offset += counter * stride * sizeof(double);
		}
}

/**
 * \brief Get received data that we own and update our copy
 *
 * Other processes must communicate us any change that they make to
 * any data owned by us, so that we can aggregate all the partial
 * results from other processes. This function copies the received
 * buffer to our local copy, which is the important one.
 *
 * \param send_buffer The receive buffer where preliminary results
 *        are located.
 * \param offset Array containing the start address, in the receive
 *        buffer, of the chunk coming from each process that used
 *        our data.
 */
void local_data_double::RecvOwnedData(char* recv_buffer, const int* offset){

	for (int i = 0; i < nparts; i++)
		if (i != proc_id){

			int curr_offset = offset[i];
			double* buffer =
				reinterpret_cast<double*>(recv_buffer + curr_offset);
			int counter = 0;
			if (stride == 1)
				for (int j = owned_offset[i]; j < owned_offset[i + 1];
				     j++, counter++)

					local_array[owned[j]] += buffer[counter];
			else
				for (int j = owned_offset[i]; j < owned_offset[i + 1];
				     j++, counter++)
					for( int k = 0 ; k < stride ; k++ )

						local_array[owned[j] * stride + k] +=
							buffer[counter * stride + k];
			curr_offset += counter * stride * sizeof(double);
		}
}


void local_data_double::InitWriteGhosts(){

	for (int i = 0; i < nparts; i++){
		if (stride == 1)
			for (int j = ghosts_offset[i]; j < ghosts_offset[i + 1]; j++)
				local_array[ghosts[j]] = 0.0;
		else
			for (int j = ghosts_offset[i]; j < ghosts_offset[i + 1]; j++)
				for (int k = 0; k < stride; k++)
					local_array[ghosts[j] * stride + k] = 0.0;
	}
}


void local_data_double::PopulateLocalGhosts
(local_data* source_array, int source){

	local_data_double* source_array_double =
		reinterpret_cast<local_data_double*>(source_array);
	int* source_owned_offset = source_array_double->owned_offset;
	int* source_owned = source_array_double->owned;
	double* source_local_array = source_array_double->local_array;

	if (stride == 1)
		for (int j = ghosts_offset[source], k = source_owned_offset[0];
		     j < ghosts_offset[source + 1]; j++, k++){

			local_array[ghosts[j]] = source_local_array[source_owned[k]];
		}
	else
		for (int j = ghosts_offset[source], k = source_owned_offset[0];
		     j < ghosts_offset[source + 1]; j++, k++){
			for (int l = 0; l < stride; l++)
				local_array[ghosts[j] * stride + l] =
					source_local_array[source_owned[k] * stride + l];
		}
}

void local_data_double::UpdateLocalOwned(local_data* source_array, int source){

	local_data_double* source_array_double =
		reinterpret_cast<local_data_double*>(source_array);
	int* source_ghosts_offset = source_array_double->ghosts_offset;
	int* source_ghosts = source_array_double->ghosts;
	double* source_local_array = source_array_double->local_array;
	if (stride == 1)
		for (int j = owned_offset[source], k = source_ghosts_offset[0];
		     j < owned_offset[source + 1]; j++, k++)
			local_array[owned[j]] += source_local_array[source_ghosts[k]];
	else
		for (int j = owned_offset[source], k = source_ghosts_offset[0];
		     j < owned_offset[source + 1]; j++, k++)
			for (int l = 0; l < stride; l++)
				local_array[owned[j] * stride + l] +=
					source_local_array[source_ghosts[k] * stride + l];
}

void local_data_double::PopulateGlobalArray(){

	if (!is_read_only && !is_constrained){
		memset((void*)orig_array,0,orig_array_size*stride*sizeof(double));
		for( int i = 0 ; i < local_array_size ; i++ ){
			int global_index = l_to_g[i];
			if( global_index != -1 )
				for( int j = 0 ; j < stride ; j++ )
					orig_array[global_index*stride+j] = local_array[i*stride+j];
		}
		MPI_Allreduce(MPI_IN_PLACE, orig_array, orig_array_size * stride,
		              MPI_DOUBLE, MPI_SUM, global_comm::team_communicator);
	}
}


int local_data_double::pipe_populate_send_buf(int iter, int proc, char* buf){

	std::vector<int> indexes = pipeSendIndexes[iter][proc];
	for (std::vector<int>::iterator idxIt = indexes.begin(),
		     idxEnd = indexes.end();
	     idxIt != idxEnd;
	     ++idxIt)

		memcpy(buf, local_array + *idxIt, sizeof(double) * stride);
}


void local_data_double::pipe_update(int iter, int proc, char* buf){

	std::vector<int> indexes = pipeRecvIndexes[iter][proc];
	for (std::vector<int>::iterator idxIt = indexes.begin(),
		     idxEnd = indexes.end();
	     idxIt != idxEnd;
	     ++idxIt)

		memcpy(local_array + *idxIt, buf, sizeof(double));
}
