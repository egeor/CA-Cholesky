#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <upc.h>
#include <upc_coll_mpi.h>
#include "cholesky.h"

typedef shared[] double* mytype_ptr;

shared mytype_ptr distr_A[THREADS];

int main(int argc, char **argv) {

	
	const int matrixDim = atoi(argv[1]);
	const int blockDim = atoi(argv[2]);
	const int c_rep = atoi(argv[3]);
	const int big_blockDim = atoi(argv[4]);
	
	const int num_blocks_dim = matrixDim/blockDim;
	const int num_pes_dim = sqrt(THREADS/c_rep);
	
	const int my_num_blocks_dim = num_blocks_dim/num_pes_dim;
	const int num_big_blocks_dim = matrixDim/big_blockDim;
	const int num_small_in_big_blk = big_blockDim/blockDim;
	const int my_num_small_in_big_blk = num_small_in_big_blk/num_pes_dim;
	
	const int threads_per_layer = THREADS/c_rep;
	
	upccoll_team_t my_row_team, my_column_team, my_peers_team, my_layer_team;
	int row_rank, column_rank, my_layer_rank, my_layer_id, my_row_rank, my_column_rank, my_peer_rank, layer_rank;
	
    
    /* Create teams per layer */
    upccoll_team_split(UPC_TEAM_ALL, MYTHREAD/threads_per_layer, MYTHREAD%threads_per_layer, &my_layer_team);
    upccoll_team_rank(my_layer_team, &my_layer_rank);
    
    /* Create teams for intra-layer broadcasts */
    upccoll_team_split(my_layer_team, my_layer_rank/num_pes_dim, my_layer_rank%num_pes_dim, &my_row_team);
    upccoll_team_rank(my_row_team, &my_row_rank);
    
    upccoll_team_split(my_layer_team, my_layer_rank%num_pes_dim, my_layer_rank/num_pes_dim, &my_column_team);
    upccoll_team_rank(my_column_team, &my_column_rank);
	
    
    /* Create teams for inter-layer broadcasts (peers per layer) */
    upccoll_team_split(UPC_TEAM_ALL, MYTHREAD%threads_per_layer, MYTHREAD/threads_per_layer, &my_peers_team);
    upccoll_team_rank(my_peers_team, &my_peer_rank);
	
    int layer_rank = my_layer_rank;
	double * priv_A;
	
	int i_big, width, blocks_to_sub, blocks_to_proceed ;
		
	/* Count how many elements each process should allocate, 2 cases: 
	 (1) If a process is above diagonal, then it owns (my_num_blocks_dim * (my_num_blocks_dim-1))/2 blocks  
	 (2) If a process is on/below diagonal then it owns (my_num_blocks_dim * (my_num_blocks_dim+1))/2 blocks  */
	
	received_factorized_block[MYTHREAD] = (mytype_ptr) upc_alloc(chunk_dim * chunk_dim * sizeof(double));
	
	int above_diag = (my_layer_rank/num_pes_dim) < (my_layer_rank%num_pes_dim);
	if (above_diag) { // case (1)
		distr_A[MYTHREAD] = (mytype_ptr) upc_alloc(((my_num_blocks_dim * (my_num_blocks_dim-1)) * blockDim * blockDim / 2)* sizeof(double));
		priv_A = (double *) distr_A[MYTHREAD];
	} else { // case(2)
		distr_A[MYTHREAD] = (mytype_ptr) upc_alloc(((my_num_blocks_dim * (my_num_blocks_dim+1)) * blockDim * blockDim / 2)* sizeof(double));
		priv_A = (double *) distr_A[MYTHREAD];
	}
	
	/* Create symmetric positive definite (Lehmer) matrix and store lower half */
	/* TODO */
	
	/* Replicate input matrix A */
	/* TODO */
	 
	for (i_big=0; i_big < num_big_blocks_dim; i_big++) {
		/* Offset private A by the nymber of big blocks we have already factorized */
		width = i_big * my_num_small_in_big_blk - 1 + above_diag ;
		blocks_to_sub = width * (1+width) / 2 ;
		blocks_to_proceed = my_num_blocks_dim * i_big * my_num_small_in_big_blk - blocks_to_sub;
		
		
		
	}
	 
	 
	
}