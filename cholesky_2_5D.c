#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <upc.h>
#include <upc_coll_mpi.h>
#include "cholesky.h"


/* Execution parameters */
#define threads 256 
#define threads_per_layer 64
#define sqrt_threads_per_layer 8
#define c 4
#define big_block_num
#define total_small_block_num

/*	
 	Derived quantities will be:
	(*) block_size = n/total_small_block_num
 	(*) panel_size = n/big_block_num
 	(*) blocks_per_panel = total_small_block_num/big_block_num
 	(*) blocks_per_thread_per_panel = blocks_per_panel/sqrt_threads_per_layer
 
	(*) THINK: Define panel size in terms of number of small blocks
 */

/* TODO: Add timers */

shared mytype_ptr A[big_block_num][big_block_num][c][sqrt_threads_per_layer][sqrt_threads_per_layer];

shared mytype_ptr global_Schur_complement[big_block_num][big_block_num][c][sqrt_threads_per_layer][sqrt_threads_per_layer];
shared mytype_ptr column_buffer[THREADS];
shared mytype_ptr received_factorized_block[THREADS];
shared mytype_ptr received_factorized_column[THREADS];
shared mytype_ptr received_from_diagonal[THREADS];
shared mytype_ptr reduction_buffer[THREADS];

int main(int argc, char **argv) {

	const int sqrt_tpl = sqrt_threads_per_layer;
	
	int i, j, k, l, layer, big_index, small_index_row, small_index_column, move_ptr;
	
	mytype_ptr local_A[big_block_num][big_block_num][c][sqrt_threads_per_layer][sqrt_threads_per_layer];
	mytype_ptr Schur_complement[big_block_num][big_block_num][c][sqrt_threads_per_layer][sqrt_threads_per_layer];
	
	upccoll_initialize(&argc, &argv);

	
	/* Read input parameters */
	/* The input matrix is n x n */
	const int n = atoi(argv[1]);
	/* The panel dim specifies the number of small blocks per "fat" panel */
	const int panel_dim = atoi(argv[2]);
	
    const int chunk_dim = n/(big_block_num*sqrt_tpl);
    const int panel_dime = chunk_dim * panel_dim;
	
	/* Create local directory */  
    
    for (i=0; i<big_block_num; i++)
        for (j=0; j<big_block_num; j++)
            for (k=0; k<sqrt_tpl; k++)
                for (l=0; l<sqrt_tpl; l++)
                    for (layer = 0; layer < c; layer++) {
                        owner_directory[i][j][layer][k][l] = (int) upc_threadof(&A[i][j][layer][k][l]);
                    }
    
    
    /* Memory allocation */
    
    for (big_index = 0; big_index < big_block_num; big_index++) {
        for (small_index_row = 0; small_index_row < sqrt_tpl; small_index_row++) {
            for (small_index_column = 0; small_index_column < sqrt_tpl; small_index_column++) {
                for (layer = 0; layer < c; layer++) {
                    /* Allocation of columns for elements on or below diagonal */
                    if (small_index_column <= small_index_row) {
                        if (MYTHREAD == owner_directory[big_index][big_index][layer][small_index_row][small_index_column]) {
                            column_buffer[MYTHREAD] = (mytype_ptr) upc_alloc(chunk_dim * chunk_dim * (big_block_num-big_index) * sizeof(double));
                            for (I=0; I<chunk_dim * chunk_dim * (big_block_num-big_index); I++)
                                *(column_buffer[MYTHREAD]+I)=0;
                            /* make appropriate pointer assignments */
                            move_ptr = 0;
                            for (column_part = big_block_num-1; column_part >= big_index; column_part--) {
                                A[column_part][big_index][layer][small_index_row][small_index_column] = column_buffer[MYTHREAD]+move_ptr*(chunk_dim * chunk_dim);
                                move_ptr++;
                            }
                        }
                    }
					
                    /* Allocation of columns for elements above diagonal */
                    if (small_index_column > small_index_row) {
                        if (MYTHREAD == owner_directory[big_index][big_index][layer][small_index_row][small_index_column]) {
                            column_buffer[MYTHREAD] = (mytype_ptr) upc_alloc(chunk_dim * chunk_dim * (big_block_num-big_index-1) * sizeof(double));
                            for (I=0; I<chunk_dim * chunk_dim * (big_block_num-big_index-1); I++)
                                *(column_buffer[MYTHREAD]+I)=0;
                            /* make appropriate pointer assignments */
                            move_ptr = 0;
                            for (column_part = big_block_num-1; column_part > big_index; column_part--) {
                                A[column_part][big_index][layer][small_index_row][small_index_column] = column_buffer[MYTHREAD]+move_ptr*(chunk_dim * chunk_dim);
                                move_ptr++;
                            }
                        }
                    }
                    
                    /* Allocation of Schur complement - meaningful for a part of the matrix */
                    if (big_index > 0) {
                        /* allocation of columns for elements on or below diagonal */
                        if (small_index_column <= small_index_row) {
                            if (MYTHREAD == owner_directory[big_index][big_index][layer][small_index_row][small_index_column]) {
                                column_buffer[MYTHREAD] = (mytype_ptr) upc_alloc(chunk_dim * chunk_dim * (big_block_num-big_index) * sizeof(double));
                                for (I=0; I<chunk_dim * chunk_dim * (big_block_num-big_index); I++)
                                    *(column_buffer[MYTHREAD]+I)=0;
                                /* make appropriate pointer assignments */
                                move_ptr = 0;
                                for (column_part = big_block_num-1; column_part >= big_index; column_part--) {
                                    global_Schur_complement[column_part][big_index][layer][small_index_row][small_index_column] = column_buffer[MYTHREAD]+move_ptr*(chunk_dim * chunk_dim);
                                    move_ptr++;
                                }
                            }
                        }
                        
                        /* allocation of columns for elements above diagonal */
                        if (small_index_column > small_index_row) {
                            if (MYTHREAD == owner_directory[big_index][big_index][layer][small_index_row][small_index_column]) {
                                column_buffer[MYTHREAD] = (mytype_ptr) upc_alloc(chunk_dim * chunk_dim * (big_block_num-big_index-1) * sizeof(double));
                                for (I=0; I<chunk_dim * chunk_dim * (big_block_num-big_index-1); I++)
                                    *(column_buffer[MYTHREAD]+I)=0;
                                /* make appropriate pointer assignments */
                                move_ptr = 0;
                                for (column_part = big_block_num-1; column_part > big_index; column_part--) {
                                    global_Schur_complement[column_part][big_index][layer][small_index_row][small_index_column] = column_buffer[MYTHREAD]+move_ptr*(chunk_dim * chunk_dim);
                                    move_ptr++;
                                }
                            }
                        }
						
                    }
                }
            }
        }
    }
    
    upc_barrier;
	
	my_layer_id  = MYTHREAD / threads_per_layer;

}