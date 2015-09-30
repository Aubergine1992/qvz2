//
//  main.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/19/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#include "em_markov.h"
#include "em_markov_temporal.h"
#include "hmm.h"
#include "codebook.h"
#include "cluster.h"





// *******************************************************************************//
void split_data(const char* path_head, clusters clust, qv_file  qv_seqs){
    uint32_t i = 0, cluster_id, t;
    FILE **fo;
    char path_buffer[1024];
    char cluster_name[256] = "/cluster";
    
    
    fo = (FILE**)calloc(clust->num_clusters, sizeof(FILE*));
    
    for (i = 0; i < clust->num_clusters; i++) {
        cluster_name[8] = (i+48);
        strcpy(path_buffer, path_head);
        strcat(path_buffer,cluster_name);
        fo[i] = fopen(path_buffer, "w");
    }
    
    uint8_t *seq = qv_seqs->file_head;
    
    for (i = 0; i < qv_seqs->lines; i++) {
        cluster_id = clust->clusters[i];
        for (t = 0; t < qv_seqs->read_length; t++) {
            fputc((*seq)+33, fo[cluster_id]), seq++;
        }
        fputc('\n', fo[cluster_id]);
        seq++;
    }
    
    
}

// *******************************************************************************//

void compute_clusters(clusters clust, double ** r, uint64_t num_seq){
    
    uint32_t seq_idx, i = 0, tmp_model;
    
    double tmp = 0;
    
    for (i = 0; i < clust->num_clusters; i++) {
        clust->cluster_sizes[i] = 0;
    }
    
    for (seq_idx = 0; seq_idx < num_seq; seq_idx++) {
        tmp = 0;
        tmp_model = 0;
        for (i = 0; i < clust->num_clusters; i++) {
            if (r[i][seq_idx] >= tmp) {
                tmp = r[i][seq_idx];
                tmp_model = i;
            }
        }
        clust->cluster_sizes[tmp_model]++;
        clust->clusters[seq_idx] = tmp_model;
    }
}


// *******************************************************************************//

int main(int argc, const char * argv[]) {
    
    uint32_t status;
    const char* path;
    
    FILE *fo;
    
    path = argv[1];
    
    struct quality_file_t qv_info;
    
    // Load input file all at once
    qv_file qv_f = load_file(path, -1);
    status = generate_qv_struct(qv_f, &qv_info, 0);
    if (status != LF_ERROR_NONE) {
        printf("load_file returned error: %d\n", status);
        exit(1);
    }
    
    
    
    fo = fopen(argv[2], "w");
    
    uint32_t K = atoi(argv[3]);
    
    uint32_t num_iters = atoi(argv[4]);
    
    //perform_bw_temporal(qv_f, K, num_iters, fo, argv[5]);
    //perform_bw(qv_f, K, num_iters, fo, argv[5]);
    //perform_em_markov(qv_f, K, num_iters, fo, argv[5]);
    perform_em_temporal_markov(qv_f, K, num_iters, fo, argv[5]);
    
    return 0;
}



