//
//  precompression_processing.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 9/29/15.
//  Copyright © 2015 Mikel Hernaez. All rights reserved.
//

#include "codebook.h"
#include "cluster.h"

void preprocess_qvs(struct quality_file_t qv_info, struct qv_options_t *opts){
    
    struct distortion_t *dist;
    struct alphabet_t *alphabet = alloc_alphabet(ALPHABET_SIZE);
    struct hrtimer_t cluster_time, stats, total;
    
    start_timer(&total);
    
    if (opts->distortion == DISTORTION_CUSTOM) {
        dist = gen_custom_distortion(ALPHABET_SIZE, opts->dist_file);
    }
    else {
        dist = generate_distortion_matrix(ALPHABET_SIZE, opts->distortion);
    }
    
    qv_info.alphabet = alphabet;
    qv_info.dist = dist;
    qv_info.cluster_count = opts->clusters;
    
    // Set up clustering data structures
    qv_info.clusters = alloc_cluster_list(&qv_info);
    qv_info.opts = opts;
    
    // Do k-means clustering
    start_timer(&cluster_time);
    do_kmeans_clustering(&qv_info);
    stop_timer(&cluster_time);
    if (opts->verbose) {
        printf("Clustering took %.4f seconds\n", get_timer_interval(&cluster_time));
    }
    
    // Then find stats and generate codebooks for each cluster
    start_timer(&stats);
    calculate_statistics(&qv_info);
    generate_codebooks(&qv_info);
    stop_timer(&stats);
    
    if (opts->verbose) {
        printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));
        // @todo expected distortion is inaccurate due to lack of pmf
        //printf("Expected distortion: %f\n", opts->e_dist);
    }
}


