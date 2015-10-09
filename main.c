//
//  main.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/19/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//


#include "codebook.h"
#include "em_markov.h"
#include "em_markov_temporal.h"
#include "qv_compressor.h"




// *******************************************************************************//
/*
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
*/
// *******************************************************************************//

/**
 * Allocate the memory used for the clusters based on the number wanted and column config
 */
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info) {
    uint8_t j;
    struct cluster_list_t *rtn = (struct cluster_list_t *) calloc(1, sizeof(struct cluster_list_t));
    
    // Allocate array of cluster structures
    rtn->count = info->cluster_count;
    rtn->clusters = (struct cluster_t *) calloc(info->cluster_count, sizeof(struct cluster_t));
    
    // Fill in each cluster
    for (j = 0; j < info->cluster_count; ++j) {
        rtn->clusters[j].id = j;
        rtn->clusters[j].count = 0;
        rtn->clusters[j].training_stats = alloc_conditional_pmf_list(info->alphabet, info->columns);
    }
    
    return rtn;
}

/**
 *
 */
void compute_clusters(struct quality_file_t *qv_info, double ** r, uint64_t num_seq){
    
    uint32_t seq_idx, i = 0, tmp_model, block_ctr = 0, line_ctr = 0;
    
    double tmp = 0;
    
    for (i = 0; i < qv_info->cluster_count; i++) {
        qv_info->clusters->clusters[i].count = 0;
        qv_info->clusters->clusters[i].id = i;
    }
    
    for (seq_idx = 0; seq_idx < num_seq; seq_idx++) {
        tmp = 0;
        tmp_model = 0;
        for (i = 0; i < qv_info->cluster_count; i++) {
            if (r[i][seq_idx] >= tmp) {
                tmp = r[i][seq_idx];
                tmp_model = i;
            }
        }
        qv_info->clusters->clusters[tmp_model].count++;
        
        qv_info->blocks[block_ctr].lines[line_ctr++].cluster = tmp_model;
        
        if (line_ctr == qv_info->blocks[block_ctr].count) {
            block_ctr++;
            line_ctr = 0;
        }
        
    }
}

/**
 *
 */
void encoding(struct quality_file_t qv_info, struct qv_options_t *opts, const char* output_name){
    
    struct distortion_t *dist;
    struct alphabet_t *alphabet = alloc_alphabet(ALPHABET_SIZE);
    struct hrtimer_t cluster_time, stats, total, encoding;
    uint64_t bytes_used;
    FILE *fout, *funcompressed = NULL;
    
    double distortion;
    
    double **clust_prob;
    
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
    qv_info.clusters = alloc_cluster_list(&qv_info);
    qv_info.opts = opts;
    
    // Do clustering
    start_timer(&cluster_time);
    //do_kmeans_clustering(&qv_info);
    
    // Set up clustering data structures
    clust_prob = perform_em_temporal_markov(qv_info.qv_f, opts->clusters, opts->num_iters);
    compute_clusters(&qv_info, clust_prob, qv_info.lines);
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
    
    // Note that we want \r\n translation in the input
    // but we do not want it in the output
    fout = fopen(output_name, "wb");
    if (!fout) {
        perror("Unable to open output file");
        exit(1);
    }
    
    if (opts->uncompressed) {
        funcompressed = fopen(opts->uncompressed_name, "w");
        if (!funcompressed) {
            perror("Unable to open uncompressed file");
            exit(1);
        }
    }
    
    start_timer(&encoding);
    write_codebooks(fout, &qv_info);
    bytes_used = start_qv_compression(&qv_info, fout, &distortion, funcompressed);
    stop_timer(&encoding);
    stop_timer(&total);
    
    printf("%llu\t%f\n",bytes_used,distortion);
}



// *******************************************************************************//

int main(int argc, const char * argv[]) {
    
    uint32_t status;
    const char* path;
    
    path = argv[1];
    
    struct quality_file_t qv_info;
    struct qv_options_t opts;
    
    // DEFAULT OPTIONS FOR THE MOMENT
    opts.verbose = 0;
    opts.stats = 0;
    opts.ratio = 1;
    opts.clusters = atoi(argv[3]);
    opts.num_iters = atoi(argv[4]);
    opts.uncompressed = 1;
    opts.uncompressed_name = argv[6];
    opts.distortion = DISTORTION_MSE;
    opts.cluster_threshold = 4;
    opts.mode = 1;
    opts.D = atof(argv[5]);
    //////////////////////////////////////
    
    // Load input file all at once
    qv_file qv_f = load_file(path, -1);
    status = generate_qv_struct(qv_f, &qv_info, 0);
    if (status != LF_ERROR_NONE) {
        printf("load_file returned error: %d\n", status);
        exit(1);
    }
    qv_info.qv_f = qv_f;
    
    //Generate the clusters and calculate the quantizers
    encoding(qv_info, &opts, argv[2]);
    
    return 0;
}



