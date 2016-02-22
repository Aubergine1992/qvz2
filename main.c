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
    
    printf("Clustering reads....\n");
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
    
    printf("Clustering finished\n");
    
    // Then find stats and generate codebooks for each cluster
    start_timer(&stats);
    calculate_statistics(&qv_info);
    generate_codebooks(&qv_info);
    stop_timer(&stats);
    printf("Codebooks computed\n");
    
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

/**
 *
 */
void decode(char *input_file, char *output_file, struct qv_options_t *opts) {
    FILE *fin, *fout;
    struct hrtimer_t timer;
    struct quality_file_t qv_info;
    struct alphabet_t *A = alloc_alphabet(ALPHABET_SIZE);
    
    qv_info.alphabet = A;
    qv_info.opts = opts;
    
    start_timer(&timer);
    
    fin = fopen(input_file, "rb");
    fout = fopen(output_file, "wt");
    if (!fin || !fout) {
        perror("Unable to open input or output files");
        exit(1);
    }
    
    read_codebooks(fin, &qv_info);
    start_qv_decompression(fout, fin, &qv_info);
    
    fclose(fout);
    fclose(fin);
    stop_timer(&timer);
    
    if (opts->verbose) {
        printf("Decoded %llu lines in %f seconds.\n", qv_info.lines, get_timer_interval(&timer));
    }
}



// *******************************************************************************//

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(const char *name) {
    printf("Usage: %s (options) [input file] [output file]\n", name);
    printf("Options are:\n");
    printf("   -q           : Store quality values in compressed file (default)\n");
    printf("   -x           : Extract quality values from compressed file\n");
    printf("   -t           : Target average distortion, measured as specified by -d or -D (default 1)\n");
    printf("   -i           : Number of iterations of the clustering algorithm (default 50)\n");
    printf("   -d [M|L|A]   : Optimize for MSE, Log(1+L1), L1 distortions, respectively (default: MSE)\n");
    printf("   -D [FILE]    : Optimize using the custom distortion matrix specified in FILE\n");
    printf("   -c [#]       : Compress using [#] clusters (default: 1)\n");
    printf("   -u [FILE]    : Write the uncompressed lossy values to FILE (default: off)\n");
    printf("   -h           : Print this help\n");
    printf("   -s           : Print summary stats\n");
    printf("   -v           : Enable verbose output\n");
    printf("\nFor custom distortion matrices, a 72x72 matrix of values must be provided as the cost of reconstructing\n");
    printf("the x-th row as the y-th column, where x and y range from 0 to 71 (inclusive) corresponding to the possible\n");
    printf("Phred scores.\n");
}

// *******************************************************************************//

int main(int argc, const char * argv[]) {
    
    uint32_t status;
    const char* path;
    
    path = argv[1];
    
    struct quality_file_t qv_info;
    struct qv_options_t opts;
    
    double A[25]={1,2,3,4,5, 2,3,4,1,5, 6,6,4,8,7, 9,4,3,6,7, 3,5,7,1,2};
    double **c;
    int i,j,n=5;
    c = malloc((n)*sizeof(double *));
    for (i=0;i<n;i++){
        c[i] = malloc((n)*sizeof(double));
        for (j=0; j<n; j++) {
            c[i][j]=A[n*i+j];
        }
    }
    
    double **b;
    b = malloc((n)*sizeof(double *));
    for (i=0;i<n;i++)
        b[i] = malloc((n)*sizeof(double));
    
    uint8_t extract = 0;
    uint8_t file_idx = 0;
    
    const char *input_name = 0;
    const char *output_name = 0;
    
    //MxM(c, c, &A, 5, 5, 25);
    //Inverse(c, b, 5);
    
    // DEFAULT OPTIONS
    opts.verbose = 0;
    opts.stats = 0;
    opts.ratio = 1;
    opts.num_iters = 50;
    opts.distortion = DISTORTION_MSE;
    opts.cluster_threshold = 4;
    opts.mode = 1;
    opts.D = 1;
    opts.clusters = 1;
    extract = 0;
    //////////////////////////////////////
    
    
    // No dependency, cross-platform command line parsing means no getopt
    // So we need to settle for less than optimal flexibility (no combining short opts, maybe that will be added later)
    i = 1;
    while (i < argc) {
        // Handle file names and reject any other untagged arguments
        if (argv[i][0] != '-') {
            switch (file_idx) {
                case 0:
                    input_name = argv[i];
                    file_idx = 1;
                    break;
                case 1:
                    output_name = argv[i];
                    file_idx = 2;
                    break;
                default:
                    printf("Garbage argument \"%s\" detected.\n", argv[i]);
                    usage(argv[0]);
                    exit(1);
            }
            i += 1;
            continue;
        }
        
        // Flags for options
        switch(argv[i][1]) {
            case 'x':
                extract = 1;
                i += 1;
                break;
            case 'q':
                extract = 0;
                i += 1;
                break;
            case 'c':
                opts.clusters = atoi(argv[i+1]);
                i += 2;
                break;
            case 'v':
                opts.verbose = 1;
                i += 1;
                break;
            case 'h':
                usage(argv[0]);
                exit(0);
            case 's':
                opts.stats = 1;
                i += 1;
                break;
            case 'u':
                opts.uncompressed = 1;
                opts.uncompressed_name = argv[i+1];
                i += 2;
                break;
            case 't':
                opts.D = atof(argv[i+1]);
                i += 2;
                break;
            case 'i':
                opts.num_iters = atoi(argv[i+1]);
                i += 2;
                break;
            case 'd':
                switch (argv[i+1][0]) {
                    case 'M':
                        opts.distortion = DISTORTION_MSE;
                        break;
                    case 'L':
                        opts.distortion = DISTORTION_LORENTZ;
                        break;
                    case 'A':
                        opts.distortion = DISTORTION_MANHATTAN;
                        break;
                    default:
                        printf("Distortion measure not supported, using MSE.\n");
                        break;
                }
                i += 2;
                break;
            case 'D':
                opts.distortion = DISTORTION_CUSTOM;
                opts.dist_file = argv[i+1];
                i += 2;
                break;
            default:
                printf("Unrecognized option -%c.\n", argv[i][1]);
                usage(argv[0]);
                exit(1);
        }
    }
    
    if (file_idx != 2) {
        printf("Missing required filenames.\n");
        usage(argv[0]);
        exit(1);
    }
    
    if (opts.verbose) {
        if (extract) {
            printf("%s will be decoded to %s.\n", input_name, output_name);
        }
        else {
            printf("%s will be encoded as %s.\n", input_name, output_name);
            if (opts.mode == MODE_RATIO)
                printf("Ratio mode selected, targeting %f compression ratio.\n", opts.ratio);
            else if (opts.mode == MODE_FIXED)
                printf("Fixed-rate mode selected, targeting %f bits per symbol.\n", opts.ratio);
            else if (opts.mode == MODE_FIXED_MSE)
                printf("Fixed-MSE mode selected, targeting %f average distortion per context.\n", opts.ratio);
            
            switch (opts.distortion) {
                case DISTORTION_MSE:
                    printf("MSE will be used as a distortion metric.\n");
                    break;
                case DISTORTION_LORENTZ:
                    printf("log(1+L1) will be used as a distortion metric.\n");
                    break;
                case DISTORTION_MANHATTAN:
                    printf("L1 will be used as a distortion metric.\n");
                    break;
                case DISTORTION_CUSTOM:
                    printf("A custom distortion metric stored in %s will be used.\n", opts.dist_file);
                    break;
            }
            
            printf("Compression will use %d clusters, with a movement threshold of %.0f.\n", opts.clusters, opts.cluster_threshold);
        }
    }

    
    if (extract) {
        decode(input_name, output_name, &opts);
    }
    else{
        // Load input file all at once
        qv_file qv_f = load_file(path, -1);
        status = generate_qv_struct(qv_f, &qv_info, 0);
        if (status != LF_ERROR_NONE) {
            printf("load_file returned error: %d\n", status);
            exit(1);
        }
        qv_info.qv_f = qv_f;
        
        printf("File loaded into memory\n");
    
        //Generate the clusters and calculate the quantizers and compress
        encoding(qv_info, &opts, output_name);
    }
    
    return 0;
}



