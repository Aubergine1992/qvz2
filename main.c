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



qv_file load_file(const char *path, uint64_t max_lines){
    
    uint32_t line_idx;
    int ch;
    char line[READ_LINEBUF_LENGTH];
    FILE *fp;
    struct stat finfo;
    uint8_t* current_qv;
    
    qv_file my_qv_file = (struct qv_file_t*)malloc(sizeof(struct qv_file_t));
    
    //printf("%s\n", path);
    
    // Open the file
    fp = fopen(path, "rt");
    if (!fp) {
        return NULL;
    }
    
    
    
    // Use the first line to figure out the read length
    fgets(line, READ_LINEBUF_LENGTH, fp);
    my_qv_file->read_length = (uint32_t)strlen(line) - 1;
    if (my_qv_file->read_length > MAX_READS_PER_LINE) {
        fclose(fp);
        return NULL;
    }
    
    rewind(fp);
    
    // Figure out how many lines we'll need depending on whether we were limited or not
    stat(path, &finfo);
    my_qv_file->lines = finfo.st_size / ((uint64_t) (my_qv_file->read_length+1));
    if (max_lines > 0 && my_qv_file->lines > max_lines) {
        my_qv_file->lines = max_lines;
    }
    
    // Right now we put a zero between sequences. This helps for the debbuging but it is not necessary and
    // is a waste of memory (not much, though).
    my_qv_file->file_head = (uint8_t*)calloc(my_qv_file->lines*(my_qv_file->read_length+1), sizeof(uint8_t));
    if (my_qv_file->file_head == NULL)
        return NULL;
    
    current_qv = my_qv_file->file_head;
    line_idx = 0;
    while (line_idx < my_qv_file->lines ) {
        while ((ch = getc(fp)) != '\n') {
            *current_qv = qv2ch(ch);
            current_qv++;
        }
        line_idx++;
        current_qv++;// TODO: remove. This introduce a zero after the sequence
    }
    
    
    
    fclose(fp);
    
    return my_qv_file;
    
}

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
    
    
    const char* path;
    
    FILE *fo;
    
    path = argv[1];
    
    qv_file qv_f = load_file(path, -1);
    
    fo = fopen(argv[2], "w");
    
    uint32_t K = atoi(argv[3]);
    
    uint32_t num_iters = atoi(argv[4]);
    
    //perform_bw_temporal(qv_f, K, num_iters, fo, argv[5]);
    //perform_bw(qv_f, K, num_iters, fo, argv[5]);
    //perform_em_markov(qv_f, K, num_iters, fo, argv[5]);
    perform_em_temporal_markov(qv_f, K, num_iters, fo, argv[5]);
    
    return 0;
}



