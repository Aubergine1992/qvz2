//
//  main.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/19/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#include "em_markov.h"
#include "hmm.h"

void compute_model_entropy(em_markov em, uint32_t num_model){
    
    
    int i = 0, j = 0, num_ctx = 1, l = 0, col_idx = 0, seq_idx = 0;
    double p = 0;
    int ctx = 0;
    long long unsigned int counts = 0;
    uint8_t *qv_ptr;
    
    int N_CTX = em->models[num_model].num_nodes;
    //      printf("%d\n",N_CTX);
    long long unsigned int **P = (long long unsigned int**) calloc(N_CTX, sizeof(long long unsigned int*));
    
    for(i = 0; i < N_CTX; ++i)
        P[i] = (long long unsigned int*)calloc(em->models->num_nodes+1, sizeof(long long unsigned int));
    
    double *H = (double*) calloc(N_CTX, sizeof(double));
    
    
    //      printf("Computing the Probability Distribution\n");
    qv_ptr = em->qv_seqs->file_head;
    
    for(seq_idx = 0; seq_idx < em->num_seq; ++seq_idx){
        
        if (em->clust->clusters[seq_idx] == num_model) {
            
            //              printf("N: %d\n", N);
            for(col_idx = 0; col_idx < em->qv_seqs->read_length; ++col_idx){
                //                      printf("qv(%d) %d\n",col, qv);
                i = 0;
                for(j = 1; j <= num_ctx; j++){
                    ctx = col_idx - j < 0 ? 0:*(qv_ptr - j);
                    //                              printf("ctx %d\n", ctx);
                    i += ((ctx & 0xff) << (8*l));
                    //                              printf("i %d\n", i);
                }
                //                      printf("%d %d\n", i, qv);
                P[i][*qv_ptr]++;
                ++P[i][em->models->num_nodes+1];
                qv_ptr++;
            }
            qv_ptr++;// TODO: remove. Necessary now because of the zero between sequences.
            
        }
        else{
            qv_ptr += em->qv_seqs->read_length;
            qv_ptr++;// TODO: remove. Necessary now because of the zero between sequences.
        }
        
    }
    
    for(i = 0; i < N_CTX; ++i){
        for(j = 0; j < em->models->num_nodes; ++j){
            if( (counts = P[i][j]) == 0)
                continue;
            if (P[i][em->models->num_nodes+1] <=0)
                printf("counts is less than 0 in context: %d\n", i);
            
            p = (double)counts/(double)P[i][em->models->num_nodes+1];
            if(p <= 0)
                printf("fooooo %d %d %llu %llu\n", i, j, P[i][42], counts);
            H[N_CTX] += -(counts*log2(p));
        }
    }
    
    printf("Size(%d) = %f Bits\n", num_ctx, H[0]);
    printf("Size(%d) = %f MB\n", num_ctx, (double)H[0]/8.0/1000000.0);

}

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





int main(int argc, const char * argv[]) {
    
    
    const char* path;
    
    FILE *fo;
    
    path = argv[1];
    
    qv_file qv_f = load_file(path, -1);
    
    fo = fopen(argv[2], "w");
    
    uint32_t K = atoi(argv[3]);
    
    uint32_t num_iters = atoi(argv[4]);
    
    //fgraph = fopen(argv[5], "w");
    
    perform_bw_temporal(qv_f, K, num_iters, fo, argv[5]);
    //perform_em_markov(qv_f, K, num_iters, fo, fgraph);
    
    return 0;
}



