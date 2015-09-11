//
//  hmm.h
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/27/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#ifndef __markov_clustering__hmm__
#define __markov_clustering__hmm__

#include "markov_clustering.h"

typedef struct hmm_t{
    double **A;
    double *pi;
    double **B;
    uint32_t num_states;
    uint32_t K; // Needed for the temporal version, it's the number of states per time unit
}*hmm;

typedef struct forwards_backwards_t{
    double **alpha;
    double **beta;
    double *c;
    double **gamma;
    double ***chi;
    uint32_t seq_length;
    uint32_t num_states;
}*forwards_backwards;

typedef struct baum_welch_t{
    struct hmm_t *model;
    forwards_backwards *fb;
    uint32_t num_states;
    uint64_t num_seq;
    uint32_t seq_length;
    qv_file qv_seqs;
}*baum_welch;


hmm alloc_hmm(uint32_t num_states);
forwards_backwards alloc_fb(uint32_t seq_length, uint32_t num_states);
baum_welch alloc_bw(qv_file qv_seqs, uint32_t num_states);
void initialize_hmm(hmm h);
void initialize_bw(baum_welch bw);
void forwards_backwards_algorithm(forwards_backwards fb, hmm model, uint8_t *seq);

void e_step_bw(baum_welch bw);
void m_step_bw(baum_welch bw);
double compute_nll_bw(baum_welch bw);

void perform_bw(qv_file qv_f, uint32_t num_states, uint32_t num_iters, FILE *fo, const char *fgraph);

void output_model(hmm h, FILE *fo);
void print_graph_hmm(hmm m, const char* graph_file_path);

void perform_bw_temporal(qv_file qv_f, uint32_t num_states_t, uint32_t num_iters, FILE *fo, const char *fgraph);

#endif /* defined(__markov_clustering__hmm__) */
