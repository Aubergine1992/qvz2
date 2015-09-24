//
//  em_markov.h
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/27/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#ifndef __markov_clustering__em_markov__
#define __markov_clustering__em_markov__

#include "markov_clustering.h"

typedef struct markov_model_t{
    double **A;
    double *pi;
    uint32_t num_nodes;
}*markov_model;

typedef struct em_markov_t{
    struct markov_model_t *models;
    double *models_prior;
    double **r;
    uint32_t num_models;
    uint64_t num_seq;
    qv_file qv_seqs;
    clusters clust;
}*em_markov;

uint32_t perform_em_markov(qv_file qv_f, uint32_t num_models, uint32_t iters, FILE *fo, const char *split_path);

#endif /* defined(__markov_clustering__em_markov__) */


