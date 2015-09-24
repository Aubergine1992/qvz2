//
//  em_markov_temporal.h
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 9/24/15.
//  Copyright © 2015 Mikel Hernaez. All rights reserved.
//

#ifndef em_markov_temporal_h
#define em_markov_temporal_h

#include "markov_clustering.h"
#include "em_markov.h"

typedef struct markov_temporal_model_t{
    double ***A;
    double *pi;
    uint32_t num_nodes;
    uint32_t num_ts;
}*markov_temporal_model;

typedef struct em_temporal_markov_t{
    struct markov_temporal_model_t *models;
    double *models_prior;
    double **r;
    uint32_t num_models;
    uint64_t num_seq;
    qv_file qv_seqs;
    clusters clust;
}*em_temporal_markov;



// ******************************************************************//


uint32_t perform_em_temporal_markov(qv_file qv_f, uint32_t num_models, uint32_t iters, FILE *fo, const char *split_path);

#endif /* em_markov_temporal_h */
