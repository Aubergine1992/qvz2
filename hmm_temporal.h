//
//  hmm_temporal.h
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 9/7/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#ifndef markov_clustering_hmm_temporal_h
#define markov_clustering_hmm_temporal_h

#include "markov_clustering.h"


void perform_bw_temporal(qv_file qv_f, uint32_t num_states, uint32_t num_iters, FILE *fo, const char *fgraph);

/* defined(__markov_clustering__hmm__) */


#endif
