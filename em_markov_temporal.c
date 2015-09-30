//
//  em_markov_temporal.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 9/24/15.
//  Copyright © 2015 Mikel Hernaez. All rights reserved.
//

#include "em_markov_temporal.h"


// *******************************************************************************//

void alloc_markov_temporal_model(markov_temporal_model mm, uint32_t seq_length){
    
    uint32_t i = 0, t = 0;
    
    mm->num_nodes = QV_ALPHABET;
    
    mm->num_ts = seq_length;
    
    mm->A = (double***)calloc(mm->num_ts, sizeof(double*));
    for (t = 0; t < mm->num_ts; t++) {
        mm->A[t] = (double**)calloc(mm->num_nodes, sizeof(double*));
        for (i = 0; i < mm->num_nodes; i++) {
            mm->A[t][i] = (double*)calloc(mm->num_nodes, sizeof(double));
        }
        
    }
    
    mm->pi = (double*)calloc(mm->num_nodes, sizeof(double));
    
}

// *******************************************************************************//

em_temporal_markov alloc_temporal_em(qv_file qv_seqs, uint32_t num_models, uint64_t num_seq){
    
    uint32_t i = 0;
    
    em_temporal_markov EM = (struct em_temporal_markov_t*)calloc(1, sizeof(struct em_temporal_markov_t));
    
    EM->num_models = num_models;
    
    EM->num_seq = num_seq;
    
    EM->qv_seqs = qv_seqs;
    
    EM->clust = (struct clusters_t*)calloc(1, sizeof(struct clusters_t));
    
    EM->clust->clusters = (uint32_t *)calloc(EM->num_seq, sizeof(uint32_t));
    
    EM->clust->cluster_sizes = (uint32_t *)calloc(EM->num_models, sizeof(uint32_t));
    
    EM->clust->num_clusters = num_models;
    
    EM->models = (struct markov_temporal_model_t*)calloc(num_models, sizeof(struct markov_temporal_model_t));
    for (i = 0; i<num_models; i++) {
        alloc_markov_temporal_model(&(EM->models[i]),qv_seqs->read_length);
    }
    
    EM->models_prior = (double*)calloc(num_models, sizeof(double));
    
    EM->r = (double**)calloc(num_models, sizeof(double*));
    for (i = 0; i<num_models; i++) {
        EM->r[i] = (double*)calloc(num_seq, sizeof(double));
    }
    
    return EM;
}

// *******************************************************************************//

void initialize_markov_temporal_model(markov_temporal_model mm){
    
    uint32_t i = 0, j = 0, t = 0;
    double s_A, s_pi = 0;
    
    for (i = 0; i < mm->num_nodes; ++i) {
        mm->pi[i] = (double)rand();
        //mm->pi[i] = 1.0/(double)mm->num_nodes;
        s_pi += mm->pi[i];
    }
    // Normalize
    for (i = 0; i < mm->num_nodes; ++i) {
        mm->pi[i] /= s_pi;
    }
    for (t = 0; t < mm->num_ts-1; t++) {
        for (i = 0; i < mm->num_nodes; ++i) {
            s_A = 0;
            for (j = 0; j < mm->num_nodes; j++) {
                mm->A[t][i][j] = (double)(rand()/(double)RAND_MAX);
                s_A += mm->A[t][i][j];
            }
            // Normalize
            for (j = 0; j < mm->num_nodes; j++) {
                mm->A[t][i][j] /= s_A;
            }
        }
    }
    
}

// *******************************************************************************//

void initialize_em_temporal_markov(em_temporal_markov em){
    
    time_t  t;
    
    srand((unsigned) time(&t));
    
    //srand(1);
    
    uint32_t model_idx;
    
    for (model_idx = 0; model_idx < em->num_models; model_idx++) {
        initialize_markov_temporal_model(&(em->models[model_idx]));
    }
    
    double s_priors = 0;
    
    for (model_idx = 0; model_idx < em->num_models; model_idx++) {
        //em->models_prior[model_idx] = (double)rand();
        em->models_prior[model_idx] = 1.0/(double)em->num_models;
        s_priors += em->models_prior[model_idx];
    }
    for (model_idx = 0; model_idx < em->num_models; model_idx++) {
        em->models_prior[model_idx] /= s_priors;
    }
}

// ************************* EM Algorithm ****************************************//
//
// *******************************************************************************//

void temporal_e_step(em_temporal_markov em){
    
    uint32_t model_idx, seq_idx, t, xt0, xt1, T, x0;
    double norm_cte = 0, trans_prob;
    double ***A;
    double *pi;
    
    uint8_t *x_head = em->qv_seqs->file_head, *x_ptr;
    
    T = em->qv_seqs->read_length;
    
    
    for (seq_idx = 0; seq_idx < em->num_seq; seq_idx++){
        
        for (model_idx = 0; model_idx < em->num_models; model_idx++){
            x_ptr = x_head;
            A = em->models[model_idx].A;
            pi = em->models[model_idx].pi;
            
            xt0 = *(x_ptr), x_ptr++;
            x0 = xt0;
            trans_prob = 1.0;
            for (t = 1; t < T; ++t) {
                xt1 = *(x_ptr);
                trans_prob *= A[t-1][xt0][xt1];
                xt0 = xt1;
                x_ptr++;
            }
            x_ptr++; // Necessary now because of the zero between sequences.
            em->r[model_idx][seq_idx] = em->models_prior[model_idx]*pi[x0]*trans_prob;
            norm_cte += em->r[model_idx][seq_idx];
        }
        assert(norm_cte!=0.0);
        // Normalize r wrt num_models
        for (model_idx = 0; model_idx < em->num_models; model_idx++){
            em->r[model_idx][seq_idx] /= norm_cte;
        }
        norm_cte = 0;
        x_head += T;
        x_head++; // Necessary now because of the zero between sequences.
        
    }
}

// *******************************************************************************//

double temporal_m_step(em_temporal_markov em){
    static double prev_nll = 0;
    uint32_t seq_idx, model_idx, node_idx, dest_node_idx, t;
    uint32_t i;
    
    double nll = 0;
    
    double temp_w,norm_cte = 0;
    
    double *rl;
    uint8_t *qvs;
    
    double temp_pi[QV_ALPHABET] = {0};
    
    const uint32_t read_length = em->qv_seqs->read_length;
    
    double ***temp_A;
    temp_A = (double***)calloc(read_length, sizeof(double**));
    for (t = 0; t < read_length; t++) {
        temp_A[t] = (double**)calloc(em->models[0].num_nodes, sizeof(double*));
        for (i = 0; i < em->models[0].num_nodes; i++) {
            temp_A[t][i] = (double*)calloc(em->models[0].num_nodes, sizeof(double));
        }
    }
    
    
    
    for (model_idx = 0; model_idx < em->num_models; ++model_idx) {
        
        // Reset the accumulators to zero
        temp_w = 0;
        memset(temp_pi, 0, em->models[model_idx].num_nodes*sizeof(double));
        for (t = 0; t < read_length; t++) {
            for (i = 0; i < em->models[0].num_nodes; i++) {
                memset(temp_A[t][i],0,em->models[0].num_nodes*sizeof(double));
            }
            
        }
        // Set the pointers to the iterators
        rl = em->r[model_idx];
        qvs = em->qv_seqs->file_head;
        
        for (seq_idx = 0; seq_idx < em->num_seq; ++seq_idx){
            
            // Update the model prior
            temp_w += *rl;
            
            // Update the models pi
            temp_pi[*qvs] += *rl;
            
            // Update the models A
            for (t = 0; t < read_length - 1; t++) {
                temp_A[t][*qvs][*(qvs+1)] += *rl;
                qvs++;
            }
            
            rl++;
            qvs++;
            qvs++; // TODO: remove. Necessary now because of the zero between sequences.
        }
        
        // Compute the nll before normalizing
        nll += temp_w*log(em->models_prior[model_idx]);
        for (node_idx = 0; node_idx < em->models[model_idx].num_nodes;node_idx++ ) {
            if (em->models[model_idx].pi[node_idx] != 0) {
                nll += temp_pi[node_idx]*log(em->models[model_idx].pi[node_idx]);
            }
        }
        for (t = 0; t < read_length-1; t++) {
            for (node_idx = 0; node_idx < em->models[model_idx].num_nodes;node_idx++ ) {
                for (dest_node_idx = 0; dest_node_idx < em->models[model_idx].num_nodes; dest_node_idx++) {
                    if (em->models[model_idx].A[t][node_idx][dest_node_idx] != 0) {
                        nll += temp_A[t][node_idx][dest_node_idx]*log(em->models[model_idx].A[t][node_idx][dest_node_idx]);
                    }
                }
            }
        }
        
        // Normalize the model Prior
        em->models_prior[model_idx] = temp_w/em->num_seq;
        
        // Normalize Pi wrt num_nodes
        norm_cte = 0;
        for (node_idx = 0; node_idx < em->models[model_idx].num_nodes; node_idx++) {
            em->models[model_idx].pi[node_idx] = temp_pi[node_idx];
            norm_cte+=temp_pi[node_idx];
        }
        for (node_idx = 0; node_idx < em->models[model_idx].num_nodes; node_idx++) {
            em->models[model_idx].pi[node_idx] /= norm_cte;
        }
        assert(norm_cte!=0);
        
        // Normalize A wrt num_dest_nodes
        // Note that here we allow for zero probability as we don't mind overfitting (for now...)
        // This is quite dangerous...
        for (t = 0; t < read_length; t++) {
            for (node_idx = 0; node_idx < em->models[model_idx].num_nodes; node_idx++) {
                norm_cte = 0;
                for (dest_node_idx = 0; dest_node_idx < em->models[model_idx].num_nodes; dest_node_idx++) {
                
                    em->models[model_idx].A[t][node_idx][dest_node_idx] = temp_A[t][node_idx][dest_node_idx];
                    norm_cte+=temp_A[t][node_idx][dest_node_idx];
                }
                if (norm_cte != 0) {
                    for (dest_node_idx = 0; dest_node_idx < em->models[model_idx].num_nodes; dest_node_idx++) {
                        em->models[model_idx].A[t][node_idx][dest_node_idx] /= norm_cte;
                    }
                }
            
            }
        }
    }
    
    //if (nll < prev_nll) {
    //  ;
    //printf("asdf/n");
    //}
    prev_nll = nll;
    for (t = 0; t < read_length; t++) {
        for (i = 0; i < em->models[0].num_nodes; i++) {
            free(temp_A[t][i]);
        }
        free(temp_A[t]);
    }
    free(temp_A);
    return nll;
}

// *******************************************************************************//

uint32_t perform_em_temporal_markov(qv_file qv_f, uint32_t num_models, uint32_t iters, FILE *fo, const char *split_path){
    
    uint32_t i = 0;
    
    double data_ll;
    
    em_temporal_markov em = alloc_temporal_em(qv_f,num_models, qv_f->lines);
    
    initialize_em_temporal_markov(em);
    
    while (i++ < iters){
        temporal_e_step(em);
        data_ll = temporal_m_step(em);
        printf("%03d: %f\n",i,data_ll);
    }
    
    compute_clusters(em->clust, em->r, em->num_seq);
    
    return 0;
    
}

// *******************************************************************************//


