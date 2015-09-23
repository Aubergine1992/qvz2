//
//  em_markov.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/27/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#include "em_markov.h"


// *******************************************************************************//

void alloc_markov_model(markov_model mm){
    
    uint32_t i = 0;
    
    mm->num_nodes = QV_ALPHABET;
    
    mm->A = (double**)calloc(mm->num_nodes, sizeof(double*));
    for (i = 0; i < mm->num_nodes; i++) {
        mm->A[i] = (double*)calloc(mm->num_nodes, sizeof(double));
    }
    
    mm->pi = (double*)calloc(mm->num_nodes, sizeof(double));
    
}

// *******************************************************************************//

em_markov alloc_em(qv_file qv_seqs, uint32_t num_models, uint64_t num_seq){
    
    uint32_t i = 0;
    
    em_markov EM = (struct em_markov_t*)calloc(1, sizeof(struct em_markov_t));
    
    EM->num_models = num_models;
    
    EM->num_seq = num_seq;
    
    EM->qv_seqs = qv_seqs;
    
    EM->clust = (struct clusters_t*)calloc(1, sizeof(struct clusters_t));
    
    EM->clust->clusters = (uint32_t *)calloc(EM->num_seq, sizeof(uint32_t));
    
    EM->clust->cluster_sizes = (uint32_t *)calloc(EM->num_models, sizeof(uint32_t));
    
    EM->models = (struct markov_model_t*)calloc(num_models, sizeof(struct markov_model_t));
    for (i = 0; i<num_models; i++) {
        alloc_markov_model(&(EM->models[i]));
    }
    
    EM->models_prior = (double*)calloc(num_models, sizeof(double));
    
    EM->r = (double**)calloc(num_models, sizeof(double*));
    for (i = 0; i<num_models; i++) {
        EM->r[i] = (double*)calloc(num_seq, sizeof(double));
    }
    
    return EM;
}

// *******************************************************************************//

void initialize_markov_model(markov_model mm){
    
    uint32_t i = 0, j = 0;
    double s_A, s_pi = 0;
    
    for (i = 0; i < mm->num_nodes; ++i) {
        //mm->pi[i] = (double)rand();
        mm->pi[i] = 1.0/(double)mm->num_nodes;
        s_pi += mm->pi[i];
    }
    // Normalize
    for (i = 0; i < mm->num_nodes; ++i) {
        mm->pi[i] /= s_pi;
    }
    
    for (i = 0; i < mm->num_nodes; ++i) {
        s_A = 0;
        for (j = 0; j < mm->num_nodes; j++) {
            mm->A[i][j] = (double)(rand()/(double)RAND_MAX);
            s_A += mm->A[i][j];
        }
        // Normalize
        for (j = 0; j < mm->num_nodes; j++) {
            mm->A[i][j] /= s_A;
        }
    }
    
}

// *******************************************************************************//

void initialize_em_markov(em_markov em){
    
    time_t  t;
    
    srand((unsigned) time(&t));
    
    //srand(1);
    
    uint32_t model_idx;
    
    for (model_idx = 0; model_idx < em->num_models; model_idx++) {
        initialize_markov_model(&(em->models[model_idx]));
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

void e_step(em_markov em){
    
    uint32_t model_idx, seq_idx, t, xt0, xt1, T;
    double norm_cte = 0, r = 1, log_trans_prob;
    double **A;
    double *pi;
    
    uint8_t *x_head = em->qv_seqs->file_head, *x_ptr;
    
    T = em->qv_seqs->read_length;
    
    
    for (seq_idx = 0; seq_idx < em->num_seq; seq_idx++){
        
        for (model_idx = 0; model_idx < em->num_models; model_idx++){
            x_ptr = x_head;
            A = em->models[model_idx].A;
            pi = em->models[model_idx].pi;
            
            xt0 = *(x_ptr), x_ptr++;
            r = em->models_prior[model_idx]*pi[xt0];
            log_trans_prob = 0.0;
            for (t = 1; t < T; ++t) {
                xt1 = *(x_ptr);
                log_trans_prob += log(A[xt0][xt1]);
                xt0 = xt1;
                x_ptr++;
            }
            x_ptr++; // TODO: remove. Necessary now because of the zero between sequences.
            em->r[model_idx][seq_idx] = r*exp(log_trans_prob);
            norm_cte += em->r[model_idx][seq_idx];
        }
        assert(norm_cte!=0.0);
        // Normalize r wrt num_models
        for (model_idx = 0; model_idx < em->num_models; model_idx++){
            em->r[model_idx][seq_idx] /= norm_cte;
        }
        norm_cte = 0;
        x_head += T;
        x_head++; // TODO: remove. Necessary now because of the zero between sequences.
        
    }
}

// *******************************************************************************//

double m_step(em_markov em){
    static double prev_nll = 0;
    uint32_t seq_idx, model_idx, node_idx, dest_node_idx, t;
    uint32_t i;
    
    double nll = 0;
    
    double temp_w,norm_cte = 0;
    
    double *rl;
    uint8_t *qvs;
    
    double temp_pi[QV_ALPHABET] = {0};
    
    double **temp_Ajk;
    temp_Ajk = (double**)calloc(em->models[0].num_nodes, sizeof(double*));
    for (i = 0; i < em->models[0].num_nodes; i++) {
        temp_Ajk[i] = (double*)calloc(em->models[0].num_nodes, sizeof(double));
    }
    
    const uint32_t read_length = em->qv_seqs->read_length;
    
    for (model_idx = 0; model_idx < em->num_models; ++model_idx) {
        
        // Reset the accumulators to zero
        temp_w = 0;
        memset(temp_pi, 0, em->models[model_idx].num_nodes*sizeof(double));
        for (i = 0; i < em->models[model_idx].num_nodes; i++) {
            memset(temp_Ajk[i],0,em->models[model_idx].num_nodes*sizeof(double));
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
                temp_Ajk[*qvs][*(qvs+1)] += *rl;
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
        for (node_idx = 0; node_idx < em->models[model_idx].num_nodes;node_idx++ ) {
            for (dest_node_idx = 0; dest_node_idx < em->models[model_idx].num_nodes; dest_node_idx++) {
                if (em->models[model_idx].A[node_idx][dest_node_idx] != 0) {
                    nll += temp_Ajk[node_idx][dest_node_idx]*log(em->models[model_idx].A[node_idx][dest_node_idx]);
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
        for (node_idx = 0; node_idx < em->models[model_idx].num_nodes; node_idx++) {
            norm_cte = 0;
            for (dest_node_idx = 0; dest_node_idx < em->models[model_idx].num_nodes; dest_node_idx++) {
                
                em->models[model_idx].A[node_idx][dest_node_idx] = temp_Ajk[node_idx][dest_node_idx];
                norm_cte+=temp_Ajk[node_idx][dest_node_idx];
            }
            if (norm_cte != 0) {
                for (dest_node_idx = 0; dest_node_idx < em->models[model_idx].num_nodes; dest_node_idx++) {
                    em->models[model_idx].A[node_idx][dest_node_idx] /= norm_cte;
                }
            }
            
        }
    }
    
    //if (nll < prev_nll) {
      //  ;
        //printf("asdf/n");
    //}
    prev_nll = nll;
    
    for (i = 0; i < em->models[0].num_nodes; i++) {
        free(temp_Ajk[i]);
    }
    free(temp_Ajk);
    return nll;
}
// *******************************************************************************//

void compute_clusters(em_markov em){
    
    uint32_t seq_idx, i = 0, tmp_model;
    
    double tmp = 0;
    
    for (i = 0; i < em->num_models; i++) {
        em->clust->cluster_sizes[i] = 0;
    }
    
    for (seq_idx = 0; seq_idx < em->num_seq; seq_idx++) {
        tmp = 0;
        tmp_model = 0;
        for (i = 0; i < em->num_models; i++) {
            if (em->r[i][seq_idx] >= tmp) {
                tmp = em->r[i][seq_idx];
                tmp_model = i;
            }
        }
        em->clust->cluster_sizes[tmp_model]++;
        em->clust->clusters[seq_idx] = tmp_model;
    }
}

// *******************************************************************************//

void print_graph(em_markov m, FILE * graph_file){
    uint32_t i = 0, j = 0, k = 0;
    
    fprintf(graph_file, "digraph {\n");
    
    for (i = 0; i < m->models->num_nodes; i++) {
        fprintf(graph_file, "%d\n",i);
    }
    
    char tmp_char[512] = "";
    char edge_color[256] = "%d->%d[label=%.3f color=";
    
    char color[256][12] = {"red","blue","green"};
    
    for (k = 0; k < m->num_models; k++) {
        for (i = 0; i < m->models[k].num_nodes; i++) {
            for (j = 0; j < m->models[k].num_nodes; j++) {
                if (m->models[k].A[i][j] > 0.1) {
                    strcpy(tmp_char, edge_color);
                    strcat(tmp_char,color[k]);
                    strcat(tmp_char,"]\n");
                    fprintf(graph_file,tmp_char,i,j,m->models[k].A[i][j]);
                    *tmp_char = 0;
                }
            }
        }
    }
    
    fprintf(graph_file, "}\n");
}

// *******************************************************************************//
void split_data(const char* path_head, em_markov em){
    uint32_t i = 0, cluster_id, t;
    FILE **fo;
    char path_buffer[1024];
    char cluster_name[256] = "/cluster";
    
    
    fo = (FILE**)calloc(em->num_models, sizeof(FILE*));
    
    for (i = 0; i < em->num_models; i++) {
        cluster_name[8] = (i+48);
        strcpy(path_buffer, path_head);
        strcat(path_buffer,cluster_name);
        fo[i] = fopen(path_buffer, "w");
    }
    
    uint8_t *seq = em->qv_seqs->file_head;
    
    for (i = 0; i < em->num_seq; i++) {
        cluster_id = em->clust->clusters[i];
        for (t = 0; t < em->qv_seqs->read_length; t++) {
            fputc((*seq)+33, fo[cluster_id]), seq++;
        }
        fputc('\n', fo[cluster_id]);
        seq++;
    }
    
    
}

// *******************************************************************************//

uint32_t perform_em_markov(qv_file qv_f, uint32_t num_models, uint32_t iters, FILE *fo, const char *split_path){
    
    uint32_t i = 0, j = 0;
    
    double data_ll;
    
    em_markov em = alloc_em(qv_f,num_models, qv_f->lines);
    
    initialize_em_markov(em);
    
    while (i++ < iters){
        e_step(em);
        data_ll = m_step(em);
        printf("%03d: %f\n",i,data_ll);
    }
    
    compute_clusters(em);
    for (j = 0; j < em->num_seq; j++) {
        fprintf(fo, "%d\n",em->clust->clusters[j]);
    }
    //print_graph(em, fgraph);
    split_data(split_path, em);
    
    return 0;
    
}

// *******************************************************************************//


