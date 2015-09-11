//
//  hmm_temporal.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 9/7/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

//
//  hmm.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/27/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#include "hmm.h"




// *******************************************************************************//

void initizlize_A_matrix(hmm h){
    
    uint32_t i,j,l;
    
    for (i = 0; i < h->num_states; i++) {
        for (j = i +(h->K - i%h->K), l = 0; l < h->K; j++, l++) {
            h->A[i][j] = 1.0/(double)h->K;
        }
    }
    
}

void initialize_pi(hmm h){
    uint32_t j;
    for (j = 0; j < h->K; j++) {
        h->pi[j] = 1.0/(double)h->K;
    }
}

void initialize_hmm_temporal(hmm h){
    
    uint32_t i = 0, j = 0;
    
    //uint64_t s_A, s_B, s_pi;
    double s_B, s_pi;
    
    s_pi = 0;
    for (i = 0; i < h->num_states; i++) {
        
        s_B = 0;
        for (j = 0; j < QV_ALPHABET; j++) {
            h->B[i][j] = rand();
            //h->B[i][j] = 1.0/(double)QV_ALPHABET;
            s_B += h->B[i][j];
        }
        //Normalize
        for (j = 0; j < QV_ALPHABET; j++) {
            h->B[i][j] /= (double)s_B;
        }
    }
    
    initizlize_A_matrix(h);
    
    initialize_pi(h);
    
}

// *******************************************************************************//

void initialize_bw_temporal(baum_welch bw){
    
    initialize_hmm_temporal(bw->model);
}

// *******************************************************************************//

void forwards_backwards_algorithm_temporal(forwards_backwards fb, hmm model, uint8_t *seq){
    
    int32_t t = 0, i = 0, j = 0, l = 0;
    
    double tmp = 0, s = 0.0;
    
    uint32_t K = model->K;
    
    // Initialize Alpha
    s = 0.0;
    for (j = 0; j < fb->num_states; j++) {
        fb->alpha[0][j] = model->pi[j]*model->B[j][seq[0]];
        s += fb->alpha[0][j];
    }
    fb->c[0] = 1/s;
    // Normalize wrt j;
    for (j = 0; j < fb->num_states; j++) {
        fb->alpha[0][j] /= s;
    }
    
    // Forward recursion [Alpha].
    for (t = 1; t < fb->seq_length; t++) {
        s = 0.0;
        for (j = 0; j < fb->num_states; j++) {
            tmp = 0.0;
            for (i = j - K - j%K, l = 0; (i >= 0 && l < K); i++, l++) {
                tmp += model->A[i][j]*fb->alpha[t-1][i];
            }
            fb->alpha[t][j] = tmp*model->B[j][seq[t]];
            s += fb->alpha[t][j];
        }
        fb->c[t] = 1/s;
        //Normalize
        for (j = 0; j < fb->num_states; j++) {
            fb->alpha[t][j] /= s;
        }
    }
    
    // Initialize Beta
    s = 0.0;
    for (j = 0; j < fb->num_states; j++) {
        //fb->beta[fb->seq_length - 1][j] = 1/fb->c[fb->seq_length - 1];
        fb->beta[fb->seq_length - 1][j] = 1;
    }
    
    // Backwards recursion [Beta].
    for (t = fb->seq_length-1; t > 0; --t) {
        for (i = 0; i < fb->num_states; i++) {
            fb->beta[t-1][i] = 0;
            for (j = i +(K - i%K), l = 0; (j < fb->num_states&&l < K); j++, l++) {
                fb->beta[t-1][i] += fb->beta[t][j]*model->A[i][j]*model->B[j][seq[t]];
            }
            fb->beta[t-1][i] *= fb->c[t];
        }
    }
}

// *******************************************************************************//




double em_step_bw_temporal(baum_welch bw){
    
    uint32_t K = bw->model->K;
    uint32_t j,k,i,t,l;
    
    uint8_t* seq = bw->qv_seqs->file_head;
    double s = 0.0, nll = 0.0;
    
    double **A, **B, *pi;
    
    A = (double**)calloc(bw->num_states, sizeof(double*));
    for (j = 0; j < bw->num_states; j++) {
        A[j] = (double*)calloc(bw->num_states, sizeof(double));
    }
    B = (double**)calloc(bw->num_states, sizeof(double*));
    for (j = 0 ; j < bw->num_states; j++) {
        B[j] = (double*)calloc(QV_ALPHABET, sizeof(double));
    }
    
    pi = (double*)calloc(bw->num_states, sizeof(double));
    
    for (i = 0; i < bw->num_seq; i++) {
        
        //printf("Seq: %d\n",i);
        
        forwards_backwards_algorithm_temporal(bw->fb[0], bw->model, seq);
        
        // Compute the new A
        for (j = 0; j < bw->num_states; j++) {
            for (k = j+(K-j%K), l = 0; (k < bw->model->num_states && l < K); k++, l++) {
                for (t = 0; t < bw->seq_length-1; t++) {
                    A[j][k] += bw->fb[0]->alpha[t][j]*bw->model->A[j][k]*bw->model->B[k][seq[t+1]]*bw->fb[0]->c[t+1]*bw->fb[0]->beta[t+1][k];
                }
            }
        }
        
        // Compute the new pi
        for (j = 0; j < bw->num_states; j++) {
            pi[j] += bw->fb[0]->alpha[0][j]*bw->fb[0]->beta[0][j];
            //pi[j] += bw->fb[0]->gamma[0][j];
        }
        
        // Compute the new B
        for (j = 0; j < bw->num_states; j++) {
            for (t = 0; t < bw->seq_length; t++) {
                B[j][seq[t]] += bw->fb[0]->alpha[t][j]*bw->fb[0]->beta[t][j];
                //B[j][seq[t]] += bw->fb[0]->gamma[t][j];
            }
        }
        
        // Compute nll
        for (t = 0; t < bw->fb[0]->seq_length; t++) {
            nll+= log(bw->fb[0]->c[t]);
        }
        
        seq += bw->seq_length;
        seq++;// TODO: remove. Necessary now because of the zero between sequences.
        
    }
    
    // Normalize the model
    
    //Normalize A
    for (i = 0; i < bw->num_states; i++) {
        s = 0.0;
        for (j = i +(K - i%K), l = 0; (j < bw->num_states&&l < K); j++, l++){
            bw->model->A[i][j] = A[i][j];
            s += bw->model->A[i][j];
        }
        
        for (j = 0; j <bw->num_states; j++)
            bw->model->A[i][j] /= s;
    }
    
    //Normalize B
    for (i = 0; i < bw->num_states; i++) {
        s = 0.0;
        for (l = 0; l <QV_ALPHABET; l++){
            bw->model->B[i][l] = B[i][l];
            s += bw->model->B[i][l];
        }
        
        for (l = 0; l <QV_ALPHABET; l++)
            bw->model->B[i][l] /= s;
    }
    
    //Normalize Pi
    s = 0.0;
    for (i = 0; i < bw->num_states; i++){
        bw->model->pi[i] = pi[i];
        s+= bw->model->pi[i];
    }
    for (i = 0; i < bw->num_states; i++)
        bw->model->pi[i] /= s;
    
    ///// Free Aux Variables ////
    for (j = 0; j < bw->num_states; j++)
        free(A[j]);
    free(A);
    for (j = 0; j < bw->num_states; j++)
        free(B[j]);
    free(B);
    
    free(pi);
    
    return nll;
}

// *******************************************************************************//

void output_model_temporal(hmm h, FILE *fo){
    
    uint32_t i = 0, j = 0;
    for (i = 0; i < QV_ALPHABET; i++) {
        for (j = 0; j < h->num_states-1; j++) {
            fprintf(fo,"%f,",h->B[j][i]);
        }
        fprintf(fo,"%f\n",h->B[j][i]);
    }
    for (j = 0; j < h->num_states-1; j++) {
        fprintf(fo,"%f,",h->pi[j]);
    }
    fprintf(fo,"%f\n",h->pi[j]);
}

// *******************************************************************************//

void print_graph_hmm_temporal(hmm m, const char* graph_file_path){
    uint32_t i = 0, j = 0;
    
    static uint32_t ctr = 0;
    
    char graph_path[1024];
    char ch_ctr[256];
    
    FILE *graph_file;
    
    sprintf(ch_ctr, "%05d",ctr);
    
    strcpy(graph_path, graph_file_path);
    strcat(graph_path, ch_ctr);
    
    
    
    graph_file = fopen(graph_path, "w");
    
    fprintf(graph_file, "digraph {\n");
    
    for (i = 0; i < m->num_states; i++) {
        fprintf(graph_file, "%d\n",i+1);
    }
    
    
    for (i = 0; i < m->num_states; i++) {
        for (j = 0; j < m->num_states; j++) {
            if (m->A[i][j] > 0.00001)
                fprintf(graph_file,"%d->%d[label=%.3f]\n",i+1,j+1,m->A[i][j]);
        }
    }
    
    
    fprintf(graph_file, "}\n");
    
    fclose(graph_file);
    ctr++;
}

// *******************************************************************************//

void perform_bw_temporal(qv_file qv_f, uint32_t num_states_t, uint32_t num_iters, FILE *fo, const char *fgraph){
    
    uint32_t i = 0;
    
    double nll;
    
    uint32_t num_states = num_states_t*qv_f->read_length;
    
    baum_welch bw = alloc_bw(qv_f, num_states);
    
    bw->model->K = num_states_t;
    
    initialize_bw_temporal(bw);
    
    while (i++ < num_iters) {
        nll = em_step_bw_temporal(bw);
        //print_graph_hmm(bw->model, fgraph);
        printf("%f\n",nll);
    }
    output_model(bw->model, fo);
}

// *******************************************************************************//

