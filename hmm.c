//
//  hmm.c
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/27/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#include "hmm.h"


hmm alloc_hmm(uint32_t num_states){
    
    uint32_t i = 0;
    
    hmm h = (struct hmm_t*)calloc(1, sizeof(struct hmm_t));
    
    h->num_states = num_states;
    
    h->pi = (double*)calloc(num_states, sizeof(double));
    
    h->A = (double**) calloc(num_states, sizeof(double*));
    for (i = 0; i < num_states; i++) {
        h->A[i] = (double*) calloc(num_states, sizeof(double));
    }
    
    h->B = (double**) calloc(num_states, sizeof(double*));
    for (i = 0; i < num_states; i++) {
        h->B[i] = (double*) calloc(QV_ALPHABET, sizeof(double));
    }
    
    return h;
    
}

forwards_backwards alloc_fb(uint32_t seq_length, uint32_t num_states){
    
    uint32_t i, j;
    
    forwards_backwards fb  = (struct forwards_backwards_t*)calloc(1, sizeof(struct forwards_backwards_t));
    
    fb->seq_length = seq_length;
    fb->num_states = num_states;
    
    fb->c = (double*)calloc(fb->seq_length, sizeof(double));
    fb->alpha = (double**)calloc(fb->seq_length, sizeof(double*));
    fb->beta = (double**)calloc(fb->seq_length, sizeof(double*));
    fb->gamma = (double**)calloc(fb->seq_length, sizeof(double*));
    fb->chi = (double***)calloc(fb->seq_length, sizeof(double**));
    for (i = 0; i < fb->seq_length; i++) {
        fb->alpha[i] = (double*)calloc(fb->num_states, sizeof(double));
        fb->beta[i] = (double*)calloc(fb->num_states, sizeof(double));
        fb->gamma[i] = (double*)calloc(fb->num_states, sizeof(double));
        fb->chi[i] = (double**)calloc(fb->num_states, sizeof(double*));
        for (j = 0; j < fb->num_states; j++) {
            fb->chi[i][j] = (double*)calloc(fb->num_states, sizeof(double));
        }
    }
    
    return fb;
}


// *******************************************************************************//

baum_welch alloc_bw(qv_file qv_seqs, uint32_t num_states)
{
    uint32_t i = 0;
    
    baum_welch bw = (struct baum_welch_t*) calloc(num_states, sizeof(struct baum_welch_t));
    
    bw->qv_seqs = qv_seqs;
    
    bw->seq_length = qv_seqs->read_length;
    
    bw->num_seq = qv_seqs->lines;
    
    bw->num_states = num_states;
    
    bw->model = alloc_hmm(num_states);
    
    bw->fb = (forwards_backwards*)calloc(bw->num_seq, sizeof(forwards_backwards));
    for (i = 0; i < 1; i++) {
        bw->fb[i] = alloc_fb(bw->seq_length, bw->num_states);
    }
    return bw;
}

// *******************************************************************************//

void initialize_hmm(hmm h){
    
    uint32_t i = 0, j = 0;
    
    //uint64_t s_A, s_B, s_pi;
    double s_A, s_B, s_pi;
    
    s_pi = 0;
    for (i = 0; i < h->num_states; i++) {
        s_A = 0;
        for (j = 0; j < h->num_states; j++) {
            //h->A[i][j] = rand();
            h->A[i][j] = 1.0/(double)h->num_states;
            s_A += h->A[i][j];
        }
        //Normalize
        for (j = 0; j < h->num_states; j++) {
            h->A[i][j] /= (double)s_A;
        }
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
        
        //h->pi[i] = rand();
        h->pi[i] = 1.0/(double)h->num_states;
        s_pi +=h->pi[i];
        
    }
    //Normalize
    for (i = 0; i < h->num_states; i++) {
        h->pi[i] /= (double)s_pi;
    }
    
}

// *******************************************************************************//

void initialize_bw(baum_welch bw){
    
    initialize_hmm(bw->model);
}

// *******************************************************************************//

void forwards_backwards_algorithm(forwards_backwards fb, hmm model, uint8_t *seq){
    
    int32_t t = 0, i = 0, j = 0;
    
    double tmp = 0, s = 0.0;
    
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
            for (i = 0; i < fb->num_states; i++) {
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
        fb->beta[fb->seq_length - 1][j] = 1;
    }
    
    // Backwards recursion [Beta].
    for (t = fb->seq_length-1; t > 0; --t) {
        for (i = 0; i < fb->num_states; i++) {
            fb->beta[t-1][i] = 0;
            for (j = 0; j < fb->num_states; j++) {
                fb->beta[t-1][i] += fb->beta[t][j]*model->A[i][j]*model->B[j][seq[t]];
            }
            fb->beta[t-1][i] *= fb->c[t];
        }
    }
}

// *******************************************************************************//




double em_step_bw(baum_welch bw){
    
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
        
        forwards_backwards_algorithm(bw->fb[0], bw->model, seq);
        
        // Compute the new A
        for (j = 0; j < bw->num_states; j++) {
            for (k = 0; k < bw->num_states; k++) {
                for (t = 0; t < bw->seq_length-1; t++) {
                    A[j][k] += bw->fb[0]->alpha[t][j]*bw->model->A[j][k]*bw->model->B[k][seq[t+1]]*bw->fb[0]->c[t+1]*bw->fb[0]->beta[t+1][k];
                    //A[j][k] += bw->fb[0]->chi[t][j][k];
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
        for (j = 0; j <bw->num_states; j++){
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

void output_model(hmm h, FILE *fo){
    
    uint32_t i = 0, j = 0;
    for (i = 0; i < QV_ALPHABET; i++) {
        for (j = 0; j < h->num_states-1; j++) {
            fprintf(fo,"%f,",h->B[j][i]);
        }
        fprintf(fo,"%f\n",h->B[j][i]);
    }
    for (i = 0; i < h->num_states; i++) {
        for (j = 0; j < h->num_states-1; j++) {
            fprintf(fo,"%f,",h->A[i][j]);
        }
        fprintf(fo,"%f\n",h->A[i][j]);
    }
    for (j = 0; j < h->num_states-1; j++) {
        fprintf(fo,"%f,",h->pi[j]);
    }
    fprintf(fo,"%f\n",h->pi[j]);
}

// *******************************************************************************//

void print_graph_hmm(hmm m, const char* graph_file_path){
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

void perform_bw(qv_file qv_f, uint32_t num_states, uint32_t num_iters, FILE *fo, const char *fgraph){
    
    uint32_t i = 0;
    
    double nll;

    baum_welch bw = alloc_bw(qv_f, num_states);
    
    initialize_bw(bw);
    
    while (i++ < num_iters) {
        nll = em_step_bw(bw);
        print_graph_hmm(bw->model, fgraph);
        printf("%f\n",nll);
    }
    output_model(bw->model, fo);
}

// *******************************************************************************//
// *******************************************************************************//
// *******************************************************************************//
// *******************************************************************************//

/* void e_step_bw(baum_welch bw){
 
 uint32_t i;
 
 uint8_t *seq;
 
 seq = bw->qv_seqs->file_head;
 
 for (i = 0; i < bw->num_seq; i++) {
 forwards_backwards_algorithm(bw->fb[i], bw->model, seq);
 seq += bw->seq_length;
 seq++;// TODO: remove. Necessary now because of the zero between sequences.
 }
 }
 
 
 
 void m_step_bw(baum_welch bw){
 
 uint32_t j,k,i,t,l;
 
 double Ajk = 0.0, pi = 0.0, Nj = 0.0;
 double s;
 
 uint8_t* seq = bw->qv_seqs->file_head;
 
 // Compute the new A
 for (j = 0; j < bw->num_states; j++) {
 
 Nj = 0.0;
 for (i = 0; i < bw->num_seq; i++) {
 for (t = 0; t < bw->seq_length-1; t++) {
 Nj += bw->fb[i]->gamma[t][j];
 }
 }
 
 s = 0.0;
 for (k = 0; k < bw->num_states; k++) {
 Ajk = 0.0;
 for (i = 0; i < bw->num_seq; i++) {
 for (t = 0; t < bw->seq_length-1; t++) {
 Ajk += bw->fb[i]->chi[t][j][k];
 }
 }
 bw->model->A[j][k] = Ajk;
 s += bw->model->A[j][k];
 }
 // Normalize wrt k
 if (s != 0) {
 for (k = 0; k < bw->num_states; k++) {
 bw->model->A[j][k] /= s;
 }
 }
 }
 
 // Compute the new pi
 s = 0.0;
 for (j = 0; j < bw->num_states; j++) {
 pi = 0.0;
 for (i = 0; i < bw->num_seq; i++) {
 pi += bw->fb[i]->gamma[0][j];
 }
 bw->model->pi[j] = pi;
 s += pi;
 }
 // Normalize wrt j
 if (s != 0) {
 for (j = 0; j < bw->num_states; j++) {
 bw->model->pi[j] /= s;
 }
 }
 
 // Compute the new B
 for (j = 0; j < bw->num_states; j++) {
 
 Nj = 0.0;
 seq = bw->qv_seqs->file_head;
 for (i = 0; i < bw->num_seq; i++) {
 for (t = 0; t < bw->seq_length; t++) {
 bw->model->B[j][seq[t]] += bw->fb[i]->gamma[t][j];
 Nj += bw->fb[i]->gamma[t][j];
 }
 seq++;// TODO: remove. Necessary now because of the zero between sequences.
 seq += bw->seq_length;
 }
 // Normalize
 for (l = 0; l < QV_ALPHABET; l++) {
 if (Nj != 0) {
 bw->model->B[j][l] /= Nj;
 }
 }
 }
 }*/
// *******************************************************************************//
/*
 double compute_nll_bw(baum_welch bw){
 
 uint32_t i,k,j,t;
 
 double nll = 0.0, nll1, nll2, nll3;
 
 double pi, Ajk;
 
 uint8_t* seq = bw->qv_seqs->file_head;
 
 uint8_t x;
 
 for (k = 0; k < bw->num_states; k++) {
 if (bw->model->pi[k] == 0) {
 continue;
 }
 pi = 0.0;
 for (i = 0; i < bw->num_seq; i++) {
 pi += bw->fb[i]->gamma[0][k];
 }
 nll += pi*log(bw->model->pi[k]);
 }
 nll1 = nll;
 //printf("nll: %f\t",nll1);
 
 for (j = 0; j < bw->num_states; j++) {
 for (k = 0; k < bw->num_states; k++) {
 Ajk = 0.0;
 for (i = 0; i < bw->num_seq; i++) {
 for (t = 0; t < bw->seq_length-1; t++) {
 Ajk += bw->fb[i]->chi[t][j][k];
 
 }
 }
 if (bw->model->A[j][k] != 0) {
 nll += Ajk*log(bw->model->A[j][k]);
 }
 
 }
 }
 nll2 = nll - nll1;
 //printf("%f\t",nll2);
 
 for (k = 0; k < bw->num_states; k++) {
 seq = bw->qv_seqs->file_head;
 for (i = 0; i < bw->num_seq; i++) {
 for (t = 0; t < bw->seq_length; t++) {
 x = *seq, seq++;
 if (bw->model->B[k][x] != 0) {
 nll += bw->fb[i]->gamma[t][k]*log(bw->model->B[k][x]);
 }
 }
 seq++;// TODO: remove. Necessary now because of the zero between sequences.
 }
 }
 nll3 = nll - nll1 - nll2;
 //printf("%f\t%f\n",nll3, nll);
 return nll;
 }
 
 
 *****************************************************
// Computing the posterior transition probability [chi].
// Computing the posterior probability [gamma].
for (t = 0; t < fb->seq_length-1; t++) {
 //s = 0.0;
 //C = 0.0;
 for (j = 0; j < fb->num_states;j++) {
 fb->gamma[t][j] = fb->alpha[t][j]*fb->beta[t][j];
 //C += fb->alpha[t][j]*fb->beta[t][j];
 }
 for (j = 0; j < fb->num_states; j++) {
 for (i = 0; i < fb->num_states; i++) {
 fb->chi[t][i][j] = fb->alpha[t][i]*model->A[i][j]*model->B[j][seq[t+1]]*fb->c[t+1]*fb->beta[t+1][j];
 //fb->chi[t][i][j] /= C;
 //s += fb->chi[t][i][j];
 }
 }
 // Chi should be already normalized
 //if (s < 0.99999998 || s> 1.000000002 || C < 0.99999998 || C> 1.000000002) {
 //printf("foooooo");
 //for (j = 0; j < fb->num_states; j++) {
 // for (i = 0; i < fb->num_states; i++) {
 //    fb->chi[t][i][j] /= s;
 // }
 //}
 //}
 }*/
//double test[200][200] = {{0}};
// Computing the posterior probability [gamma].
//for (t = 0; t < fb->seq_length; t++) {
//for (j = 0; j < fb->num_states; j++) {
//fb->gamma[t][j] = fb->alpha[t][j]*fb->beta[t][j];
//}
//}
// Computing the posterior probability [gamma].
/*for (t = 0; t < fb->seq_length; t++) {
 for (j = 0; j < fb->num_states; j++) {
 for (i = 0; i < fb->num_states; i++) {
 test[t][j] += fb->chi[t][j][i];
 }
 }
 }
 
 for (t = 0; t < fb->seq_length-1; t++) {
 s = 0.0;
 for (j = 0; j < fb->num_states; j++) {
 if (test[t][j] - fb->gamma[t][j] > 0.00000001 || test[t][j] - fb->gamma[t][j] < -0.00000001 ) {
 printf("asdfasdf");
 }
 }
 }
 */
// *******************************************************************************//
