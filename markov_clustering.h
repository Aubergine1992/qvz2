//
//  markov_clustering.h
//  markov_clustering
//
//  Created by Mikel Hernaez-Arrazola on 8/27/15.
//  Copyright (c) 2015 Mikel Hernaez. All rights reserved.
//

#ifndef markov_clustering_markov_clustering_h
#define markov_clustering_markov_clustering_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

//#define log2(a) log(a)/log(2.0)

#define qv2ch(a) a-33


#define QV_ALPHABET 42
// This limits us to chunks that aren't too big to fit into a modest amount of memory at a time
#define MAX_LINES_PER_BLOCK			1000000
#define MAX_READS_PER_LINE			1022
#define READ_LINEBUF_LENGTH			(MAX_READS_PER_LINE+2)

typedef struct qv_file_t{
    uint8_t* file_head;
    uint64_t lines;
    uint32_t read_length;
} *qv_file;

typedef struct clusters_t{
    uint32_t* clusters;
    uint32_t* cluster_sizes;
}*clusters;

#endif
