#include "codebook.h"
//#include "cluster.h"

#include <stdio.h>
#include <assert.h>

#if defined(LINUX) || defined(__APPLE__)
	#include <arpa/inet.h>
#endif

/**
 * To compute stats for the training data, we will need a set of conditional PMFs, one
 * per column
 * @param alphabet The symbol alphabet for each column
 * @param columns The number of columns to allocate conditional PMFs for
 */
struct cond_pmf_list_t *alloc_conditional_pmf_list(const struct alphabet_t *alphabet, uint32_t columns) {
	uint32_t count = 1 + alphabet->size*(columns-1);
	uint32_t i;
	struct cond_pmf_list_t *list = (struct cond_pmf_list_t *) calloc(1, sizeof(struct cond_pmf_list_t));

	// We need one array of PMF pointers that will index into the buffer allocated above, for the columns
	list->columns = columns;
	list->alphabet = alphabet;
	list->pmfs = (struct pmf_t **) calloc(count, sizeof(struct pmf_t *));

	// All PMFs are stored in a flat array, the accessor function will resolve a PMF's address
	for (i = 0; i < count; ++i) {
		list->pmfs[i] = alloc_pmf(alphabet);
	}

	return list;
}

/**
 * Deallocate the PMF list given and unallocate the two allocated memory blocks
 * @param list The conditional pmf list to deallocate
 */
void free_conditional_pmf_list(struct cond_pmf_list_t *list) {
	uint32_t count = 1 + list->alphabet->size * (list->columns - 1);
	uint32_t i;

	for (i = 0; i < count; ++i) {
		free_pmf(list->pmfs[i]);
	}

	free_pmf_list(list->marginal_pmfs);
    
    free(list);
}

/**
 * Allocate the quantizer list structure and the first level of array based on columns
 * @param columns The number of columns for which we have quantizers
 * @return Pointer to conditional quantizer list structure
 */
struct cond_quantizer_list_t *alloc_conditional_quantizer_list(uint32_t columns) {
	struct cond_quantizer_list_t *rtn = (struct cond_quantizer_list_t *) calloc(1, sizeof(struct cond_quantizer_list_t));
	rtn->columns = columns;
	rtn->input_alphabets = (struct alphabet_t **) calloc(columns, sizeof(struct alphabet_t *));
	rtn->q = (struct quantizer_t ***) calloc(columns, sizeof(struct quantizer_t **));
    rtn->ratio = (double **) calloc(columns, sizeof(double *));
	rtn->qratio = (uint8_t **) calloc(columns, sizeof(uint8_t *));
	return rtn;
}

/**
 * Deallocate the quantizer list as well as any alphabets or pmfs that are stored
 * @param list The conditional quantizer list to deallocate
 */
void free_cond_quantizer_list(struct cond_quantizer_list_t *list) {
	uint32_t i, j;

	for (i = 0; i < list->columns; ++i) {
		if (list->q[i]) {
			for (j = 0; j < list->input_alphabets[i]->size; ++j) {
				if (list->q[i][j])
					free_quantizer(list->q[i][j]);
			}
			free_alphabet(list->input_alphabets[i]);
			free(list->q[i]);
			free(list->ratio[i]);
			free(list->qratio[i]);
		}
	}

	free(list->qratio);
	free(list->ratio);
	free(list->q);
	free(list->input_alphabets);
	free(list);
}

/**
 * Initialize the information within a quantizer for the given column. This can't be done
 * at allocation time because we don't know everything about this column until we get here
 * during the optimization process
 * @param list The conditional quantizer list to update
 * @param column The column to initialize
 * @param input_union The alphabet of all possible left context symbols
 */
void cond_quantizer_init_column(struct cond_quantizer_list_t *list, uint32_t column, const struct alphabet_t *input_union) {
	list->input_alphabets[column] = duplicate_alphabet(input_union);

	// Low and high quantizer per element of the input union
	list->q[column] = (struct quantizer_t **) calloc(input_union->size*2, sizeof(struct quantizer_t *));

	// One ratio per element of input union
	list->ratio[column] = (double *) calloc(input_union->size, sizeof(double));
	list->qratio[column] = (uint8_t *) calloc(input_union->size, sizeof(uint8_t));
}

/**
 * Find a PMF for a specific column with the specific previous value
 */
struct pmf_t *get_cond_pmf(struct cond_pmf_list_t *list, uint32_t column, symbol_t prev) {
	if (column == 0)
		return list->pmfs[0];
	return list->pmfs[1 + (column-1)*list->alphabet->size + prev];
}

/**
 * Get a quantizer by its indexed location within the quantizer list for a column
 */
struct quantizer_t *get_cond_quantizer_indexed(struct cond_quantizer_list_t *list, uint32_t column, uint32_t index) {
	return list->q[column][index];
}

/**
 * Get a quantizer by its left context symbol
 */
struct quantizer_t *get_cond_quantizer(struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev) {
	uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
	if (idx != ALPHABET_SYMBOL_NOT_FOUND)
		return get_cond_quantizer_indexed(list, column, idx);
	return NULL;
}

/**
 * Stores the given quantizers at the appropriate index corresponding to the left context symbol given
 * for the specific column
 */
void store_cond_quantizers(struct quantizer_t *restrict lo, struct quantizer_t *restrict hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev) {
	uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
	store_cond_quantizers_indexed(lo, hi, ratio, list, column, idx);
}

/**
 * Stores the given quantizers directly at the given index based on the previous context symbol. Faster when
 * we know what the previous index was in addition to what the previous symbol was
 */
void store_cond_quantizers_indexed(struct quantizer_t *restrict lo, struct quantizer_t *restrict hi, double ratio, struct cond_quantizer_list_t *list, uint32_t column, uint32_t idx) {
    list->q[column][2*idx] = lo;
	list->q[column][2*idx + 1] = hi;
    list->ratio[column][idx] = ratio;
	list->qratio[column][idx] = (uint8_t) (ratio * 128.);
}

/**
 * Selects a quantizer for the given column from the quantizer list with the appropriate ratio
 */
struct quantizer_t *choose_quantizer(struct cond_quantizer_list_t *list, struct well_state_t *well, uint32_t column, symbol_t prev, uint32_t *q_idx) {
	uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
	assert(idx != ALPHABET_SYMBOL_NOT_FOUND);
	if (well_1024a_bits(well, 7) >= list->qratio[column][idx]) {
        *q_idx = 2*idx+1;
		return list->q[column][2*idx+1];
	}
    *q_idx = 2*idx;
	return list->q[column][2*idx];
}

/**
 * Converts a quality score into a state encoded value, which is the same as doing a symbol index lookup
 * in the output alphabet. This needs to be inlined.
 */
uint32_t find_state_encoding(struct quantizer_t *q, symbol_t value) {
	return get_symbol_index(q->output_alphabet, value);
}

/**
 * Calculates the statistics, producing a conditional pmf list per cluster and storing
 * it directly inside the cluster in question
 */
void calculate_statistics(struct quality_file_t *info) {
	uint32_t block, line_idx, column;
	uint32_t j;
	uint8_t c;
	struct line_t *line;
	struct cluster_t *cluster;
	struct cond_pmf_list_t *pmf_list;
    uint32_t read_length = 0;

	for (block = 0; block < info->block_count; ++block) {
		for (line_idx = 0; line_idx < info->blocks[block].count; ++line_idx) {
			line = &info->blocks[block].lines[line_idx];
            read_length = info->blocks[block].lines[line_idx].line_length;
			cluster = &info->clusters->clusters[line->cluster];
			pmf_list = cluster->training_stats;

			// First, find conditional PMFs
			pmf_increment(get_cond_pmf(pmf_list, 0, 0), qv2ch_mmap(line->m_data[0]));
			for (column = 1; column < read_length; ++column) {
				pmf_increment(get_cond_pmf(pmf_list, column, qv2ch_mmap(line->m_data[column-1])), qv2ch_mmap(line->m_data[column]) );
			}
		}
	}

	// Then find unconditional PMFs for each cluster once the full conditional ones are ready
	for (c = 0; c < info->cluster_count; ++c) {
		cluster = &info->clusters->clusters[c];
		pmf_list = cluster->training_stats;

		pmf_list->marginal_pmfs = alloc_pmf_list(info->columns, pmf_list->alphabet);
		combine_pmfs(get_cond_pmf(pmf_list, 0, 0), pmf_list->marginal_pmfs->pmfs[0], 1.0, 0.0, pmf_list->marginal_pmfs->pmfs[0]);
		for (column = 1; column < info->columns; ++column) {
			for (j = 0; j < pmf_list->alphabet->size; ++j) {
				combine_pmfs(pmf_list->marginal_pmfs->pmfs[column], get_cond_pmf(pmf_list, column, j), 1.0, get_probability(pmf_list->marginal_pmfs->pmfs[column-1], j), pmf_list->marginal_pmfs->pmfs[column]);
			}
		}
	}
}

/**
 * Compute the slope of the D(R)
 */
double compute_DR_slope(struct pmf_t *pmf, struct distortion_t *dist){
    
    struct pmf_t *pmf_temp = alloc_pmf(pmf->alphabet);
    struct quantizer_t *q_temp;
    int K = ALPHABET_SIZE+2,N = 2;
    
    int i, n = 0,j;
    double sum_x, sum_x2, sum_y, sum_xy;
    double prev_yi = 0, yi=1.0, prev_xi, xi;
    //double *x = (double*)calloc(K*N,sizeof(double));
    //double *y = (double*)calloc(K*N,sizeof(double));
    
    double *x = (double*)malloc(K*N*sizeof(double));
    double *y = (double*)malloc(K*N*sizeof(double));
    
    q_temp = generate_quantizer(pmf, dist, 1);
    prev_yi = q_temp->mse;
    if (prev_yi == 0) {
        // In this case we are in a zero pmf or in a one-point pmf.
        // On state will suffice in the computation of the quantizer
        free(x);
        free(y);
        free_quantizer(q_temp);
        free_pmf(pmf_temp);
        return 1.0;
    }
    prev_xi = get_entropy(apply_quantizer(q_temp, pmf, pmf_temp));
    
    //assert(prev_xi ==0.0);

    for (i = 2; i<K ; i++) {
        q_temp = generate_quantizer(pmf, dist, i);
        yi = q_temp->mse;
        if (prev_yi == yi) {
            // This should not happend, however the discrete Lloyd-max
            // as implemented is suboptimal.
            continue;
        }
        xi = get_entropy(apply_quantizer(q_temp, pmf, pmf_temp));
        if (yi == 0) {
            // We are at the last point of the D(R).
            // We need to treat it differently, as log(0)=-Inf
            // The last point is also misleading a it has different behaviour
            y[n] = log(prev_yi);
            n++;
            break;
        }
        //if (prev_yi > yi) {
            // This souldnt happend but again, the quantizer is suboptimal
            // Assuming that it does not matter, thus, do nothing
        //}
        
        // Fill the vectors with piece-wise linear functions
        for (j=0; j<N; j++) {
            //y[n] = log(prev_yi)*(1.0 - (double)j/(double)N) + log(yi)*((double)j/(double)N);
            y[n] = log( prev_yi*(1.0 - (double)j/(double)N) + yi*((double)j/(double)N) );
            x[n] = prev_xi*(1.0 - (double)j/(double)N) + xi*((double)j/(double)N);
            assert(!isnan(x[n]));
            n++;
        }
        
        //y[n] = log(yi);
        //x[n] = get_entropy(apply_quantizer(q_temp, pmf, pmf_temp));
        prev_yi = yi;
        prev_xi = xi;
        free_quantizer(q_temp);
    }
    if (n==1) {
        // We are in the case of a 2-point pmf [(0,D) (R,0)]
        // Treat this case differently
        free(x);
        free(y);
        free_quantizer(q_temp);
        free_pmf(pmf_temp);
        return 1.0;
    }
    assert(n!=0);
    
    sum_x=0.0;
    sum_x2=0.0;
    sum_xy=0.0;
    sum_y=0.0;
    for (i=0; i<n; i++) {
        sum_x += x[i];
        sum_x2 += x[i]*x[i];
        sum_xy += x[i]*y[i];
        sum_y += y[i];
    }
    free_pmf(pmf_temp);
    double b=(sum_xy-(1.0/(double)n)*sum_x*sum_y)/(sum_x2-(1.0/(double)n)*sum_x*sum_x);
    
    assert(!isnan(b));

    free(x);
    free(y);
    free_quantizer(q_temp);
    return b;
    
}

/**
 * Searches (linearly) for the pair of quantizers that surround the target Distortion by guessing and checking the number of states
 * @param pmf The pmf that is to be quantized
 * @param dist The distortion metric to quantize against
 * @param lo Place to store the pointer to the low quantizer
 * @param hi Place to store the pointer to the high quantizer
 * @return double The ratio necessary to combine these two quantizers to achieve the target
 */
double optimize_for_distortion(struct pmf_t *pmf, struct distortion_t *dist, double target, struct quantizer_t **lo, struct quantizer_t **hi) {
    
    struct quantizer_t *q_temp;
    double lo_D, hi_D;
    struct pmf_t *pmf_temp = alloc_pmf(pmf->alphabet);
    uint32_t states = 1;
    
    q_temp = generate_quantizer(pmf, dist, states);
    hi_D = q_temp->mse;
    
    if (hi_D < target) { // Target distortion too high, we are in the zero-rate case;
        *lo = generate_quantizer(pmf, dist, 1);
        *hi = generate_quantizer(pmf, dist, 1);
        
        //if (verbose) {
        //	printf("Optimization targeted zero entropy, 1 state quantizers produced automatically.\n");
        //}
        
        free_pmf(pmf_temp);
        return 1.0;
    }
    
    *hi = q_temp;
    *lo = alloc_quantizer(pmf->alphabet);
    states = 0;
    do {
        free_quantizer(*lo);
        *lo = *hi;
        lo_D = hi_D;
        
        states += 1;
        q_temp = generate_quantizer(pmf, dist, states);
        hi_D = q_temp->mse;
        *hi = q_temp;
    } while (hi_D > target && states < pmf->alphabet->size);
    
    free_pmf(pmf_temp);
    
    //hi_entropy = get_entropy(apply_quantizer(*hi, pmf, pmf_temp));
    //lo_entropy = get_entropy(apply_quantizer(*lo, pmf, pmf_temp));
    
    // Assign ratio based on how we did against our distortion target
    if (hi_D == target || hi_D == lo_D)
        return 0.0;
    else
        return (target - hi_D) / (lo_D - hi_D);
}

/**
 * Searches (linearly) for the pair of quantizers that surround the target entropy by guessing and checking the number of states
 * @param pmf The pmf that is to be quantized
 * @param dist The distortion metric to quantize against
 * @param lo Place to store the pointer to the low quantizer
 * @param hi Place to store the pointer to the high quantizer
 * @return double The ratio necessary to combine these two quantizers to achieve the target
 */
double optimize_for_entropy(struct pmf_t *pmf, struct distortion_t *dist, double target, struct quantizer_t **lo, struct quantizer_t **hi) {
	struct quantizer_t *q_temp;
	double lo_entropy, hi_entropy;
	struct pmf_t *pmf_temp = alloc_pmf(pmf->alphabet);
	uint32_t states = 1;

	if (target == 0.0) {
		*lo = generate_quantizer(pmf, dist, 1);
		*hi = generate_quantizer(pmf, dist, 1);
		
		free_pmf(pmf_temp);
		return 1.0;
	}
	
	q_temp = generate_quantizer(pmf, dist, states);
	hi_entropy = get_entropy(apply_quantizer(q_temp, pmf, pmf_temp));
	*hi = q_temp;
	*lo = alloc_quantizer(pmf->alphabet);

	do {
		free_quantizer(*lo);
		*lo = *hi;
		lo_entropy = hi_entropy;
		
		states += 1;
		q_temp = generate_quantizer(pmf, dist, states);
		hi_entropy = get_entropy(apply_quantizer(q_temp, pmf, pmf_temp));
		*hi = q_temp;
	} while (hi_entropy < target && states < pmf->alphabet->size);

	free_pmf(pmf_temp);

	// Assign ratio based on how we did against our entropy target
	if (hi_entropy < target)
		return 0.0;
	else if (lo_entropy >= target || hi_entropy == lo_entropy)
		return 1.0;
	else
		return (target - hi_entropy) / (lo_entropy - hi_entropy);
}

/**
 * 
 */
void compute_qpmf_quan_list(struct quantizer_t *q_lo, struct quantizer_t *q_hi, struct pmf_list_t *q_x_pmf, double ratio, struct alphabet_t *q_output_union) {
    symbol_t x;
    uint32_t q_symbol, idx;
    
    for (x = 0; x < q_lo->alphabet->size; x++) {
        for (idx = 0; idx < q_output_union->size; idx++) {
            q_symbol = q_output_union->symbols[idx];
            
            if (q_lo->q[x] == q_symbol)
                q_x_pmf->pmfs[x]->pmf[idx] += ratio;
            
            if (q_hi->q[x] == q_symbol)
                q_x_pmf->pmfs[x]->pmf[idx] += (1-ratio);
        }
    }
}

void compute_qpmf_list(struct pmf_list_t *qpmf_list, struct cond_pmf_list_t *in_pmfs, uint32_t column, struct pmf_list_t *prev_qpmf_list, struct alphabet_t * q_alphabet_union, struct alphabet_t * prev_q_alphabet_union, struct cond_quantizer_list_t *q_list) {
    symbol_t x;
    double p_q_xq = 0.0, p_temp = 0.0;
    uint32_t q_symbol, idx, k, j;
    struct quantizer_t *q_hi, *q_lo;
    
    // compute P(Q_i | X_i)
    for (k = 0; k < qpmf_list->size; k++) {
        // compute P(Q_i | X_i = k)
        for (idx = 0; idx < q_alphabet_union->size; idx++) {
            q_symbol = q_alphabet_union->symbols[idx];
            
            // compute P(Q_i = q_symbol | X_i = k)
            for (j = 0; j < prev_q_alphabet_union->size; j++) {
                p_q_xq = 0.0;

                // extract the jth quantizers of X_i;
                q_lo = get_cond_quantizer_indexed(q_list, column-1, 2*j);
                q_hi = get_cond_quantizer_indexed(q_list, column-1, (2*j)+1);
                
                // Given the quantizers q_lo and q_hi, compute P(Q_i = q_symbol|X_i = k ,Q_{i-1} chooses the jth quantizer of X_i)
                if (q_lo->q[k] == q_symbol)
                    p_q_xq += q_lo->ratio;
                
                if (q_hi->q[k] == q_symbol)
                    p_q_xq += q_hi->ratio;
                
                p_temp = 0;
                for (x = 0; x < prev_qpmf_list->size; ++x) {
                    p_temp += get_probability(prev_qpmf_list->pmfs[x], j) * get_probability(get_cond_pmf(in_pmfs, column-1, x), k) * get_probability(in_pmfs->marginal_pmfs->pmfs[column-2], x);
                }
                qpmf_list->pmfs[k]->pmf[idx] += p_q_xq * p_temp;
            }
        }
        
        // Normilize P(Q_i | X_i = k)
        qpmf_list->pmfs[k]->pmf_ready = 1;
        renormalize_pmf(qpmf_list->pmfs[k]);
    }
}

void compute_xpmf_list(struct pmf_list_t *qpmf_list, struct cond_pmf_list_t *in_pmfs, uint32_t column, struct pmf_list_t *xpmf_list, struct alphabet_t * q_alphabet_union){
    symbol_t x;
    uint32_t idx, k;
    
    // compute P(X_{i+1} | Q_i)
    for (idx = 0; idx < q_alphabet_union->size; idx++) {
        // compute P(X_{i+1} | Q_i = q)
        for (k = 0; k < qpmf_list->size; k++) {
            // compute P(X_{i+1} = k | Q_i = q)
            for (x = 0; x < qpmf_list->size; ++x) {
                xpmf_list->pmfs[idx]->pmf[k] += get_probability(qpmf_list->pmfs[x], idx) * get_probability(get_cond_pmf(in_pmfs, column, x), k) * get_probability(in_pmfs->marginal_pmfs->pmfs[column-1], x);
            }
        }
        // Normilize P(X_{i+1} | Q_i = q)
        xpmf_list->pmfs[idx]->pmf_ready = 1;
        renormalize_pmf(xpmf_list->pmfs[idx]);
    }
}

/**
 * For a set of already clustered data, generate codebooks for each cluster and
 * store them inside the cluster data structure
 */
void generate_codebooks(struct quality_file_t *info) {
	// Stuff for state allocation and mixing
	double ratio;
    
	// Miscellaneous variables
	uint32_t column, j;
	double total_mse;
    double hi = 0,h1 = 0,target_dist;
    
	// Output list of conditional quantizers
	struct cond_quantizer_list_t *q_list;
    
	// Constant alphabet of all possible input symbols
	const struct alphabet_t *A = info->alphabet;
    
	// Temporary/extra pointers
	struct quantizer_t *q_lo;
    struct quantizer_t *q_hi;
    
	// List of conditionally quantized PMFs after quantizer has been added out
	struct pmf_list_t *xpmf_list;
    
	// List of conditionally quantized PMFs after the next quantizer was applied
	struct pmf_list_t *qpmf_list;
	struct pmf_list_t *prev_qpmf_list;
    
	// Alphabet of all possible quantizer outputs from the previous column
	struct alphabet_t *q_output_union;
    struct alphabet_t *q_prev_output_union;

	uint8_t cluster_id;
	struct cond_pmf_list_t *in_pmfs;
	struct qv_options_t *opts = info->opts;
	struct distortion_t *dist = info->dist;

    target_dist = opts->D;
	for (cluster_id = 0; cluster_id < info->cluster_count; ++cluster_id) {
		q_list = alloc_conditional_quantizer_list(info->columns);
		info->clusters->clusters[cluster_id].qlist = q_list;
		in_pmfs = info->clusters->clusters[cluster_id].training_stats;
    
    	// For the column 0 the quantizers aren't conditional, so find them directly
    	q_output_union = alloc_alphabet(1);
    	cond_quantizer_init_column(q_list, 0, q_output_union);
		q_list->options = opts;
    
    	// Initialize the new pmfs (dummy)
    	qpmf_list = alloc_pmf_list(A->size, q_output_union);
    
    	// Handle column zero specially
		// @todo handle fixed mse target
		if (opts->mode == MODE_RATIO)
			ratio = optimize_for_entropy(get_cond_pmf(in_pmfs, 0, 0), dist, get_entropy(get_cond_pmf(in_pmfs, 0, 0))*opts->ratio, &q_lo, &q_hi);
        else{
            //h1 = compute_DR_slope(get_cond_pmf(in_pmfs, 0, 0), dist);
			ratio = optimize_for_distortion(get_cond_pmf(in_pmfs, 0, 0), dist, opts->D, &q_lo, &q_hi);
        }
		q_lo->ratio = ratio;
		q_hi->ratio = 1-ratio;
		total_mse = ratio*q_lo->mse + (1-ratio)*q_hi->mse;
    	store_cond_quantizers(q_lo, q_hi, ratio, q_list, 0, 0);
    
    	// free the used pmfs and alphabet
    	// (do not free q_prev_output_union and prev_qpmf_output as it's the first assignment).
    	q_prev_output_union = q_output_union;
    	prev_qpmf_list = qpmf_list;
    
    	// Start computing the quantizers of the rest of the columns
    	for (column = 1; column < info->columns; column++) {
        	// Compute the next output alphabet union over all quantizers for this column
			q_output_union = duplicate_alphabet(get_cond_quantizer_indexed(q_list, column-1, 0)->output_alphabet);
			for (j = 1; j < 2*q_prev_output_union->size; ++j) {
				alphabet_union(q_output_union, get_cond_quantizer_indexed(q_list, column-1, j)->output_alphabet, q_output_union);
			}
        	cond_quantizer_init_column(q_list, column, q_output_union);
        	
        	// Initialize the new pmfs
        	qpmf_list = alloc_pmf_list(A->size, q_output_union);
        	xpmf_list = alloc_pmf_list(q_output_union->size, A);
        
        	// Compute P(Q_i|X_i)
        	if (column == 1)
        	    compute_qpmf_quan_list(q_lo, q_hi, qpmf_list, ratio, q_output_union);
        	else
        	    compute_qpmf_list(qpmf_list, in_pmfs, column, prev_qpmf_list, q_output_union, q_prev_output_union, q_list);
        	
        	// Compute P(X_{i+1}|Q_i)
        	compute_xpmf_list(qpmf_list, in_pmfs, column, xpmf_list, q_output_union);
        
        	// for each previous value Q_i compute the quantizers
        	for (j = 0; j < q_output_union->size; ++j) {
        	    // Find and save quantizers
				// @todo handle fixed mse target
				if (opts->mode == MODE_RATIO)
					ratio = optimize_for_entropy(xpmf_list->pmfs[j], dist, get_entropy(xpmf_list->pmfs[j])*opts->ratio, &q_lo, &q_hi);
                else{
                    // For version 2
                    //hi = compute_DR_slope(xpmf_list->pmfs[j], dist);
                    //if (hi>0) {
                    //    hi = h1;
                    //}
                    //target_dist = (hi/h1)*opts->D;
                    ///////
                    target_dist = opts->D;
					ratio = optimize_for_distortion(xpmf_list->pmfs[j], dist,target_dist , &q_lo, &q_hi);
                }
				q_lo->ratio = ratio;
				q_hi->ratio = 1-ratio;
        	    store_cond_quantizers_indexed(q_lo, q_hi, ratio, q_list, column, j);

				// This actually needs to be scaled by the probability of this quantizer pair being used to be accurate, uniform assumption is an approximation
				total_mse = (ratio*q_lo->mse + (1-ratio)*q_hi->mse);
        	}
        
        	// deallocated the memory of the used pmfs and alphabet
        	free(q_prev_output_union);
        	q_prev_output_union = q_output_union;
	        free_pmf_list(prev_qpmf_list);
			prev_qpmf_list = qpmf_list;
	        free_pmf_list(xpmf_list);
    	}
    	
		// Final cleanup, things we saved at the end of the final iteration that aren't needed
		free_pmf_list(qpmf_list);
    	free(q_output_union);
	}
}
/**
 * Writes a codebook to a file that will be used by the decoder to initialize the arithmetic decoder
 * identically to how it was set up during encoding. The format for the file is:
 * Line 1: 1 byte ratio offset by 33 to be human readable
 * Line 2: 1 byte per quantizer symbol for column 0, low
 * Line 3: 1 byte per quantizer symbol for column 0, high
 * Lines (1+3j, 2+3j, 3+3j):
 *  1: 1 byte per ratio per unique output of previous column
 *  2: 1 byte per quantizer symbol for each low quantizer in order of symbols in previous column
 *  3: 1 byte per quantizer symbol for each high quantizer in order of symbols in previous column
 */
void write_codebook(FILE *fp, struct cond_quantizer_list_t *quantizers) {
	uint32_t i, j, k;
	uint32_t columns = quantizers->columns;

	struct quantizer_t *q_temp = get_cond_quantizer_indexed(quantizers, 0, 0);
	uint32_t size = q_temp->alphabet->size;
	uint32_t buflen = columns > size ? columns : size;
	char *eol = "\n";
	char *linebuf = (char *) _alloca(sizeof(char)*buflen);

	// First line, ratio for zero context quantizer
	linebuf[0] = quantizers->qratio[0][0] + 33;
	linebuf[1] = eol[0];
	fwrite(linebuf, sizeof(char), 2, fp);

	// Second line is low quantizer
	COPY_Q_TO_LINE(linebuf, q_temp->q, i, size);
	fwrite(linebuf, sizeof(char), size, fp);
	fwrite(eol, sizeof(char), 1, fp);

	// Third line is high quantizer
	q_temp = get_cond_quantizer_indexed(quantizers, 0, 1);
	COPY_Q_TO_LINE(linebuf, q_temp->q, i, size);
	fwrite(linebuf, sizeof(char), size, fp);
	fwrite(eol, sizeof(char), 1, fp);

	// Now for the rest of the columns, use the same format
	for (i = 1; i < columns; ++i) {
		// First a line containing ratios for each previous context
		for (j = 0; j < quantizers->input_alphabets[i]->size; ++j) {
			linebuf[j] = quantizers->qratio[i][j] + 33;
		}
		fwrite(linebuf, sizeof(char), quantizers->input_alphabets[i]->size, fp);
		fwrite(eol, sizeof(char), 1, fp);

		// Next, the low quantizers in index order
		for (j = 0; j < quantizers->input_alphabets[i]->size; ++j) {
			q_temp = get_cond_quantizer_indexed(quantizers, i, 2*j);
			COPY_Q_TO_LINE(linebuf, q_temp->q, k, size);
			fwrite(linebuf, sizeof(char), size, fp);
		}
		fwrite(eol, sizeof(char), 1, fp);

		// Finally, the high quantizers in index order
		for (j = 0; j < quantizers->input_alphabets[i]->size; ++j) {
			q_temp = get_cond_quantizer_indexed(quantizers, i, 2*j+1);
			COPY_Q_TO_LINE(linebuf, q_temp->q, k, size);
			fwrite(linebuf, sizeof(char), size, fp);
		}
		fwrite(eol, sizeof(char), 1, fp);
	}
}

/**
 * Writes all of the codebooks for the set of quantizers given, along with necessary
 * metadata (columns, lines, cluster counts) first
 */
void write_codebooks(FILE *fp, struct quality_file_t *info) {
    uint32_t columns, lines;
    uint32_t j;
    char linebuf[1];
    
    // Header line is number of clusters (1 byte)
    // number of columns (4), total number of lines (4), then a newline
    columns = htonl(info->columns);
    lines = htonl((uint32_t)info->lines);
    linebuf[0] = info->cluster_count;
    fwrite(linebuf, sizeof(char), 1, fp);
    fwrite(&columns, sizeof(uint32_t), 1, fp);
    fwrite(&lines, sizeof(uint32_t), 1, fp);
    // Now, write each cluster's codebook in order
    for (j = 0; j < info->cluster_count; ++j) {
        write_codebook(fp, info->clusters->clusters[j].qlist);
    }
}

/**
 * Reads in all of the codebooks for the clusters from the given file
 */
void read_codebooks(FILE *fp, struct quality_file_t *info) {
	uint8_t j;
	char line[9];

	// Figure out how many clusters we have to set up cluster sizes
	fread(line, sizeof(char), 9, fp);
	info->cluster_count = line[0];

	// Recover columns and lines as 32 bit integers
	info->columns = (line[1] & 0xff) | ((line[2] << 8) & 0xff00) | ((line[3] << 16) & 0xff0000) | ((line[4] << 24) & 0xff000000);
	info->columns = ntohl(info->columns);
	info->lines = (line[5] & 0xff) | ((line[6] << 8) & 0xff00) | ((line[7] << 16) & 0xff0000) | ((line[8] << 24) & 0xff000000);
	info->lines = ntohl(info->lines);
	
	// Can't allocate clusters until we know how many columns there are
	info->clusters = alloc_cluster_list(info);

	// Read codebooks in order
	for (j = 0; j < info->cluster_count; ++j) {
		info->clusters->clusters[j].qlist = read_codebook(fp, info);
	}
}

/**
 * Reads a single codebook and sets up the quantizer list
 */
struct cond_quantizer_list_t *read_codebook(FILE *fp, struct quality_file_t *info) {
	uint32_t column, size;
	uint32_t i, j;
	struct quantizer_t *q_lo, *q_hi;
	struct cond_quantizer_list_t *qlist;
	struct alphabet_t *uniques;
	char line[MAX_CODEBOOK_LINE_LENGTH];
	uint8_t qratio;
	struct alphabet_t *A = info->alphabet;
	uint32_t columns = info->columns;

	uniques = alloc_alphabet(1);
	qlist = alloc_conditional_quantizer_list(info->columns);
	cond_quantizer_init_column(qlist, 0, uniques);
	free_alphabet(uniques);

	// Next line is qratio for zero quantizer offset by 33
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	qratio = line[0] - 33;

	// Allocate some quantizers and copy the tables from lines 3 and 4
	q_lo = alloc_quantizer(A);
	q_hi = alloc_quantizer(A);
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	COPY_Q_FROM_LINE(line, q_lo->q, j, A->size);
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	COPY_Q_FROM_LINE(line, q_hi->q, j, A->size);

	// Fill in missing uniques information and store
	find_output_alphabet(q_lo);
	find_output_alphabet(q_hi);
	uniques = alloc_alphabet(0);
	alphabet_union(q_lo->output_alphabet, q_hi->output_alphabet, uniques);
	store_cond_quantizers_indexed(q_lo, q_hi, 0.0, qlist, 0, 0);
	qlist->qratio[0][0] = qratio;

	// Now handle the remaining columns uniformly
	for (column = 1; column < columns; ++column) {
		// Initialize the column information so we can write to it directly
		cond_quantizer_init_column(qlist, column, uniques);
		size = uniques->size;
		free_alphabet(uniques);
		uniques = alloc_alphabet(0);
		
		// First line is the ratios
		fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
		for (i = 0; i < size; ++i) {
			qlist->qratio[column][i] = line[i] - 33;
		}

		// Next line is a number of low quantizers
		for (i = 0; i < size; ++i) {
			q_lo = alloc_quantizer(A);
			fread(line, A->size*sizeof(symbol_t), 1, fp);
			COPY_Q_FROM_LINE(line, q_lo->q, j, A->size);
			
			find_output_alphabet(q_lo);
			qlist->q[column][2*i] = q_lo;
			alphabet_union(uniques, q_lo->output_alphabet, uniques);
		}

		// Kill the line with fgets to handle \n or \r\n automatically
		(void) fgets(line, 2, fp);

		// Next line is a number of high quantizers
		for (i = 0; i < size; ++i) {
			q_hi = alloc_quantizer(A);
			fread(line, A->size*sizeof(symbol_t), 1, fp);
			COPY_Q_FROM_LINE(line, q_hi->q, j, A->size);

			find_output_alphabet(q_hi);
			qlist->q[column][2*i+1] = q_hi;
			alphabet_union(uniques, q_hi->output_alphabet, uniques);
		}

		// Kill the line with fgets again
		(void) fgets(line, 2, fp);
	}

	// We don't use the uniques from the last column
	free_alphabet(uniques);

	return qlist;
}

/**
 * Print out a codebook by printing all of the quantizers
 */
void print_codebook(struct cond_quantizer_list_t *q) {
	struct alphabet_t *A;
	uint32_t j;
	uint32_t column;

	for (column = 0; column < q->columns; ++column) {
		A = q->input_alphabets[column];
		for (j = 0; j < 2*A->size; ++j) {
			print_quantizer(q->q[column][j]);
		}
	}
}

