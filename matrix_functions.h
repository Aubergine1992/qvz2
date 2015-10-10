//
//  matrix_functions.h
//  qvz2
//
//  Created by Mikel Hernaez-Arrazola on 10/9/15.
//  Copyright © 2015 Mikel Hernaez. All rights reserved.
//

#ifndef matrix_functions_h
#define matrix_functions_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Inverse(double** a, double **inv, int n);
void Transpose(double **a,int n);
void MxM(double** a, double** c, double**results, int n, int m, int k);

#endif /* matrix_functions_h */
