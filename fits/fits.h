#include <iostream>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multimin.h>


void fit1dgaus_neldermead (gsl_vector *m, double *fit);
void fit1dgaus( gsl_vector * m, double *fit);

void fit2dgaus( gsl_matrix * m, double *fit);
void fit2dgaus_neldermead (gsl_matrix *m, double *fit);
