#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

#include "/lab/software/apparatus3/cpp/qini/qini.h"

#include "/usr/local/CCfits/config.h"

#include <CCfits/CCfits>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>
#include "/lab/software/apparatus3/cpp/utils/tiffio.h"

#define ROW 779
#define COL 1034

#include <math.h>

using namespace std;

int  makeShotPaths_Basler (char *shot, string & shotnum, string & report, string & atoms);
int  makeShotPaths (char *shot, string & shotnum, string & report, string & atoms, string & noatoms, string & atomsref, string & noatomsref);

int  NLines( string datafile);

bool ReadFluorImg(string & datafile, double img[ROW][COL]);
gsl_matrix * ReadFluorImg_gsl_matrix (string & datafile);
bool ReadFitsImg (string & datafile, valarray < unsigned long >&imgdata);
gsl_matrix * ReadFitsImg_gsl_matrix (string & datafile); 

bool save_gsl_matrix_ASCII(gsl_matrix *m, string & file);
void toTiffImage (gsl_matrix * m, string & filename, bool invert = false);

void to_dat_file( gsl_vector *vecs[], unsigned int N, string shot, string datfile ); 


void getmaxRowCol(gsl_matrix *m, gsl_vector * max_row, gsl_vector * max_col); 
void findpeak( gsl_matrix *m, unsigned int * i_max_ptr, unsigned int * j_max_ptr, double * max_ptr); 
void findpeak_running_avg( gsl_matrix *m, unsigned int * i_max_ptr, unsigned int * j_max_ptr, double * max_ptr, unsigned int ravg);
void findmoments(gsl_matrix *m, unsigned int *ci, unsigned int *cj,  double *peak, unsigned int *wi1e, unsigned int *wj1e);
void findcenter( gsl_matrix *m, unsigned int * i_max_ptr, unsigned int * j_max_ptr, double * max_ptr); 
void findFWHM ( gsl_matrix * m, unsigned int * FWHM_i, unsigned int * FWHM_j);

gsl_matrix *mask ( gsl_matrix *m, double factor=5);
gsl_matrix * smooth(gsl_matrix * raw, unsigned int bins);  
gsl_matrix * autocropImage(gsl_matrix * raw, double nFWHM); 
gsl_matrix * cropImage_ROI( unsigned int roi[], gsl_matrix * raw);
gsl_matrix * cropImage( string & reportfile, gsl_matrix * raw);
gsl_matrix * subtract( gsl_matrix* m1, gsl_matrix* m2);

unsigned int coerce_matrix_index( unsigned int i, unsigned int size); 

double img_counts( gsl_matrix *m);
double img_peak ( gsl_matrix *m, double *pos);

double counts2atoms(string &reportfile);
