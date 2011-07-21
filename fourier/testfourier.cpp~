/*
 * Project:  Testing fitting routines
 *
 * Author:   Pedro M Duarte 2011-02
 * 
 */

#include "/lab/software/apparatus3/cpp/utils/utils.h"
#include "fourier.h"

using namespace std;

bool VERBOSE;

int
main (int argc, char **argv)
{

  VERBOSE = false;

  string file ("4504atoms.fits");
  gsl_matrix *img = ReadFitsImg_gsl_matrix (file);

  gsl_matrix *FT = fft2d (img);


  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string ft_path ("");
  ft_path += path;
  ft_path += "/";
  ft_path += "fft2d_";
  ft_path += "4504";
  ft_path += ".TIFF";
  toTiffImage (FT, ft_path);

  gsl_matrix_free (FT);

  cout << "testing ..." << endl;

  return 0;
}
