/*
 * Project:  App3 Data analysis 
 *
 * File:     analyze.cpp
 *
 * Author:   Pedro M Duarte 2010-11
 * 
 */

#include "/lab/software/apparatus3/cpp/utils/utils.h"
#include "/lab/software/apparatus3/cpp/fits/fits.h"
#include "/lab/software/apparatus3/cpp/vt100_macros.h"
#include <getopt.h>
#include <sstream>

#include "Basler.h"

bool VERBOSE;


int processArgsBasler (int argc, char **argv, struct params &p);

int
main (int argc, char **argv)
{
  struct params p;
  processArgsBasler (argc, argv, p);
  VERBOSE = p.verbose;

  Basler *b = new Basler (&p);
  b->LoadFLUOR ();
  b->ComputeColumnDensity ();



  b->NAtoms ();
//  setINI_num (p.reportfile, "FLUORCPP", "number", f->number);
//  setINI_num (p.reportfile, "FLUORCPP", "nsp", f->nsp);
//  setINI_num (p.reportfile, "FLUORCPP", "maxod", f->maxOD);

  b->Fit2DGauss ();
  return EXIT_SUCCESS;
}

/*
  setINI_num (p.reportfile, "FLUORCPP", "offset",
	      f->gaus2dfit[5] * f->GetNPixels ());

  double nfit = f->gaus2dfit[4] * 3.14159 * f->gaus2dfit[1] * f->gaus2dfit[3];
  double peakd = f->gaus2dfit[4] / (pow (M_PI, 0.5) * f->gaus2dfit[1]) / pow (p.magnif * 1e-4, 3);	// cm^-3

  setINI_num (p.reportfile, "FLUORCPP", "nfit", nfit);
  setINI_num (p.reportfile, "FLUORCPP", "peakd", peakd);
  setINI_num (p.reportfile, "FLUORCPP", "ax0w", f->gaus2dfit[1]);
  setINI_num (p.reportfile, "FLUORCPP", "ax1w", f->gaus2dfit[3]);
  setINI_num (p.reportfile, "FLUORCPP", "ax0c", f->abs_ci);
  setINI_num (p.reportfile, "FLUORCPP", "ax1c", f->abs_cj);
  cout << p.
    shotnum << "  nfit = " << nfit << setprecision (4) << " , peakd = " <<
    peakd << " , ax0w = " << f->gaus2dfit[1] *
    p.
    magnif << " , center = (" << setiosflags (ios::
					      fixed) << setprecision (0) <<
    f->abs_ci << "," << f->abs_cj << ")" << endl;
  cout << resetiosflags (ios::fixed) << setprecision (6);

  f->ComputeRadialAxialDensity ();

//  cout << p.shotnum << "  ncounts = " << f->number << endl;



//  double pos[2];

//  double cts = img_counts (signal);
//  double peak = img_peak (signal, pos);

//  setINI_num (p.reportfile, "FLUORCPP", "peak", peak);
//  setINI_num (p.reportfile, "FLUORCPP", "ipeak", pos[0]);
//  setINI_num (p.reportfile, "FLUORCPP", "jpeak", pos[1]);
//  cout << "#" << p.
//    shotnum << " Counts=" << img_counts (signal) << " Peak=" << peak << endl;

  return EXIT_SUCCESS;
}

*/
int
processArgsBasler (int argc, char **argv, struct params &p)
{
/*  Read command line arguments */


  if (argc == 1)
    {
      cout << endl;
      cout << "  usage:  " << argv[0] << " SHOTNUM [OPTIONS] " << endl;
      cout << endl;
      cout << set_bold (true) << "  OPTIONS " << set_bold (false) << endl <<
	endl;

      cout << set_bold (true) << "\t-C, --center [c0,c1]" << set_bold (false)
	<< endl <<
	"\t\tmanuallly specify the initial guess for the cloud center" << endl
	<< "\t\t(not implemented yet, doesn't do anything)" << endl << endl;

      cout << set_bold (true) << "\t-c, --crop" << set_bold (false) << endl <<
	"\t\tcrop images before any further analysis" << endl << endl;

      cout << set_bold (true) << "\t-f, --force" << set_bold (false) << endl
	<< "\t\tforce reanalysis of the shot" << endl << endl;

      cout << set_bold (true) << "\t-p, --plots" << set_bold (false) << endl
	<< "\t\tplots profiles when fitting" << endl << endl;

      cout << set_bold (true) << "\t-r [PATH], --ref [PATH]" <<
	set_bold (false) << endl << "\t\tindicates path of reference image" <<
	endl << endl;

      cout << set_bold (true) << "\t-R, --roi [ax0pos,ax1pos,ax0size,ax1size]"
	<< set_bold (false) << endl << "\t\tsets the atoms region of interest"
	<< endl << endl;

      cout << set_bold (true) << "\t-S, --roisize [ax0size,ax1size]"
	<< set_bold (false) << endl <<
	"\t\tsets the size for the atoms region of interest" << endl <<
	"\t\tprogram attempts to center ROI around cloud" << endl << endl;

      cout << set_bold (true) << "\t-v, --verbose" << set_bold (false) << endl
	<< "\t\tshow messages to explain what is being done" << endl << endl;

      exit (2);
    }

  makeShotPaths_Basler (argv[1], p.shotnum, p.reportfile, p.atomsfile);

  p.shot = atoi (p.shotnum.c_str ());
  p.verbose = false;
  p.reanalyze = false;
  p.center = false;
  p.crop = false;
  p.plots = false;
  p.roi_user = false;
  p.roisize_user = false;

  int c;
  string temp;
  stringstream ss (stringstream::in | stringstream::out);

  while (1)
    {
      static struct option long_options[] = {
	{"center", required_argument, 0, 'C'},
	{"crop", no_argument, 0, 'c'},
	{"force", no_argument, 0, 'f'},
	{"plots", no_argument, 0, 'p'},
	{"roi", required_argument, 0, 'R'},
	{"roisize", required_argument, 0, 'S'},
	{"verbose", no_argument, 0, 'v'},
	{0, 0, 0, 0}
      };


      int option_index = 0;
      c = getopt_long (argc, argv, "CcfpR:S:v", long_options, &option_index);
      if (c == -1)
	break;

      switch (c)
	{
	case 'C':
	  p.center = 1;
	  temp = optarg;
	  replace (temp.begin (), temp.end (), ',', '\n');
	  ss << temp;
	  ss >> (p.centerpt)[0];
	  ss >> (p.centerpt)[1];
	  break;

	case 'c':
	  p.crop = 1;
	  break;

	case 'f':
	  p.reanalyze = 1;
	  break;

	case 'p':
	  p.plots = true;
	  break;

	case 'R':
	  p.roi_user = true;
	  temp = optarg;
	  replace (temp.begin (), temp.end (), ',', '\n');
	  ss << temp;
	  ss >> (p.roi)[0];
	  ss >> (p.roi)[1];
	  ss >> (p.roi)[2];
	  ss >> (p.roi)[3];
	  break;

	case 'S':
	  p.roisize_user = true;
	  temp = optarg;
	  replace (temp.begin (), temp.end (), ',', '\n');
	  ss << temp;
	  ss >> (p.roisize)[0];
	  ss >> (p.roisize)[1];
	  break;

	case 'v':
	  p.verbose = 1;
	  break;

	case '?':
	  break;

	default:
	  abort ();
	}
    }

  if (p.verbose)
    {
      cout << "---- OPTIONS LIST ----" << endl
	<< " verbose      = " << p.verbose << endl
	<< " reanalyze    = " << p.reanalyze << endl
	<< " center       = " << p.center << endl
	<< " crop         = " << p.crop << endl
	<< " plots        = " << p.plots << endl
	<< " roi_user     = " << p.roi_user << endl
	<< " roisize_user = " << p.roisize_user << endl;
    }


  if (sectionExists (p.reportfile, "FLUORCPP") && !p.reanalyze)
    {
      cout << endl;
      cout << " Shot " << p.shotnum << " has already been analyzed." << endl;
      cout << " Use option -f  to force analysis." << endl << endl;
      exit (2);
    }
//  cout << p.shotnum << endl << p.reportfile << endl << p.
//  atomsfile << endl << p.noatomsfile << endl;


  /************* PARAMETERS OBTAINED FROM REPORT *******/
  p.texp = (double) getINI_num (p.reportfile, "BASLER", "exp");	// Exposure in milliseconds
  p.det = (double) getINI_num (p.reportfile, "BASLER", "imgdettrap");	// Image detuning in MHz

  double imgpowrep = (double) getINI_num (p.reportfile, "BASLER", "imgpowrep");	//Repump power for image
  double imgpowtrap = (double) getINI_num (p.reportfile, "BASLER", "imgpowtrap");	//Trap power for image

  if (p.det != 0.0 && p.det != -30.0)
    {
      cout << endl <<
	" Imaging is not calibrated  with the detuning used in this shot" <<
	endl << " imgdettrap = " << p.det << endl;
    }
  if (imgpowrep != 3.68 || imgpowtrap != 15.9)
    {
      cout << endl <<
	"   Imaging is not calibrated with the repump and trap powers used in this shot"
	<< endl <<
	" imgpowrep = " << imgpowrep << "   ; imgpowtrap = " << imgpowtrap <<
	endl;
      exit (2);
    }



  /************* HARD-CODED PARAMETER VALUES ***********/


  p.lambda = 670.977e-7;	// cm
  p.hc = 1.98644521e-25;	// J*m 
  p.isat = 5.1;			// mW/cm^2
  p.gamma = 5.9e6;		// Hertz
  p.tau = 27.02e-6;		// ms 


  p.quantumeff = 13.06;		//photons per count
  p.solidangle = 5.55e-4;	//spheres
  p.magnif = 0.0025;		//cm/pixel
  return EXIT_SUCCESS;
}
