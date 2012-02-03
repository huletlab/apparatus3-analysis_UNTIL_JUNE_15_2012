/*
 *
 *  Evaluator class for absorption images of Fermions taken with Apparatus 3
 *
 *
 *  Author : Pedro M Duarte  Feb 2011
 *
 * 
 */


// This structure contains all parameters relevant for analysis

// Convention is  size1 : 0 : radial : i
//                size2 : 1 : axial  : j 


#include <math.h>

extern bool VERBOSE;


struct params
{
  unsigned int shot;

  string shotnum, reportfile, atomsfile, noatomsfile, atomsreffile,
    noatomsreffile;

  unsigned int roi[4], roisize[2], centerpt[2];

  bool verbose, center, crop, plots, reanalyze, roi_user, roisize_user,
    fitradialaxial, fermi2d, fermiazimuth, showfermi, show_B;
  double azimuth_maxr, azimuth_chop;

  double lambda, hc, gamma, magnif, kbm, decay;

  // Image parameters
  double texp, odttof, det;

  // Evap parameters
  double free, p1, t1, tau, beta, p0, offset, t2, tau2, image, U0;
  double trapdepth, TF_factor;

  // Trap frequencies and geometry factors
  bool w_user;
  double w;
  double wx0, wy0, wz0, a, b, B1, B2;

  // Camera efficiency calibration
  double andoreff_10MHz_14bit_x1_Electron_Mult,
    andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset, eff;

  double mJ_per_photon, photon_per_mJ;

  // Absorption imaging cailbration
  double alphastar;

  // Phase-Contrast analysis
  bool phc;
  double phcsign;


};

double
b_i (double w, double t)
{
  return w * w / (1 + w * w * t * t);
};

void
init_params (struct params *p)
{
  /************* HARD-CODED PARAMETER VALUES ***********/


  p->lambda = 670.977e-7;	// cm
  p->hc = 1.98644521e-25;	// J*m 
  p->gamma = 5.9e6;		// Hertz
  p->decay = 27.2e-9;		// Seconds

  p->mJ_per_photon = 1.e3 * (p->hc / (p->lambda * 1.e-2));
  p->photon_per_mJ = 1. / p->mJ_per_photon;


  p->kbm = 1385.;		// um^2/ms^2/uK This is kb/m

  double emp = 0.88;
  p->andoreff_10MHz_14bit_x1_Electron_Mult =
    p->mJ_per_photon * (67.9 / 0.95 / emp);

  // mJ/count = (mJ per photon) * (electrons per A/D count  / QuantumEff / EmpiricalCorrectionFactor) 
  emp = 0.95;
  p->andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset =
    p->mJ_per_photon * (12.9 / 0.95 / emp);

  //p->eff = p->andoreff_10MHz_14bit_x1_Electron_Mult;
  p->eff = p->andoreff_10MHz_14bit_x5_Electron_Mult_300_BaselineOffset;

  p->magnif = 16. / 5.;		// 16um/pixel for the Andor, using the 5x obj with the telephoto at f=200

  p->U0 = 280.;

  // Calculate trap depth 
  char py[MAXPATHLEN];
  sprintf (py,
	   "evapramp2.py %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e > out.temp",
	   (p->free + p->image) / 1000., p->free, p->p1, p->t1, p->tau,
	   p->beta, p->p0, p->offset, p->t2, p->tau2, p->U0);
  std::system (py);

  ifstream f ("out.temp");
  f >> p->trapdepth;
  f >> p->TF_factor;
  f.close ();

  //DEBUG: 
  //printf("trap depth = %f\n", p->trapdepth); 
  //printf("TF_factor  = %f\n", p->TF_factor); 

  // Trap frequencies in kilohertz
  if (!p->w_user)
    p->w = 3.800;

  p->wx0 = p->w * 2 * M_PI;
  p->wy0 = p->w * 2 * M_PI;
  p->wz0 = p->w * 2 * M_PI / 8.;

  p->a = cos (52.5 * M_PI / 180.);
  p->b = sin (52.5 * M_PI / 180.);


  double ws = sqrt (p->trapdepth / p->U0);
  double wx = ws * p->wx0;
  double wy = ws * p->wy0;
  double wz = ws * p->wz0;

  //  b(w,t) = w^2/(1+w^2*t^2) 

  double bx = b_i (wx, p->odttof);
  double by = b_i (wy, p->odttof);
  double bz = b_i (wz, p->odttof);

  double a = p->a * p->a;
  double b = p->b * p->b;

  p->B1 =
    -1 * (a * b * (bz - bx) * (bz - bx) / (b * bx + a * bz) - a * bx -
	  b * bz);
  p->B2 = by;

  //DEBUG 
  if (p->show_B)
    {
      //printf("%s: w = %.3f , 1/t^2 = %.5f , B1 = %.5f, B2 = %.5f\n", (p->shotnum).c_str(), p->w, 1/(p->odttof*p->odttof), p->B1, p->B2); 
      printf ("%.3f\t%.5f\t%.5f\t%.5f\n", p->w, 1 / (p->odttof * p->odttof),
	      p->B1, p->B2);
      exit (0);
    }
};



class Fermions
{
public:

  Fermions (struct params *params)
  {
    p = params;
  }
   ~Fermions ()
  {
    return;
  }

  void LoadFITS ();
  void ComputeColumnDensity ();
  void NAtoms ();
  double GetNPixels ()
  {
    return columndensity->size1 * columndensity->size2;
  }
  void Fit1DGauss (bool col_row);
  void Fit2DGauss ();
  void FitScatt2DGauss ();
  void FitProbe2DGauss ();
  void Fit2DFermi ();
  void Get2DCuts (bool gauss, bool fermi);
  void ComputeRadialAxialDensity ();
  void GetAzimuthalAverage ();
  void GetAzimuthalAverageEllipse ();
  void FitAzimuthalFermi ();
  void MakeAzimuthalPlots ();
  void Fit1DFermi (bool radial_axial);

  struct params *p;

  double number, number_fit, maxI, maxOD, maxCD, maxPHI, maxSP, maxNSP, maxPEE, maxCA,
    minCA, maxCN, minCN, Tp0, Tsp, Nsp, gaus2dfit[6], gaus2dfit_err[6],
    scatt2dfit[6], probe2dfit[5], fermi2dfit[7], fermi_azimuth_fit[5];
  double abs_ci, abs_cj;	// centers of cloud in the uncropped pict

  double nfit, nfit_err, peakd;

  double TF;

  double TF_rd, TF_ax, TF_2d, TF_az;	// T/TF for various Fermi fits obtained from BetaMu
  double T_rd, T_ax, T_2d_rd, T_2d_ax, T_az;	// T for various Fermi fits obtained from cloud size


  double offset, peak, ci_, cj_, wi_1e, wj_1e;	// 2D Gauss parameters
  double n0, BetaMu, cj_Fermi, rj_Fermi, ci_Fermi, ri_Fermi, B_Fermi, Fugacity_Fermi, f_Fermi;	// 2D Fermi parameters 
  double n0_az, BetaMu_az, r_az, B_az, mx_az, Fugacity_az, f_az;	// Azimuthal Fermi parameters

private:
  gsl_matrix * atoms;
  gsl_matrix *noatoms;
  gsl_matrix *atomsref;
  gsl_matrix *noatomsref;

  gsl_matrix *probe;

  gsl_matrix *columndensity;
  gsl_matrix *residuals;

  gsl_matrix *scattered_ph;

  gsl_vector *cutIdata;
  gsl_vector *cutJdata;
  gsl_vector *cutIgauss;
  gsl_vector *cutJgauss;
  gsl_vector *cutIfermi;
  gsl_vector *cutJfermi;

  gsl_vector *axialdensity;
  gsl_vector *radialdensity;

  unsigned int nbins, usedbins;
  double binsize;
  gsl_vector *azimuthal_r;
  gsl_vector *azimuthal_dat;
  gsl_vector *azimuthal_gauss;
  gsl_vector *azimuthal_fermi2d;
  gsl_vector *azimuthal_fermi;



  double norm_noat;		//normalization constant for noatoms pict 
  unsigned int ci, cj, FWHMi, FWHMj;

};


void
Fermions::LoadFITS ()
{
  /********************************************
  This function loads the FITS data from disk and computes the column density
  matrices.   
 
  It takes care of any necessary cropping. 

  It calls the method ComputeColumnDensity()  to compute the column density.

  ********************************************/

  norm_noat = 1.0;
  abs_ci = 0.;
  abs_cj = 0.;

  atoms = ReadFitsImg_gsl_matrix (p->atomsfile);
  noatoms = ReadFitsImg_gsl_matrix (p->noatomsfile);
  atomsref = ReadFitsImg_gsl_matrix (p->atomsreffile);
  noatomsref = ReadFitsImg_gsl_matrix (p->noatomsreffile);

  if (p->roisize_user)
    {

      if (p->roisize[0] > atoms->size1 || p->roisize[1] > atoms->size2)
	{
	  cout << "\tERROR:  Size of ROI is larger than image" << endl;
	  exit (2);
	}

      p->roi[2] = p->roisize[0];
      p->roi[3] = p->roisize[1];
      ComputeColumnDensity ();
      if (VERBOSE)
	{
	  cout << endl << "------------ Finding ROI box ------------";
	}
      findpeak (columndensity, &ci, &cj, &peak);
      gsl_matrix_free (columndensity);

      if (VERBOSE)
	{
	  cout << "\tROI Size = " << p->roi[2] << ", " << p->roi[3] << endl;
	  cout << "\tci=" << ci << ",  cj=" << cj << endl;
	  cout << "\tSTART BOX = [ " << (int) ci -
	    (int) p->roisize[0] / 2 << " : " << (int) ci +
	    (int) p->roisize[0] / 2 << " , " << (int) cj -
	    (int) p->roisize[1] / 2 << " : " << (int) cj +
	    (int) p->roisize[1] / 2 << " ]" << endl;
	}

      while ((int) ci - (int) p->roi[2] / 2 > 0
	     && (int) ci + (int) p->roi[2] / 2 > (int) atoms->size1)
	{
	  ci--;
	}

      while ((int) cj - (int) p->roi[3] / 2 > 0
	     && (int) cj + (int) p->roi[3] / 2 > (int) atoms->size2)
	{
	  cj--;
	}

      while ((int) ci - (int) p->roi[2] / 2 < 0
	     && (int) ci + (int) p->roi[2] / 2 < (int) atoms->size1)
	{
	  ci++;
	}

      while ((int) cj - (int) p->roi[3] / 2 < 0
	     && (int) cj + (int) p->roi[3] / 2 < (int) atoms->size2)
	{
	  cj++;
	}

      if (VERBOSE)
	{
	  cout << "\tFINAL BOX = [ " << (int) ci -
	    (int) p->roisize[0] / 2 << " : " << (int) ci +
	    (int) p->roisize[0] / 2 << " , " << (int) cj -
	    (int) p->roisize[1] / 2 << " : " << (int) cj +
	    (int) p->roisize[1] / 2 << " ]" << endl;
	}

      p->roi[0] = ci - p->roisize[0] / 2;
      p->roi[1] = cj - p->roisize[1] / 2;

      p->roi_user = true;
      if (VERBOSE)
	{
	  cout << "\tROI found at " << p->roi[0] << "," << p->
	    roi[1] << "," << p->roi[2] << "," << p->roi[3] << endl;
	}

    }


  if (p->crop)
    {
      abs_ci += p->roi[0];
      abs_cj += p->roi[1];

      if (VERBOSE)
	cout << endl << "------------ Cropping Images ------------";
      gsl_matrix *catoms = p->roi_user ? cropImage_ROI (p->roi,
							atoms) :
	cropImage (p->reportfile,
		   atoms);
      gsl_matrix *cnoatoms = p->roi_user ? cropImage_ROI (p->roi,
							  noatoms) :
	cropImage (p->reportfile,
		   noatoms);
      gsl_matrix *catomsref = p->roi_user ? cropImage_ROI (p->roi,
							   atomsref) :
	cropImage (p->reportfile,
		   atomsref);
      gsl_matrix *cnoatomsref = p->roi_user ? cropImage_ROI (p->roi,
							     noatomsref) :
	cropImage (p->reportfile,
		   noatomsref);
      gsl_matrix_free (atoms);
      gsl_matrix_free (noatoms);
      gsl_matrix_free (atomsref);
      gsl_matrix_free (noatomsref);
      atoms = catoms;
      noatoms = cnoatoms;
      atomsref = catomsref;
      noatomsref = cnoatomsref;
    }


  if (VERBOSE)
    cout << endl;

  unsigned int s1 = atoms->size1;
  unsigned int s2 = atoms->size2;

  if (s1 != noatoms->size1 || s2 != noatoms->size2)
    {
      cout << "Atoms and NoAtoms matrix dimensions differ!!" << endl;
      exit (EXIT_FAILURE);
    }

  if (s1 != atomsref->size1 || s2 != atomsref->size2)
    {
      cout << "Atoms and Reference matrix dimensions differ!!" << endl;
      exit (EXIT_FAILURE);
    }

  if (s1 != noatomsref->size1 || s2 != noatomsref->size2)
    {
      cout << "Atoms and noatoms Reference matrix dimensions differ!!" <<
	endl;
      exit (EXIT_FAILURE);
    }

  //Override this values to define a different normalization region
  unsigned int imin, imax = s1;
  unsigned int jmin, jmax = s2;
  imin = imax - imax / 10;
  jmin = jmax - jmax / 10;

  gsl_matrix *norm = gsl_matrix_alloc (imax - imin, jmax - jmin);
  double an = 0.0;		// sum_{i} atoms_{i} * noatoms_{i}
  double nn = 0.0;		// sum_{i} noatoms_{i}^2 

  double a_ = 0.0;		// sum_{i} atoms_{i}
  double n_ = 0.0;		// sum_{i} noatoms_{i} 


  double at, noat, atref, noatref, c0, c1;

  for (unsigned int i = imin; i < imax; i++)
    {
      for (unsigned int j = jmin; j < jmax; j++)
	{

	  gsl_matrix_set (norm, i - imin, j - jmin,
			  gsl_matrix_get (noatoms, i, j));

	  at = gsl_matrix_get (atoms, i, j);
	  noat = gsl_matrix_get (noatoms, i, j);
	  atref = gsl_matrix_get (atomsref, i, j);
	  noatref = gsl_matrix_get (noatomsref, i, j);

	  c0 = noat - noatref;
	  c1 = at - atref;

	  an += c1 * c0;
	  nn += c0 * c0;

	  a_ += c1;
	  n_ += c0;



	}
    }

//  cout << "Norm (ave) = " << a_/n_ << endl; 
//  cout << "Norm (min) = " << an/nn << endl; 

  // Bad idea to use normalization at the moment because noatoms shot still has light 
  // coupled into the fiber for 45 ms when the probe AOM is off
  // norm_noat = an / nn;
  norm_noat = 1.0;

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string norm_path ("");
  norm_path += path;
  norm_path += "/";
  norm_path += "normpatch_";
  norm_path += p->shotnum;
  norm_path += ".TIFF";
  toTiffImage (norm, norm_path);
  gsl_matrix_free (norm);

  ComputeColumnDensity ();

  return;
}


void
Fermions::ComputeColumnDensity ()
{
  /********************************************
  Three matrices are calculted here  

  columndensity = the number of atoms in each pixel
  scattered_ph = the number of photons that were scattered per pixel
  probe = the number of photons that were incident on each pixel 

  ********************************************/

  unsigned int s1 = atoms->size1;
  unsigned int s2 = atoms->size2;


  columndensity = gsl_matrix_alloc (s1, s2);
  scattered_ph = gsl_matrix_alloc (s1, s2);
  probe = gsl_matrix_alloc (s1, s2);

  /************ ABSORPTION IMAGING VARIABLES ************/

  double at, noat, atref, noatref, c0, c1, cd, OD, i0, i0eff, p0, term1,
    term2;
  double maxT1, maxT2;		//maxima for term1 and term2 in the column density
  maxT1 = -1e6;
  maxT2 = -1e6;

  double sigma0 = 3 * M_PI * pow (p->lambda / 2 / M_PI, 2);	//This corresponds to half the maximal cross-section, because of our polarization.  
  double isat = 10.2;		// This corresponds to twice the minimal Isat because of our polarization

  //should be equal to  gamma h v / sigma0 .... test this.  


  double eff = p->eff;

  double det = 0.;		//for now all  absorption imaging is on resonance 
  maxI = 0.;

  maxCD = 0.;			//maximum column density
  maxPHI = 0.;			//maximum phase shift in phase constrast
  maxNSP = 0.;
  maxSP = 0.;
  Tp0 = 0.;			//total photons in noatoms frame 
  double sp = 0.;		//scattered photons in pixel 
  Tsp = 0.;			//total number of scattered photons 
  double nsp = 0.;		//number from scattered photons in pixel
  Nsp = 0.;			//total number from scattered photons
  double pee = 0.;
  maxPEE = 0.;

  maxOD = 0.;
  double OD_err = 5;		//any higher OD will be cutoff and trigger a warning
  bool OD_errmsg = false;
  bool high_phase_shift_flag = false;
  bool sanity_check_flag = false; 

  minCA = 1.e6;			// minimum counts in atoms frame
  minCN = 1.e6;			// minimum counts in noatoms frame
  maxCA = 0.;			// maximum counts in atoms frame
  maxCN = 0.;			// maximum counts in noatoms frame



  /************ PHASE-CONTRAST IMAGING VARIABLES ************/
  double signal, phase_shift, sig1, sig2;
  if (p->phc)
    {
      sigma0 = 6 * M_PI * pow (p->lambda / 2. / M_PI, 2);	//This is the maximal cross section
      isat = 5.1;		// This is the standard saturation intensity
      det = p->det * 1.e6 / p->gamma;	//detuning in units of gamma
      cout << endl; 
      cout << "   Andor magnification     = " << p->magnif << " um/pixel" << endl;
      cout << "   Phase-contrast detuning = " << det << " Gamma  = " << det * p->gamma/1.e6 << " MHz" << endl;
    }

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{

	  at = gsl_matrix_get (atoms, i, j);
	  noat = gsl_matrix_get (noatoms, i, j);
	  atref = gsl_matrix_get (atomsref, i, j);
	  noatref = gsl_matrix_get (noatomsref, i, j);

	  c0 = norm_noat * (noat - noatref);	// counts in no atoms frame
	  c1 = at - atref;	// counts in atoms frame

	  if (c0 < minCN)
	    minCN = c0;
	  if (c0 > maxCN)
	    maxCN = c0;

	  if (c1 < minCA)
	    minCA = c1;
	  if (c1 > maxCA)
	    maxCA = c1;

	  // eff is in mJ/count
	  p0 = c0 * p->eff * p->photon_per_mJ;	// number of photons in pixel in noatoms frame
	  Tp0 += p0;		// total number of photons in noatoms frame 
	  i0 = (c0 * eff / (pow (p->magnif * 1e-4, 2)) / p->texp / (isat));	// intensity in noatoms frame in units of isat
	  i0eff = (c0 * eff / (pow (p->magnif * 1e-4, 2)) / p->texp / (isat * p->alphastar));	// effective intensity in noatoms frame
	  if (i0 > maxI)
	    maxI = i0;

	  /*********** ABSORPTION IMAGING ***************/

	  if (!p->phc)
	    {
	      OD = log (fabs (c0 / c1));

	      if (OD > OD_err)
		{
		  OD = OD_err;

		  if (!OD_errmsg)
		    cout << endl <<
		      "  ******   WARNING: Optical density too high !!! ******"
		      << endl;
		  OD_errmsg = true;
		}

	      else if (OD > maxOD)
		maxOD = OD;

	      term1 = p->alphastar * (1. + 4. * det * det) * OD;


	      pee = i0eff / (1 + 2 * i0eff);
	      if (pee > maxPEE)
		maxPEE = pee;
	      //sp = scattered photons in pixel
	      sp = (c0 - c1) * eff * p->photon_per_mJ;
	      if (sp > maxSP)
		maxSP = sp;
	      //nsp = number of atoms in pixel, calculated from number of scattered photons
	      nsp = sp / (pee * p->texp / p->decay);
	      if (nsp > maxNSP)
		maxNSP = nsp;

	      Tsp += sp;	// total number of scattered photons
	      Nsp += nsp;	// total number of atoms from scattered photons 

	      term2 = 2 * i0 * (1. - c1 / c0);

	      if (term1 > maxT1)
		maxT1 = term1;
	      if (term2 > maxT2)
		maxT2 = term2;

	      cd = (term1 + term2) / sigma0 * pow (p->magnif * 1e-4, 2);

	    }

	  /*********** PHASE - CONTRAST IMAGING ***************/

	  else
	    {
	      signal = p->phcsign * (c1 / c0 - 1);
              phase_shift = -1.0 * det * signal / (0.5 - det);

              //cout << "signal = " << signal << " phi = " << phase_shift << endl; 
                
              if ( fabs(phase_shift) > maxPHI ) maxPHI = phase_shift;  

	      //if (signal / (0.5 - det) > 0.5 || signal > 0.6)
	      if ( phase_shift > M_PI / 5.)
		{
		  if (!high_phase_shift_flag)
		    cout << endl <<
		      "  ******   WARNING: Maximum phase shift exceeds pi/5  !!! ******" <<
		      endl;
		  high_phase_shift_flag = true;
		}

	      i0 = 2 * i0;	// The polarizer cuts it down by a factor of 2.  

	      sig1 = (0.5 - det) * sigma0 / pow(p->magnif * 1e-4,2) / (1 + 4 * det * det + 2 * i0);
	      sig2 =
		(0.5 - det - 13.3) * sigma0 / pow(p->magnif * 1e-4,2) / (1 + 4 * det * det + 2 * i0);

	      cd = 2 * signal / (sig1 + sig2);	//columdensity 
	      nsp = 0.;		//this is not relevant for phase contrast so just make it zero.

	      // sanity check
              /*if ( !sanity_check_flag){
                    printf("det    = %.3e * Gamma\n", det );
                    printf("det^2  = %.3e * Gamma^2\n", det*det );
                    printf("i0     = %.3e * Isat\n", i0 );
                    printf("sigma0 = %.3e\n", sigma0 );
                    printf("sig1   = %.3e\n", sig1);
                    printf("sig2   = %.3e\n", sig2);
                    printf("cd     = %.3e\n", cd);
                    printf("signal = %.3e\n",signal);
              }*/

	    }
	  if (cd > maxCD)
	    maxCD = cd;

	  gsl_matrix_set (scattered_ph, i, j, nsp);
	  gsl_matrix_set (columndensity, i, j, cd);
	  gsl_matrix_set (probe, i, j, c0);
	}
    }

  maxT1 = maxT1 / sigma0 * pow (p->magnif * 1e-4, 2);
  maxT2 = maxT2 / sigma0 * pow (p->magnif * 1e-4, 2);


  if (VERBOSE)
    {
      cout << endl << "------------ Column Density Stats ------------" <<
	endl;
      if (!p->phc)
	{

	  cout << endl << "----> Method used : absorption" << endl;
	  cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
	  cout << "\tmax probe intensity = " << maxI << " Isat " << endl;
	  cout << "\tmax PHI = " << maxPHI << endl;
	  cout << "\tmax CD = " << maxCD << endl;
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n", minCA,
		  maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n", minCN,
		  maxCN);
	  cout << "\ntotal number of photons in probe pulse = " << Tp0 <<
	    endl;
	  cout << "First few elements of column density matrix: " << endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (columndensity, i, j));
		}
	      cout << "... " << endl;
	    }


	}
      else
	{
	  cout << endl << "----> Method used : phase contrast" << endl;
	  cout << "\ts1=" << s1 << ", s2=" << s2 << endl;
	  cout << "\talphastar = " << p->alphastar << endl;
	  cout << "\tmax term1 = " << maxT1 << endl;
	  cout << "\tmax term2 = " << maxT2 << endl;
	  cout << "\tmax OD = " << maxOD << endl;
	  cout << "\tmax probe intensity = " << maxI << " Isat " << endl;
	  cout << "\tmax CD = " << maxCD << endl;
	  cout << "\tmax scattered photons = " << maxSP << endl;
	  cout << "\tmax number from scattered photons = " << maxNSP << endl;
	  cout << "\tmax rho_{ee} = " << maxPEE << endl;
	  printf ("\tcounts in atoms picture:     ( %f to %f )\n", minCA,
		  maxCA);
	  printf ("\tcounts in no atoms picture:  ( %f to %f )\n", minCN,
		  maxCN);
	  cout << "\ttotal scattered photons = " << Tsp << endl;
	  cout << "\ttotal number from scattered photons = " << Nsp << endl;
	  cout << "\ntotal number of photons in probe pulse = " << Tp0 <<
	    endl;
	  cout << "First few elements of column density matrix: " << endl;
	  for (int i = 0; i < 5; i++)
	    {
	      for (int j = 0; j < 5; j++)
		{
		  printf ("%.5f\t ", gsl_matrix_get (columndensity, i, j));
		}
	      cout << "... " << endl;
	    }
	}

    }

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string column_path ("");
  column_path += path;
  column_path += "/";

  string scatt_path ("");
  scatt_path = column_path;

  column_path += "column_";
  column_path += p->shotnum;
  string column_ascii_path ("");
  column_ascii_path = column_path;
  column_path += ".TIFF";
  column_ascii_path += ".ascii";

  toTiffImage (columndensity, column_path);
  save_gsl_matrix_ASCII (columndensity, column_ascii_path);

  scatt_path += "scatt_";
  scatt_path += p->shotnum;
  string scatt_path_ASCII ("");
  scatt_path_ASCII = scatt_path;
  scatt_path += ".TIFF";
  scatt_path_ASCII += ".ascii";

  toTiffImage (scattered_ph, scatt_path);
  save_gsl_matrix_ASCII (scattered_ph, scatt_path_ASCII);



  return;
}

void
Fermions::NAtoms ()
{
  unsigned int s1 = columndensity->size1;
  unsigned int s2 = columndensity->size2;

  double num = 0.;
  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  num += gsl_matrix_get (columndensity, i, j);
    }}
  if (VERBOSE)
    {
      cout << endl << "------------ NAtoms ------------" << endl;
      cout << "\tnumber = " << num << endl;
    }
  number = num;

  return;

}

void
Fermions::Fit1DGauss (bool col_row)
{


  //True fits a column , False fits a row
  unsigned int ii = col_row ? cj : ci;
  if (VERBOSE)
    {
      cout << endl << "------------ Fitting 1D Gaussian ------------" << endl;
      cout << "Column density ( " << columndensity->
	size1 << ", " << columndensity->size2 << ")  ";
      if (col_row)
	cout << "Fitting COLUMN " << ii << endl;
      if (!col_row)
	cout << "Fitting ROW " << ii << endl;
    }

  if (col_row && ii > columndensity->size2 || ii < 0)
    {
      cout << "Cannot fit column with index " << ii << " : out of bounds!" <<
	endl;
      return;
    }
  if (!col_row && ii > columndensity->size1 || ii < 0)
    {
      cout << "Cannot fit row with index " << ii << " : out of bounds!" <<
	endl;
      return;
    }

  unsigned int s = col_row ? columndensity->size1 : columndensity->size2;

  gsl_vector *profile = gsl_vector_alloc (s);
  for (unsigned int index = 0; index < s; index++)
    {
      unsigned int i = col_row ? index : ii;
      unsigned int j = col_row ? ii : index;
      gsl_vector_set (profile, index, gsl_matrix_get (columndensity, i, j));
    }

  //findpeak (columndensity, &ci, &cj, &peak);
  //findFWHM (columndensity, &FWHMi, &FWHMj);

  double center = (double) (col_row ? ci : cj);
  double FWHM = (double) (col_row ? FWHMi : FWHMj);


  if (VERBOSE)
    cout << endl;
  double gausfit1d[4] = { center, FWHM, peak, 0.1 };
  fit1dgaus (profile, gausfit1d);

  gausfit1d[1] = fabs (gausfit1d[1]);

  if (col_row)
    {
      ci_ = gausfit1d[0];
      wi_1e = gausfit1d[1];
    }
  if (!col_row)
    {
      cj_ = gausfit1d[0];
      wj_1e = gausfit1d[1];
    }

  gsl_vector_free (profile);

  return;
}



void
Fermions::Fit2DGauss ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT COLUMN DENSITY WITH 2D GAUSISAN ------------" << endl;
  //findpeak (columndensity, &ci, &cj, &peak);
  findcenter (columndensity, &ci, &cj, &peak);
  findFWHM (columndensity, &FWHMi, &FWHMj);

  if (VERBOSE)
    cout << endl;
  Fit1DGauss (false);		//False is to fit row
  if (VERBOSE)
    cout << endl;
  Fit1DGauss (true);		//True is to fit column

  gaus2dfit[0] = ci_;
  gaus2dfit[1] = wi_1e;
  gaus2dfit[2] = cj_;
  gaus2dfit[3] = wj_1e;
  gaus2dfit[4] = peak;
  gaus2dfit[5] = 0.1;

  gaus2dfit_err[0] = 1e15;
  gaus2dfit_err[1] = 1e15;
  gaus2dfit_err[2] = 1e15;
  gaus2dfit_err[3] = 1e15;
  gaus2dfit_err[4] = 1e15;
  gaus2dfit_err[5] = 1e15;


  if (VERBOSE)
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;
  fit2dgaus_err (columndensity, gaus2dfit, gaus2dfit_err);
  //Override for debugging
  /* gaus2dfit[0] = 278.8;
     gaus2dfit[1] = 79.5;
     gaus2dfit[2] = 291.5;
     gaus2dfit[3] = 92.0;
     gaus2dfit[4] = 6.876e+01;
     gaus2dfit[5] = 0.58; */
  if (VERBOSE)
    cout << endl;

  gaus2dfit[1] = fabs (gaus2dfit[1]);
  gaus2dfit[3] = fabs (gaus2dfit[3]);

  ci_ = gaus2dfit[0];
  wi_1e = gaus2dfit[1];
  cj_ = gaus2dfit[2];
  wj_1e = gaus2dfit[3];
  peak = gaus2dfit[4];
  offset = gaus2dfit[5];

  nfit = gaus2dfit[4] * M_PI * gaus2dfit[1] * gaus2dfit[3];

  TF = p->TF_factor * pow (6 * nfit, 1. / 3.);

  nfit_err =
    pow (pow (gaus2dfit_err[4] * nfit / gaus2dfit[4], 2) +
	 pow (gaus2dfit_err[1] * nfit / gaus2dfit[1],
	      2) + pow (gaus2dfit_err[3] * nfit / gaus2dfit[3], 2), 0.5);

  peakd = gaus2dfit[4] / (pow (M_PI, 0.5) * gaus2dfit[1]) / pow (p->magnif * 1e-4, 3);	// cm^-3

  abs_ci += ci_;
  abs_cj += cj_;

  if (VERBOSE)
    {
      printf ("..............  GAUSSIAN 2D FIT RESULTS ..............\n");
      printf ("ci     = %.1f pixels\n", gaus2dfit[0]);
      printf ("wi     = %.1f pixels\n", gaus2dfit[1]);
      printf ("cj     = %.1f pixels\n", gaus2dfit[2]);
      printf ("wj     = %.1f pixels\n", gaus2dfit[3]);
      printf ("peak   = %.3e \n", gaus2dfit[4]);
      printf ("offset = %.3e \n", gaus2dfit[5]);
      printf ("N from fit = %.3e \n", nfit);
    }

  residuals = gsl_matrix_alloc (columndensity->size1, columndensity->size2);

  ci = (unsigned int) floor (ci_) - 1;
  cj = (unsigned int) floor (cj_) - 1;

  ci = coerce_matrix_index (ci, columndensity->size1);
  cj = coerce_matrix_index (cj, columndensity->size2);

  if (VERBOSE)
    {
      cout << "ci = " << ci << " ; cj = " << cj << endl;
    }

  number_fit = 0.0;

  double model, data;
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;
	  double xj = (double) j;
	  model =
	    offset +
	    peak * exp (-1. *
			(pow ((xi - ci_) / wi_1e, 2.) +
			 pow ((xj - cj_) / wj_1e, 2.)));
	  data = gsl_matrix_get (columndensity, i, j);
	  gsl_matrix_set (residuals, i, j, data - model);
	  number_fit += model;

	}
    }
  return;
}

void
Fermions::FitScatt2DGauss ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT SCATTERED PHOTONS WITH 2D GAUSISAN ------------" <<
      endl;
  if (VERBOSE)
    cout << endl;

  scatt2dfit[0] = gaus2dfit[0];
  scatt2dfit[1] = gaus2dfit[1];
  scatt2dfit[2] = gaus2dfit[2];
  scatt2dfit[3] = gaus2dfit[3];
  scatt2dfit[4] = maxNSP;
  scatt2dfit[5] = 0.1;

  if (VERBOSE)
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;
  fit2dgaus (scattered_ph, scatt2dfit);
  //Override for debugging
  /* gaus2dfit[0] = 278.8;
     gaus2dfit[1] = 79.5;
     gaus2dfit[2] = 291.5;
     gaus2dfit[3] = 92.0;
     gaus2dfit[4] = 6.876e+01;
     gaus2dfit[5] = 0.58; */
  if (VERBOSE)
    cout << endl;

  scatt2dfit[1] = fabs (scatt2dfit[1]);
  scatt2dfit[3] = fabs (scatt2dfit[3]);

  double nscattph = scatt2dfit[4] * M_PI * scatt2dfit[1] * scatt2dfit[3];

  if (VERBOSE)
    {
      printf
	("..............  SCATTERED PHOTONS 2D FIT RESULTS ..............\n");
      printf ("ci     = %.1f pixels\n", scatt2dfit[0]);
      printf ("wi     = %.1f pixels\n", scatt2dfit[1]);
      printf ("cj     = %.1f pixels\n", scatt2dfit[2]);
      printf ("wj     = %.1f pixels\n", scatt2dfit[3]);
      printf ("peak   = %.3e \n photons", scatt2dfit[4]);
      printf ("offset = %.3e \n", scatt2dfit[5]);
      printf ("N from scattered photons = %.3e \n", nscattph);	// / (p->hc * 1e3 / (p->lambda * 1e-2)) / ( ( maxI/ (1+2*maxI) ) * 2 * M_PI * 5.9e6 * p->texp ) );
    }

  return;
}

void
Fermions::FitProbe2DGauss ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT PROBE WITH 2D GAUSISAN ------------" << endl;
  if (VERBOSE)
    cout << endl;

  probe2dfit[0] = gaus2dfit[0];
  probe2dfit[1] = gaus2dfit[1];
  probe2dfit[2] = gaus2dfit[2];
  probe2dfit[3] = gaus2dfit[3];
  probe2dfit[4] = maxCN;

  if (VERBOSE)
    cout << endl << "------------ Fitting with 2D Gaussian ------------" <<
      endl;
  fit2dgaus_no_offset (probe, probe2dfit);
  //Override for debugging
  /* gaus2dfit[0] = 278.8;
     gaus2dfit[1] = 79.5;
     gaus2dfit[2] = 291.5;
     gaus2dfit[3] = 92.0;
     gaus2dfit[4] = 6.876e+01;
     gaus2dfit[5] = 0.58; */
  if (VERBOSE)
    cout << endl;

  probe2dfit[1] = fabs (probe2dfit[1]);
  probe2dfit[3] = fabs (probe2dfit[3]);


  if (VERBOSE)
    {
      printf ("..............  PROBE 2D FIT RESULTS ..............\n");
      printf ("ci     = %.1f pixels\n", probe2dfit[0]);
      printf ("wi     = %.1f pixels\n", probe2dfit[1]);
      printf ("cj     = %.1f pixels\n", probe2dfit[2]);
      printf ("wj     = %.1f pixels\n", probe2dfit[3]);
      printf ("peak   = %.3e \n photons", probe2dfit[4]);
    }

  return;
}

void
Fermions::Get2DCuts (bool gauss, bool fermi)
{

  if (VERBOSE)
    cout << endl <<
      "------------ PRODUCE 1D CUTS FROM RESULTS OF 2D FITS ------------" <<
      endl;

  cutIdata = gsl_vector_alloc (columndensity->size1);
  cutJdata = gsl_vector_alloc (columndensity->size2);

  if (gauss)
    {
      cutIgauss = gsl_vector_alloc (columndensity->size1);
      cutJgauss = gsl_vector_alloc (columndensity->size2);
    }

  if (fermi)
    {
      cutIfermi = gsl_vector_alloc (columndensity->size1);
      cutJfermi = gsl_vector_alloc (columndensity->size2);
    }


  double model_gauss, model_fermi, data;
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;
	  double xj = (double) j;
	  if (gauss)
	    {
	      model_gauss =
		offset +
		peak * exp (-1. *
			    (pow ((xi - ci_) / wi_1e, 2.) +
			     pow ((xj - cj_) / wj_1e, 2.)));
	    }

	  if (fermi)
	    {
	      model_fermi =
		n0 / f1 (BetaMu) * f1 (BetaMu -
				       fq (BetaMu) *
				       (pow ((xj - cj_Fermi) / rj_Fermi, 2) +
					pow ((xi - ci_Fermi) / ri_Fermi,
					     2))) + B_Fermi;
	    }

	  data = gsl_matrix_get (columndensity, i, j);
	  if (i == ci)
	    {
	      gsl_vector_set (cutJdata, j, data);
	      if (gauss)
		gsl_vector_set (cutJgauss, j, model_gauss);
	      if (fermi)
		gsl_vector_set (cutJfermi, j, model_fermi);
	    }
	  if (j == cj)
	    {
	      gsl_vector_set (cutIdata, i, data);
	      if (gauss)
		gsl_vector_set (cutIgauss, i, model_gauss);
	      if (fermi)
		gsl_vector_set (cutIfermi, i, model_fermi);
	    }


	}
    }

  if (gauss)
    {
      gsl_vector *cuts[4] = { cutIdata, cutJdata, cutIgauss, cutJgauss };
      to_dat_file (cuts, 4, p->shotnum, "cuts.dat");
    }

  if (fermi)
    {
      gsl_vector *cuts[6] =
	{ cutIdata, cutJdata, cutIgauss, cutJgauss, cutIfermi, cutJfermi };
      to_dat_file (cuts, 6, p->shotnum, "cuts.dat");
    }

  //--- On screen plot
  std::ofstream gpl;
  gpl.open ("temp.gpl");
  gpl << "set size 1.0,0.45" << endl;
  gpl << "set multiplot" << endl;
  gpl << "set origin 0.0,0.0" << endl;

  gpl << "plot \"";
  gpl << p->shotnum << "_cuts.dat\" u 0:1 pt 7 ps 1 notit";
  if (gauss)
    gpl << ",\\" << endl << "\"\" u 0:3 w lines title \"gauss FIT I\"";
  if (fermi)
    gpl << ",\\" << endl << "\"\" u 0:5 w lines title \"fermi FIT I\"";
  gpl << endl;

  gpl << "set origin 0.0,0.5" << endl;

  gpl << "plot \"";
  gpl << p->shotnum << "_cuts.dat\" u 0:2 pt 7 ps 1 notit";
  if (gauss)
    gpl << ",\\" << endl << "\"\" u 0:4 w lines title \"gauss FIT J\"";
  if (fermi)
    gpl << ",\\" << endl << "\"\" u 0:6 w lines title \"fermi FIT J\"";
  gpl << endl;

  gpl << "unset multiplot" << endl;
  if (p->plots)
    std::system ("gnuplot -persist temp.gpl");
  gpl.close ();
  remove ("temp.gpl");

  //--- PNG plot
  gpl.open ("temp.gpl");
  gpl << "set terminal png medium" << endl;
  gpl << "set output \"cuts_" << p->shotnum << ".png\"" << endl;
  gpl << "set size 1.0,1.0" << endl;
  gpl << "set multiplot" << endl;
  gpl << "set size 1.0,0.45" << endl;
  gpl << "set origin 0.0,0.0" << endl;

  gpl << "plot \"";
  gpl << p->shotnum << "_cuts.dat\" u 0:1 pt 7 ps 1 notit";
  if (gauss)
    gpl << ",\\" << endl << "\"\" u 0:3 w lines title \"gauss FIT I\"";
  if (fermi)
    gpl << ",\\" << endl << "\"\" u 0:5 w lines title \"fermi FIT I\"";
  gpl << endl;

  gpl << "set origin 0.0,0.5" << endl;

  gpl << "plot \"";
  gpl << p->shotnum << "_cuts.dat\" u 0:2 pt 7 ps 1 notit";
  if (gauss)
    gpl << ",\\" << endl << "\"\" u 0:4 w lines title \"gauss FIT J\"";
  if (fermi)
    gpl << ",\\" << endl << "\"\" u 0:6 w lines title \"fermi FIT J\"";
  gpl << endl;

  gpl << "unset multiplot" << endl;
  gpl.close ();
  std::system ("gnuplot temp.gpl");
  remove ("temp.gpl");

  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);
  string residuals_path ("");
  residuals_path += path;
  residuals_path += "/";
  residuals_path += "res_";
  residuals_path += p->shotnum;
  residuals_path += ".TIFF";
  toTiffImage (residuals, residuals_path);

  return;
}


void
Fermions::Fit2DFermi ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT COLUMN DENSITY WITH 2D FERMI-DIRAC ------------" <<
      endl;

  fermi2dfit[0] = peak;
  fermi2dfit[1] = -5.0;
  fermi2dfit[2] = wi_1e;
  fermi2dfit[3] = wj_1e;
  fermi2dfit[4] = ci_;
  fermi2dfit[5] = cj_;
  fermi2dfit[6] = offset;

  //override so that it doesn't take forever
/*  fermi2dfit[0] = 1.63e+02;
  fermi2dfit[1] = 1.76e+00;
  fermi2dfit[2] = 68.4;
  fermi2dfit[3] = 87.9;
  fermi2dfit[4] = 265.0;
  fermi2dfit[5] = 297.4;
  fermi2dfit[7] = 1.38; */
  fit2dfermi_neldermead (columndensity, fermi2dfit);
  if (VERBOSE)
    cout << endl;

  fermi2dfit[2] = fabs (fermi2dfit[2]);
  fermi2dfit[3] = fabs (fermi2dfit[3]);

  BetaMu = fermi2dfit[1];

  TF_2d = pow (6 * f2 (BetaMu), -0.333);
  Fugacity_Fermi = exp (BetaMu);
  f_Fermi = fq (BetaMu);




  ri_Fermi = fermi2dfit[2];
  rj_Fermi = fermi2dfit[3];

  T_2d_rd =
    p->B1 / (2 * p->kbm) * ri_Fermi * ri_Fermi * p->magnif * p->magnif /
    f_Fermi;
  T_2d_ax =
    p->B2 / (2 * p->kbm) * rj_Fermi * rj_Fermi * p->magnif * p->magnif /
    f_Fermi;

  if (VERBOSE || p->showfermi)
    {
      printf ("..............  FERMI 2D FIT RESULTS ..............\n");
      printf ("n0     = %.2e\n", fermi2dfit[0]);
      printf ("BetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	      fermi2dfit[1], f2 (fermi2dfit[1]), TF_2d);
      printf ("ri     = %.1f pixels\n", fermi2dfit[2]);
      printf ("rj     = %.1f pixels\n", fermi2dfit[3]);
      printf ("ci     = %.1f pixels\n", fermi2dfit[4]);
      printf ("cj     = %.1f pixels\n", fermi2dfit[5]);
      printf ("B      = %.2e \n", fermi2dfit[6]);
      printf ("fugacity   = %.2e\n", Fugacity_Fermi);
      printf ("f(z)       = %.2e\n", f_Fermi);
      printf ("T_2d_rd   = %.2f uK\n", T_2d_rd);
      printf ("T_2d_ax   = %.2f uK\n", T_2d_ax);
      printf ("TF     = %.2f uK\n", TF);
    }

  n0 = fermi2dfit[0];
  ci_Fermi = fermi2dfit[4];
  cj_Fermi = fermi2dfit[5];
  B_Fermi = fermi2dfit[6];

  return;
}

void
Fermions::GetAzimuthalAverage ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ GET AZIMUTHAL AVERAGE OF CLOUD ------------" << endl;


  nbins = 1024;
  //ci_Fermi and cj_Fermi are the centers obtained by the 2D Fermi fit
  // 0 is sum of values
  // 1 is N values
  gsl_vector *azim_histo_0 = gsl_vector_alloc (nbins);
  gsl_vector *azim_histo_1 = gsl_vector_alloc (nbins);

  for (unsigned int j = 0; j < nbins; j++)
    {
      gsl_vector_set (azim_histo_0, j, 0.0);
      gsl_vector_set (azim_histo_1, j, 0.0);
    }

  binsize = 1.2;

  //find all possible distances 
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;
	  double xj = (double) j;
	  unsigned int dist =
	    (unsigned int)
	    floor (pow (pow (xi - ci_, 2) + pow (xj - cj_, 2), 0.5) /
		   binsize);
	  //DEBUG: printf( " dist(%d,%d) = %d\n", i, j, dist); 

	  gsl_vector_set (azim_histo_0, dist,
			  gsl_vector_get (azim_histo_0,
					  dist) +
			  gsl_matrix_get (columndensity, i, j));

	  gsl_vector_set (azim_histo_1, dist,
			  gsl_vector_get (azim_histo_1, dist) + 1);

	}
    }
  // After computing the sums, divide by the number to get the azimuthal average
  // 0 is the distance in pixels
  // 1 is the data
  // 2 is the 2D Gauss Fit
  // 3 is the 2D Fermi Fit
  usedbins = 1;
  while (gsl_vector_get (azim_histo_1, usedbins) != 0)
    usedbins++;

  if (usedbins * binsize > p->azimuth_maxr)
    {
      usedbins = (unsigned int) ceil (p->azimuth_maxr / binsize);
    }
  else if (p->azimuth_maxr == 512.)
    {
      usedbins -= (unsigned int) ceil (p->azimuth_chop / binsize);
    }

  azimuthal_r = gsl_vector_alloc (usedbins);
  azimuthal_dat = gsl_vector_alloc (usedbins);

  double r;
  for (unsigned int index = 0; index < usedbins; index++)
    {
      r = index * binsize;
      double hsum = gsl_vector_get (azim_histo_0, index);
      double hnum = gsl_vector_get (azim_histo_1, index);
      //DEBUG: if (hnum == 0) { printf("ERROR: no entries in azimuthal histogram\n");}
      gsl_vector_set (azimuthal_r, index, r);
      gsl_vector_set (azimuthal_dat, index, hsum / hnum);
    }

  return;
}

void
Fermions::GetAzimuthalAverageEllipse ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ GET AZIMUTHAL AVERAGE OF CLOUD (ELLIPSE) ------------" <<
      endl;


  nbins = 1024;
  //ci_Fermi and cj_Fermi are the centers obtained by the 2D Fermi fit
  // 0 is sum of values
  // 1 is N values
  gsl_vector *azim_histo_0 = gsl_vector_alloc (nbins);
  gsl_vector *azim_histo_1 = gsl_vector_alloc (nbins);

  for (unsigned int j = 0; j < nbins; j++)
    {
      gsl_vector_set (azim_histo_0, j, 0.0);
      gsl_vector_set (azim_histo_1, j, 0.0);
    }

  binsize = 1.2;

  //find all possible distances 
  for (unsigned int i = 0; i < columndensity->size1; i++)
    {
      for (unsigned int j = 0; j < columndensity->size2; j++)
	{
	  double xi = (double) i;	// this is y  radial
	  double xj = (double) j;	// this is g  axial
	  unsigned int dist =
	    (unsigned int)
	    floor (pow
		   ((p->B2 / p->B1) * pow (xi - ci_, 2) + pow (xj - cj_, 2),
		    0.5) / binsize);
	  //DEBUG: printf( " dist(%d,%d) = %d\n", i, j, dist); 

	  gsl_vector_set (azim_histo_0, dist,
			  gsl_vector_get (azim_histo_0,
					  dist) +
			  gsl_matrix_get (columndensity, i, j));

	  gsl_vector_set (azim_histo_1, dist,
			  gsl_vector_get (azim_histo_1, dist) + 1);

	}
    }
  // After computing the sums, divide by the number to get the azimuthal average
  // 0 is the distance in pixels
  // 1 is the data
  // 2 is the 2D Gauss Fit
  // 3 is the 2D Fermi Fit
  usedbins = 1;
  while (gsl_vector_get (azim_histo_1, usedbins) != 0)
    usedbins++;

  if (usedbins * binsize > p->azimuth_maxr)
    {
      usedbins = (unsigned int) ceil (p->azimuth_maxr / binsize);
    }
  else if (p->azimuth_maxr == 512.)
    {
      usedbins -= (unsigned int) ceil (p->azimuth_chop / binsize);
    }

  azimuthal_r = gsl_vector_alloc (usedbins);
  azimuthal_dat = gsl_vector_alloc (usedbins);

  double r;
  for (unsigned int index = 0; index < usedbins; index++)
    {
      r = index * binsize;
      double hsum = gsl_vector_get (azim_histo_0, index);
      double hnum = gsl_vector_get (azim_histo_1, index);
      //DEBUG: if (hnum == 0) { printf("ERROR: no entries in azimuthal histogram\n");}
      gsl_vector_set (azimuthal_r, index, r);
      gsl_vector_set (azimuthal_dat, index, hsum / hnum);
    }

  return;
}

void
Fermions::FitAzimuthalFermi ()
{
  if (VERBOSE)
    cout << endl <<
      "------------ FIT AZIMUTHAL AVERAGE WITH FERMI-DIRAC ------------" <<
      endl;


  fermi_azimuth_fit[0] = peak;
  fermi_azimuth_fit[1] = -5.0;
  //fermi_azimuth_fit[2] = pow(wi_1e*wj_1e,0.5); 
  fermi_azimuth_fit[2] = wj_1e;
  fermi_azimuth_fit[3] = offset;
  fermi_azimuth_fit[4] = 0.1;

  gsl_vector *azimuthal_[2] = { azimuthal_r, azimuthal_dat };


  fit1dfermi_azimuthal_neldermead (azimuthal_, fermi_azimuth_fit);


  if (VERBOSE)
    cout << endl;

  fermi_azimuth_fit[2] = fabs (fermi_azimuth_fit[2]);
  TF_az = pow (6 * f2 (fermi_azimuth_fit[1]), -0.333);
  BetaMu_az = fermi_azimuth_fit[1];

  Fugacity_az = exp (BetaMu_az);
  f_az = fq (BetaMu_az);

  r_az = fermi_azimuth_fit[2];

  T_az = p->B1 / (2 * p->kbm) * r_az * r_az * p->magnif * p->magnif / f_az;

  if (VERBOSE || p->showfermi)
    {
      printf ("..............  FERMI AZIMUTHAL FIT RESULTS ..............\n");
      printf ("n0         = %.2e\n", fermi_azimuth_fit[0]);
      printf ("BetaMu     = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	      fermi_azimuth_fit[1], f2 (fermi_azimuth_fit[1]), TF_az);
      printf ("r          = %.1f pixels\n", fermi_azimuth_fit[2]);
      printf ("offset     = %.1f pixels\n", fermi_azimuth_fit[3]);
      printf ("mx         = %.1f pixels\n", fermi_azimuth_fit[4]);
      printf ("fugacity   = %.2e\n", Fugacity_az);
      printf ("f(z) = %.2e\n", f_az);
      printf ("T_az = %.2f uK\n", T_az);
      printf ("TF   = %.2f uK\n", TF);
    }


  n0_az = fermi_azimuth_fit[0];
  B_az = fermi_azimuth_fit[3];
  mx_az = fermi_azimuth_fit[4];

  return;
}


void
Fermions::MakeAzimuthalPlots ()
{
  if (VERBOSE)
    cout << endl << "------------ MAKE AZIMUTHAL PLOTS ------------" << endl;

  azimuthal_gauss = gsl_vector_alloc (nbins);
  azimuthal_fermi2d = gsl_vector_alloc (nbins);
  azimuthal_fermi = gsl_vector_alloc (nbins);

  double r, model_gauss, model_fermi2d, model_fermi;
  for (unsigned int index = 0; index < usedbins; index++)
    {
      r = index * binsize;

      if (p->fermi2d)
	{
	  model_fermi2d =
	    n0 / f1 (BetaMu) * f1 (BetaMu -
				   fq (BetaMu) * pow (r,
						      2) / rj_Fermi /
				   ri_Fermi) + B_Fermi;
	}
      else
	{
	  model_fermi2d = 0.;
	}

      if (p->fermiazimuth)
	{
	  model_fermi =
	    n0_az / f1 (BetaMu_az) * f1 (BetaMu_az -
					 fq (BetaMu_az) * pow (r / r_az,
							       2)) +
	    B_az + mx_az * r;
	}
      else
	{
	  model_fermi = 0.;
	}

      model_gauss = offset + peak * exp (-1. * pow (r, 2) / wi_1e / wj_1e);


      gsl_vector_set (azimuthal_gauss, index, model_gauss);
      gsl_vector_set (azimuthal_fermi2d, index, model_fermi2d);
      gsl_vector_set (azimuthal_fermi, index, model_fermi);
    }

  gsl_vector *azimuthal_avg[5] =
    { azimuthal_r, azimuthal_dat, azimuthal_gauss, azimuthal_fermi2d,
    azimuthal_fermi
  };

  to_dat_file (azimuthal_avg, 5, p->shotnum, "azimuthal.dat");


  //--- On screen plot
  std::ofstream gpl;
  gpl.open ("temp.gpl");
  gpl << "set xr [0:512]" << endl;
  gpl << "plot \"" << p->shotnum << "_azimuthal.dat\" ";
  gpl << "u 1:2 pt 7 ps 0.5 title 'data' ";
  gpl << ",\\" << endl << "\"\" u 1:3 w lines title 'gauss' ";
  gpl << ",\\" << endl << "\"\" u 1:4 w lines title 'fermi' ";
  gpl << ",\\" << endl << "\"\" u 1:5 w lines title 'fermi-azimuthal' ";

  if (p->plots)
    std::system ("gnuplot -persist temp.gpl");
  gpl.close ();
  remove ("temp.gpl");

  //--- PNG plot
  gpl.open ("temp2.gpl");
  gpl << "set xr [0:512]" << endl;
  gpl << "set terminal png medium" << endl;
  gpl << "set output \"azimuthal_" << p->shotnum << ".png\"" << endl;
  gpl << "plot \"" << p->shotnum << "_azimuthal.dat\" ";
  gpl << "u 1:2 pt 7 ps 0.5 title 'data' ";
  gpl << ",\\" << endl << "\"\" u 1:3 w lines title 'gauss' ";
  gpl << ",\\" << endl << "\"\" u 1:4 w lines title 'fermi' ";
  gpl << ",\\" << endl << "\"\" u 1:5 w lines title 'fermi-azimuthal' ";
  gpl.close ();
  std::system ("gnuplot temp2.gpl");
  return;
}


void
Fermions::ComputeRadialAxialDensity ()
{

  unsigned int s1 = columndensity->size1;
  unsigned int s2 = columndensity->size2;

  radialdensity = gsl_vector_alloc (s1);
  axialdensity = gsl_vector_alloc (s2);

  gsl_vector_set_all (radialdensity, 0.0);
  gsl_vector_set_all (axialdensity, 0.0);

  double cd_ij = 0.0;

  for (unsigned int i = 0; i < s1; i++)
    {
      for (unsigned int j = 0; j < s2; j++)
	{
	  cd_ij = gsl_matrix_get (columndensity, i, j);

	  gsl_vector_set (radialdensity, i,
			  gsl_vector_get (radialdensity, i) + cd_ij);
	  gsl_vector_set (axialdensity, j,
			  gsl_vector_get (axialdensity, j) + cd_ij);
	}
    }

  if (VERBOSE)
    cout << endl;

  double radialfit[4] = { ci, wi_1e, gsl_vector_max (radialdensity), 0.1 };
  fit1dgaus (radialdensity, radialfit);

  double axialfit[4] = { cj, wj_1e, gsl_vector_max (axialdensity), 0.1 };
  fit1dgaus (axialdensity, axialfit);

  double radialfit_fermi[5] =
    { gsl_vector_max (radialdensity), -5.0, wi_1e, ci, radialfit[3] };
  fit1dfermi_neldermead (radialdensity, radialfit_fermi);

  double axialfit_fermi[5] =
    { gsl_vector_max (axialdensity), -5.0, wj_1e, cj, axialfit[3] };
  fit1dfermi_neldermead (axialdensity, axialfit_fermi);

  TF_rd = pow (6 * f2 (radialfit_fermi[1]), -0.333);
  TF_ax = pow (6 * f2 (axialfit_fermi[1]), -0.333);

  if (VERBOSE)
    {
      cout << endl <<
	"------------ RADIAL AND AXIAL DENSITY GAUS FIT RESULTS ------------"
	<< endl;
      printf ("\tc_rd = %.1f\n", radialfit[0]);
      printf ("\tw_rd = %.2f\n", radialfit[1]);
      printf ("\tA_rd = %.3e\n", radialfit[2]);
      printf ("\tB_rd = %.3e\n", radialfit[3]);
      printf ("\n");

      printf ("\tc_ax = %.1f\n", axialfit[0]);
      printf ("\tw_ax = %.2f\n", axialfit[1]);
      printf ("\tA_ax = %.3e\n", axialfit[2]);
      printf ("\tB_ax = %.3e\n", axialfit[3]);

      cout << endl <<
	"------------ RADIAL AND AXIAL DENSITY FERMI FIT RESULTS ------------"
	<< endl;
      printf ("\tn0_rd     = %.1f\n", radialfit_fermi[0]);
      printf ("BetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	      radialfit_fermi[1], f2 (radialfit_fermi[1]), TF_rd);
      printf ("\tr_rd      = %.3e\n", radialfit_fermi[2]);
      printf ("\tc_rd      = %.3e\n", radialfit_fermi[3]);
      printf ("\tB_rd      = %.3e\n", radialfit_fermi[4]);
      printf ("\n");

      printf ("\tn0_ax     = %.1f\n", axialfit_fermi[0]);
      printf ("BetaMu = %.2e -->  f2(BetaMu) = %.4f, T/TF = %.2f\n",
	      axialfit_fermi[1], f2 (axialfit_fermi[1]), TF_ax);
      printf ("\tr_ax      = %.3e\n", axialfit_fermi[2]);
      printf ("\tc_ax      = %.3e\n", axialfit_fermi[3]);
      printf ("\tB_ax      = %.3e\n", axialfit_fermi[4]);
    }



  char path[MAXPATHLEN];
  getcwd (path, MAXPATHLEN);

  FILE *radial_F, *axial_F;
  string radial_path (path);
  radial_path += "/";
  radial_path += p->shotnum;
  radial_path += "_radial.ndat";
  string axial_path (path);
  axial_path += "/";
  axial_path += p->shotnum;
  axial_path += "_axial.ndat";

  radial_F = fopen (radial_path.c_str (), "w+");
  axial_F = fopen (axial_path.c_str (), "w+");

  fprintf (radial_F, "#pixel\tdat\tmodel\n");
  fprintf (axial_F, "#pixel\tdat\tmodel\n");

  double gaus_model, fermi_model, data;

  for (unsigned int i = 0; i < s1; i++)
    {
      double xi = (double) i;
      gaus_model =
	radialfit[3] +
	radialfit[2] * exp (-1.0 *
			    pow ((xi - radialfit[0]) / radialfit[1], 2.));
      fermi_model =
	radialfit_fermi[0] / f32 (radialfit_fermi[1]) *
	f32 (radialfit_fermi[1] -
	     fq (radialfit_fermi[1]) * pow ((i - radialfit_fermi[3]) /
					    radialfit_fermi[2],
					    2)) + radialfit_fermi[4];
      data = gsl_vector_get (radialdensity, i);
      fprintf (radial_F, "%e\t%e\t%e\t%e\n", xi, data, gaus_model,
	       fermi_model);
    }

  for (unsigned int j = 0; j < s2; j++)
    {
      double xj = (double) j;
      gaus_model =
	axialfit[3] +
	axialfit[2] * exp (-1.0 * pow ((xj - axialfit[0]) / axialfit[1], 2.));
      fermi_model =
	axialfit_fermi[0] / f32 (axialfit_fermi[1]) *
	f32 (axialfit_fermi[1] -
	     fq (axialfit_fermi[1]) * pow ((j - axialfit_fermi[3]) /
					   axialfit_fermi[2],
					   2)) + axialfit_fermi[4];
      data = gsl_vector_get (axialdensity, j);
      fprintf (axial_F, "%e\t%e\t%e\t%e\n", xj, data, gaus_model,
	       fermi_model);
    }

  fclose (radial_F);
  fclose (axial_F);

  // Save plots to png file
  std::ofstream gpl;
  gpl.open ("temp.gpl");
  gpl << "set terminal png medium" << endl;
  gpl << "set output \"radial_" << p->shotnum << ".png\"" << endl;
  gpl << "set size 1.0,0.45" << endl;
  gpl << "plot \"" << p->
    shotnum << "_radial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_radial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << " radial gaus\" ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_radial.ndat\" u 1:4 w lines title \"" << p->
    shotnum << " radial fermi\" " << endl;
  gpl << "set output \"axial_" << p->shotnum << ".png\"" << endl;
  gpl << "plot \"" << p->
    shotnum << "_axial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_axial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << " axial gaus\" ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_axial.ndat\" u 1:4 w lines title \"" << p->
    shotnum << " axial fermi\" " << endl;
  std::system ("gnuplot -persist temp.gpl");
  gpl.close ();

  // Show plots on screen if the plots option is selected 
  gpl.open ("temp.gpl");
  gpl << "set size 1.0,0.45" << endl;
  gpl << "set multiplot" << endl;
  gpl << "set origin 0.0,0.0" << endl;
  gpl << "plot \"" << p->
    shotnum << "_radial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_radial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << "radial gaus\" ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_radial.ndat\" u 1:4 w lines title \"" << p->
    shotnum << "radial fermi\" " << endl;
  gpl << "set origin 0.0,0.5" << endl;
  gpl << "plot \"" << p->
    shotnum << "_axial.ndat\" u 1:2 pt 7 ps 1 notit ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_axial.ndat\" u 1:3 w lines title \"" << p->
    shotnum << "axial gaus\" ,\\" << endl;
  gpl << "\"" << p->
    shotnum << "_axial.ndat\" u 1:4 w lines title \"" << p->
    shotnum << "axial fermi\" " << endl;
  gpl << "unset multiplot" << endl;
  gpl.close ();
  if (p->plots)
    std::system ("gnuplot -persist temp.gpl");

  return;
}
