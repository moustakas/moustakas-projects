/* John Moustakas, 2000 September 1, UofA */

/* generate a Sakai, Madore & Freedman (SMF) (gaussian smoothed)
   luminosity function */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/deep1/ioannis/lib/c/export.h"
#include "/deep1/ioannis/lib/c/nrutil.h"

IDL_LONG smf_idl     /* IDL entry point */

(int argc, void *argv[])

{
  void smf_lf_();

  /* declare the IDL variables as pointers */

  int i, j;
  int argct = 0;
  float *mags, *merr, *marray, coef, exparg, *total;
  IDL_LONG *nbins, *nstars, *phi;
  double pi = 4.*atan(1.);

  /* assign names to the IDL variables (order matters!) */

  mags = (float *)argv[argct++];
  merr = (float *)argv[argct++];
  marray = (float *)argv[argct++];
  nbins = (IDL_LONG *)argv[argct++];
  nstars = (IDL_LONG *)argv[argct++];
  phi = (IDL_LONG *)argv[argct++];

  /* allocate memory to the new variables */
  
  total = vector(0,*nbins-1);

  for (i=0; i < *nbins; i++)   /* loop on each bin */
    {

      total[i] = 0.;
      
      /* gaussian-smoothed luminosity function (see SMF 1996 Appendix A) */
      
      for (j=0; j < *nstars; j++) /* loop on each star */
	{

	  coef = 1./(sqrt(2.*pi)*merr[j]);      /* Gaussian normalization */
	  exparg = (mags[j]-marray[i])/merr[j]; /* exponent argument */

	  total[i] = total[i] + coef*exp(-0.5*(pow(2.,exparg)));

	}
      
      printf("%f %f %f \n",coef,exparg,total[i]);
      phi[i] = total[i];

    };

  return 0;

}
