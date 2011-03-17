/* John Moustakas, 2000 September, UofA */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/deep1/ioannis/lib/c/export.h"
#include "/deep1/ioannis/lib/c/nrutil.h"

IDL_LONG mlike_idl     /* IDL entry point */

(int argc, void *argv[])

{
  /* declare the IDL variables as pointers */

  float *marray, *cutarray, *alpha, *beta, *gamma, *lnlike;
  IDL_LONG *nbins, *nmbins, *nstars, *ncut, *indx;

  /* declare internal variables */

  int i, j, k;
  float cut, temp, *total, *logtot, **g, **x;
  int argct = 0;   /* argument counter */

  /* assign names to the IDL variables (order matters!) */
  
  nbins = (IDL_LONG *)argv[argct++];
  marray = (float *)argv[argct++];
  cutarray = (float *)argv[argct++];

  alpha = (float *)argv[argct++];
  beta = (float *)argv[argct++];
  gamma = (float *)argv[argct++];

  nstars = (IDL_LONG *)argv[argct++];
  ncut = (IDL_LONG *)argv[argct++];
  indx = (IDL_LONG *)argv[argct++];
  nmbins = (IDL_LONG *)argv[argct++];
  lnlike = (float *)argv[argct++];

  /* allocate memory to the internal variables */
  
  g = matrix(0,*nbins-1,0,*nbins-1);  /* three-piece stellar distribution function */
  x = matrix(0,*nbins-1,0,*nbins-1);  /* 2D magnitude array */
  total = vector(0,*nbins-1);         /* array sum */
  logtot = vector(0,*nbins-1);        /* log_10(array) sum */

  /* generate the 2D magnitude array */

  for (i=0; i < *nbins; i++)
    {
      for (j=0; j < *nbins; j++)
	{
	  x[j][i] = marray[j]-marray[i];  /* shifts the TRGB origin */
	  /* x[j][i] = *(marray+i*(*nbins)+j); */
	};
    };

  /* main part of the program */
  
  for (i=0; i < *ncut; i++)  /* loop on the discontinuity width */
    {
      cut = cutarray[i];

      for (j=0; j < *nbins; j++) /* loop on each magnitude array */
	{
	  total[j] = 0.; 
	  logtot[j] = 0.;

	  for (k=0; k < *nbins; k++) /* for a fixed magnitude array, loop on each magnitude */
	    {
	      /* evaluate the three-piece distribution function */

	      if (x[k][j] >= 0) g[k][j] = pow(10.,(*alpha*x[k][j]));
	      if ((x[k][j] < 0) && (x[k][j] >= -cut/(*gamma))) g[k][j] = pow(10.,(*gamma*x[k][j]));
	      if (x[k][j] < -cut/(*gamma)) g[k][j] = pow(10.,-cut)*pow(10.,(*beta*(x[k][j]+cut/(*gamma))));
	      
	      total[j] = total[j] + g[k][j];
	    };

	  for (k=0; k < *nmbins; k++) /* sum over the log of the binned magnitudes */
	    {
	      logtot[j] = logtot[j] + log10(g[indx[k]][j]);
	    };
	
	  *(lnlike+j+i*(*nbins)) = - (*nstars) * log10(total[j]) + logtot[j];
	};
    };
  
  /*  for (k=0; k < *nbins; k++)
      {printf("%f %f %f \n",-(*nstars)*log10(total[k]),logtot[k],(-(*nstars)*log10(total[k]))+logtot[k]);}*/

  free_matrix(g, 0, *nbins-1, 0, *nbins-1);
  free_matrix(x, 0, *nbins-1, 0, *nbins-1);
  free_vector(total, 0, *nbins-1);
  free_vector(logtot, 0, *nbins-1);
	
  return 0;

}
