/* John Moustakas, 2000 September 1, UofA */
/* generate SMF LFs and response functions for bootstrap resampled starlists */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "export.h"

IDL_LONG bootsmf_idl

(int argc,
 void * argv[])

{
  /* declare variables */

  IDL_LONG *iter;
  float *bootmags, *booterrs, *marray;
  IDL_LONG *nbins, *smfhist;
  float *boothist, *bootresplin, *bootresplog, *trgbmaglin, *trgbmaglog;
  float *resplin, *resplog;
  float *mags, merr;

  int i, j;
  int argct = 0;

  iter = (IDL_LONG *)argv[argct++];
  bootmags = (float *)argv[argct++];
  booterrs = (float *)argv[argct++];
  marray = (float *)argv[argct++];
  nbins = (IDL_LONG *)argv[argct++];
  boothist = (float *)argv[argct++];
  bootresplin = (float *)argv[argct++];
  bootresplog = (float *)argv[argct++];
  trgbmaglin = (float *)argv[argct++];
  trgbmaglog = (float *)argv[argct++];
  
  void smf_lf(float *mags, float *merr, IDL_LONG *marray, IDL_LONG *nbins, 
	      IDL_LONG *smfhist, float *resplin, float *resplog);

  mags = (float *) malloc (sizeof (float) * (size_t)nbins);
  merr = (float *) malloc (sizeof (float) * (size_t)nbins);

  for (i=0; i < *iter; i++)
    {

      for (j=0; j < *nbins; j++)
	{
	  mags[j] = bootmags[j,i]; /* bootstrapped starlist */
	  merr[j] = booterrs[j,i]; /* bootstrapped magnitude errors */
	}

      /* generate the new LF and response functions */

      smf_lf_(mags,merr,marray,nbins,smfhist,resplin,resplog);

      for (j=0; j < *nbins; j++)
	{
	  boothist[j,i] = smfhist[j];
      
	  bootresplin[j,i] = resplin[j];
	  bootresplog[j,i] = resplog[j];
	  
	  /* constrain the endpoints of the response */

	  if ((j = 0) || (j = 1)) bootresplin[j,i] = 0.;
	  if ((j = *nbins-1) || (j = *nbins-1)) bootresplin[j,i] = 0.;

	  if ((j = 0) || (j = 1)) bootresplog[j,i] = 0.;
	  if ((j = *nbins-1) || (j = *nbins-1)) bootresplog[j,i] = 0.;

	  /* search for the maximum between minmag and the peak of the LF */
	  
	  if (resplin[j] = max(resplin)) trgbmaglin[i] = marray[j];
	  if (resplog[j] = max(resplog)) trgbmaglog[i] = marray[j];

	}
    }

  return 0;

}
