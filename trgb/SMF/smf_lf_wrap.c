/* John Moustakas, 2000 September 1, UofA */
/* C-wrapper for calling SMF_LF; see SMF_LF for more documentation */

#include <stdio.h>
#include "/deep1/ioannis/lib/c/export.h"

IDL_LONG smf_wrap_idl     /* IDL entry point */

(int argc, void *argv[])

{
  void smf_lf_();

  /* declare the IDL variables as pointers */

  float *mags, *merr, *marray, *resplin, *resplog;
  IDL_LONG *nbins, *smfhist;

  /* declare internal variables */

  int i;
  int argct = 0;

  /* assign names to the IDL variables (order matters!) */

  mags = (float *)argv[argct++];
  merr = (float *)argv[argct++];
  marray = (float *)argv[argct++];
  nbins = (IDL_LONG *)argv[argct++];
  nstars = (IDL_LONG *)argv[argct++];
  smfhist = (IDL_LONG *)argv[argct++];
  resplin = (float *)argv[argct++];
  resplog = (float *)argv[argct++];

  /* generate the luminosity function */

  smf_lf_(mags,merr,marray,nbins,nstars,smfhist,resplin,resplog);

}




