/* John Moustakas, 2000 September, UofA */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/deep1/ioannis/lib/c/export.h"

IDL_LONG mlike_idl     /* IDL entry point */

(int argc, void *argv[])

{
  /* declare the IDL variables as pointers */

  float *g, *total, *logtot;
  IDL_LONG *indx, *nbins, *nmbins;

  /* declare internal variables */

  int i, j;
  int argct = 0;   /* argument counter */

  /* assign names to the IDL variables (order matters!) */
  
  g = (float *)argv[argct++];
  indx = (IDL_LONG *)argv[argct++];
  nbins = (IDL_LONG *)argv[argct++];
  nmbins = (IDL_LONG *)argv[argct++];
  total = (float *)argv[argct++];
  logtot = (float *)argv[argct++];
  
  for (i=0; i < *nbins; i++) {

    for (j=0; j < *nbins; j++) {
      *(total+i) = *(total+i) + *(g+j+i*(*nbins));}

    for (j=0; j < *nmbins; j++) {
      *(logtot+i) = *(logtot+i) + log10(*(g+*(indx+j)+i*(*nbins)));}

  }
	
  return 0;

}
