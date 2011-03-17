pro trgb_lfunction, mags, merr, lf, binsize=binsize, minmag=minmag, maxmag=maxmag
;+
; NAME:
;	TRGB_LFUNCTION
;
; PURPOSE:
;	Generate a luminosity function.
;
; INPUTS:
;	mags	: magnitude array
;	merr	: magnitude error array (for completeness)
;
; OPTIONAL INPUTS:
;	binsize	: histogram bin size
;	minmag	: starting luminosity function magnitude
;	maxmag	: ending luminosity function magnitude
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	lf	: luminosity function structure
;	  .mags    : magnitude array
;	  .merr    : magnitude error array
;	  .binsize : luminosity function bin size
;	  .minmag  : starting luminosity function magnitude
;	  .maxmag  : ending luminosity function magnitude
;	  .hist    : luminosity function (histogram)
;	  .nhist   : number of luminosity function bins
;	  .mag     : luminosity function ordinate (magnitudes)
;
; COMMON BLOCKS:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	John Moustakas, 4 August 2000, UCB
;-

	if not keyword_set(binsize) then binsize = 0.5
	if not keyword_set(minmag) then minmag = min(mags)
	if not keyword_set(maxmag) then maxmag = max(mags)
        
        hist = histogram(mags,binsize=binsize,min=minmag,max=maxmag)
        nhist = n_elements(hist)

        xnorm = max(findgen(nhist))/(maxmag-minmag)
        xmag = findgen(nhist)/xnorm+minmag   ; magnitude array

; pack everything into a structure

        lf = {mags:	float(mags), $
              merr:	float(merr), $
              binsize:	float(binsize), $
              minmag:	float(minmag), $
              maxmag:	float(maxmag), $
              hist:	long(hist), $
              nhist:	long(nhist), $
              mag:	float(xmag)}
        
return
end
