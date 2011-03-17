pro boot_arrays, iter, lf, bootmags, booterrs
;+
; NAME: 
;	BOOT_ARRAYS
;
; PURPOSE:
;	Generate a uniformly resampled list of magnitudes and
;	magnitude errors.
;
; INPUTS:
;	iter : number of times to resample the starlist
;	lf   : stellar luminosity function (see TRGB_LFUNCTION)
;
; OUTPUTS:
;	bootmags : resampled magnitudes with dimensions (nstars,iter)
;	booterrs : resampled magnitude errors with dimensions (nstars,iter)
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 September 1
;-

	nstars = n_elements(lf.mags) ; total number of stars

        bootindx = floor(randomu(seed,nstars,iter)*nstars) ; uniform resampling

; perturb the new magnitudes by a Gaussian whose width is the magnitude error

        bootmags = lf.mags[bootindx] + lf.merr[bootindx]*randomn(nseed,nstars,iter)
        booterrs = lf.merr[bootindx]

return
end
