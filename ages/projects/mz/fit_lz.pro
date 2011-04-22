function fit_lz, mag, oh, weight, oh_err=oh_err, binsize=binsize, $
  brightmag=brightmag, faintmag=faintmag, mingal=mingal, verbose=verbose
; jm09mar19nyu - fit a generic luminosity-metallicity relation;
;   BRIGHTMAG should be the least negative number (the brightest) 
; jm10oct13ucsd - updated
    
    ngal = n_elements(mag)

; compute the binned LZ relation; let IM_MEDXBIN() deal with WEIGHT
; undefined
    lzbin = im_medxbin(mag,oh,binsize,weight=weight,$
      minx=brightmag,maxx=faintmag,minpts=mingal,$
      verbose=verbose)
    
; now fit an ordinary least-squares bisector
    pivotmag = lz_pivotmag()
    if (n_elements(brightmag) eq 0L) then brightmag = min(mag)
    if (n_elements(faintmag) eq 0L) then faintmag = max(mag)
    if (n_elements(weight) eq 0L) then weight = mag*0.0+1.0

    sixlin, mag-pivotmag, oh, a, siga, b, sigb, weight=weight;/oh_err[these]^2
;   these = where((mag lt faintmag) and (mag gt brightmag),number)
;   sixlin, mag[these]-pivotmag, oh[these], a, $
;     siga, b, sigb, weight=weight[these];/oh_err[these]^2

    lindx = 2L ; ordinary least-squares bisector
    coeff = [a[lindx],b[lindx]]
    coeff_err = [siga[lindx],sigb[lindx]]
    scatter = djsig(oh-poly(mag-pivotmag,coeff),sigrej=4.0) ; 1-sigma
    
    fit = {$
      ngal:       ngal, $
      coeff:      coeff, $
      coeff_err:  coeff_err, $
      scatter:    scatter,$
      brightmag:  brightmag,$
      faintmag:   faintmag,$
      mingal:     mingal,$
      bin_mag:    lzbin.xbin,$
      bin_oh:     lzbin.medy,$
      bin_oh_err: lzbin.sigymean}

return, fit
end

