function fit_mz, mass, oh, weight, pivotmass=pivotmass, $
  binsize=binsize, minmass=minmass, maxmass=maxmass, npoly=npoly, $
  mingal=mingal, verbose=verbose
; jm09mar19nyu - fit a generic mass-metallicity relation
    
    ngal = n_elements(mass)
    if (n_elements(pivotmass) eq 0) then pivotmass = 0.0 ; median(mass)
    if (n_elements(binsize) eq 0) then binsize = 0.1
    if (n_elements(minmass) eq 0) then minmass = min(mass)
    if (n_elements(maxmass) eq 0) then maxmass = max(mass)
    if (n_elements(npoly) eq 0) then npoly = 2
    if (n_elements(mingal) eq 0) then mingal = 100

; let IM_MEDXBIN() handle the unweighted case    
    run = im_medxbin(mass,oh,binsize,weight=weight,$
      minx=minmass,maxx=maxmass,minpts=mingal,$
      verbose=verbose)
    coeff = poly_fit(run.xbin-pivotmass,run.medy,npoly,sigma=coeff_err)
;   im_poly_iter, run.xbin-pivotmass, run.medy, npoly, coeff=coeff, $
;     e_coeff=coeff_err, nsig=3.0
;   coeff = robust_poly_fit(run.xbin-pivotmass,run.medy,npoly)
;   coeff_err = coeff*0.0
    scatter = djsig(oh-poly(mass-pivotmass,coeff),sigrej=4.0) ; 1-sigma

;   massaxis = im_array(min(run.xbin),max(run.xbin),0.05)
;   massfit = poly(massaxis,coeff)

    fit = {$
      ngal:         ngal, $
      coeff:        reform(coeff), $
      coeff_err:    reform(coeff_err), $
      scatter:      scatter, $
;     massaxis:     massaxis, $
;     massfit:      massfit, $
      minmass:      minmass, $
      maxmass:      maxmass, $
      pivotmass:    pivotmass, $
      binsize:      binsize, $
      bin_mass:     run.xbin, $
      bin_oh:       run.medy}
    
return, fit
end

