function sings_lzfit, absmag, oh, pivot=pivot, $
  minabsmag=minabsmag, maxabsmag=maxabsmag
; fit the luminosity-metallicity relation
    
    lindx = 2 ; SIXLIN least-squares bisector    
    if (n_elements(pivot) eq 0L) then pivot = -18.0

    good = where((oh gt -900.0) and (absmag gt -900.0),ngood)
    if (n_elements(minabsmag) eq 0) then minabsmag = min(absmag[good])
    if (n_elements(maxabsmag) eq 0) then maxabsmag = max(absmag[good])
    
    these = where((absmag[good] le maxabsmag) and (absmag[good] ge minabsmag),number)
    sixlin, absmag[good[these]]-pivot, oh[good[these]], a, siga, b, sigb
    coeff = [a[lindx],b[lindx]]
    coeff_err = [siga[lindx],sigb[lindx]]

    scatter = djsig(oh[good[these]]-poly(absmag[good[these]]-pivot,coeff),sigrej=3.0)

    new_coeff = [coeff[0]-coeff[1]*pivot,coeff[1]]
    new_coeff_err = coeff_err
;   new_coeff_err = [sqrt(coeff_err[0]^2+(coeff_err[1]*pivot)^2),coeff_err[1]]
    
    fit = {$
      ngal:          ngood,$
      number:        number,$
      coeff:         coeff, $
      coeff_err:     coeff_err, $
      new_coeff:     new_coeff, $
      new_coeff_err: new_coeff_err, $
      pivot:         pivot, $
      scatter:       scatter}

return, fit
end

