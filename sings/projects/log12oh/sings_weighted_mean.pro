function sings_weighted_mean, values, errors, wmean_err=wmean_err, $
  wsigma=wsigma, wquantile=wquantile, quant=quant
; jm10jul18ucsd - 
    
    npts = n_elements(values)
    if (npts eq 0L) or (n_elements(errors) eq 0L) then begin
       doc_library, 'sings_weighted_mean'
       return, -1.0
    endif

    if (n_elements(errors) ne npts) then begin
       splog, 'Dimensions of VALUES and ERRORS do not agree.'
       return, -1.0
    endif

    bad = where(errors le 0,nbad)
    if (nbad ne 0L) then message, 'Bad errors'
    weights = 1.0/errors^2

; compute the weighted mean and weighted standard deviation; also
; calculate the error on mean from the distribution itself (talk to
; E. Sheldon's about where this comes from)
    wtot = total(weights,double=double)
    wmean = total(weights*values,double=double)/wtot 

    if arg_present(wmean_err) then wmean_err = sqrt(1.0/wtot)
    if arg_present(wsigma) then wsigma = $
      sqrt(total(weights*(values-wmean)^2,double=double)/wtot)

; finally, optionally compute the weighted median
    if arg_present(wquantile) then wquantile = $
      weighted_quantile(values,weights,quant=quant)
    
return, wmean
end
