function read_bcgmstar_sersic, cluster, results=results, radius=radius, $
  model=model, firstsersic=firstsersic, secondsersic=secondsersic, $
  phot=phot, band=band, covar=covar
; jm14aug12siena - read the appropriate Sersic-fitting results, which
; varies from cluster to cluster
    
; note that the individual structures cannot be 'stacked' because of a
; varying number of bandpasses (the covariance matrices are all
; different sizes)    
    
    if n_elements(cluster) eq 0 then begin
       splog, 'CLUSTER input required'
       return, -1
    endif

    sersicpath = bcgmstar_path(/sersic)
    
; originally I used the wavelength-independent model for every cluster
; except macs2129, macs1149, rxj2129, and macs1720 (and obviously a209
; and a2261); however, further tests showed very little improvement in
; the reduced chi2 or in the visual quality of the fits, so I moved to
; the wavelength-independent model for every cluster    
    suffix = '-alphabetazero'
    
;    case cluster of
;; double-Sersic 
;       'a209': suffix = '-alphabetazero'
;       'a2261': suffix = '-alphabetazero'
;; single-Sersic 
;       'macs2129': suffix = '-alphabetazero'
;       'macs1149': suffix = '-alphabetazero'
;       'rxj2129': suffix = '-alphabetazero'
;       'macs1720': suffix = '-alphabetazero'
;       else: suffix = '-alphabetazero'
;    endcase

    sersic = mrdfits(sersicpath+cluster+'-allsersic'+suffix+'.fits.gz',1,/silent)
    if arg_present(results) or arg_present(covar) then results = mrdfits(sersicpath+cluster+$
      '-allsersic'+suffix+'-results.fits.gz',1,/silent)

; get the photometry if requested
    if arg_present(phot) then phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
    
; chose one or more bandpasses
    if n_elements(band) ne 0 then begin
       match, strlowcase(strtrim(sersic.band,2)), strlowcase(band), m1, m2
       if n_elements(m2) ne n_elements(band) or m1[0] eq -1 then begin
          splog, 'Unknown bandpass(es): '+strjoin(band,',')
          return, -1
       endif
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       sersic = sersic[m1]
       if arg_present(phot) then phot = phot[m1]
    endif

; extract the covariance matrix, if requested; need to eliminate the
; columns that refer to the alpha1 and beta1 parameters (which are
; held fixed at zero) 
    if arg_present(covar) then covar = (results.sersic_covar[[0,1,lindgen(12)+4],*])[*,[0,1,lindgen(12)+4]]    
    
; evaluate the model if requested
    if n_elements(radius) ne 0 then begin
       nfilt = n_elements(sersic)
       model = fltarr(n_elements(radius),nfilt)
       firstsersic = model*0
       secondsersic = model*0
       for ib = 0, nfilt-1 do begin
; double-Sersic
          if cluster eq 'a209' or cluster eq 'a2261' or cluster eq 'rxj2248' then begin
             model[*,ib] = bcgmstar_sersic2_func(radius,params=sersic[ib],/allbands)
             firstsersic[*,ib] = bcgmstar_sersic_func(radius,[sersic[ib].sersic2_all_sbe1,$
               sersic[ib].sersic2_all_re1,sersic[ib].sersic2_all_n1])
             secondsersic[*,ib] = bcgmstar_sersic_func(radius,[sersic[ib].sersic2_all_sbe2,$
               sersic[ib].sersic2_all_re2,sersic[ib].sersic2_all_n2])
          endif else begin
; single-Sersic
             model[*,ib] = bcgmstar_sersic_func(radius,params=sersic[ib],/allbands)
          endelse 
       endfor 
    endif 

return, sersic
end
    
