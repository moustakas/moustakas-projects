;
; NAME:
;   SG1120_TWEAK_ZEROPOINTS
;
; PURPOSE:
;   Tweak the SG1120 photometric zeropoints by iteratively fitting
;   K-correct templates and finding the average zeropoint tweaks that
;   minimize chi2 for all objects.
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Aug 03, UCSD - written
;-

function tweak_mpfunc, maggies, zpt_tweak, ngal=ngal
    zpt_tweak1 = cmreplicate(zpt_tweak,ngal)
return, maggies + zpt_tweak1
end
    
pro sg1120_tweak_zeropoints

; start with the zeropoints derived by comparing with the SDSS
; photometry
    filterlist = sg1120_filterlist()
    nband = n_elements(filterlist)

;   zpt = [-0.298,-0.125,+0.085,+0.345,+0.235,0.0]
    zpt = [-0.275,-0.116,+0.097,+0.364,+0.243,0.0]
    use_zpt = zpt

    for ii = 0, 10 do begin

; compute K-corrections
       sg1120_kcorrect, sg1120, rr, use_zpt=use_zpt, /tweak;, /nowrite
       ngal = n_elements(rr)

       if (ii eq 0) then chi2 = fltarr(ngal,10)
       chi2[*,ii] = rr.chi2
       
; for each filter, derive the linear tweak to the flux that minimizes
; chi2 over the full set of galaxies
       parinfo = replicate({value: 1E-9, fixed: 0},nband)
       flux_tweak = mpfitfun('tweak_mpfunc',rr.maggies,rr.bestmaggies,$
         weights=rr.ivarmaggies,parinfo=parinfo,functargs={ngal: ngal},$
         quiet=1)
       
;      old = total(rr.ivarmaggies*(rr.maggies-rr.bestmaggies)^2,1,/double)
;      new = total(rr.ivarmaggies*(rr.maggies-rr.bestmaggies+cmreplicate(flux_tweak,ngal))^2,1,/double)

; convert the flux tweaks into the median zeropoint correction in
; *magnitudes*
       zpt_tweak = use_zpt*0.0
       for jj = 0, nband-1 do zpt_tweak[jj] = median(2.5*alog10(1.0+flux_tweak[jj]/rr.maggies[jj]))
;      for jj = 0, nband-1 do zpt_tweak[jj] = weighted_quantile($
;        flux_tweak[jj]/rr.maggies[jj],rr.ivarmaggies[jj])
       
       use_zpt = zpt + zpt_tweak ; note sign!
       splog, 'Iteration ', ii
       splog, '####################'
       niceprint, filterlist, zpt, zpt_tweak, use_zpt

;; compute the weighted average       
;       zpt_tweak = use_zpt*0.0
;       for jj = 0, nband-1 do zpt_tweak[jj] = $
;         total(rr.ivarmaggies[jj,*]*rr.maggies[jj,*]*rr.bestmaggies[jj,*],/double)/$
;         total(rr.ivarmaggies[jj,*]*rr.bestmaggies[jj,*]^2.0,/double)
;;      for jj = 0, nband-1 do zpt_tweak[jj] = $
;;        weighted_quantile(rr.maggies[jj]/rr.bestmaggies[jj],rr.ivarmaggies[jj])

;;     use_zpt = zpt + 2.5*alog10(zpt_tweak) ; note sign!
;      splog, 'Iteration ', ii
;      splog, '####################'
;      niceprint, filterlist, zpt, 2.5*alog10(zpt_tweak), use_zpt
stop
    endfor

stop    
    
return
end
    
