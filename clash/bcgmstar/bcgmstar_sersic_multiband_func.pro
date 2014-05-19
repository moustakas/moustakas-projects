function bcgmstar_sersic_multiband_func, rr, pp, parinfo=parinfo, wave=wave
; fit multiple bands simultaneously 
; 
; pp = [sbe,re,n]
; sbe - surface brightness at re
; re - effective (half-light) radius
; n - Sersic index

; figure out how many bands we're fitting
    uwave = wave[uniq(wave,reverse(sort(wave)))]
    nband = n_elements(uwave)
    
; make sure the parameters don't go outside the boundaries,
; since mpfit does not 
    use_pp = pp
    if n_elements(parinfo) ne 0 then begin
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[0] then $
         use_pp[ii] = use_pp[ii] > parinfo[ii].limits[0]
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[1] then $
         use_pp[ii] = use_pp[ii] < parinfo[ii].limits[1]
    endif

    use_re = use_pp[0] ; = 1.0/Re (see Graham & Driver 2005)
    use_n = use_pp[1]
    use_sbe = use_pp[2:2+nband-1]

; see equation 6 in Graham & Driver 2005; this expression is
; -2.5*alog10() of the Sersic model
    for ib = 0, nband-1 do begin
       ww = where(uwave[ib] eq wave)
       model1 = use_sbe[ib] + (2.5/alog(10))*get_sersicb(use_n)*((rr[ww]/use_re)^(1D/use_n)-1D)
       if ib eq 0 then model = model1 else model = [model,model1]
    endfor

return, model
end

