function bcgmstar_sersic2_multiband_func, rr, pp, parinfo=parinfo, wave=wave
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

    use_re1 = use_pp[0] ; = 1.0/Re (see Graham & Driver 2005)
    use_n1 = use_pp[1]
    use_re2 = use_pp[2] ; = 1.0/Re (see Graham & Driver 2005)
    use_n2 = use_pp[3]

    use_sbe1 = use_pp[4:4+nband-1]
    use_sbe2 = use_pp[4+nband:4+2*nband-1]

; fit the sum of two Sersic models     
    for ib = 0, nband-1 do begin
       ww = where(uwave[ib] eq wave)
       model1 = use_sbe1[ib]*exp(-get_sersicb(use_n1)*((rr[ww]/use_re1)^(1D/use_n1)-1D)) + $
         use_sbe2[ib]*exp(-get_sersicb(use_n2)*((rr[ww]/use_re2)^(1D/use_n2)-1D))
       if ib eq 0 then model = model1 else model = [model,model1]
    endfor

return, model
end

