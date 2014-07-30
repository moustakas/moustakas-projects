function bcgmstar_sersic2_allbands_func, rr, pp, parinfo=parinfo, wave=wave
; fit all the bands simultaneously 
; 
; pp = [n1_ref,n2_ref,re1_ref,re2_ref,alpha1,alpha2,beta1,beta2,sbe]
; sbe - surface brightness at re
; re - effective (half-light) radius
; n - Sersic index

; n1(wave) = n1_ref(wave/wave_ref)^alpha1
; n2(wave) = n2_ref(wave/wave_ref)^alpha2
; re1(wave) = re1_ref(wave/wave_ref)^beta1
; re2(wave) = re2_ref(wave/wave_ref)^beta2

; figure out how many bands we're fitting
    uwave = wave[uniq(wave,sort(wave))]
;   uwave = wave[uniq(wave,reverse(sort(wave)))]
    nband = n_elements(uwave)

    get_element, uwave, 15326.0, this
    wave_ref = uwave[this] ; =F160W

; make sure the parameters don't go outside the boundaries,
; since mpfit does not 
    use_pp = pp
    if n_elements(parinfo) ne 0 then begin
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[0] then $
         use_pp[ii] = use_pp[ii] > parinfo[ii].limits[0]
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[1] then $
         use_pp[ii] = use_pp[ii] < parinfo[ii].limits[1]
    endif

    use_n1_ref = 10D^use_pp[0]
    use_n2_ref = 10D^use_pp[1]
    use_re1_ref = 10D^use_pp[2]
    use_re2_ref = 10D^use_pp[3]

    use_alpha1 = use_pp[4]
    use_alpha2 = use_pp[5]
    use_beta1 = use_pp[6]
    use_beta2 = use_pp[7]

    use_n1 = use_n1_ref*(uwave/wave_ref)^use_alpha1
    use_n2 = use_n2_ref*(uwave/wave_ref)^use_alpha2
    use_re1 = use_re1_ref*(uwave/wave_ref)^use_beta1
    use_re2 = use_re2_ref*(uwave/wave_ref)^use_beta2

    use_sbe1 = use_pp[8:8+nband-1]
    use_sbe2 = use_pp[8+nband:8+2*nband-1]

; fit the sum of two Sersic models     
    for ib = 0, nband-1 do begin
       ww = where(uwave[ib] eq wave)
       model1 = use_sbe1[ib]*exp(-get_sersicb(use_n1[ib])*((rr[ww]/use_re1[ib])^(1D/use_n1[ib])-1D)) + $
         use_sbe2[ib]*exp(-get_sersicb(use_n2[ib])*((rr[ww]/use_re2[ib])^(1D/use_n2[ib])-1D))
       if ib eq 0 then model = model1 else model = [model,model1]
    endfor

return, model
end

