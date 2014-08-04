function bcgmstar_sersic2_allbands_func, rr, pp, parinfo=parinfo, wave=wave
; fit all the bands simultaneously 

    common com_sersic2_allbands, nlookup, blookup
    
; pp = [n1_ref,n2_ref,re1_ref,re2_ref,alpha1,alpha2,beta1,beta2,sbe]
; sbe - surface brightness at re
; re - effective (half-light) radius
; n - Sersic index

; n1(wave) = n1_ref(wave/wave_ref)^alpha1
; n2(wave) = n2_ref(wave/wave_ref)^alpha2
; re1(wave) = re1_ref(wave/wave_ref)^beta1
; re2(wave) = re2_ref(wave/wave_ref)^beta2

; figure out how many bands we're fitting
;   uwave = wave[uniq(wave,sort(wave))]
    uwave = wave[uniq(wave,reverse(sort(wave)))] 
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

; provide a simple look-up table of the Sersic b parameter on a grid
; of Sersic n (see Graham & Driver (2005), equations 1 and 4 and my
; get_sersicb() function); here we use the analytic approximations of
; Ciotti & Bertin for n>0.36 and MacArthur et al. 2003 for smaller n
    if n_elements(blookup) eq 0 then begin
       nlookup = range(0.1D,20D,500,/log)
       blookup = nlookup*0.0
       hi = where(nlookup gt 0.36,comp=lo)
       blookup[hi] = 2*nlookup[hi] - 1D/3D + 4D/(405*nlookup[hi]) + 46D/(25515D*nlookup[hi]^2) + $
         131D/(1148175D*nlookup[hi]^3) - 2194697D/(30690717750D*nlookup[hi]^4) ; +order(n^-5)
       blookup[lo] = poly(nlookup[lo],[0.01945D,-0.8902D,10.95D,-19.67D,13.43D])
;      plot, nlookup, blookup, /xlog, /ylog
    endif

    use_b1 = interpol(blookup,nlookup,use_n1)
    use_b2 = interpol(blookup,nlookup,use_n2)

; fit the sum of two Sersic models     
    for ib = 0, nband-1 do begin
       ww = where(uwave[ib] eq wave)
       model1 = use_sbe1[ib]*exp(-use_b1[ib]*((rr[ww]/use_re1[ib])^(1D/use_n1[ib])-1D)) + $
         use_sbe2[ib]*exp(-use_b2[ib]*((rr[ww]/use_re2[ib])^(1D/use_n2[ib])-1D))
;      model1 = use_sbe1[ib]*exp(-get_sersicb(use_n1[ib])*((rr[ww]/use_re1[ib])^(1D/use_n1[ib])-1D)) + $
;        use_sbe2[ib]*exp(-get_sersicb(use_n2[ib])*((rr[ww]/use_re2[ib])^(1D/use_n2[ib])-1D))
       if ib eq 0 then model = model1 else model = [model,model1]
    endfor

return, model
end

