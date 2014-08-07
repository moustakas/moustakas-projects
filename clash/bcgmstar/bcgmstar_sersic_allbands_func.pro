function bcgmstar_sersic_allbands_func, rr, pp, parinfo=parinfo, wave=wave
; fit all the bands simultaneously 

    common com_sersic_allbands, nlookup, blookup
    
; pp = [n_ref,re_ref,alpha1,beta1,sbe]
; sbe - surface brightness at re
; re - effective (half-light) radius
; n - Sersic index

; n(wave) = n_ref(wave/wave_ref)^alpha1
; re(wave) = re_ref(wave/wave_ref)^beta1

; figure out how many bands we're fitting
    uwave = wave[uniq(wave,reverse(sort(wave)))] 
;   uwave = wave[uniq(wave,sort(wave))]
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

    use_n_ref = use_pp[0]
    use_re_ref = use_pp[1]
;   use_n_ref = 10D^use_pp[0]
;   use_re_ref = 10D^use_pp[1]

    use_alpha1 = use_pp[2]
    use_beta1 = use_pp[3]

    use_n = use_n_ref*(uwave/wave_ref)^use_alpha1
    use_re = use_re_ref*(uwave/wave_ref)^use_beta1
    use_sbe = use_pp[4:4+nband-1]

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

    use_b = interpol(blookup,nlookup,use_n)>0

; see equation 6 in Graham & Driver 2005; this expression is
; -2.5*alog10() of the Sersic model
    for ib = 0, nband-1 do begin
       ww = where(uwave[ib] eq wave)

;      model1 = use_sbe[ib]*exp(-use_b[ib]*((rr[ww]/use_re[ib])^(1D/use_n[ib])-1D))
       model1 = use_sbe[ib] + (2.5/alog(10))*use_b[ib]*((rr[ww]/use_re[ib])^(1D/use_n[ib])-1D)

       if ib eq 0 then model = model1 else model = [model,model1]
       if total(finite(model) eq 0) ne 0 then stop
    endfor

return, model
end

