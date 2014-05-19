function bcgmstar_sersic_func, rr, pp, params=params, parinfo=parinfo
; pp = [sbe,re,n]
; sbe - surface brightness at re
; re - effective (half-light) radius
; n - Sersic index

    if n_elements(params) ne 0 then pp = $
      [params.sersic_sbe,params.sersic_re,params.sersic_n]

; make sure the parameters don't go outside the boundaries,
; since mpfit does not 
    use_pp = pp
    if n_elements(parinfo) ne 0 then begin
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[0] then $
         use_pp[ii] = use_pp[ii] > parinfo[ii].limits[0]
       for ii = 0, n_elements(pp)-1 do if parinfo[ii].limited[1] then $
         use_pp[ii] = use_pp[ii] < parinfo[ii].limits[1]
    endif

    use_sbe = use_pp[0]
    use_re = use_pp[1] ; = 1.0/Re (see Graham & Driver 2005)
    use_n = use_pp[2]

; see equation 6 in Graham & Driver 2005; this expression is
; -2.5*alog10() of the Sersic model
    model = use_sbe + (2.5/alog(10))*get_sersicb(use_n)*((rr/use_re)^(1D/use_n)-1D)
    
;   model = alog(use_sb0)-use_k*rr^(1D/use_n)
;   model = use_sb0*exp(-(use_k*rr)^(1D/use_n))
;   model = use_sb0*exp(-get_sersicb(use_n)*((rr/use_k)^(1D/use_n)-1))
;   model = alog(use_sb0)-get_sersicb(use_n)*((rr/use_k)^(1D/use_n)-1D)

return, model
end

