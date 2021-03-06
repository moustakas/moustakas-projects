function bcgmstar_sersic_func, rr, pp, params=params, parinfo=parinfo, allbands=allbands
; pp = [sbe,re,n]
; sbe - surface brightness at re
; re - effective (half-light) radius
; n - Sersic index

    common com_sersic, nlookup, blookup

; support input from the bcgmstar_sersic_allbands_func function; see
; bcgmstar_sersicfit, /qaplot_sbprofiles for proper usage
    if keyword_set(allbands) then begin
       if n_elements(params) ne 0 then pp = $
         [params.sersic_all_sbe,params.sersic_all_re,params.sersic_all_n]
    endif else begin
       if n_elements(params) ne 0 then pp = $
         [params.sersic_sbe,params.sersic_re,params.sersic_n]
    endelse

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

    use_b = interpol(blookup,nlookup,use_n)

; see equation 6 in Graham & Driver 2005; this expression is
; -2.5*alog10() of the Sersic model
    model = use_sbe + (2.5/alog(10))*use_b*((rr/use_re)^(1D/use_n)-1D)
;   model = use_sbe + (2.5/alog(10))*get_sersicb(use_n)*((rr/use_re)^(1D/use_n)-1D)
    
;   model = alog(use_sb0)-use_k*rr^(1D/use_n)
;   model = use_sb0*exp(-(use_k*rr)^(1D/use_n))
;   model = use_sb0*exp(-get_sersicb(use_n)*((rr/use_k)^(1D/use_n)-1))
;   model = alog(use_sb0)-get_sersicb(use_n)*((rr/use_k)^(1D/use_n)-1D)

return, model
end

