function bcgmstar_sersic2_func, rr, pp, params=params, parinfo=parinfo, allbands=allbands
; pp = [sbe1,re1,n1,sbe2,re2,n2]

    common com_sersic2, nlookup, blookup

; support input from the bcgmstar_sersic2_allbands_func function; see
; bcgmstar_sersicfit, /qaplot_sbprofiles for proper usage
    if keyword_set(allbands) then begin
       if n_elements(params) ne 0 then pp = $
         [params.sersic2_all_sbe1,params.sersic2_all_re1,params.sersic2_all_n1,$
         params.sersic2_all_sbe2,params.sersic2_all_re2,params.sersic2_all_n2]
    endif else begin
       if n_elements(params) ne 0 then pp = $
         [params.sersic2_sbe1,params.sersic2_re1,params.sersic2_n1,$
         params.sersic2_sbe2,params.sersic2_re2,params.sersic2_n2]
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

    use_sbe1 = 10D^(-0.4*pp[0])
    use_re1 = pp[1]
    use_n1 = pp[2]
    use_sbe2 = 10D^(-0.4*pp[3])
    use_re2 = pp[4]
    use_n2 = pp[5]

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
    model = use_sbe1*exp(-use_b1*((rr/use_re1)^(1D/use_n1)-1D)) + $
      use_sbe2*exp(-use_b2*((rr/use_re2)^(1D/use_n2)-1D))
;   model = use_sbe1*exp(-get_sersicb(use_n1)*((rr/use_re1)^(1D/use_n1)-1D)) + $
;     use_sbe2*exp(-get_sersicb(use_n2)*((rr/use_re2)^(1D/use_n2)-1D))
    model = -2.5*alog10(model)

return, model
end

